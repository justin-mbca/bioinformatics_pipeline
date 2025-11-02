#!/usr/bin/env python3
"""Local integrator script.
Usage: integrate_genomics_clinical.py genes.csv clinical_summary.csv combined_report.csv
"""
import sys
import pandas as pd
from difflib import SequenceMatcher


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Integrate gene results with clinical summary and optional scRNA markers.')
    parser.add_argument('genes', help='path to genes CSV (DE results)')
    parser.add_argument('clinical', help='path to clinical summary CSV')
    parser.add_argument('out', help='output combined report CSV')
    parser.add_argument('--markers', help='optional markers CSV (scRNA markers)')
    parser.add_argument('--mapping', help='optional mapping CSV with columns source_id,canonical_symbol')
    args = parser.parse_args()

    genes = pd.read_csv(args.genes)
    clinical = pd.read_csv(args.clinical)
    out = args.out
    markers_path = args.markers
    mapping_path = args.mapping
    # Build standardized actionable report schema
    # Preferred columns: gene_id, symbol, log2FoldChange, padj, clinical_association, cohort_frequency, note
    report = pd.DataFrame()

    # gene_id - prefer human-readable name columns before numeric index-like columns
    gene_id_col = None
    for candidate in ['gene', 'Unnamed: 0', 'annotation', 'symbol', 'gene_id']:
        if candidate in genes.columns:
            gene_id_col = candidate
            break
    if gene_id_col is None:
        # pick the first non-numeric column if possible
        for c in genes.columns:
            if genes[c].dtype == object:
                gene_id_col = c
                break
    if gene_id_col is None:
        gene_id_col = genes.columns[0]

    report['gene_id'] = genes[gene_id_col].astype(str)

    # symbol - try genes, then annotation, then external mapping file
    if 'symbol' in genes.columns:
        report['symbol'] = genes['symbol']
    elif 'annotation' in genes.columns:
        report['symbol'] = genes['annotation']
    else:
        # attempt to look up gene symbol from local annotation file
        try:
            ann = pd.read_csv('workflows_local/gene_annotation.csv')
            # map columns 'gene' or 'gene_id' -> symbol/annotation
            key_col = None
            if 'gene' in ann.columns:
                key_col = 'gene'
            elif 'gene_id' in ann.columns:
                key_col = 'gene_id'
            if key_col:
                mapping = dict(zip(ann[key_col].astype(str), ann.get('annotation', ann.get('symbol', ann[key_col])).astype(str)))
                report['symbol'] = report['gene_id'].map(mapping).fillna(report['gene_id'])
            else:
                report['symbol'] = report['gene_id']
        except Exception:
            report['symbol'] = report['gene_id']

    # effect size and p-values
    report['log2FoldChange'] = genes.get('log2FoldChange', pd.NA)
    report['padj'] = genes.get('padj', pd.NA)

    # clinical association: try to pick a representative clinical column name
    clinical_assoc = None
    if clinical.shape[1] > 0:
        clinical_assoc = clinical.columns[0]
    report['clinical_association'] = clinical_assoc

    # cohort_frequency: if clinical provides a count column, use it; else placeholder
    if 'clinical_count' in genes.columns:
        # clinical_count might already be present after earlier merge
        report['cohort_frequency'] = genes['clinical_count']
    elif 'count' in clinical.columns:
        total = clinical['count'].sum() if clinical['count'].dtype.kind in 'biufc' else None
        if total and total > 0:
            report['cohort_frequency'] = clinical['count'] / total
        else:
            report['cohort_frequency'] = pd.NA
    else:
        report['cohort_frequency'] = pd.NA

    # note column for actionable recommendation (heuristics)
    report['note'] = pd.NA
    try:
        padj_vals = pd.to_numeric(report.get('padj', pd.Series([pd.NA]*len(report))), errors='coerce')
        l2fc_vals = pd.to_numeric(report.get('log2FoldChange', pd.Series([pd.NA]*len(report))), errors='coerce')
        notes = []
        for p, l in zip(padj_vals, l2fc_vals):
            if pd.notna(p) and p < 0.05:
                notes.append('validate by qPCR')
            elif pd.notna(l) and abs(l) > 1.0:
                notes.append('high priority')
            else:
                notes.append(pd.NA)
        report['note'] = notes
    except Exception:
        report['note'] = pd.NA

    report.to_csv(out, index=False)
    print(f"Wrote integrated report to {out}")

    # If markers provided, append marker info (cluster, avg_log2FC) when gene names match
    if markers_path:
        try:
            markers = pd.read_csv(markers_path)
            # markers likely have columns 'gene' and 'avg_log2FC' or 'avg_log2FC'
            mcols = markers.columns
            gene_col = 'gene' if 'gene' in mcols else markers.columns[-1]
            score_col = None
            for c in ['avg_log2FC', 'avg_log2FC', 'avg_logFC', 'avg_log2FC']:
                if c in mcols:
                    score_col = c
                    break
            cluster_col = 'cluster' if 'cluster' in mcols else None

            # Map markers to report by matching gene symbol or gene_id
            marker_info = markers[[gene_col] + ([cluster_col] if cluster_col else []) + ([score_col] if score_col else [])].copy()
            marker_info = marker_info.rename(columns={gene_col: 'gene_match'})
            # Try matching on symbol first, then gene_id
            report['marker_cluster'] = pd.NA
            report['marker_avg_log2FC'] = pd.NA
            # If a mapping file is provided, load it and create canonical symbols
            mapping = None
            if mapping_path:
                try:
                    mapping = pd.read_csv(mapping_path)
                    # validate mapping columns
                    if not ({'source_id', 'canonical_symbol'} <= set(mapping.columns)):
                        raise ValueError('Mapping file must contain columns: source_id, canonical_symbol')
                    map_dict = dict(zip(mapping['source_id'].astype(str), mapping['canonical_symbol'].astype(str)))
                    # normalize report symbol/gene_id to canonical
                    report['canonical'] = report['symbol'].fillna(report['gene_id']).astype(str).map(map_dict).fillna(report['symbol'].fillna(report['gene_id']).astype(str))
                except Exception as e:
                    print(f'Warning: failed to load/validate mapping file {mapping_path}: {e}')
                    mapping = None
                    report['canonical'] = report['symbol'].fillna(report['gene_id']).astype(str)
            if mapping is None:
                report['canonical'] = report['symbol'].fillna(report['gene_id']).astype(str)
            # build lower-case lookup series for conservative fuzzy matching
            sym_series = report['canonical'].astype(str).fillna('').str.lower()
            gid_series = report['gene_id'].astype(str).fillna('').str.lower()

            for i, row in marker_info.iterrows():
                gm_raw = row['gene_match']
                # if mapping provided, map marker gene to canonical before matching
                if mapping is not None and 'source_id' in mapping.columns and 'canonical_symbol' in mapping.columns:
                    gm = mapping.set_index('source_id')['canonical_symbol'].to_dict().get(str(gm_raw), str(gm_raw))
                else:
                    gm = str(gm_raw)
                gm = gm.lower()

                # 1) exact case-insensitive match on symbol
                matches = sym_series == gm
                # 2) fallback: exact on gene_id
                if not matches.any():
                    matches = gid_series == gm
                # 3) fallback: substring match on symbol
                if not matches.any() and 'symbol' in report.columns:
                    matches = sym_series.str.contains(gm, na=False)
                # 4) fallback: substring match on gene_id
                if not matches.any():
                    matches = gid_series.str.contains(gm, na=False)

                idx = []
                # If we have any straightforward matches, use them
                if matches.any():
                    idx = report.index[matches].tolist()
                else:
                    # aggressive fallback: try best fuzzy match using SequenceMatcher
                    combined = pd.concat([sym_series, gid_series])
                    best_score = 0.0
                    best_pos = None
                    for pos, cand in enumerate(combined):
                        try:
                            score = SequenceMatcher(None, gm, str(cand)).ratio()
                        except Exception:
                            score = 0.0
                        if score > best_score:
                            best_score = score
                            best_pos = pos
                    # conservative threshold to avoid false positives
                    if best_score >= 0.70 and best_pos is not None:
                        if best_pos < len(sym_series):
                            idx = [sym_series.index[best_pos]]
                        else:
                            idx = [gid_series.index[best_pos - len(sym_series)]]

                # apply matched info to report rows
                for j in idx:
                    if cluster_col:
                        report.at[j, 'marker_cluster'] = row.get(cluster_col, pd.NA)
                    if score_col:
                        report.at[j, 'marker_avg_log2FC'] = row.get(score_col, pd.NA)
            report.to_csv(out, index=False)
            print(f"Appended marker info from {markers_path} and updated {out}")
        except Exception as e:
            print(f"Failed to append markers: {e}")


if __name__ == '__main__':
    main()
