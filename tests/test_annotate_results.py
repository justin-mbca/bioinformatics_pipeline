import pandas as pd
import tempfile
import os
from rna_seq.annotate_results import annotate


def write_csv(df, path, index=True):
    df.to_csv(path, index=index)


def test_annotate_with_index_gene_ids(tmp_path):
    # create a DE results with gene IDs in the index
    de = pd.DataFrame({'log2FoldChange': [1.2, -0.5], 'padj': [0.01, 0.2]}, index=['G1', 'G2'])
    ann = pd.DataFrame({'gene_id': ['G1', 'G2'], 'symbol': ['GeneOne', 'GeneTwo']})

    de_fp = tmp_path / 'de.csv'
    ann_fp = tmp_path / 'ann.csv'
    out_fp = tmp_path / 'out.csv'

    write_csv(de, de_fp, index=True)
    write_csv(ann, ann_fp, index=False)

    annotate(str(de_fp), str(ann_fp), str(out_fp))

    out = pd.read_csv(out_fp)
    assert 'symbol' in out.columns
    assert 'G1' in out['gene_id'].values


def test_annotate_with_gene_header(tmp_path):
    # create DE results where first column is 'gene'
    de = pd.DataFrame({'gene': ['G1', 'G2'], 'log2FoldChange': [1.2, -0.5], 'padj': [0.01, 0.2]})
    ann = pd.DataFrame({'gene': ['G1', 'G2'], 'symbol': ['GeneOne', 'GeneTwo']})

    de_fp = tmp_path / 'de2.csv'
    ann_fp = tmp_path / 'ann2.csv'
    out_fp = tmp_path / 'out2.csv'

    write_csv(de, de_fp, index=False)
    write_csv(ann, ann_fp, index=False)

    annotate(str(de_fp), str(ann_fp), str(out_fp))

    out = pd.read_csv(out_fp)
    # ensure symbol column present
    assert 'symbol' in out.columns


def test_annotate_with_explicit_merge_key(tmp_path):
    # create DE results and annotation that use 'entrez' as the common key
    de = pd.DataFrame({'entrez': ['E1', 'E2'], 'log2FoldChange': [0.5, -1.1], 'padj': [0.05, 0.2]})
    ann = pd.DataFrame({'entrez': ['E1', 'E2'], 'symbol': ['One', 'Two']})

    de_fp = tmp_path / 'de3.csv'
    ann_fp = tmp_path / 'ann3.csv'
    out_fp = tmp_path / 'out3.csv'

    de.to_csv(de_fp, index=False)
    ann.to_csv(ann_fp, index=False)

    # call annotate with explicit merge key
    annotate(str(de_fp), str(ann_fp), str(out_fp), merge_key='entrez')

    out = pd.read_csv(out_fp)
    assert 'symbol' in out.columns
    assert 'E1' in out['entrez'].values
