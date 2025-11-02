import sys
import subprocess
import pandas as pd
import os

script = os.path.join(os.path.dirname(__file__), '..', 'workflows_local', 'integrate_genomics_clinical.py')


def run_script(python, genes, clinical, out, markers=None, mapping=None):
    cmd = [python, script, str(genes), str(clinical), str(out)]
    if markers:
        cmd += ['--markers', str(markers)]
    if mapping:
        cmd += ['--mapping', str(mapping)]
    subprocess.run(cmd, check=True)


def test_mapping_valid(tmp_path):
    python = sys.executable
    genes = tmp_path / "genes.csv"
    clinical = tmp_path / "clinical.csv"
    markers = tmp_path / "markers.csv"
    mapping = tmp_path / "mapping.csv"
    out = tmp_path / "out.csv"

    pd.DataFrame({'gene_id':['GeneA','GeneB'],'log2FoldChange':[0.1,-0.2],'padj':[0.1,0.2]}).to_csv(genes,index=False)
    pd.DataFrame({'count':[10]}).to_csv(clinical,index=False)
    pd.DataFrame({'gene':['Gene145','Gene383'],'avg_log2FC':[7.02,6.01],'cluster':[0,0]}).to_csv(markers,index=False)
    pd.DataFrame({'source_id':['Gene145','Gene383'],'canonical_symbol':['GeneA','GeneB']}).to_csv(mapping,index=False)

    run_script(python, genes, clinical, out, markers=markers, mapping=mapping)

    df = pd.read_csv(out)
    assert 'canonical' in df.columns
    assert df['marker_cluster'].notna().sum() >= 1


def test_mapping_invalid_columns(tmp_path):
    python = sys.executable
    genes = tmp_path / "genes.csv"
    clinical = tmp_path / "clinical.csv"
    markers = tmp_path / "markers.csv"
    mapping = tmp_path / "bad_mapping.csv"
    out = tmp_path / "out2.csv"

    pd.DataFrame({'gene_id':['G1','G2'],'log2FoldChange':[0.1,-0.2],'padj':[0.1,0.2]}).to_csv(genes,index=False)
    pd.DataFrame({'count':[10]}).to_csv(clinical,index=False)
    pd.DataFrame({'gene':['M1'],'avg_log2FC':[1.2],'cluster':[1]}).to_csv(markers,index=False)
    # bad mapping missing required columns
    pd.DataFrame({'foo':['M1'],'bar':['G1']}).to_csv(mapping,index=False)

    run_script(python, genes, clinical, out, markers=markers, mapping=mapping)

    df = pd.read_csv(out)
    # mapping invalid -> canonical column still present but mapping ignored
    assert 'canonical' in df.columns


def test_no_mapping(tmp_path):
    python = sys.executable
    genes = tmp_path / "genes.csv"
    clinical = tmp_path / "clinical.csv"
    markers = tmp_path / "markers.csv"
    out = tmp_path / "out3.csv"

    pd.DataFrame({'gene_id':['G1','G2'],'log2FoldChange':[0.1,-0.2],'padj':[0.1,0.2],'symbol':['G1sym','G2sym']}).to_csv(genes,index=False)
    pd.DataFrame({'count':[10]}).to_csv(clinical,index=False)
    pd.DataFrame({'gene':['M1'],'avg_log2FC':[1.2],'cluster':[1]}).to_csv(markers,index=False)

    run_script(python, genes, clinical, out, markers=markers, mapping=None)

    df = pd.read_csv(out)
    assert 'canonical' in df.columns
