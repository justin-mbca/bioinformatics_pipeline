import os
import pandas as pd
import importlib
import sys
import os

# Ensure repo root is on sys.path so 'workflows' package can be imported
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)


def test_integration_tmp(tmp_path):
    genes = tmp_path / "genes.csv"
    clinical = tmp_path / "clinical.csv"
    out = tmp_path / "out.csv"

    pd.DataFrame({'gene':['G1','G2'],'score':[1.2,3.4]}).to_csv(genes,index=False)
    pd.DataFrame({'cohort_metric':[0.1]}).to_csv(clinical,index=False)

    mod = importlib.import_module('workflows.integrate_genomics_clinical')
    mod.main_args = [str(genes), str(clinical), str(out)]

    # call main by patching sys.argv
    import sys
    old = sys.argv
    sys.argv = ['integrate_genomics_clinical.py', str(genes), str(clinical), str(out)]
    try:
        mod.main()
    finally:
        sys.argv = old

    assert os.path.exists(out)
    df = pd.read_csv(out)
    assert 'gene' in df.columns
    assert 'clinical_cohort_metric' in df.columns
