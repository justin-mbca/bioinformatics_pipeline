import os
import pandas as pd
import importlib
import sys
import os

# Ensure repo root is on sys.path so 'workflows' package can be imported
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)



def test_clinical_pipeline_tmp(tmp_path, capsys):
    # create a small raw clinical CSV
    raw = tmp_path / "raw_clinical.csv"
    df = pd.DataFrame({
        'sample': ['s1', 's2'],
        'age': [50, 60],
        'outcome': [1, 0]
    })
    df.to_csv(raw, index=False)

    cleaned = tmp_path / "cleaned.csv"
    analysis = tmp_path / "analysis.csv"

    clinical_main_args = [str(raw), str(cleaned), str(analysis)]

    # run main by importing and calling
    clinical_module = importlib.import_module('workflows.clinical_pipeline')
    import sys
    old_argv = sys.argv
    sys.argv = ["clinical_pipeline.py"] + clinical_main_args
    try:
        clinical_module.main()
    finally:
        sys.argv = old_argv

    assert os.path.exists(cleaned)
    assert os.path.exists(analysis)

    df_clean = pd.read_csv(cleaned)
    assert list(df_clean['sample']) == ['s1', 's2']
