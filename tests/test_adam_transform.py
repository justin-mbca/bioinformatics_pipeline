import pandas as pd
import sys
import os

# Ensure repo root is on sys.path for imports
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)

from workflows.clinical.adam_transform import create_adam_demographics


def test_create_adam_demographics_basic():
    rows = [
        {'USUBJID':'s1','AGE':30,'HEIGHT':170,'WEIGHT':70},
        {'USUBJID':'s2','AGE':65,'HEIGHT':160,'WEIGHT':80}
    ]
    df = create_adam_demographics(rows)
    assert 'BMI' in df.columns
    assert 'AGEGRP' in df.columns
    assert df.loc[0,'BMI'] is not None
    assert df.loc[1,'AGEGRP'] == '>60'
