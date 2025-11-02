"""Prefect flow for clinical ingestion and ADaM transforms.
This example shows how to orchestrate steps: ingest -> validate -> transform -> save.
"""
from prefect import flow, task
import pandas as pd
import json
from typing import List, Dict

@task
def ingest_sdtm(path: str) -> List[Dict]:
    df = pd.read_csv(path)
    return df.to_dict(orient='records')

@task
def validate_sdtm(rows: List[Dict]) -> List[Dict]:
    from workflows.clinical.schemas import validate_sdtm_dataset
    return validate_sdtm(rows)

@task
def transform_adam(rows: List[Dict]) -> str:
    from workflows.clinical.adam_transform import create_adam_demographics
    adsl = create_adam_demographics(rows)
    out = 'clinical/adsl.csv'
    adsl.to_csv(out, index=False)
    return out

@flow
def clinical_workflow(sdtm_path: str):
    rows = ingest_sdtm(sdtm_path)
    valid = validate_sdtm(rows)
    adm_out = transform_adam(valid)
    print(f"ADaM demographics written to {adm_out}")

if __name__ == '__main__':
    clinical_workflow('clinical/raw_clinical.csv')
