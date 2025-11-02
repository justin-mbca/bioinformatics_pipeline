"""Simple ADaM-like transformation helpers: create analysis datasets from SDTM-like inputs.
This is intentionally small and illustrative; production ADaM creation is complex and study-specific.
"""
import pandas as pd
from typing import List, Dict


def create_adam_demographics(sdtm_dem: List[Dict]) -> pd.DataFrame:
    """Create a minimal ADaM-like demographics table from SDTM DM dataset.
    - compute BMI if HEIGHT and WEIGHT are present
    - categorize age groups
    """
    df = pd.DataFrame(sdtm_dem)
    if 'HEIGHT' in df.columns and 'WEIGHT' in df.columns:
        # HEIGHT in cm, WEIGHT in kg -> BMI
        df['BMI'] = df.apply(lambda r: r['WEIGHT'] / ((r['HEIGHT']/100.0)**2) if pd.notnull(r['HEIGHT']) and pd.notnull(r['WEIGHT']) else None, axis=1)
    df['AGEGRP'] = pd.cut(df['AGE'], bins=[0,40,60,120], labels=['<=40','41-60','>60'])
    # Ensure ADaM-style column names
    df = df.rename(columns={'USUBJID':'SUBJID'})
    return df


def create_adam_survival(adsl: pd.DataFrame, events: List[Dict]) -> pd.DataFrame:
    # Basic example: join enrollment to event table and compute time-to-event
    events_df = pd.DataFrame(events)
    merged = events_df.merge(adsl, left_on='USUBJID', right_on='SUBJID', how='left')
    # assume event has EVENT_DATE and ENROLL_DATE
    merged['TTE'] = (pd.to_datetime(merged['EVENT_DATE']) - pd.to_datetime(merged['ENROLL_DATE'])).dt.days
    return merged
