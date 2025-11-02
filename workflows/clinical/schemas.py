"""Lightweight SDTM-like schema validation with optional pydantic support.
This module avoids importing heavy dependencies at import time so tests can run
in minimal environments. If pydantic is available it will be used for stricter
validation; otherwise a permissive runtime check is performed.
"""
from typing import Dict, Any, List


def _has_pydantic():
    try:
        import pydantic  # type: ignore
        return True
    except Exception:
        return False


def validate_subject_record(rec: Dict[str, Any]) -> Dict[str, Any]:
    """Validate a single SDTM subject record. Returns cleaned record or raises ValueError.
    Expected minimal fields: USUBJID (subject id), AGE, SEX
    """
    if _has_pydantic():
        from pydantic import BaseModel, ValidationError

        class Subject(BaseModel):
            USUBJID: str
            AGE: float
            SEX: str

        try:
            subj = Subject(**rec)
            return subj.dict()
        except ValidationError as ve:
            raise ValueError(str(ve))
    else:
        # Light validation
        if 'USUBJID' not in rec:
            raise ValueError('Missing USUBJID')
        if 'AGE' not in rec:
            raise ValueError('Missing AGE')
        if 'SEX' not in rec:
            rec['SEX'] = 'U'
        # Coerce types simply
        try:
            rec['AGE'] = float(rec['AGE'])
        except Exception:
            raise ValueError('AGE must be numeric')
        rec['USUBJID'] = str(rec['USUBJID'])
        return rec


def validate_sdtm_dataset(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    cleaned = []
    for r in rows:
        cleaned.append(validate_subject_record(r))
    return cleaned
