"""Airflow DAG skeleton for clinical processing.
This module avoids importing Airflow at top-level to keep test environment light. To enable
in a real Airflow environment, replace the import guard with normal DAG/operators.
"""
from datetime import datetime

def make_dag(dag_id='clinical_pipeline', schedule='@daily'):
    try:
        from airflow import DAG
        from airflow.operators.python import PythonOperator
    except Exception:
        # Airflow not installed; return a lightweight description instead
        return {
            'dag_id': dag_id,
            'schedule': schedule,
            'tasks': ['ingest', 'validate', 'transform']
        }

    def ingest():
        print('ingest')

    def validate():
        print('validate')

    def transform():
        print('transform')

    dag = DAG(dag_id=dag_id, start_date=datetime(2025,1,1), schedule_interval=schedule)

    t1 = PythonOperator(task_id='ingest', python_callable=ingest, dag=dag)
    t2 = PythonOperator(task_id='validate', python_callable=validate, dag=dag)
    t3 = PythonOperator(task_id='transform', python_callable=transform, dag=dag)

    t1 >> t2 >> t3
    return dag
