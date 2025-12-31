"""
Kafka Producer for Medical Device Data
Generates and sends medical device data to Kafka topics
"""
import json
import random
import time
from datetime import datetime
from kafka import KafkaProducer


def generate_vital_signs(patient_id):
    """
    Generate simulated vital signs data

    Args:
        patient_id: Patient identifier

    Returns:
        Dictionary with vital signs data
    """
    return {
        "patient_id": patient_id,
        "timestamp": datetime.now().isoformat(),
        "device_id": f"DEV-{random.randint(1000, 9999)}",
        "vital_signs": {
            "heart_rate": random.randint(60, 100),
            "blood_pressure_systolic": random.randint(110, 140),
            "blood_pressure_diastolic": random.randint(70, 90),
            "temperature": round(random.uniform(36.5, 37.5), 1),
            "oxygen_saturation": random.randint(95, 100)
        }
    }


def send_to_kafka(producer, topic, message):
    """
    Send message to Kafka topic

    Args:
        producer: KafkaProducer instance
        topic: Topic name
        message: Message to send
    """
    producer.send(topic, value=message)


if __name__ == "__main__":
    # Example usage
    producer = KafkaProducer(
        bootstrap_servers=['localhost:9092'],
        value_serializer=lambda v: json.dumps(v).encode('utf-8')
    )

    for i in range(10):
        message = generate_vital_signs(f"patient-{i}")
        send_to_kafka(producer, 'medical-devices-vitals', message)
        print(f"Sent message for patient-{i}")
        time.sleep(1)

    producer.flush()
    producer.close()
