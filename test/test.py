import requests
import pandas as pd
import json

SERVER_URL = "http://localhost.ucsd.edu:4857"

def test_api():
    url = "{}/api/smart3/search".format(SERVER_URL)

    # Preparing request
    peaks_df = pd.read_csv("./swinolide.csv", sep=None)
    data = {}
    data["peaks"] = json.dumps(peaks_df.to_dict(orient="records"))
    r = requests.post(url, data=data)
    
    r.raise_for_status()

def test_api2():
    url = "{}/api/smart3/search".format(SERVER_URL)

    # Preparing request
    peaks_df = pd.read_csv("./Calcd_42086.csv", sep=None)
    data = {}
    data["peaks"] = json.dumps(peaks_df.to_dict(orient="records"))
    r = requests.post(url, data=data)
    
    r.raise_for_status()