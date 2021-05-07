import sys
import json
import requests
import numpy as np
import pandas as pd

sys.path.insert(0, "Classifier")
import SMART_3
import get_highlight


# Loading database into memory
DB, index_super = SMART_3.load_db(db_folder="./Classifier")

def search_smart3(nmr_data_df, output_image=None, channel=1):
    nmr_mat = SMART_3._convert_data(nmr_data_df, channel)

    # Running prediction here
    query_dict = {}
    query_dict["input_1"] = nmr_mat.tolist()

    # # Handling SUPERCLASS
    if channel == 1:
        pred_url = "http://smart3-tf-server:8501/v1/models/CHANNEL1:predict"
    else:
        pred_url = "http://smart3-tf-server:8501/v1/models/CHANNEL2:predict"
    
    payload = json.dumps({"instances": [ query_dict ]})

    headers = {"content-type": "application/json"}
    json_response = requests.post(pred_url, data=payload, headers=headers)
    
    prediction_dict = json.loads(json_response.text)['predictions'][0]

    # TODO: Will update names when hyunwoo updates model
    fingerprint_prediction = np.asarray(prediction_dict["dense"])
    pred_MW = np.asarray(prediction_dict["dense_1"])
    pred_class = np.asarray(prediction_dict["class"])
    pred_gly = np.asarray(prediction_dict["glycoside"])

    fingerprint_prediction, pred_MW, pred_class_index, pred_class_prob, pred_gly = SMART_3.refine_predict_nmr([fingerprint_prediction], [pred_MW], [pred_class], [pred_gly])

    experimental_mw = None
    top_candidates = 10

    topK = None
    topK = SMART_3.search_database(fingerprint_prediction, pred_MW, DB, mw=experimental_mw, top_candidates=top_candidates)

    if output_image is not None:
        SMART_3._draw_search_image(output_image, nmr_mat, pred_class_index, pred_class_prob, pred_gly, index_super=index_super)

    return topK