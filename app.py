# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
import plotly.express as px
from dash.dependencies import Input, Output
import os
from zipfile import ZipFile
import urllib.parse
from flask import Flask, request
import json
import pandas as pd
import requests
import numpy as np
import urllib.parse

import sys
sys.path.insert(0, "Classifier")
import SMART_3
import get_highlight


server = Flask(__name__)
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server


NAVBAR = dbc.Navbar(
    children=[
        dbc.NavbarBrand(
            html.Img(src="https://gnps-cytoscape.ucsd.edu/static/img/GNPS_logo.png", width="120px"),
            href="https://gnps.ucsd.edu"
        ),
        dbc.Nav(
            [
                dbc.NavItem(dbc.NavLink("SMART NMR 3", href="#")),
            ],
        navbar=True)
    ],
    color="light",
    dark=False,
    sticky="top",
)

DASHBOARD = [
    dcc.Location(id='url', refresh=False),
    dbc.CardHeader(html.H5("SMART NMR 3")),
    dbc.CardBody(
        [
            html.Div(id='version', children="Version - 0.1"),
            html.Br(),
            dbc.Textarea(className="mb-3", id='query_text', placeholder="NMR Spectrum", rows="20"),
            dcc.Loading(
                id="structure",
                children=[html.Div([html.Div(id="loading-output-5")])],
                type="default",
            ),
            dcc.Loading(
                id="classification_table",
                children=[html.Div([html.Div(id="loading-output-3")])],
                type="default",
            )
        ]
    )
]

BODY = dbc.Container(
    [
        dbc.Row([dbc.Col(dbc.Card(DASHBOARD)),], style={"marginTop": 30}),
    ],
    className="mt-12",
)

app.layout = html.Div(children=[NAVBAR, BODY])


# This enables parsing the URL to shove a task into the qemistree id
@app.callback(Output('query_text', 'value'),
              [Input('url', 'pathname')])
def display_page(pathname):
    # Otherwise, lets use the url
    if len(pathname) > 1:
        return pathname[1:]
    else:
        return dash.no_update

# This function will rerun at any 
@app.callback(
    [Output('classification_table', 'children'), Output('structure', 'children')],
    [Input('query_text', 'value')],
)
def handle_query(query_text):
    # Saving input to a file to be read
    from io import StringIO
    data = StringIO(query_text)
    nmr_data_df = pd.read_csv(data, sep=",")
    nmr_mat = SMART_3._convert_data(nmr_data_df, 1)

    print(nmr_mat, file=sys.stderr)

    # Running prediction here
    query_dict = {}
    query_dict["input_1"] = nmr_mat.tolist()

    # # Handling SUPERCLASS
    pred_url = "http://smart3-tf-server:8501/v1/models/CHANNEL1:predict"
    payload = json.dumps({"instances": [ query_dict ]})

    headers = {"content-type": "application/json"}
    json_response = requests.post(pred_url, data=payload, headers=headers)

    prediction_dict = json.loads(json_response.text)['predictions'][0]

    # TODO: Will update names when hyunwoo updates model
    fingerprint_prediction = np.asarray(prediction_dict["dense"])
    pred_MW = np.asarray(prediction_dict["dense"])
    pred_class = np.asarray(prediction_dict["dense_1"])
    pred_gly = np.asarray(prediction_dict["glycoside"])

    return [[str(pred_MW), str(pred_class), str(pred_gly)], dash.no_update]

    isglycoside, class_results, superclass_results, pathway_results, path_from_class, path_from_superclass, n_path, fp1, fp2 = classify_structure(smiles_string)

    output_list = []

    for result in pathway_results:
        output_dict = {}
        output_dict["type"] = "pathway"
        output_dict["entry"] = result
        output_list.append(output_dict)


    for result in superclass_results:
        output_dict = {}
        output_dict["type"] = "superclass"
        output_dict["entry"] = result
        output_list.append(output_dict)


    for result in class_results:
        output_dict = {}
        output_dict["type"] = "class"
        output_dict["entry"] = result
        output_list.append(output_dict)

    if isglycoside:
        output_dict = {}
        output_dict["type"] = "is glycoside"
        output_dict["entry"] = "glycoside"
        output_list.append(output_dict)

    #Creating Table
    white_list_columns = ["type", "entry"]
    table_fig = dash_table.DataTable(
        columns=[
            {"name": i, "id": i, "deletable": False, "selectable": True} for i in white_list_columns
        ],
        data=output_list,
        editable=False,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_deletable=False,
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        page_current= 0,
        page_size= 10,
    )

    # Creating Structure Image
    img_obj = html.Img(id='image', src="https://gnps-structure.ucsd.edu/structureimg?smiles={}".format(urllib.parse.quote(smiles_string)))

    return [table_fig, img_obj]


def classify_structure(smiles):
    isglycoside = fingerprint_handler._isglycoside(smiles)

    fp = fingerprint_handler.calculate_fingerprint(smiles, 2)

    fp1 = fp[0].tolist()[0]
    fp2 = fp[1].tolist()[0]

    query_dict = {}
    query_dict["input_3"] = fp1
    query_dict["input_4"] = fp2

    # Handling SUPERCLASS
    fp_pred_url = "http://npclassifier-tf-server:8501/v1/models/SUPERCLASS:predict"
    payload = json.dumps({"instances": [ query_dict ]})

    headers = {"content-type": "application/json"}
    json_response = requests.post(fp_pred_url, data=payload, headers=headers)

    pred_super = np.asarray(json.loads(json_response.text)['predictions'])[0]
    n_super = list(np.where(pred_super>=0.3)[0])

    path_from_superclass = []
    for j in n_super:
        path_from_superclass += ontology_dictionary['Super_hierarchy'][str(j)]['Pathway']
    path_from_superclass = list(set(path_from_superclass))

    query_dict = {}
    query_dict["input_3"] = fp1
    query_dict["input_4"] = fp2

    # Handling CLASS
    fp_pred_url = "http://npclassifier-tf-server:8501/v1/models/CLASS:predict"
    payload = json.dumps({"instances": [ query_dict ]})

    headers = {"content-type": "application/json"}
    json_response = requests.post(fp_pred_url, data=payload, headers=headers)

    pred_class = np.asarray(json.loads(json_response.text)['predictions'])[0]
    n_class = list(np.where(pred_class>=0.1)[0])

    path_from_class = []
    for j in n_class:
        path_from_class += ontology_dictionary['Class_hierarchy'][str(j)]['Pathway']
    path_from_class = list(set(path_from_class))

    query_dict = {}
    query_dict["input_1"] = fp1
    query_dict["input_2"] = fp2

    # Handling PATHWAY
    fp_pred_url = "http://npclassifier-tf-server:8501/v1/models/PATHWAY:predict"
    payload = json.dumps({"instances": [ query_dict ]})

    headers = {"content-type": "application/json"}
    json_response = requests.post(fp_pred_url, data=payload, headers=headers)

    pred_path = np.asarray(json.loads(json_response.text)['predictions'])[0]
    n_path = list(np.where(pred_path>=0.5)[0])

    class_result = []
    superclass_result = []
    pathway_result = []

    # Voting on Answer
    pathway_result, superclass_result, class_result, isglycoside = prediction_voting.vote_classification(n_path, 
                                                                                                        n_class, 
                                                                                                        n_super, 
                                                                                                        pred_class,
                                                                                                        pred_super, 
                                                                                                        path_from_class, 
                                                                                                        path_from_superclass, 
                                                                                                        isglycoside, 
                                                                                                        ontology_dictionary)
    
    return isglycoside, class_result, superclass_result, pathway_result, path_from_class, path_from_superclass, n_path, fp1, fp2

# This gets you the model metadata
@server.route("/model/metadata")
def metadata():
    """Serve a file from the upload directory."""
    all_metadata = {}
    pathway_metadata = json.loads(requests.get("http://smart3-tf-server:8501/v1/models/CHANNEL1/metadata").text)
    class_metadata = json.loads(requests.get("http://smart3-tf-server:8501/v1/models/CHANNEL2/metadata").text)

    all_metadata["pathway"] = pathway_metadata
    all_metadata["class"] = class_metadata

    return json.dumps(all_metadata)

if __name__ == "__main__":
    app.run_server(debug=True, port=5000, host="0.0.0.0")
