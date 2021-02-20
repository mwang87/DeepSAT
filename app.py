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
import uuid
from flask import send_from_directory

import sys
sys.path.insert(0, "Classifier")
import SMART_3
import get_highlight


server = Flask(__name__)
app = dash.Dash(__name__, server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

# Loading database into memory
DB, index_super = SMART_3.load_db(db_folder="./Classifier")

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
    pred_MW = np.asarray(prediction_dict["dense_1"])
    pred_class = np.asarray(prediction_dict["class"])
    pred_gly = np.asarray(prediction_dict["glycoside"])

    fingerprint_prediction, pred_MW, pred_class_index, pred_class_prob, pred_gly = SMART_3.refine_predict_nmr([fingerprint_prediction], [pred_MW], [pred_class], [pred_gly])

    experimental_mw = None
    top_candidates = 10

    topK = None
    topK = SMART_3.search_database(fingerprint_prediction, pred_MW, DB, mw=experimental_mw, top_candidates=top_candidates)

    table_fig = dash_table.DataTable(
        columns=[
            {"name": i, "id": i, "deletable": False, "selectable": True} for i in topK.columns
        ],
        data=topK.to_dict(orient="records"),
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


    # Drawing image
    output_nmr_image = os.path.join("output", str(uuid.uuid4()) + ".png")
    SMART_3._draw_search_image(output_nmr_image, nmr_mat, pred_class_index, pred_class_prob, pred_gly, index_super=index_super)

    return [table_fig, [html.Img(src="/plot/{}".format(os.path.basename(output_nmr_image)), width="1200px")]]

@server.route("/plot/<uuid_save>")
def download(uuid_save):
    """Serve a file from the upload directory."""
    return send_from_directory("./output", uuid_save)

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
