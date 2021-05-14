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
import urllib.parse

import smart3wrapper


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
            html.H5("NMR Peaks Entry"),
            dbc.Textarea(className="mb-3", id='query_text', placeholder="NMR Spectrum", rows="20"),
            html.Br(),
            dbc.InputGroup(
                [
                    dbc.InputGroupAddon("Channels", addon_type="prepend"),
                    dbc.Select(
                        id="channel",
                        options=[
                            {"label": "Normal HSQC", "value": "1"},
                            {"label": "Edited HSQC", "value": "2"},
                        ],
                        value="1"
                    )
                ],
                className="mb-3",
            ),
            dbc.InputGroup(
                [
                    dbc.InputGroupAddon("Molecular Weight (Optional)", addon_type="prepend"),
                    dbc.Input(id="mw_filter", placeholder="MW Filter", type="number"),
                ],
                className="mb-3",
            ),
            html.Hr(),
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
    [Input('query_text', 'value'), Input('channel', 'value'), Input("mw_filter", "value")],
)
def handle_query(query_text, channel, mw_filter):
    # Saving input to a file to be read
    from io import StringIO
    data = StringIO(query_text)
    nmr_data_df = pd.read_csv(data, sep=None)

    # Determing if a molecular weight as entered
    try:
        mw_filter = float(mw_filter)
    except:
        mw_filter = None

    # Drawing image
    output_nmr_image = os.path.join("output", str(uuid.uuid4()) + ".png")
    top_search_results_df = smart3wrapper.search_smart3(nmr_data_df, output_image=output_nmr_image, channel=int(channel), mw_filter=mw_filter)

    # Reformatting the results
    results_list = top_search_results_df.to_dict(orient="records")
    for result_dict in results_list:
        all_names_list = [name[:20] for name in result_dict["Name"]]
        result_dict["Name"] = "\n".join(all_names_list)

    if len(results_list) == 0:
        return ["No Results", [html.Img(src="/plot/{}".format(os.path.basename(output_nmr_image)), width="1200px")]]

    top_search_results_df = pd.DataFrame(results_list)
    top_search_results_df['structure'] = top_search_results_df["SMILES"].apply(lambda x: '![SMILES](https://gnps-structure.ucsd.edu/structureimg?smiles={})'.format(urllib.parse.quote(x)))
    top_search_results_df = top_search_results_df[["Name", "MW", "Cosine score", "structure"]]

    table_fig = dash_table.DataTable(
        columns=[
            {"name": i, "id": i, "deletable": False, "selectable": True, 'presentation':'markdown'} for i in top_search_results_df.columns
        ],
        data=top_search_results_df.to_dict(orient="records"),
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
        style_cell={
            'whiteSpace': 'pre-line' #to allow breaks to new line, https://community.plotly.com/t/creating-new-line-within-datatable-cell/44145/2
        }
    )

    return [table_fig, [html.Img(src="/plot/{}".format(os.path.basename(output_nmr_image)), width="1200px")]]

@server.route("/plot/<uuid_save>")
def download(uuid_save):
    """Serve a file from the upload directory."""
    return send_from_directory("./output", uuid_save)

@server.route('/api/smart3/search', methods=['POST', 'GET'])
def apismart3search():
    nmr_data_df = pd.DataFrame(json.loads(request.values["peaks"]))
    channel = int(request.values.get("channel", 1))
    topK = smart3wrapper.search_smart3(nmr_data_df, channel=channel)

    return json.dumps(topK.to_dict(orient="records"))

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
