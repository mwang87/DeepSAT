# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_table
import plotly.express as px
from dash.dependencies import Input, Output, State

from urllib.parse import urlencode, quote, parse_qs

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

INPUT_DASHBOARD = [
    dbc.CardHeader(dbc.Row([
            dbc.Col(
                html.H5("Input Data"),
            ),
            dbc.Col(
                html.A(
                    dbc.Button("Link to Plot", 
                        color="primary", size="sm", 
                        className="mr-1", 
                        style={
                            "float" : "right"
                        }
                    ),
                    id="plot_link", 
                )
            )
    ])),
    dbc.CardBody([
        html.Div(id='version', children="Version - 0.1"),
        html.Br(),
        html.H5("SMART NMR 3 Peaks Entry"),
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
    ])
]

DRAWING_DASHBOARD = [
    dbc.CardHeader(dbc.Row([
            dbc.Col(
                html.H5("Drawing Results")
            ),
    ])),
    dbc.CardBody([
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
    ])
]

EXAMPLE_DASHBOARD = [
    dbc.CardHeader(dbc.Row([
            dbc.Col(
                html.H5("Examples")
            ),
    ])),
    dbc.CardBody([
        html.A("Swinolide A", href='/?peaks=1H%2C13C%0A5.79%2C113.2%0A7.58%2C153.3%0A1.88%2C12.3%0A6.08%2C142.2%0A2.46%2C37.4%0A2.18%2C37.4%0A4.14%2C66.6%0A1.58%2C40.8%0A1.73%2C40.8%0A4.51%2C65.7%0A5.69%2C129.8%0A5.78%2C123.2%0A1.82%2C29.9%0A2.27%2C29.9%0A3.86%2C65.8%0A1.46%2C33.8%0A2.14%2C33.8%0A4.01%2C75.1%0A3.35%2C57.4%0A1.68%2C41%0A0.81%2C9.4%0A3.83%2C73.8%0A1.62%2C38.4%0A3.98%2C71.3%0A1.75%2C41.3%0A0.97%2C9.2%0A5.36%2C74.3%0A1.95%2C37.6%0A0.84%2C9.1%0A3.12%2C76%0A1.65%2C33.2%0A0.99%2C17.7%0A1.27%2C23.9%0A1.38%2C23.9%0A1.3%2C29.3%0A1.9%2C29.3%0A4.02%2C71.4%0A1.6%2C34.8%0A1.82%2C34.8%0A3.53%2C73.2%0A3.33%2C55.2%0A1.18%2C38.8%0A1.96%2C38.8%0A3.69%2C64.5%0A1.2%2C21.7'),
        html.Br(),
        html.A("Homonojirimycin 7-glycoside with Edited HSQC", href='/?peaks=13C%2C1H%2CIntensity%0A102.16999816894501%2C4.8899998664856%2C8042.07861328125%0A76.4599990844727%2C3.43415474891663%2C10728.3544921875%0A75.69000244140629%2C3.4524886608123797%2C8217.021484375%0A75.5500030517578%2C3.2774889469146697%2C9897.5654296875%0A71.02000427246091%2C3.41148996353149%2C4376.31005859375%0A70.4700012207031%2C3.36314058303833%2C5323.1181640625%0A70.16000366210939%2C3.41975212097168%2C8310.3291015625%0A70.09999847412111%2C3.48369264602661%2C5734.6123046875%0A68.3000030517578%2C3.9418740272522004%2C-6017.26416015625%0A68.3000030517578%2C3.7912034988403303%2C-7367.166015625%0A62.6100006103516%2C3.9221179485321005%2C-5700.578125%0A62.6100006103516%2C3.6975493431091295%2C-5259.1328125%0A62.0099983215332%2C3.83098268508911%2C-9722.73828125%0A62.0099983215332%2C3.65098285675049%2C-8217.021484375%0A59.6099967956543%2C3.20346140861511%2C6186.80517578125%0A55.0900001525879%2C2.6643035411834703%2C8042.27197265625%0A&channel=2')
    ])
]


BODY = dbc.Container(
    [   
        dcc.Location(id='url', refresh=False),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Card(INPUT_DASHBOARD),
                        html.Br(),
                        dbc.Card(EXAMPLE_DASHBOARD)
                    ],
                    className="col-4",
                    id="left_panel_col",
                ),
                dbc.Col(
                    [
                        dbc.Card(DRAWING_DASHBOARD),
                    ],
                    className="col-8",
                    id="right_panel_col",
                ),
            ],
            style={"marginTop": 30},
        ),
    ],
    fluid=True,
    className="",
)

app.layout = html.Div(children=[NAVBAR, BODY])


# This enables parsing the URL to shove a task into the qemistree id
@app.callback([
                  Output('query_text', 'value'),
                  Output('channel', 'value')
              ],
              [
                  Input('url', 'pathname')
              ],
              [
                  State('url', 'search')
              ])
def url_params(pathname, search):
    try:
        params = parse_qs(search[1:])
    except:
        params = {}
    
    return [
        params.get("peaks", dash.no_update)[0],
        params.get("channel", dash.no_update)[0]
    ]

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
    top_search_results_df['Cosine score'] = top_search_results_df['Cosine score'].round(decimals=2)


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



@app.callback(
    [
        Output("plot_link", "href")
    ],
    [
        Input('query_text', 'value'), 
        Input('channel', 'value')
    ],
)
def update_link(query_text, channel):
    return ["/?peaks={}&channel={}".format(quote(query_text), channel)]

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
