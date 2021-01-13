import os
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px

import pandas as pd

base_path = 'D:/data/Naba_deep_matrisome/01102021/'
folders = [f for f in os.listdir(base_path) if '.' not in f]
pep_path_list = [base_path+each+'/peptide.tsv' for each in folders]
file_name_number_dict = {i:pep_tsv.split('/')[-2] for pep_tsv, i in zip(pep_path_list,range(len(pep_path_list)))}

df = pd.read_csv('D:/data/Naba_deep_matrisome/01102021/dash_info.csv')

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    dcc.Graph(id='graph-with-slider'),
    dcc.Slider(
        id='different time points',
        min=df['file_number'].min(),
        max=df['file_number'].max(),
        value=df['file_number'].min(),
        marks=file_name_number_dict,
        step=None
    )
])


@app.callback(
    Output('graph-with-slider', 'figure'),
    Input('different time points', 'value'))
def update_figure(selected_year):
    filtered_df = df[df.file_number == selected_year]

    fig = px.scatter(filtered_df, x="coverage", y="length",
                     size="spec_count",color='ecm_class', hover_name="gene",
                     log_x=False, size_max=55)

    fig.update_layout(transition_duration=500)

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)