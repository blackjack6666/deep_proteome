import os
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import random
import pandas as pd
from protein_coverage import fasta_reader
from calculations_and_plot import protein_mass
from tsv_reader import protein_info_from_fasta


# base_path = 'D:/data/Naba_deep_matrisome/01102021/'
# folders = [f for f in os.listdir(base_path) if '.' not in f]
# pep_path_list = [base_path+each+'/peptide.tsv' for each in folders]
# file_number_name_dict = {i:pep_tsv.split('/')[-2] for pep_tsv, i in zip(pep_path_list,range(len(pep_path_list)))}
#
# # read dataframe and generate color scheme
# df = pd.read_csv('D:/data/Naba_deep_matrisome/01102021/dash_info_new_1_20.csv')
# ecm_class = df.ecm_class.unique()
# length_of_color_scheme = len(ecm_class)
# color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
#              for i in range(length_of_color_scheme)]
#
# color_map = {each:c for each,c in zip(ecm_class,color)}
# # print (color_map)
#
# # dash implementation
# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
#
# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
#
# app.layout = html.Div([
#     dcc.Graph(id='graph-with-slider'),
#     html.Div(id='my-output'),
#     dcc.Slider(
#         id='different time points',
#         min=df['file_number'].min(),
#         max=df['file_number'].max(),
#         value=df['file_number'].min(),
#         marks=file_number_name_dict,
#         step=None
#     )
# ])
#
#
# @app.callback(
#     Output('graph-with-slider', 'figure'),
#     Output('my-output', 'children'),
#     Input('different time points', 'value'))
# def update_figure(selected_file_number):
#     filtered_df = df[df.file_number == selected_file_number]
#
#     fig = px.scatter(filtered_df, x="uniprot_num", y="length",
#                      size="coverage",color='ecm_class', hover_name="gene",
#                      log_x=False, size_max=55, color_discrete_map=color_map,
#                      title='{}'.format(file_number_name_dict[selected_file_number]))
#
#     fig.update_layout(transition_duration=500)
#
#     return fig, "{} ECM proteins identified".format(filtered_df.shape[0])

# read dataframe and generate color scheme
# protein_dict = fasta_reader('D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta')
# protein_info_dict = protein_info_from_fasta('D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta')

# df_raw = pd.read_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist_centroid_mean.xlsx', index_col=0)
# df_raw = df_raw.dropna()
# prot_list = df_raw.index
# df_raw['gene'] = [protein_info_dict[each][0] for each in prot_list]
# df_raw['description'] = [protein_info_dict[each][1] for each in prot_list]
# df_raw['length'] = [len(protein_dict[each]) for each in prot_list]
# df_raw['mass'] = [protein_mass(protein_dict[each]) for each in prot_list]
# df_raw.to_excel('D:/data/native_protein_digestion/10282021/h20_cov_dist_centroid_mean_nadrop.xlsx')

df_raw = pd.read_excel('D:/uic/lab/native_protein_digestion/h20_cov_dist_centroid_mean_nadrop.xlsx',index_col=0)
df = pd.DataFrame(dict(distance=df_raw['01h_h2o'].tolist()+df_raw['02h_h2o'].tolist()+
                                df_raw['04h_h2o'].tolist()+df_raw['20h_h2o'].tolist(),
                       time_point=[1]*df_raw.shape[0]+[2]*df_raw.shape[0]+[4]*df_raw.shape[0]+[20]*df_raw.shape[0],
                       gene=df_raw['gene'].tolist()*4,
                       description=df_raw['description'].tolist()*4,len=df_raw['length'].tolist()*4,
                       mass=df_raw['mass'].tolist()*4),index=df_raw.index.tolist()*4)



# print (color_map)

# dash implementation
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    dcc.Graph(id='graph-with-slider'),
    html.Div(id='my-output'),
    dcc.Slider(
        id='different time points',
        min=1,
        max=20,
        value=1,
        marks={
        1: {'label': '1-hour', 'style': {'color': '#FF3F33'}},
        2: {'label': '2-hour', 'style': {'color': '#36FF33'}},
        4: {'label': '4-hour', 'style': {'color': '#F033FF'}},
        20: {'label': '20-hour', 'style': {'color': '#33B8FF'}}
        },
        step=None
    )
])


@app.callback(
    Output('graph-with-slider', 'figure'),
    Output('my-output','children'),
    Input('different time points', 'value'))

def update_figure(time_point):
    color_map = {1: '#FF3F33', 2: '#36FF33', 4: '#F033FF', 20: '#33B8FF'}
    filtered_df = df[df.time_point == time_point]

    fig = px.scatter(filtered_df, x="mass", y="distance",
                      hover_name="gene",
                     log_x=False, size_max=30, trendline="ols",
                     title='Time point: %s h' % str(time_point))
    fig.data[1].line.color = 'black'
    fig.update_traces(marker=dict(
        color=color_map[time_point]))
    fig.update_yaxes(range=[0, 180])
    fig.update_layout(transition_duration=500)

    return fig, "Average distance: %.2f" % filtered_df.distance.mean()




if __name__ == '__main__':
    app.run_server(debug=True)

# import re
# df = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx',index_col=0)
# gene_id = df['gene_id'].tolist()
# gene_id_filter = []
#
# for each in gene_id:
#     if re.search(r'[a-z]',each) == None:
#         gene_id_filter.append(each)
#     else:
#         continue
#
# new_df = pd.DataFrame(columns=['human matrisome genes'])
# new_df['human matrisome genes'] = gene_id_filter
# new_df.to_excel('D:/data/Naba_deep_matrisome/matrisome_human_genes.xlsx')
