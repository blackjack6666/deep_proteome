from tsv_reader import protein_info_from_fasta,map_psm_file
from protein_coverage import fasta_reader
import aho_corasick
import multiprocessing_naive_algorithym
from calculations_and_plot import identified_proteome_cov, creat_ID_pep_dict
from collections import Counter
import pandas as pd
import os
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import random
import plotly.graph_objects as go

# construct dash dataframe for plotting
# fasta_path = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640.fasta'
# protein_dict = fasta_reader(fasta_path)
#
# ID_list, seq_list = multiprocessing_naive_algorithym.extract_UNID_and_seq(protein_dict)
# seq_line = multiprocessing_naive_algorithym.creat_total_seq_line(seq_list)
# protein_info_dict = protein_info_from_fasta(fasta_path)
#
# pos_id_dict = multiprocessing_naive_algorithym.read_position_ID_into_dict(ID_list,seq_list,seq_line)
#
# info_list = []
# psm_path = 'D:/data/deep_proteome/20210114_chymo/psm.tsv'
# file_psm_dict = map_psm_file(psm_path)
# del file_psm_dict['CT_37C_0min']
# del file_psm_dict['CT_50C_0min']

# get overlapped proteins from all files
# identified_protein_list = [p for p in protein_dict]
# for f in file_psm_dict:
#     psm_list = file_psm_dict[f]
#     psm_dict = Counter(psm_list)
#     peptide_list = [p for p in psm_dict]
#     automaton = aho_corasick.automaton_trie(peptide_list)
#     aho_result = aho_corasick.automaton_matching(automaton, seq_line)
#     id_pep_dict = creat_ID_pep_dict(aho_result, pos_id_dict)
#     identified_protein_list = [p for p in id_pep_dict if p in identified_protein_list]

# print (len(identified_protein_list))

# for ind,val in enumerate(file_psm_dict):
#     print (val)
#     psm_list = file_psm_dict[val]
#     psm_dict = Counter(psm_list)
#     peptide_list = [p for p in psm_dict]
#     automaton = aho_corasick.automaton_trie(peptide_list)
#     aho_result = aho_corasick.automaton_matching(automaton,seq_line)
#     cov_dict = identified_proteome_cov(aho_result,protein_dict)[1]
#     id_pep_dict = creat_ID_pep_dict(aho_result,pos_id_dict)
#
#     spec_dict = {}
#     for id in id_pep_dict:
#         spec_count = 0
#         for pep in id_pep_dict[id]:
#             spec_count+=psm_dict[pep]
#         spec_dict[id] = spec_count
#     file_info_list = [[p,cov_dict[p],
#                        len(protein_dict[p]),
#                        spec_dict[p],
#                        protein_info_dict[p][0],
#                        protein_info_dict[p][-1],
#                        val,
#                        ind] for p in id_pep_dict]
#     info_list+=file_info_list
#
# df = pd.DataFrame(info_list, columns=['protein id','coverage','protein len','spectra count','gene name','class','file name','file number'])
# df.to_csv('D:/data/deep_proteome/20210114_chymo/dash_info.csv')


df = pd.read_csv('D:/uic/lab/data/Deep_proteome/dash_info_1_19.csv')
file_number_name_map = {f_number:f_name for f_name, f_number in zip(df['file name'],df['file number'])}
protein_class = df['class'].unique()
length_of_color_scheme = len(protein_class)
color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(length_of_color_scheme)]

color_map = {each:c for each,c in zip(protein_class,color)}

# dash implementation
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    dcc.Graph(id='graph-with-slider'),
    html.Div(id='my-output'),
    dcc.Slider(
        id='different time points',
        min=df['file number'].min(),
        max=df['file number'].max(),
        value=df['file number'].min(),
        marks=file_number_name_map,
        step=None
    )
])


@app.callback(
    Output('graph-with-slider', 'figure'),
    Output('my-output', 'children'),
    Input('different time points', 'value'))
def update_figure(selected_file_number):
    filtered_df = df[df['file number'] == selected_file_number]

    fig = px.scatter(filtered_df, x="uniprot_num", y="protein len",
                     color='class', hover_name="protein id",
                     log_x=False, color_discrete_map=color_map, size= 'coverage', size_max=20,
                     title='{}'.format(file_number_name_map[selected_file_number]),opacity=0.5)

    # fig = px.scatter_3d(filtered_df, x='uniprot_num', y='protein len', z='coverage',
    #                     color='class', size_max=5, color_discrete_map=color_map,log_x=True,
    #                     hover_name="protein id",opacity=0.7)

    fig.update_layout(transition_duration=500)

    return fig, "{} identified proteins, {} sp, {} tr, average coverage {:.2f}%, {} proteins coverage>80%".format(filtered_df.shape[0],
                                                                         filtered_df[filtered_df['class']=='sp'].shape[0],
                                                                         filtered_df[filtered_df['class']=='tr'].shape[0],
                                                                         filtered_df['coverage'].mean(),
                                                                         filtered_df[filtered_df['coverage']>80].shape[0])


if __name__ == '__main__':
    app.run_server(debug=True)
