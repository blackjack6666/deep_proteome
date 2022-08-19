import dash
from dash import html
from dash import dcc
from dash.dependencies import Input, Output
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from smart_crawler import ptm_map
from protein_coverage import fasta_reader
from tsv_reader import modified_peptide_from_psm, protein_info_from_fasta


def ptm_dist_dash(id_ptm_idx_dict, protein_dict, uniprot_gene_dict, protein_entry):
    """
    generate ptm distribution in html file with plotly and dash app
    :param id_ptm_idx_dict: first return from ptm_map in smart_crawler.py
    :param protein_dict: protein seq dictionary
    :param queried protein, uniprot ID
    :return:
    """

    ptm_set = list(set([ptm for each in id_ptm_idx_dict for ptm in id_ptm_idx_dict[each]]))
    print(ptm_set)

    df = pd.DataFrame(columns=['PTM', 'freq', 'protein_id','len','gene'])
    tp_list = [(ptm, len(id_ptm_idx_dict[each][ptm]), each, len(protein_dict[each]), uniprot_gene_dict[each])
               for each in id_ptm_idx_dict for ptm in ptm_set]
    ptm_list, freq_list, protein_id_list, length_list, gene_list = zip(*tp_list)
    df['PTM'] = ptm_list
    df['freq'] = freq_list
    df['protein_id'] = protein_id_list
    df['len'] = length_list
    df['gene'] = gene_list
    df = df.sort_values(by=['len'])

    # dash app initiation
    app = dash.Dash()
    app_options = [{'label':each.replace('\\',''),'value':each} for each in ptm_set]
    app.layout = html.Div(id='parent', children=[
        html.H1(id='H1', children='PTM distributions', style={'textAlign': 'center',
                                                                          'marginTop': 40, 'marginBottom': 40}),

        dcc.Dropdown(id='dropdown',
                     options=app_options,
                     value=ptm_set[0]),
        dcc.Graph(id='line_plot')
    ])

    @app.callback(Output(component_id='line_plot', component_property='figure'),
                  [Input(component_id='dropdown', component_property='value')])
    def graph_update(ptm):

        df_select = df.loc[df['PTM'] == ptm]

        # line plot showing ptm distribution across proteins
        x = range(df_select.shape[0])
        y = df_select['freq']
        df_plot_line = pd.DataFrame(dict(x=x, freq=y,
                                         prot=df_select['protein_id'], len=df_select['len'], gene=df_select['gene']))
        fig1 = px.line(data_frame=df_plot_line, x='x', y='freq', custom_data=['prot', 'gene', 'len'])
        fig1.update_traces(line_color='#808080', line_width=3,
                           hovertemplate='Uniprot_id: %{customdata[0]}<br>gene: %{customdata[1]}<br>length: %{customdata[2]}aa')

        # dot showing ptm freq for one protein
        df_plot_scatter = df_plot_line.loc[df_plot_line['prot'] == protein_entry]
        fig2 = px.scatter(df_plot_scatter, x='x', y='freq', custom_data=['prot','gene', 'len'])
        fig2.update_traces(marker=dict(
            color='LightSkyBlue',
            size=15,
            line=dict(
                color='black',
                width=2
            )
        ), hovertemplate='Uniprot_id: %{customdata[0]}<br>gene: %{customdata[1]}<br>length: %{customdata[2]}aa')
        fig3 = go.Figure(data=fig1.data + fig2.data)
        fig3.update_layout(xaxis_title="ECM proteins", yaxis_title="frequency")

        return fig3

    app.run_server(debug=True)


if __name__=="__main__":

    # use subset of protein dict
    protein_dict = fasta_reader('D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta')
    ecm_df = pd.read_excel('F:/matrisomedb2.0/AnnotatedDeDuplicated.xlsx',index_col=0)
    ecm_list = ecm_df[ecm_df['Division']!='Non-matrisome'].index.tolist()
    ecm_protein_dict = {prot:protein_dict[prot] for prot in ecm_list}
    ecm_prot_gene_dict = {prot:ecm_df.at[prot,'gene'] for prot in ecm_list}

    # ptm mapping
    psm_tsv = 'D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_1080D/psm.tsv'
    psm_list = modified_peptide_from_psm(psm_tsv)
    ptm_map_result = ptm_map(psm_list,ecm_protein_dict)

    # plot dash app
    ptm_dist_dash(ptm_map_result[0],ecm_protein_dict,ecm_prot_gene_dict,'P11276')