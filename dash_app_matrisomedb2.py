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
from bokeh.embed import components
from bokeh.io import show
from bokeh.models import CustomJS, Dropdown, ColumnDataSource, HoverTool,IndexFilter,CDSView, Select
# from bokeh.models.widgets import Select
from bokeh.plotting import figure
from bokeh.layouts import column,row


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
        html.H1(id='H1', children='PTM distributions in Matrisome', style={'textAlign': 'center',
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


def ptm_dist_bokeh(id_ptm_idx_dict, protein_dict, uniprot_gene_dict, protein_entry):

    """
    use bokeh, callback requires JS code, still in development
    :param id_ptm_idx_dict:
    :param protein_dict:
    :param uniprot_gene_dict:
    :param protein_entry:
    :return:
    """

    ptm_set = list(set([ptm for each in id_ptm_idx_dict for ptm in id_ptm_idx_dict[each]]))
    print(ptm_set)

    df = pd.DataFrame(columns=['PTM', 'freq', 'protein_id', 'len', 'gene'])
    tp_list = [(ptm, len(id_ptm_idx_dict[each][ptm]), each, len(protein_dict[each]), uniprot_gene_dict[each])
               for each in id_ptm_idx_dict for ptm in ptm_set]
    ptm_list, freq_list, protein_id_list, length_list, gene_list = zip(*tp_list)
    df['PTM'] = ptm_list
    df['freq'] = freq_list
    df['protein_id'] = protein_id_list
    df['len'] = length_list
    df['gene'] = gene_list
    df = df.sort_values(by=['len'])
    df_select = df.loc[df['PTM'] == ptm_set[0]].copy()
    df_select.loc[:,'ind'] = range(0,df_select.shape[0])

    source = ColumnDataSource(data={'PTM':df_select.PTM,'freq':df_select.freq,
                                    'protein_id':df_select.protein_id,'ind':df_select.ind,
                                    'len':df_select.len,'gene':df_select.gene})
    hover = HoverTool(names=['line','scatter'],tooltips=[('Uniprot_id', '@protein_id'), ('gene', '@gene'),('length','@len'),('frequency','@freq'),])

    p = figure(x_range=(0,df_select.shape[0]),
               # y_range=(0,2),
               tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                      'reset', hover],
               plot_height=500, plot_width=1600,
               toolbar_location='right',
               title='',
               x_axis_label='ECM proteins')

    p.line(x='ind',y='freq',source=source,line_width=2,color='#000000',name='line')

    df_scatter = df_select.loc[df_select['protein_id']==protein_entry]
    scatter_source = ColumnDataSource(data={'PTM':df_scatter.PTM,'freq':df_scatter.freq,
                                    'protein_id':df_scatter.protein_id,'ind':df_scatter.ind,
                                    'len':df_scatter.len,'gene':df_scatter.gene})
    p.scatter(x='ind',y='freq',color='red',source=scatter_source,name='scatter',size=10)

    # Define a callback function
    filter_ = IndexFilter(indices=[])
    callback = CustomJS(args=dict(src=source,filter=filter_),code='''
    const indices = []
    for (var i = 0; i < src.get_length(); i++) {
    console.log(i, src.data['PTM'][i], cb_obj.value)
    if (src.data['PTM'][i] == cb_obj.value) {
      indices.push(i)
        }
    }
    filter.indices = indices
    src.change.emit()
  ''')

    # create dropdown select widget
    select = Select(title="PTMs", options=ptm_set, value=ptm_set[0])
    # Attach the update_plot callback to the 'value' property of select
    select.js_on_change('value', callback)
    view = CDSView(source=source,filters=[filter_])
    layout = column(select, p)
    show(layout)

if __name__ == "__main__":

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

    # plot using bokeh app,updates 8/21
    # ptm_dist_bokeh(ptm_map_result[0],ecm_protein_dict,ecm_prot_gene_dict,'P11276')