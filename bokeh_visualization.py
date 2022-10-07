"""
For plotting and creating database for native digestion data, visualizing cleavage to center distance,
proximal atoms, sasa, etc..
"""

import numpy as np
from bokeh.plotting import figure, show, output_file,save
from bokeh.models import ColumnDataSource, Slider, CustomJS, Range1d, RadioButtonGroup, HoverTool, DataTable, TableColumn, HTMLTemplateFormatter, LegendItem, Legend
from bokeh.layouts import column
import pandas as pd
import time

## an example to use customJS callbacks to control graph
# slider_slope = Slider(start = 0, end = 1, value = 0.5, step = 0.1)
# slider_intercept = Slider(start = 0, end = 20, value = 10, step = 1)
#
# s = slider_slope.value  # slope.
# i = slider_intercept.value  # intercept.
#
# x = np.linspace(-40, 20, 200)
# y = [(s * xx + i) for xx in x]
#
# p = figure(title = "simple line example", plot_height = 500, plot_width = 600, y_range = Range1d(start = -80, end = 40), background_fill_color = '#efefef')
# r = p.line(x, y, color = "red", line_width = 1.5, alpha = 0.8)  # 1st line. This line can be controlled by sliders.
# q = p.line(x, 2 * x + 1.2, color = "blue", line_width = 1.9, alpha = 0.2)  # 2nd line.
#
# slider_code = '''   var i = slider_intercept.value
#
#                     var s = slider_slope.value
#                     var x = r.data_source.data['x'];
#                     var y = [];
#
#                     for (var index = 0; index < x.length; index ++)
#                         y.push((s * x[index]) + i);
#
#                     r.data_source.data['y'] = y
#                     r.data_source.change.emit(); '''
#
# slider_callback = CustomJS(args = dict(slider_slope = slider_slope,
#                                 slider_intercept = slider_intercept,
#                                 r = r), code = slider_code)
#
# slider_slope.js_on_change('value',slider_callback)
# slider_intercept.js_on_change('value',slider_callback)
#
# layout = column(p, slider_slope, slider_intercept)
# show(layout, notebook_handle = True)


def line_bokeh(dataframes, dataframe_source, protein_entry, save_path=None):
    """

    :param dataframes: a list of pandas dataframes, such as ditance, density, sasa, etc.
    :param dataframe_source: a list of string to name each df
    :param protein_entry: target protein
    :return:
    """
    time_start = time.time()
    values,times,source = [], [], []

    for each_df, name in zip(dataframes,dataframe_source):

        values += each_df.loc[protein_entry,:].tolist()
        times += list(range(1,each_df.shape[1]+1))
        source += [name]*each_df.shape[1]

    df_one_protein = pd.DataFrame(dict(value=values,time=times,source=source))

    labels = df_one_protein['source'].unique().tolist()

    # plot lines
    # hover = HoverTool(tooltips=[('domain', '@domain'), ('start position', '@start'),('end position','@end'),('domain coverage','@coverage{:.1%}'),])

    p = figure(plot_width=1000, plot_height=400,x_range=(0,8),tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                      'reset'],x_axis_label='time point')

    for each in labels:
        x = df_one_protein[df_one_protein['source']==each]['time'].tolist()
        y = df_one_protein[df_one_protein['source']==each]['value'].tolist()
        if each == 'distance':
            distance = p.line(x=x,y=y,line_width=3, line_color='#4875bd')
        elif each == 'density':
            density = p.line(x=x, y=y, line_width=3, line_color='#4875bd')
            density.visible = False # hide it at first
        elif each == 'sasa':
            sasa = p.line(x=x, y=y, line_width=3, line_color='#4875bd')
            sasa.visible = False # hide it at first
        else:
            sphere = p.line(x=x, y=y, line_width=3, line_color='#4875bd')
            sphere.visible = False # hide at first

    p.yaxis.visible = False
    p.xaxis.ticker = list(range(1,8))
    p.xaxis.major_label_overrides = {1: '10min', 2: '30min', 3: '60min',
                                     4: '120min', 5: '240min', 6: '1320min',7:'1740min'}
    color_map_dict = {'Alphafold':'#4875bd','PDB':'#2bb55e'}
    # dummy glyphs to help draw legend
    legend_gly = [p.line(x=[1, 1], y=[1, 1], line_width=3, color=c, name='dummy_for_legend')
                  for c in [v for v in color_map_dict.values()]]

    legend = Legend(title='', background_fill_color='white',
                    border_line_color='black', border_line_width=1,
                    border_line_alpha=0.7,
                    items=[LegendItem(label=lab, renderers=[gly])
                           for lab, gly in zip([d for d in color_map_dict.keys()], legend_gly)])
    p.add_layout(legend)
    # labels for showing
    button_labels = ['Cleavage to center distance', 'Proximal atoms', 'Solvent accessibility surface area', 'Structuremap sphere index']
    # create a radiobutton, active is index of labels
    bg = RadioButtonGroup(labels=button_labels, active=0)
    js_code = """
    var line_list = [line0,line1,line2,line3];

    var lab_len=cb_obj.labels.length;

    for (var i=0;i<lab_len;i++) {
    if (cb_obj.active == i) {
    line_list[i].visible = true;
    } else {
    line_list[i].visible = false;
    }
    }
    """

    bg_callback = CustomJS(args=dict(line0=distance,line1=density,line2=sasa,line3=sphere),code=js_code)
    bg.js_on_change('active', bg_callback)

    layout = column(bg,p)
    # show(layout)
    if save_path:
        output_file(save_path,mode='inline')
    save(layout)
    print (f'time used:{time.time()-time_start}')


def line_bokeh_alphafold_pdb(dataframes_af, dataframe_source, dataframes_pdb, uniport_pdb_dict, protein_entry, save_path=None):
    """
    plot both line trend from alphafold pdb or pdb from pdb.org after seqmappdb mapping
    :param dataframes_af: a list of pandas dataframes, such as ditance, density, sasa, etc., from alphafold pdb
    :param dataframe_source: df names
    :param dataframes_pdb: a list of pandas dataframes, such as ditance, density, sasa, etc., from mapped pdb files
    :param uniport_pdb_dict: mapping of uniport id to pdb id and chain name, {'uniprotid':'pdb_chain'}
    :param protein_entry:
    :param save_path:
    :return:
    """
    time_start = time.time()
    values_af, values_pdb, times, source = [], [], [], []
    df_shape1 = dataframes_af[0].shape[1]

    # getting one protein data from alphafold dfs
    for each_df, name in zip(dataframes_af, dataframe_source):
        values_af += each_df.loc[protein_entry, :].tolist()
        times += list(range(1, each_df.shape[1] + 1))
        source += [name] * each_df.shape[1]

    uniport_pdb = protein_entry+ '_' + uniport_pdb_dict[protein_entry] if protein_entry in uniprot_pdb_dict else None # uniport_pdb_chain
    map_pdb = True if uniport_pdb in dataframes_pdb[0].index.tolist() else False  # decide if any uniport-mapped pdb
    # print (map_pdb)
    # one protein df for alphafold data
    df_one_protein = pd.DataFrame(dict(value=values_af, time=times, source=source))
    # one protein df for pdb data
    if map_pdb:
        for each_df_pdb, name in zip(dataframes_pdb,dataframe_source):
            values_pdb += each_df_pdb.loc[uniport_pdb, :].tolist()
        df_one_protein_pdb = pd.DataFrame(dict(value=values_pdb,time=times,source=source))
        # print (df_one_protein_pdb)

    # legend label and colors
    af_label = 'Alphafold: '+protein_entry
    pdb_label = 'PDB: '+uniprot_pdb_dict[protein_entry] if map_pdb else 'PDB: not available'
    color_map_dict = {af_label: '#4875bd', pdb_label: '#2bb55e'}

    labels = df_one_protein['source'].unique().tolist()

    # plot lines
    # hover = HoverTool(tooltips=[('domain', '@domain'), ('start position', '@start'),('end position','@end'),('domain coverage','@coverage{:.1%}'),])

    p = figure(plot_width=1000, plot_height=400, x_range=(0, 8), tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                                                                        'reset'], x_axis_label='time point')
    # plot from alphafold data
    for each in labels:
        x = df_one_protein[df_one_protein['source'] == each]['time'].tolist()
        y = df_one_protein[df_one_protein['source'] == each]['value'].tolist()
        if each == 'distance':
            distance = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[af_label])
        elif each == 'density':
            density = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[af_label])
            density.visible = False  # hide it at first
        elif each == 'sasa':
            sasa = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[af_label])
            sasa.visible = False  # hide it at first
        else:
            sphere = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[af_label])
            sphere.visible = False  # hide at first

    # plot from mapped pdb data
    for each in labels:
        if map_pdb:
            x = df_one_protein_pdb[df_one_protein_pdb['source']==each]['time'].tolist()
            y = df_one_protein_pdb[df_one_protein_pdb['source']==each]['value'].tolist()
        else:
            x, y = [0]*df_shape1, [0]*df_shape1
        if each == 'distance':
            distance_pdb = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[pdb_label])
        elif each == 'density':
            density_pdb = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[pdb_label])
            density_pdb.visible = False  # hide it at first
        elif each == 'sasa':
            sasa_pdb = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[pdb_label])
            sasa_pdb.visible = False  # hide it at first
        else:
            sphere_pdb = p.line(x=x, y=y, line_width=3, line_color=color_map_dict[pdb_label])
            sphere_pdb.visible = False  # hide at first

    p.yaxis.visible = True
    p.xaxis.ticker = list(range(1, 8))
    p.xaxis.major_label_overrides = {1: '10min', 2: '30min', 3: '60min',
                                     4: '120min', 5: '240min', 6: '1320min', 7: '1740min'}

    # dummy glyphs to help draw legend
    legend_gly = [p.line(x=[1, 1], y=[1, 1], line_width=3, color=c, name='dummy_for_legend')
                  for c in [v for v in color_map_dict.values()]]

    legend = Legend(title='', background_fill_color='white',
                    border_line_color='black', border_line_width=1,
                    border_line_alpha=0.7,
                    items=[LegendItem(label=lab, renderers=[gly])
                           for lab, gly in zip([d for d in color_map_dict.keys()], legend_gly)])
    p.add_layout(legend)
    # labels for showing
    button_labels = ['Cleavage to center distance', 'Proximal atoms', 'Solvent accessibility surface area',
                     'Structuremap sphere index']
    # create a radiobutton, active is index of labels
    bg = RadioButtonGroup(labels=button_labels, active=0)
    js_code = """
        var line_list = [line0,line1,line2,line3];
        var line_list_pdb = [line4,line5,line6,line7];
        var lab_len=cb_obj.labels.length;

        for (var i=0;i<lab_len;i++) {
        if (cb_obj.active == i) {
        line_list[i].visible = true;
        line_list_pdb[i].visible = true;
        } else {
        line_list[i].visible = false;
        line_list_pdb[i].visible = false;
        }
        }
        """

    bg_callback = CustomJS(args=dict(line0=distance, line1=density, line2=sasa, line3=sphere,
                                     line4=distance_pdb, line5=density_pdb,line6=sasa_pdb,line7=sphere_pdb), code=js_code)
    bg.js_on_change('active', bg_callback)

    layout = column(bg, p)
    # show(layout)
    if save_path:
        output_file(save_path, mode='inline')
    save(layout)
    print(f'time used:{time.time() - time_start}')


def protein_table(df, save_path=None):
    """
    output a table with bokeh
    :param df:
    :return:
    """
    df['uniprot'] = df.index
    df['url'] = [each+'.html' for each in df.index]
    # print (df.columns)
    source = ColumnDataSource(df)
    columns = [TableColumn(field='uniprot',title='Uniprot ID',
                               formatter=HTMLTemplateFormatter(template='<a href="<%= url %>"><%= value %></a>'))]
    columns += [TableColumn(field=each,title=each)
                for each in df.columns[:-2]]
    table = DataTable(source=source,columns=columns, width=1000, height=600, editable=True)
    if save_path:
        output_file(save_path,mode='inline')
    save(table)
    # show(table)


if __name__ == "__main__":

    df_distance = pd.read_excel('D:/data/native_protein_digestion/12072021/control/KRtocenter_dist_unique.xlsx',index_col=0)
    df_density = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cov_KR_density_15A.xlsx', index_col=0)
    # df_sasa = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa.xlsx', index_col=0)
    df_sasa = pd.read_excel('F:/native_digestion/sasa_area_15Atotal.xlsx',index_col=0)
    df_sphere = pd.read_excel('D:/data/native_protein_digestion/12072021/control/aa_exposure_structuremap.xlsx',index_col=0)
    df_names = ['distance', 'density', 'sasa', 'sphere']
    df_list = [df_distance,df_density,df_sasa,df_sphere]

    df_distance_pdb = pd.read_excel('D:/data/native_protein_digestion/12072021/control/mappdb_distance.xlsx',index_col=0)
    df_density_pdb = pd.read_excel('D:/data/native_protein_digestion/12072021/control/mappdb_KR_density_15A.xlsx',index_col=0)
    df_sasa_pdb = pd.read_excel('D:/data/native_protein_digestion/12072021/control/mappdb_KR_sasa_15A.xlsx',index_col=0)
    df_sphere_pdb = df_sasa_pdb
    df_list_pdb = [df_distance_pdb,df_density_pdb,df_sasa_pdb,df_sphere_pdb]

    df_prot_pdb = pd.read_csv('C:/tools/seqmappdb/human/fully_covered_unique_PDB.csv')
    uniprot_pdb_dict = {prot.split('>')[1]: '_'.join(pdb.split('>')[1].split('_')[:2])
                        for prot, pdb in zip(df_prot_pdb['queryID'], df_prot_pdb['pdbchainID'])}
    # for prot in df_distance.index:
    #     protein_entry = prot
    #     line_bokeh(df_list,df_names,protein_entry=protein_entry,save_path='D:/data/native_protein_digestion/12072021/control_html/'+protein_entry+'.html')
    # line_bokeh(df_list, df_names, protein_entry='P39019')

    line_bokeh_alphafold_pdb(df_list,df_names,df_list_pdb,uniprot_pdb_dict,protein_entry='Q9Y281',save_path='D:/data/native_protein_digestion/12072021/control_html/test_Q9Y281.html')

    # df_coverage = pd.read_excel('D:/data/native_protein_digestion/12072021/control/aggre_cov.xlsx',index_col=0)
    # df_coverage = df_coverage.round(decimals=2)
    #
    # protein_table(df_coverage, save_path='D:/data/native_protein_digestion/12072021/control_html/table.html')