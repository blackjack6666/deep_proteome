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
            sphere.visible = False

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
    df_sasa = pd.read_excel('D:/data/native_protein_digestion/12072021/control/sasa.xlsx', index_col=0)
    df_sphere = pd.read_excel('D:/data/native_protein_digestion/12072021/control/aa_exposure_structuremap.xlsx',index_col=0)
    df_names = ['distance', 'density', 'sasa', 'sphere']
    df_list = [df_distance,df_density,df_sasa,df_sphere]

    for prot in df_distance.index:
        protein_entry = prot
        line_bokeh(df_list,df_names,protein_entry=protein_entry,save_path='D:/data/native_protein_digestion/12072021/control_html/'+protein_entry+'.html')
    # line_bokeh(df_list, df_names, protein_entry='P39019')

    # df_coverage = pd.read_excel('D:/data/native_protein_digestion/12072021/control/aggre_cov.xlsx',index_col=0)
    # protein_table(df_coverage, save_path='D:/data/native_protein_digestion/12072021/control_html/table.html')