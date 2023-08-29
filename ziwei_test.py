import pandas as pd
from ast import literal_eval
import requests
import re
import time
import pickle as pp

### map pubmedID with journal name
def pubmd_journal_name():
    path = r'C:\Users\gao lab computer\Downloads/Sub-task 5.csv'
    df = pd.read_csv(path, converters={'Pub Med': literal_eval})
    ### PUBMED esummary api
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id='
    data_dict = pp.load(open('pride_journal_pats_dict.p', 'rb'))

    for project, pub_med in zip(df['Project'], df['Pub Med']):
        if pub_med != [] and project not in data_dict:
            pubid = pub_med[0]  # only use first pub med id

            try:
                req = requests.get(base_url + str(pubid))
                data = req.text
                journal = data.split('<Item Name="FullJournalName" Type="String">')[1].split('</Item>')[0]
                data_dict[project] = journal
                print (project,journal)
                time.sleep(1)
            except:
                data_dict[project] = 'no Esummary'
                print (pubid, 'no esummay')
    pp.dump(data_dict, open('pride_journal_pats_dict.p','wb'))


def journal_name_stats():
    from collections import Counter
    # get journal name from all projets and stats
    data_dict = pp.load(open('pride_journal_pats_dict.p', 'rb'))
    journal_name_list = [v for v in data_dict.values()]
    journal_freq = Counter(journal_name_list)
    df = pd.DataFrame(columns=['journal_name','frequency'])
    df['journal_name'] = [k for k in journal_freq.keys()]
    df['frequency'] = [v for v in journal_freq.values()]

    df.to_csv('PRIDE_journal_freq.csv')


def filter_by_impact_factor():
    # project_journal_name = pp.load(open('pride_journal_pats_dict.p', 'rb'))
    f_1 = r'C:\Users\gao lab computer\Downloads/scimagojr 2021  Subject Area - Biochemistry, Genetics and Molecular Biology.csv'
    df1 = pd.read_csv(f_1,sep=';')
    journal_impact_rank = {journal.lower():rank for journal, rank in zip(df1['Title'],df1['Rank'])}
    mannual_journal = ['nature','analytical chemistry','molecular and cellular proteomics : mcp','proteomics',
                       "npj parkinson's disease",'chemical science','journal of proteomics','journal of proteome research',
                       'proceedings of the national academy of sciences of the united states of america','science advances',
                       'british journal of anaesthesia']
    for each in mannual_journal:
        journal_impact_rank[each] = 'manual selection'
    # size_filtered_df = pd.read_csv(r'C:\Users\gao lab computer\Downloads/filtered.csv',index_col=0)
    # journal_name_list = [project_journal_name[each].replace('&amp;','and').lower() if each in project_journal_name else 'None'
    #                      for each in size_filtered_df['PXD_id']]
    #
    # size_filtered_df['journal_name'] = journal_name_list
    # size_filtered_df['impact_factor_rank'] = [journal_impact_rank[name] if name in journal_impact_rank else 'None'
    #                                           for name in journal_name_list]
    # size_filtered_df.to_csv(r'C:\Users\gao lab computer\Downloads/filter_by_size_IF.csv')

    size_filtered_df = pd.read_csv(r'C:\Users\gao lab computer\Downloads/filter_by_size_IF.csv')
    size_filtered_df['impact_rank'] = [journal_impact_rank[each] if each in journal_impact_rank else 'None' for each in size_filtered_df['journal_name']]
    size_filtered_df.to_csv(r'C:\Users\gao lab computer\Downloads/filter_by_size_IF.csv')

# filter_by_impact_factor()
# journal_name_stats()

# df = pd.read_csv(r'D:\data\pats\PRIDE_filter/Sub-task 5.csv', index_col=0)
# print (df.shape)
# print (df[df['Project']=='PXD028979'])

import numpy as np

from bokeh.plotting import figure, show, curdoc

from bokeh.models import ColumnDataSource, ColorBar, HoverTool, CustomJS

from bokeh.models.widgets import CheckboxGroup

from bokeh.transform import linear_cmap

from bokeh.palettes import Cividis10

from bokeh.models.mappers import LinearColorMapper

from bokeh.layouts import row

gene_list = ['SAMD11', 'HES4', 'ISG15']

file_paths = [
    r'C:\Users\gao lab computer\Downloads/12be772d-49bd-4abc-b08e-76684c2ef1f2_KPMPAtlas_PREMIER_062921.csv',
    r'C:\Users\gao lab computer\Downloads/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.csv']

checkbox_group = CheckboxGroup(labels=['File 1', 'File 2'], active=[0])

def checkbox_change(attr, old, new):

    if new:

        selected_file_index = new[0]

        if selected_file_index >= len(file_paths):

            print("Selected file index is out of range.")

            return

        selected_file_path = file_paths[selected_file_index]

        # df = pd.read_feather(selected_file_path, columns=gene_list + ['umap1', 'umap2', 'index'])
        df = pd.read_csv(selected_file_path,index_col=0)

        df = df.replace({'index' : 'cell_type'})

        median_score = []

        for idx, r in df.iterrows():

            score = np.sum(r[gene_list])

            median_score.append(score)

        df['score'] = median_score

        source.data = df.to_dict(orient='list')


checkbox_group.on_change('active', checkbox_change)


# Create the initial plot

# df = pd.read_feather(file_paths[0], columns=gene_list + ['umap1', 'umap2', 'index'])
df = pd.read_csv(r'C:\Users\gao lab computer\Downloads/12be772d-49bd-4abc-b08e-76684c2ef1f2_KPMPAtlas_PREMIER_062921.csv', index_col=0)

median_score = []

for idx, r in df.iterrows():
    score = np.sum(r[gene_list])

    median_score.append(score)

df['score'] = median_score

source = ColumnDataSource(df)

mapper = linear_cmap(field_name='score', palette=Cividis10, low=df['score'].min(), high=df['score'].max())

p = figure(title='UMAP Visualization', x_axis_label='umap1', y_axis_label='umap2', sizing_mode='stretch_width',
           height=1500, toolbar_location='above')

hover = HoverTool(tooltips=[('Cell Name', '@index')], mode='mouse')

p.add_tools(hover)

p.scatter("umap1", "umap2", color=mapper, source=source)

color_mapper = LinearColorMapper(palette=Cividis10, low=df['score'].min(), high=df['score'].max())

color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12, location=(0, 0), title='Score')

p.add_layout(color_bar, 'right')

layout = row([checkbox_group, p])

show(layout)
