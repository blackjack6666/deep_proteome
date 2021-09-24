import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from statannot import add_stat_annotation
from collections import defaultdict
import pylab
from brokenaxes import brokenaxes
from matplotlib.gridspec import GridSpec
import matplotlib
from protein_coverage import fasta_reader
from tsv_reader import peptide_counting

font = {'family' : 'Arial',
        'weight' : 'medium',
        'size'   : 8}

matplotlib.rc('font', **font)
### filter sp and tr entries into different
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
sp_tr_dict = defaultdict(set)
with open(fasta_path,'r') as f:

    f_split = f.read().split('\n>')
    for each in f_split[1:]:
        sp_tr_dict[each.split('|')[0]].add(each.split('|')[1])
    sp_tr_dict[f_split[0].split('|')[0].split('>')[1]].add(f_split[0].split('|')[1])

df_ecm_aggre = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
df_ecm_aggre = df_ecm_aggre.copy()

df_ecm_aggre = df_ecm_aggre[df_ecm_aggre.index.isin(sp_tr_dict['sp'])]

category_list = df_ecm_aggre['category']

ecm_class_color_dict = {"Collagens": '#0584B7', 'ECM-affiliated Proteins':'#F4511E',
                        'ECM Regulators':"#F9A287","Secreted Factors":"#FFE188",
                        "ECM Glycoproteins":"#133463", "Proteoglycans":"#59D8E6"}

sort_ecm_order = {"ECM Glycoproteins":0,"Collagens":1,"Proteoglycans":2,"ECM-affiliated Proteins":3,"ECM Regulators":4,
                  "Secreted Factors":5}
sort_category = ["ECM Glycoproteins","Collagens","Proteoglycans","ECM-affiliated Proteins","ECM Regulators",
                  "Secreted Factors"]
color_map = {prot:ecm_class_color_dict[ecm_class] for prot,ecm_class in zip(df_ecm_aggre.index,df_ecm_aggre['category'])}

df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_ecm_aggre.index]
normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_ecm_aggre.index]
normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_ecm_aggre.index]
normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_ecm_aggre.index]


df_plot = pd.DataFrame(dict(gene=df_ecm_aggre['gene'],
                            category=category_list,
                            mw=df_ecm_aggre['MW_kDa'],
                            gfp_18_agg=df_ecm_aggre['GFP_seq_1080_ave_aggre_cov'],
                            gfp_18_standard=normal18GFP_cov,
                            sned_18_agg=df_ecm_aggre['SNED1_seq_1080_ave_aggre_cov'],
                            sned_18_standard=normal18SNED_cov),index=df_ecm_aggre.index)
### scatter plot with sizes
"""
def rand_jitter(arr):
    stdev = .01 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev

df_sort = df_plot.sort_values(by=['category'], key= lambda x:x.map(sort_ecm_order))
# df_sort = df_sort[df_sort['category'].isin(['ECM Glycoproteins','Collagens','Proteoglycans'])]
df_sort['gfp_agg_standard_ratio']=(df_plot['gfp_18_agg']-df_plot['gfp_18_standard'])/df_plot['gfp_18_standard']*100
df_sort['scatter_x'] = rand_jitter(df_sort['mw'])

# df_sort.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/scatter_plot_table.xlsx')
df_sort=df_sort.replace([np.inf, -np.inf], np.nan).dropna(axis=0)


fig,ax = plt.subplots(1,1,figsize=(10,6))
g = sns.scatterplot(data=df_sort,ax=ax,x='mw',y='gfp_18_agg',hue='category',size='gfp_agg_standard_ratio',sizes=(10, 500), palette=ecm_class_color_dict,alpha=0.9)
# sns.regplot(data=df_sort,ax=ax,x='scatter_x',y='gfp_18_agg', scatter=False, line_kws={'color':'black','ls':'dashed'},ci=None,truncate=False)
# g = sns.scatterplot(data=df_sort,ax=ax,x='scatter_x',y='gfp_agg_standard_ratio',hue='category',size='gfp_18_agg',sizes=(10, 500), palette=ecm_class_color_dict)
# sns.regplot(data=df_sort,ax=ax,x='scatter_x',y='gfp_agg_standard_ratio', scatter=False, line_kws={'color':'black','ls':'dashed'},ci=None,truncate=False)
h,l = g.get_legend_handles_labels()
plt.xlabel('MW in kDa')
plt.ylabel('18h aggre. coverage%')
# plt.xticks([])
# plt.legend([])
# plt.legend(h,l,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=13)
ax.get_legend().remove()
# figlegend = pylab.figure(figsize=(3,6))
# pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left')
# figlegend.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/legend.png',dpi=300)
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/GFP_scatter_x_MW_alpha.png', dpi=300)
plt.show()
"""

### violin plot for oveall
"""
fig,ax = plt.subplots(1,1)
df_plot_sub = pd.DataFrame(dict(agg_or_standard=['GFP_18h_standard']*df_plot.shape[0]+['GFP_18h_agg']*df_plot.shape[0],
                                coverage=df_plot['gfp_18_standard'].tolist()+df_plot['gfp_18_agg'].tolist()))

sns.violinplot(data=df_plot_sub,x='agg_or_standard',y='coverage',ax =ax)
add_stat_annotation(ax, data=df_plot_sub,x='agg_or_standard',y='coverage',box_pairs=[("GFP_18h_standard", "GFP_18h_agg")],
                    test='t-test_paired', text_format='star',loc='inside', verbose=2)
ax.set_xlabel('')
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/GFP_18h_agg_standard_violin_all.png', dpi=300)
plt.show()
"""

### violin plot for each category
"""
fig,axs = plt.subplots(2,3,figsize=(10,5))

for each, ax in zip(sort_category,[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]]):
    color = ecm_class_color_dict[each]
    sub_df = df_plot[df_plot['category']==each]
    sub_df_plot = pd.DataFrame(dict(agg_or_standard=['SNED_18h']*sub_df.shape[0]+['SNED_18h_agg']*sub_df.shape[0],
                                    coverage=sub_df['sned_18_standard'].tolist()+sub_df['sned_18_agg'].tolist()))
    # axs[ax[0],ax[1]].set_ylim([0,100])
    sns.violinplot(data=sub_df_plot, x='agg_or_standard', y='coverage', ax=axs[ax[0],ax[1]], color=color)
    add_stat_annotation(ax=axs[ax[0],ax[1]], data=sub_df_plot, x='agg_or_standard', y='coverage',
                        box_pairs=[("SNED_18h", "SNED_18h_agg")],
                        test='t-test_paired', text_format='star',loc='inside', verbose=2)
    axs[ax[0],ax[1]].set_xlabel('')
    axs[ax[0], ax[1]].set_ylabel('')

plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/SNED_violin_category.png', dpi=300)
plt.show()
"""

### heatmap
# fig,ax = plt.subplots(1,1,figsize=(8,15))
# df_plot = pd.DataFrame(dict(gene=df_ecm_aggre['gene'],
#                             category=category_list,
#                             gfp_30_agg=df_ecm_aggre['GFP_seq_30_ave_aggre_cov'],
#                             gfp_120_agg=df_ecm_aggre['GFP_seq_120_ave_aggre_cov'],
#                             gfp_240_agg=df_ecm_aggre['GFP_seq_240_ave_aggre_cov'],
#                             gfp_1080_agg=df_ecm_aggre['GFP_seq_1080_ave_aggre_cov'],
#                             ),index=df_ecm_aggre.index)
"""
df_sub = df_plot.sort_values(by=['category'])

df_heatmap = pd.DataFrame([np.array(row[2:]) for row in df_sub.itertuples(index=False)], columns=list(df_sub)[2:],index=df_sub.index)
g = sns.heatmap(data=df_heatmap, ax=ax,cbar_kws={'label': 'Coverage in %','shrink': 0.5},cmap='BuPu', yticklabels=True)
# ax.axes.yaxis.set_visible(False)
ax.set_xticklabels(labels=['30min','2h','4h','18h'], rotation=30, fontsize=10, ha='right')
for ytick,prot in zip(plt.gca().get_yticklabels(),df_heatmap.index):

    ytick.set_color(color_map[prot])


plt.setp(ax.get_yticklabels(), fontsize=5)
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figures/SNED1_agg_heatmap.png', dpi=300)
plt.show()
"""
### plot coverage difference between time points in categories
"""
fig,axs = plt.subplots(6,1, figsize=(8,8))
x=[0.5,2,4,8]
for each_cat,ax in zip(sort_category,axs):
    sub_df = df_ecm_aggre[df_ecm_aggre['category']==each_cat]
    for each in sub_df.itertuples():
        y = [0]+[each[4+i+1]-each[4+i] for i in range(3)]
        # bax = brokenaxes(xlims=((0,8),(15,20)),hspace=1.5)
        # y_log10 = np.log10(y)
        ax.plot(x, y, color=color_map[each[0]], linestyle='-')
        # if y[-1]-y[0] == 0: # if line is flat, make it transparent
        #     ax.plot(x,y_log10, color='lightcoral', linestyle='-',alpha=0.1)
        # else:
        #     ax.plot(x,y_log10, color='lightcoral', linestyle='-')
    ax.set_xticks([0.5, 2, 4, 8])
    ax.set_xticklabels(['0.5h','2h','4h','18h'], fontsize=8)
    ax.set_yticks([0,10,20])
    ax.set_yticklabels(['0','10','20'])
    # ax.set_xlim(0.5,8)
    # ax.set_title(each_cat, fontsize=12)
    # ax.set_xlabel('time point', fontsize=12)
    # ax.set_ylabel('%coverage',fontsize=12)
    ax.set(yticklabels=[])
plt.xlabel('time point')
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/GFP_agg_lineplot_yaxis_cov_diff.png', dpi=300)
# plt.tight_layout()
plt.show()
"""
### line plot showing aggreated coverage for each category
"""
fig,axs = plt.subplots(6,1, figsize=(8,10))
x=[0.5,2,4,8]
for each_cat,ax in zip(sort_category,axs):
    sub_df = df_ecm_aggre[df_ecm_aggre['category']==each_cat]
    for each in sub_df.itertuples():
        y = each[8:12]
        # bax = brokenaxes(xlims=((0,8),(15,20)),hspace=1.5)
        # y_log10 = np.log10(y)
        # ax.plot(x, y, color=color_map[each[0]], linestyle='-')
        if y[-1]-y[0] == 0: # if line is flat, make it transparent
            ax.plot(x, y, color=color_map[each[0]], linestyle='-',alpha=0.2)
        else:
            ax.plot(x, y, color=color_map[each[0]], linestyle='-')
    ax.set_xticks([0.5, 2, 4, 8])
    ax.set_xticklabels(['0.5h','2h','4h','18h'], fontsize=8)
    # ax.set_yticks([0,10,20])
    # ax.set_yticklabels(['0','10','20'])
    # ax.set_xlim(0.5,8)
    # ax.set_title(each_cat, fontsize=12)
    # ax.set_xlabel('time point', fontsize=12)
    # ax.set_ylabel('%coverage',fontsize=12)
    # ax.set(yticklabels=[])
plt.xlabel('time point')
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/SNED1_agg_lineplot.png', dpi=300)
plt.tight_layout()
plt.show()
"""

### scatter plot with density
"""
df_agg_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_aggregated_D_F.xlsx',index_col=0)
all_05_aggre = np.mean([np.array(df_agg_summary['GFP_seq_D_30_aggre_coverage']), np.array(df_agg_summary['GFP_seq_F_30_aggre_coverage'])],axis=0)
all_2_aggre = np.mean([np.array(df_agg_summary['GFP_seq_D_120_aggre_coverage']), np.array(df_agg_summary['GFP_seq_F_120_aggre_coverage'])],axis=0)
all_4_aggre = np.mean([np.array(df_agg_summary['GFP_seq_D_240_aggre_coverage']), np.array(df_agg_summary['GFP_seq_F_240_aggre_coverage'])],axis=0)
all_18_aggre = np.mean([np.array(df_agg_summary['GFP_seq_D_1080_aggre_coverage']), np.array(df_agg_summary['GFP_seq_F_1080_aggre_coverage'])],axis=0)

all_2_normal = np.mean([np.array(df_summary['GFP_120D_coverage']), np.array(df_summary['GFP_120F_coverage'])], axis=0)
all_18_normal = np.mean([np.array(df_summary['GFP_1080D_coverage']), np.array(df_summary['GFP_1080F_coverage'])], axis=0)


df_plot = pd.DataFrame(dict(gene=df_ecm_aggre['gene'],

                            gfp_1080_parallel=df_ecm_aggre['GFP_seq_30_ave_aggre_cov'],
                            gfp_1080_agg=df_ecm_aggre['GFP_seq_1080_ave_aggre_cov'],
                            ),index=df_ecm_aggre.index)
# df_plot = df_plot[(df_plot['gfp_1080_parallel']!=0)&(df_plot['gfp_1080_agg']!=0)]

fig,ax = plt.subplots()
# x,y = df_plot['gfp_1080_parallel'], df_plot['gfp_1080_agg']
# xx, yy = np.meshgrid(x,y)  # convert 1d array to 2d grid
# z = np.sqrt(xx**2+yy**2)
# plt.contour(xx,yy,z,cmap='RdGy')
x = np.log2(df_plot['gfp_1080_parallel']+1)
y = np.log2(df_plot['gfp_1080_agg']+1)
ax.plot(x, y, 'k.', markersize=3)
# plt.show()

import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter

def myplot(x, y, s, bins=1000):
    # scatter density plot
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent

img, extent = myplot(x, y, 16)
ax.imshow(img, extent=extent, origin='lower', cmap=cm.binary)
ax.set_xlabel('log2 (0.5h aggre. cov%)', fontsize=8)
ax.set_ylabel('log2 (18h aggre. cov%)',fontsize=8)
ax.set_xticks(range(7))
ax.set_xticklabels([str(i) for i in range(7)])
ax.set_yticks(range(7))
ax.set_yticklabels([str(i) for i in range(7)])
plt.xlim(0,6.5)
plt.ylim(0,6.5)
ax.set_aspect('equal', adjustable='box')
ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/05h_aggre_18h_aggre_ECM.png', dpi=300)
plt.show()
"""

### jointplot
"""
g = sns.jointplot(data=df_plot, ax=ax, x='gfp_1080_parallel', y='gfp_1080_agg', hue='cat', palette=ecm_class_color_dict)

g.ax_joint.legend_._visible=False
g.ax_joint.set_xlabel('18h parallel cov%')
g.ax_joint.set_ylabel('18h aggre. cov%')
g.ax_marg_x.set_xlim(0, 100)
g.ax_marg_y.set_ylim(0,100)
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figures/ECM_GFP_18_normal_aggre.png', dpi=300)
plt.show()
"""
### matrisome db coverage
"""
matrisome_df = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx')

fig,axs = plt.subplots(2,3,figsize=(10,5))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
for each, ax in zip(sort_category,[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]]):
    color = ecm_class_color_dict[each]
    sub_df = matrisome_df[matrisome_df['category']==each]
    ave_cov = sub_df['cov'].mean()*100
    text = 'Ave Cov:\n%.2f%%' % ave_cov
    axs[ax[0],ax[1]].hist(sub_df['cov']*100,50,color=color)
    axs[ax[0],ax[1]].text(0.75, 0.95, text, transform=axs[ax[0],ax[1]].transAxes, fontsize=8,
        verticalalignment='top', bbox=props)
    axs[ax[0],ax[1]].set_xlabel(each+' coverage')
    axs[ax[0],ax[1]].set_ylabel('Frequency')

plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figure_update/matrisome_ave_cov.png', dpi=300)
plt.show()
"""
### generate 1d coverage bar
from calculations_and_plot import one_d_covearge_bar
from tsv_reader import modified_peptide_from_psm
protein_dict = fasta_reader(fasta_path)
time_points = ['30','120','240','1080']
peptide_list = [peptide_counting('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'D/peptide.tsv')+
                peptide_counting('D:/data/Naba_deep_matrisome/07232021_secondsearch/SNED1_seq_'+each+'F/peptide.tsv') for each in time_points]


# peptide_list = [df_summary.at['P28301','SNED1_seq_'+each+'D_total peptides identified'].split(', ')+\
#                df_summary.at['P28301','SNED1_seq_'+each+'F_total peptides identified'].split(', ') for each in time_points]

html_template = 'D:/data/Naba_deep_matrisome/html_template.html'
new_str = one_d_covearge_bar(html_template,peptide_list,protein_dict['Q8TER0'],
                             output_html_path='D:/data/Naba_deep_matrisome/07232021_secondsearch/coverage_1d/Q8TER0_SNED_beta.html',
                             screenshot='Q8TER0_SNED_beta.png')
