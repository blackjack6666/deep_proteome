import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from statannot import add_stat_annotation

df_ecm_aggre = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
df_ecm_aggre = df_ecm_aggre.copy()
category_list = df_ecm_aggre['category']

ecm_class_color_dict = {"Collagens": '#3DA8FB', 'ECM-affiliated Proteins':'#F3A343',
                        'ECM Regulators':"#FADDBE","Secreted Factors":"#F3BEDD",
                        "ECM Glycoproteins":"#C47CFD", "Proteoglycans":"#77F9F1"}
color_map = {prot:ecm_class_color_dict[ecm_class] for prot,ecm_class in zip(df_ecm_aggre.index,df_ecm_aggre['category'])}
df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_ecm_aggre.index]
normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_ecm_aggre.index]
normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_ecm_aggre.index]
normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_ecm_aggre.index]

##
df_plot = pd.DataFrame(dict(gene=df_ecm_aggre['gene'],
                            category=category_list,
                            gfp_18_agg=df_ecm_aggre['GFP_seq_1080_ave_aggre_cov'],
                            gfp_18_standard=normal18GFP_cov,
                            sned_18_agg=df_ecm_aggre['SNED1_seq_1080_ave_aggre_cov'],
                            sned_18_standard=normal18SNED_cov),index=df_ecm_aggre.index)

df_sort = df_plot.sort_values(by=['category'])
df_sort['gfp_agg_standard_ratio']=(df_plot['gfp_18_agg']-df_plot['gfp_18_standard'])/df_plot['gfp_18_standard']*100
df_sort['scatter_x'] = [i for i in range(df_sort.shape[0])]
df_sort=df_sort.replace([np.inf, -np.inf], np.nan).dropna(axis=0)


fig,ax = plt.subplots(1,1)
g = sns.scatterplot(data=df_sort,ax=ax,x='scatter_x',y='gfp_18_agg',hue='category',size='gfp_agg_standard_ratio',sizes=(10, 500), palette=ecm_class_color_dict)
h,l = g.get_legend_handles_labels()
plt.xlabel('ECM proteins')
plt.ylabel('18h aggregated coverage%')
plt.xticks([])
plt.legend([])
# plt.legend(h[:7],l[:7],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=13)
ax.get_legend().remove()
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figures/GFP_18h_agg_std_scatter_size_changeratio.png', dpi=300)
plt.show()


### violin plot for oveall
"""
fig,ax = plt.subplots(1,1)
df_plot_sub = pd.DataFrame(dict(agg_or_standard=['sned_18h_standard']*df_plot.shape[0]+['sned_18h_agg']*df_plot.shape[0],
                                coverage=df_plot['sned_18_standard'].tolist()+df_plot['sned_18_agg'].tolist()))

sns.violinplot(data=df_plot_sub,x='agg_or_standard',y='coverage',ax =ax)
add_stat_annotation(ax, data=df_plot_sub,x='agg_or_standard',y='coverage',box_pairs=[("sned_18h_standard", "sned_18h_agg")],
                    test='t-test_paired', text_format='star',loc='inside', verbose=2)
ax.set_xlabel('')
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figures/SNED_18h_agg_standard_violin_all.png', dpi=300)
plt.show()
"""
### violin plot for each category
"""
fig,axs = plt.subplots(2,3,figsize=(10,5))

for each, ax in zip(df_plot.category.unique(),[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]]):
    color = ecm_class_color_dict[each]
    sub_df = df_plot[df_plot['category']==each]
    sub_df_plot = pd.DataFrame(dict(agg_or_standard=['sned_18h_standard']*sub_df.shape[0]+['sned_18h_agg']*sub_df.shape[0],
                                    coverage=sub_df['sned_18_standard'].tolist()+sub_df['sned_18_agg'].tolist()))
    # axs[ax[0],ax[1]].set_ylim([0,100])
    sns.violinplot(data=sub_df_plot, x='agg_or_standard', y='coverage', ax=axs[ax[0],ax[1]], color=color)
    add_stat_annotation(ax=axs[ax[0],ax[1]], data=sub_df_plot, x='agg_or_standard', y='coverage',
                        box_pairs=[("sned_18h_standard", "sned_18h_agg")],
                        test='t-test_paired', text_format='star',loc='inside', verbose=2)
    axs[ax[0],ax[1]].set_xlabel('')
    axs[ax[0], ax[1]].set_ylabel('')

plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figures/SNED_18h_agg_standard_violin_category.png', dpi=300)
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
fig,axs = plt.subplots(6,1, figsize=(10,10))
x=[0.5,2,4,18]
for each_cat,ax in zip(ecm_class_color_dict,axs):
    sub_df = df_ecm_aggre[df_ecm_aggre['category']==each_cat]
    for each in sub_df.itertuples():
        y = [0]+[each[8+i+1]-each[8+i] for i in range(3)]
        # y_log10 = np.log10(y)
        ax.plot(x, y, color=color_map[each[0]], linestyle='-')
        # if y[-1]-y[0] == 0: # if line is flat, make it transparent
        #     ax.plot(x,y_log10, color='lightcoral', linestyle='-',alpha=0.1)
        # else:
        #     ax.plot(x,y_log10, color='lightcoral', linestyle='-')
    ax.set_xticks([0.5, 2, 4, 18])
    ax.set_xticklabels(['0.5h','2h','4h','18h'], fontsize=10)
    ax.set_yticks([0,20])
    ax.set_yticklabels(['0','20'])
    # ax.set_title(each_cat, fontsize=12)
    # ax.set_xlabel('time point', fontsize=12)
    # ax.set_ylabel('%coverage',fontsize=12)
    ax.set(yticklabels=[])
plt.xlabel('time point')
plt.savefig('D:/data/Naba_deep_matrisome/07232021_secondsearch/figures/SNED_agg_lineplot_yaxis_cov_diff.png', dpi=300)
plt.show()
"""