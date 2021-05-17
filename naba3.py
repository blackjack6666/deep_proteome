import pandas as pd
from tsv_reader import venn_diagram_gen
import numpy as np
import pickle as p
from scipy.stats import ttest_ind
from collections import defaultdict
import matplotlib.pyplot as plt
import pickle as p
from matplotlib.lines import Line2D
from scipy.stats import ttest_rel


ecm_protein_color_dict = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/ecm_protein_color_map_dict.p','rb'))
ecm_class_color_dict = {"Collagens": '#F23F51', 'ECM-affiliated Proteins':'#23AECA',
                        'ECM Regulators':"#23CA66","Secreted Factors":"#E3EB09",
                        "ECM Glycoproteins":"#EBA709", "Proteoglycans":"#EB09DC"}
color_ecm_class_dict = {ecm_class_color_dict[cat]:cat for cat in ecm_class_color_dict}

df = pd.read_csv('D:/data/Naba_deep_matrisome/02152021_1/dash_info.csv')

protein_18_2A_list = df[df['file_number']==7]['protein_id'].tolist()
protein_18_2A_cov = df[df['file_number']==7]['coverage'].tolist()
protein_18_2A_cov_dict = {prot:cov for prot,cov in zip(protein_18_2A_list,protein_18_2A_cov)}

protein_163_3A_list = df[df['file_number']==0]['protein_id'].tolist()
protein_163_3A_cov = df[df['file_number']==0]['coverage'].tolist()
protein_163_3A_cov_dict = {prot:cov for prot,cov in zip(protein_163_3A_list,protein_163_3A_cov)}

overlapped_protein_set = set(protein_18_2A_list).intersection(set(protein_163_3A_list))
# df_overlap_protein_cov = pd.DataFrame(columns=list(overlapped_protein_set),index=['163_X','163_Y','182_X','182_Y','color', 'ECM catgory'])
# df_overlap_protein_cov.iloc[0,:] = np.random.normal(1, 0.04, size=len(overlapped_protein_set))
# df_overlap_protein_cov.iloc[1,:] = [protein_163_3A_cov_dict[prot] for prot in overlapped_protein_set]
# df_overlap_protein_cov.iloc[2,:] = np.random.normal(2, 0.04, size=len(overlapped_protein_set))
# df_overlap_protein_cov.iloc[3,:] = [protein_18_2A_cov_dict[prot] for prot in overlapped_protein_set]
# df_overlap_protein_cov.iloc[4,:] = [ecm_protein_color_dict[prot] for prot in overlapped_protein_set]
# df_overlap_protein_cov.iloc[5,:] = [color_ecm_class_dict[color] for color in [ecm_protein_color_dict[prot] for prot in overlapped_protein_set]]

df_overlap_protein_cov = p.load(open('D:/data/Naba_deep_matrisome/02152021_1/df_overlap_protein_cov_standard_digestion.p','rb'))

df_transpose = df_overlap_protein_cov.transpose()
df_transpose.columns = ['163_X','163_Y','182_X','182_Y','color', 'ECM catgory']
print (df_transpose)
fig,ax = plt.subplots(1,1, figsize=(10,15))

x_163 = df_transpose.iloc[:,0]
y_163 = df_transpose.iloc[:,1]
x_182 = df_transpose.iloc[:,2]
y_182 = df_transpose.iloc[:,3]
color_list = df_transpose.iloc[:,4]


sub_df = df_transpose[df_transpose['ECM catgory']=='ECM Glycoproteins']
print (sub_df.shape)
sub_x_163 = sub_df.iloc[:,0]
sub_y_163 = sub_df.iloc[:,1]
sub_x_182 = sub_df.iloc[:,2]
sub_y_182 = sub_df.iloc[:,3]

t, p = ttest_rel(y_163, y_182)
print (y_163.mean(), y_182)
print (t,p)
paired_t, paired_p = ttest_rel(sub_y_163,sub_y_182)
print (paired_t,paired_p)
colors = [v for v in ecm_class_color_dict.values()]
labels = [k for k in ecm_class_color_dict.keys()]
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='--') for c in colors]

#colors, labels, lines= ['#EBA709'], ['ECM Glycoproteins'], [Line2D([0], [0], color='#EBA709', linewidth=3, linestyle='--')]
#
ax.legend(lines, labels, framealpha=0.5,loc='upper right',fontsize=15)

ax.boxplot([y_163,y_182], labels=['SNED1 OE', 'SNED1 KO'],showfliers=False)
# plot points onto boxplot
for x,cov_list in zip([x_163,x_182],[y_163,y_182]):
    ax.plot(x,cov_list,'ko',markersize=10,alpha=0.5)

#connect dots with different color
for x1,y1,x2,y2,color in zip(x_163,y_163,x_182,y_182,color_list):
    ax.plot([x1,x2],[y1,y2],color=color, linestyle='--', alpha=0.8)

# for x1,y1,x2,y2 in zip(sub_x_163,sub_y_163,sub_x_182,sub_y_182):
#     ax.plot([x1, x2], [y1, y2], color=ecm_class_color_dict['ECM Glycoproteins'], linestyle='--', alpha=1)
ax.set_title('ECM protein coverage in +/- SNED1 standard digestion', fontsize=22)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_ylabel('%coverage',fontsize=20)
plt.show()