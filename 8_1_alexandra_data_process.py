import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from protein_coverage import fasta_reader
from scipy.stats import ttest_1samp

# xls = pd.ExcelFile('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_aggregated_D_F_Received_210727_AN.xlsx')
# df_matrisome_AN = pd.read_excel(xls,'Matrisome_Coverage',index_col=0)

def cosine_sim_calculating(v1, v2):
    """
    calculate the cosine similarity beweeen two b/y ions binned vectors
    :param v1:
    :param v2:
    :return:
    """
    from scipy import spatial
    return 1-spatial.distance.cosine(v1,v2)
"""
fasta_path = 'D:/data/Naba_deep_matrisome/uniprot-proteome_UP000000589_mouse_human_SNED1.fasta'
protein_dict = fasta_reader(fasta_path)

matrisomedb_excel = 'D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx'
df_matrisome = pd.read_excel(matrisomedb_excel)
ecm_category_dict = {i:j for i,j in zip(df_matrisome['protein_id'].tolist(),df_matrisome['category'].tolist())}
print (len(ecm_category_dict))

df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx',index_col=0)
ecm_category_dict_2  = {i:j for i,j in zip(df_aggre_coverage.index.tolist(),df_aggre_coverage['category'].tolist())}
print (len(ecm_category_dict_2))
ecm_category_dict.update(ecm_category_dict_2)

df_summary = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx',index_col=0)
"""
# df_aggre_coverage = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx', index_col=0)
# normal18GFP_cov = [np.mean([df_summary.at[prot,'GFP_1080D_coverage'], df_summary.at[prot,'GFP_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2GFP_cov = [np.mean([df_summary.at[prot,'GFP_120D_coverage'], df_summary.at[prot,'GFP_120F_coverage']]) for prot in df_aggre_coverage.index]
# normal18SNED_cov = [np.mean([df_summary.at[prot,'SNED1_1080D_coverage'], df_summary.at[prot,'SNED1_1080F_coverage']]) for prot in df_aggre_coverage.index]
# normal2SNED_cov = [np.mean([df_summary.at[prot,'SNED1_120D_coverage'], df_summary.at[prot,'SNED1_120F_coverage']]) for prot in df_aggre_coverage.index]



# for each_tp in df_aggre_coverage.itertuples():
#     print (each_tp[4:8],each_tp[-4:], cosine_sim_calculating(each_tp[4:8],each_tp[-4:]), np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[-4:])))
#     df_aggre_coverage.at[each_tp[0],'cosine_sim'] = cosine_sim_calculating(each_tp[4:8],each_tp[-4:])
#     df_aggre_coverage.at[each_tp[0],'euclidean_distance with +-'] = np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[8:12])) \
#         if each_tp[7]-each_tp[11]>=0 else -np.linalg.norm(np.array(each_tp[4:8])-np.array(each_tp[8:12]))
# df_aggre_coverage.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/8_1_matrisome_average_aggre.xlsx')


# for each in df_aggre_coverage['category'].unique():
#     sub_df = df_aggre_coverage[df_aggre_coverage['category']==each]
#     print (each, ttest_1samp(sub_df['euclidean_distance with +-'],0.0,alternative='greater'))

# summary_df = pd.read_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/7_24_summary_D_F_squential_standard.xlsx', index_col=0)
# summary_df_slice = summary_df[summary_df.index.isin(ecm_protein_list)]

# print (summary_df_slice['GFP_1080D_total peptides identified'].unique())
# peptide_set = set()
# for column in ['SNED1_1080D_total peptides identified', 'SNED1_1080F_total peptides identified']:
#     for each in summary_df_slice[column].unique():
#         print (each)
#         if each != 0:
#             if ',' in each:
#                 for pep in each.split(', '):
#                     peptide_set.add(pep)
#             else:
#                 peptide_set.add(each)
#         else:
#             continue
# print (len(peptide_set))
"""
bcdf_combined_agg = pd.read_excel('D:/data/Naba_deep_matrisome/BCDF_combined/BCDF_ecm_aggregated.xlsx',index_col=0)
samples, times = ['GFP_seq_', 'SNED1_seq_'],['30','120','240','1080']

df_aggre_coverage = pd.DataFrame()
for prot in bcdf_combined_agg.index:
    gene = bcdf_combined_agg.at[prot,'gene']
    category = bcdf_combined_agg.at[prot,'ECM_category']
    df_aggre_coverage.at[prot,'gene'] = gene
    df_aggre_coverage.at[prot, 'ECM category'] = category
    for sample in samples:
        for time in times:
            # if time != '30':
            #     aver_cov = np.mean([bcdf_combined_agg.at[prot,sample+'D_'+time+'_aggre_coverage'],bcdf_combined_agg.at[prot,sample+'F_'+time+'_aggre_coverage']])
            #     df_aggre_coverage.at[prot,sample+time+'_ave_aggre_cov']=aver_cov
            # else:
            #     # add 30min data BC into average
            #     aver_cov = np.mean([bcdf_combined_agg.at[prot,sample+'D_'+time+'_aggre_coverage'],bcdf_combined_agg.at[prot,sample+'F_'+time+'_aggre_coverage'],
            #                         bcdf_combined_agg.at[prot, sample.split('_')[0]+'_parallel_B_30_aggre_coverage'], bcdf_combined_agg.at[prot, sample.split('_')[0]+'_parallel_C_30_aggre_coverage']])
            #     df_aggre_coverage.at[prot,sample+time+'_ave_aggre_cov']=aver_cov

            aver_cov = np.mean([bcdf_combined_agg.at[prot, sample + 'D_' + time + '_aggre_coverage'],
                                bcdf_combined_agg.at[prot, sample + 'F_' + time + '_aggre_coverage'],
                                bcdf_combined_agg.at[prot, sample.split('_')[0] + '_parallel_B_'+time+'_aggre_coverage'],
                                bcdf_combined_agg.at[prot, sample.split('_')[0] + '_parallel_C_'+time+'_aggre_coverage']])
            df_aggre_coverage.at[prot, sample + time + '_ave_aggre_cov'] = aver_cov
"""
### get matrisome non-aggregated coverage data
"""
columns = ['gene','MW','category','GFP_30_ave','GFP_120_ave','GFP_240_ave','GFP_1080_ave',
           'SNED1_30_ave','SNED1_120_ave','SNED1_240_ave','SNED1_1080_ave']
df_ecm_original = pd.DataFrame(columns=columns,index=df_aggre_coverage.index)
for each_prot in df_aggre_coverage.index:
    df_ecm_original.at[each_prot,'gene'] = df_aggre_coverage.at[each_prot,'gene']
    df_ecm_original.at[each_prot,'MW'] = df_aggre_coverage.at[each_prot,'MW_kDa']
    df_ecm_original.at[each_prot,'category'] = df_aggre_coverage.at[each_prot,'category']
    for sample in ['GFP','SNED1']:
        for time in ['_30','_120','_240','_1080']:
            df_ecm_original.at[each_prot,sample+time+'_ave'] = np.mean([df_summary.at[each_prot,sample+'_seq'+time+'D_coverage'], df_summary.at[each_prot,sample+'_seq'+time+'F_coverage']])


df_ecm_original.to_excel('D:/data/Naba_deep_matrisome/07232021_secondsearch/ecm_original_cov_ave.xlsx')
"""

tryps30_thromly30 = [[-0.25427827],
       [-0.24061336],
       [-0.11392727],
       [-0.21880266],
       [ 0.30268943],
       [ 0.09464207],
       [ 0.05344539],
       [-0.15810011],
       [-0.11902265],
       [ 0.13648006],
       [-0.04272203],
       [-0.24495608],
       [ 0.15547037],
       [-0.34659567],
       [ 0.17824247],
       [ 0.19811662],
       [ 0.04841936],
       [ 0.13629705],
       [ 0.04992884],
       [-0.12792891],
       [-0.08630499],
       [-0.13559616],
       [ 0.09711207],
       [ 0.03342704],
       [ 0.20300269],
       [ 0.20651591],
       [-0.2663266 ],
       [-0.22025561],
       [ 0.11921229],
       [-0.13724782],
       [ 0.24163309],
       [-0.28753683],
       [ 0.15460275],
       [ 0.00044304],
       [-0.26987636],
       [ 0.00049196],
       [ 0.14718263],
       [ 0.17868924],
       [ 0.22994253],
       [ 0.06783383],
       [ 0.08188912],
       [ 0.27047578],
       [-0.2310915 ],
       [ 0.34059188],
       [ 0.19695014],
       [-0.12823421],
       [ 0.13531542],
       [ 0.22831725],
       [ 0.12338001],
       [-0.06771325],
       [ 0.07726893],
       [ 0.3413931 ],
       [ 0.10312411],
       [-0.2004154 ],
       [-0.07169063],
       [ 0.15126434],
       [-0.1215061 ],
       [-0.17293437],
       [ 0.26751548],
       [ 0.17763986],
       [-0.02756455],
       [-0.19060913],
       [ 0.06746905],
       [-0.19573194]]

trypsin_4h_gluc_ON = [[ 0.06185372],
       [-0.0489403 ],
       [ 0.25576034],
       [ 0.1754957 ],
       [ 0.22090636],
       [-0.1939765 ],
       [ 0.16566238],
       [-0.25880274],
       [-0.20459655],
       [ 0.17626667],
       [-0.18189463],
       [ 0.12720956],
       [ 0.05260014],
       [ 0.23518966],
       [ 0.19258352],
       [-0.00854223],
       [-0.08335569],
       [-0.1933515 ],
       [ 0.13378626],
       [-0.06323651],
       [-0.13052419],
       [ 0.12216603],
       [ 0.17204541],
       [ 0.36867517],
       [-0.00590208],
       [-0.20140709],
       [-0.09511861],
       [-0.13541119],
       [-0.24469283],
       [-0.18348819],
       [ 0.2916169 ],
       [-0.24226987],
       [ 0.09453987],
       [ 0.05342854],
       [ 0.19493693],
       [ 0.16669123],
       [-0.32423455],
       [ 0.15493229],
       [-0.2915958 ],
       [ 0.26431865],
       [-0.07924783],
       [ 0.2878523 ],
       [ 0.01921086],
       [ 0.03266589],
       [-0.20596127],
       [-0.11908434],
       [ 0.03257464],
       [-0.10742263],
       [-0.12468215],
       [-0.22366905],
       [ 0.06384924],
       [-0.06732922],
       [-0.0970291 ],
       [-0.07221196],
       [-0.03813521],
       [-0.1681362 ],
       [-0.15765716],
       [ 0.22126731],
       [-0.20611255],
       [ 0.23542571],
       [ 0.08433133],
       [ 0.32636055],
       [ 0.28940657],
       [ 0.0812989 ]]

gluc_4h_tryps_ON = [[ 0.14569467],
       [ 0.01051134],
       [-0.03288183],
       [-0.07269703],
       [-0.01644952],
       [-0.27643105],
       [ 0.15643439],
       [ 0.21276599],
       [-0.15138891],
       [-0.08456214],
       [-0.2311876 ],
       [ 0.006546  ],
       [ 0.07059338],
       [ 0.03060225],
       [-0.00471049],
       [-0.17650129],
       [-0.10978928],
       [ 0.06132037],
       [ 0.17789531],
       [ 0.04732037],
       [-0.19272608],
       [ 0.02033071],
       [ 0.13906789],
       [ 0.18231894],
       [ 0.09107751],
       [-0.04693274],
       [ 0.03712795],
       [ 0.04850171],
       [ 0.1308713 ],
       [ 0.1485826 ],
       [ 0.11506075],
       [-0.22889495],
       [ 0.12708054],
       [-0.01148272],
       [-0.34480333],
       [-0.16133043],
       [-0.07939263],
       [-0.13405783],
       [ 0.23350139],
       [-0.01947309],
       [-0.0062887 ],
       [ 0.03890374],
       [ 0.33134216],
       [ 0.11586735],
       [-0.01691898],
       [-0.18016803],
       [-0.09909447],
       [-0.08338004],
       [ 0.114049  ],
       [-0.11543171],
       [ 0.03300363],
       [-0.08969098],
       [-0.22553116],
       [ 0.01866113],
       [ 0.1621578 ],
       [ 0.11679652],
       [-0.07492209],
       [ 0.21565065],
       [-0.13340582],
       [-0.19700155],
       [-0.12919985],
       [ 0.29137436],
       [-0.01204232],
       [ 0.21496218]]
t_ct_thermo = [[-0.2621706 ],
       [-0.20743947],
       [ 0.25795475],
       [ 0.08462656],
       [ 0.12564605],
       [ 0.23706013],
       [-0.01702294],
       [ 0.13222612],
       [ 0.05541026],
       [ 0.23680282],
       [ 0.05817756],
       [ 0.19201976],
       [ 0.17800626],
       [ 0.17772101],
       [-0.10187836],
       [-0.23260564],
       [-0.17429242],
       [-0.13190171],
       [-0.021524  ],
       [-0.11291677],
       [ 0.12829103],
       [-0.20710167],
       [ 0.18104695],
       [ 0.07732225],
       [-0.14089204],
       [ 0.02850429],
       [-0.0560386 ],
       [ 0.28708136],
       [ 0.04700142],
       [-0.26203647],
       [ 0.19817802],
       [-0.07405377],
       [-0.02122786],
       [ 0.03560288],
       [-0.24678911],
       [-0.18853126],
       [ 0.29557267],
       [-0.2256456 ],
       [-0.1853836 ],
       [-0.04216659],
       [-0.32943222],
       [-0.1361543 ],
       [-0.23742528],
       [ 0.05102573],
       [-0.16966018],
       [ 0.21199822],
       [-0.03109927],
       [-0.27259436],
       [-0.05247786],
       [-0.04079267],
       [ 0.2749393 ],
       [ 0.0089283 ],
       [-0.02832486],
       [-0.18702717],
       [ 0.10101982],
       [ 0.21861266],
       [-0.00534885],
       [ 0.01843625],
       [ 0.26579696],
       [-0.05935792],
       [ 0.09264661],
       [ 0.02370168],
       [ 0.18367328],
       [-0.10578227]]
tryps = [[-0.18721616],
       [-0.10885647],
       [-0.1757808 ],
       [ 0.36322546],
       [ 0.01640958],
       [-0.12486789],
       [-0.0843868 ],
       [-0.21095544],
       [ 0.00902342],
       [ 0.05794371],
       [ 0.16894588],
       [ 0.2292309 ],
       [ 0.1917495 ],
       [ 0.07993867],
       [ 0.04469424],
       [ 0.07576063],
       [ 0.43254188],
       [-0.03541202],
       [ 0.2943967 ],
       [-0.02264838],
       [-0.0802823 ],
       [-0.00207912],
       [-0.01223384],
       [-0.01357396],
       [ 0.13623296],
       [ 0.06126703],
       [-0.19710775],
       [ 0.09469207],
       [-0.14864469],
       [-0.01297503],
       [ 0.01077363],
       [-0.04478023],
       [ 0.02530787],
       [ 0.15275042],
       [-0.03442352],
       [ 0.02739703],
       [-0.00597449],
       [ 0.07921843],
       [ 0.01774445],
       [-0.16323234],
       [-0.0307117 ],
       [-0.02051728],
       [ 0.02376936],
       [-0.19213337],
       [-0.19385347],
       [-0.11493175],
       [-0.21694778],
       [ 0.18189111],
       [-0.22179678],
       [-0.18884704],
       [-0.12789741],
       [-0.08048175],
       [ 0.03035766],
       [ 0.00864127],
       [-0.06431194],
       [ 0.10519992],
       [ 0.08169205],
       [ 0.114754  ],
       [-0.05600975],
       [ 0.03498459],
       [-0.2952882 ],
       [ 0.1771558 ],
       [ 0.07461192],
       [-0.06698076]]
theromlysin = [[ 0.03501371],
       [-0.04993909],
       [-0.0737742 ],
       [ 0.2719467 ],
       [ 0.3139099 ],
       [ 0.18670912],
       [-0.28963792],
       [ 0.14643902],
       [-0.16983579],
       [ 0.16424125],
       [-0.14577188],
       [-0.11260509],
       [-0.04800555],
       [-0.09922897],
       [-0.20814109],
       [-0.09937082],
       [-0.2192867 ],
       [ 0.14642969],
       [-0.26771215],
       [-0.00995268],
       [-0.1402195 ],
       [-0.11223646],
       [ 0.15982765],
       [ 0.23787327],
       [-0.24026497],
       [ 0.05435216],
       [ 0.21948588],
       [ 0.08124904],
       [-0.0960879 ],
       [ 0.22037372],
       [-0.12205642],
       [-0.23647167],
       [-0.16669251],
       [-0.09587744],
       [-0.06309437],
       [-0.31133908],
       [ 0.23469837],
       [ 0.21414195],
       [-0.18830472],
       [-0.07327062],
       [-0.04732347],
       [-0.04472854],
       [-0.26582575],
       [ 0.07738394],
       [ 0.1328111 ],
       [ 0.17445298],
       [-0.22222412],
       [ 0.2776087 ],
       [ 0.09341001],
       [-0.01530282],
       [-0.13621919],
       [ 0.00113079],
       [-0.15217619],
       [ 0.26178318],
       [-0.18347812],
       [ 0.21090843],
       [-0.28332907],
       [-0.04516615],
       [-0.1855555 ],
       [-0.20949367],
       [-0.3372728 ],
       [ 0.33131206],
       [ 0.24309897],
       [ 0.00124967]]

chymotrypsin = [[ 0.26093936],
       [-0.1489989 ],
       [-0.1206305 ],
       [ 0.10191304],
       [-0.01649055],
       [ 0.08949374],
       [ 0.39009655],
       [-0.00191734],
       [ 0.00365547],
       [ 0.31394473],
       [ 0.09328517],
       [-0.14510603],
       [ 0.1749651 ],
       [-0.22823824],
       [ 0.3651351 ],
       [-0.14931786],
       [ 0.22219737],
       [ 0.19605397],
       [ 0.00721515],
       [ 0.22193977],
       [-0.18060632],
       [ 0.23905379],
       [-0.1732659 ],
       [-0.34623   ],
       [-0.04206531],
       [ 0.17990729],
       [-0.20326123],
       [ 0.2608619 ],
       [-0.22318634],
       [-0.05806597],
       [-0.24934703],
       [ 0.20400962],
       [ 0.1803197 ],
       [ 0.09376952],
       [ 0.20621479],
       [-0.22828202],
       [ 0.07505614],
       [-0.23199056],
       [ 0.10674316],
       [-0.07044643],
       [-0.13382387],
       [ 0.0644046 ],
       [ 0.08943709],
       [-0.24374098],
       [ 0.03250491],
       [ 0.18755771],
       [-0.01288175],
       [ 0.22447959],
       [ 0.24800983],
       [ 0.00170066],
       [-0.1455547 ],
       [-0.2531482 ],
       [ 0.06282921],
       [-0.25245628],
       [-0.15641087],
       [-0.13832758],
       [ 0.1614755 ],
       [-0.09253884],
       [ 0.03382566],
       [-0.03550822],
       [-0.3113521 ],
       [-0.26998234],
       [ 0.2722473 ],
       [ 0.13766457]]

trypsin_ct = [[ 2.77665108e-01],
       [-2.38220662e-01],
       [ 1.99107528e-01],
       [-6.34177029e-02],
       [ 2.67960113e-02],
       [ 3.93956376e-04],
       [ 1.65133730e-01],
       [-6.59641102e-02],
       [ 2.72693038e-01],
       [ 1.38325423e-01],
       [-1.22674644e-01],
       [ 1.50862291e-01],
       [ 4.25768569e-02],
       [-1.29140317e-01],
       [-1.78196862e-01],
       [ 2.15546280e-01],
       [ 2.58420706e-01],
       [ 2.28143111e-01],
       [ 1.62787125e-01],
       [-2.49872044e-01],
       [-8.81738737e-02],
       [-1.72886308e-02],
       [-1.74428135e-01],
       [-2.11699903e-01],
       [ 3.41351897e-01],
       [ 1.57609746e-01],
       [ 2.46945441e-01],
       [-1.87780052e-01],
       [-1.16182379e-01],
       [ 1.86055005e-01],
       [ 9.38303098e-02],
       [-8.97483379e-02],
       [-2.78487772e-01],
       [ 1.98103830e-01],
       [-7.81366825e-02],
       [ 1.59778982e-01],
       [ 1.17913596e-01],
       [ 2.33592167e-01],
       [ 3.79535854e-01],
       [ 2.24054784e-01],
       [-1.08986735e-01],
       [ 7.36222416e-02],
       [ 3.97412367e-02],
       [-2.28271663e-01],
       [ 2.69133478e-01],
       [-4.54633199e-02],
       [ 4.07320410e-01],
       [ 3.44230644e-02],
       [-2.34727696e-01],
       [ 3.07734907e-02],
       [ 1.58548936e-01],
       [ 1.71585947e-01],
       [-2.32557148e-01],
       [ 1.35077089e-01],
       [-1.03561506e-01],
       [-1.15663186e-01],
       [-5.61843067e-02],
       [ 8.90924037e-02],
       [ 2.49663576e-01],
       [ 1.73542246e-01],
       [ 7.75029063e-02],
       [ 2.68693000e-01],
       [-1.59885913e-01],
       [-2.55536973e-01]]

trypsin_test = [[-0.28439245],
       [ 0.12382675],
       [ 0.18408537],
       [-0.16984332],
       [-0.26014367],
       [-0.20787323],
       [ 0.14726478],
       [ 0.30891734],
       [ 0.10136802],
       [ 0.02889386],
       [-0.19193599],
       [ 0.154975  ],
       [ 0.23523453],
       [ 0.26107606],
       [-0.09103351],
       [ 0.05374506],
       [ 0.1967715 ],
       [ 0.00039589],
       [ 0.07123954],
       [-0.35233745],
       [ 0.18416175],
       [ 0.2985797 ],
       [-0.04285941],
       [ 0.38506365],
       [-0.19863191],
       [ 0.28397724],
       [-0.19845302],
       [ 0.08306841],
       [-0.2820916 ],
       [ 0.063764  ],
       [ 0.09093735],
       [-0.22576705],
       [ 0.33180264],
       [ 0.20558862],
       [-0.13498731],
       [ 0.2786666 ],
       [ 0.255882  ],
       [-0.15717134],
       [ 0.37603143],
       [ 0.33959106],
       [-0.22634883],
       [ 0.26674923],
       [-0.21687336],
       [ 0.19923054],
       [-0.26734075],
       [ 0.12094574],
       [-0.19510731],
       [ 0.12543595],
       [-0.16770254],
       [ 0.24199992],
       [-0.31275517],
       [-0.13233018],
       [ 0.13737296],
       [-0.20145896],
       [-0.00913237],
       [ 0.12110209],
       [-0.24076714],
       [-0.00266612],
       [ 0.27964082],
       [-0.21979076],
       [-0.124397  ],
       [ 0.0617185 ],
       [-0.01911728],
       [ 0.2077897 ]]

tryps30_thromly30 = np.reshape(np.array(tryps30_thromly30),64)
trypsin_4h_gluc_ON = np.reshape(np.array(trypsin_4h_gluc_ON),64)
gluc_4h_tryps_ON = np.reshape(np.array(gluc_4h_tryps_ON),64)
t_ct_thermo = np.reshape(np.array(t_ct_thermo),64)
tryps = np.reshape(np.array(tryps),64)
theromlysin = np.reshape(np.array(theromlysin),64)
trypsin_thermo_combined = np.max([tryps,theromlysin],axis=0)
chymotrypsin = np.reshape(np.array(chymotrypsin),64)
trypsin_ct = np.reshape(np.array(trypsin_ct),64)
trypsin_ct_combined = np.max([tryps,chymotrypsin],axis=0)

trypsin_test = np.reshape(np.array(trypsin_test),64)
for i in [tryps30_thromly30,trypsin_4h_gluc_ON,gluc_4h_tryps_ON,t_ct_thermo,tryps,theromlysin,trypsin_thermo_combined,chymotrypsin,trypsin_ct,trypsin_ct_combined]:
    cosine = cosine_sim_calculating(chymotrypsin,i)
    euclidean = np.linalg.norm(chymotrypsin-i)
    print (cosine,euclidean, euclidean/cosine)


