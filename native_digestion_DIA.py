"""data processing and analysis of native limited proteolysis, DIA-NN output, written by blackjack, 01/27/2023
"""
import pandas as pd
import numpy as np
import pickle as pk
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

def peptide_getter():
    """
    output peptide list for each time point in a dictionary
    :return:
    """
    file_peptide_dict = defaultdict(set)
    df = pd.read_csv('F:/native_digestion/01242023/time_points/report.tsv',sep='\t')
    file_list, psm_list = df['File.Name'], df['Stripped.Sequence']
    for f, psm in zip(file_list,psm_list):
        file_peptide_dict[f.split('\\')[-1].split('_rerun.mzML')[0]].add(psm)
    print ([k for k in file_peptide_dict])
    return pk.dump(file_peptide_dict, open('F:/native_digestion/01242023/time_points/f_peptides_dict.p','wb'))


def get_unique_peptides():
    """
    get unique peptides at each time point
    :return:
    """
    unique_peptides_dict = {}
    peptide_list = []
    f_list = ['tryps_0005min','tryps_0010min','tryps_0015min','tryps_0020min','tryps_0030min','tryps_0040min',
              'tryps_0050min','tryps_0060min','tryps_0120min','tryps_0180min','tryps_0240min','tryps_1440min',
              'tryps_leftover'] # order by time
    f_peptide_dict = pk.load(open('F:/native_digestion/01242023/time_points/f_peptides_dict.p','rb'))
    for f in f_list:
        unique_pep_list = [each for each in f_peptide_dict[f] if each not in peptide_list]
        peptide_list += unique_pep_list
        unique_peptides_dict[f] = unique_pep_list
    print ([(each, len(unique_peptides_dict[each])) for each in unique_peptides_dict])
    # pk.dump(unique_peptides_dict,open('F:/native_digestion/01242023/time_points/f_unique_peptides_dict.p','wb'))


def get_protein_set():
    """
    get total protein list
    :return:
    """
    protein_set = set()
    df = pd.read_csv('F:/native_digestion/01242023/time_points/report.pg_matrix.tsv',sep='\t')
    protein_list = df['Protein.Ids']
    for prot_str in protein_list:
        if ';' in prot_str:
            for prot in prot_str.split(';'):
                protein_set.add(prot)
        else:
            protein_set.add(prot_str)
    print (len(protein_set), list(protein_set)[:10])
    # pk.dump(protein_set,open('F:/native_digestion/01242023/time_points/proteinid_set.p','wb'))


def protein_intensity():
    """
    get peptide intensity vector from pg report
    :return:
    """
    protein_int_dict = {}
    df = pd.read_csv('F:/native_digestion/01242023/time_points/report.pg_matrix.tsv', sep='\t')
    for row in df.itertuples(index=False):
        protein_id_list = [row[1]] if ';' not in row[1] else row[1].split(';')
        int_array = np.array(row[-13:])
        for each in protein_id_list:
            protein_int_dict[each] = int_array
    pk.dump(protein_int_dict,open('F:/native_digestion/01242023/time_points/prot_intensity_dict.p','wb'))


def combine_distance_intensity():
    """
    combine cleavage to center distance and protein intensity by multiplying
    :return:
    """
    # protein_int_dict = pk.load(open('F:/native_digestion/01242023/time_points/prot_intensity_dict.p','rb'))
    prot_int_norm_df = pd.read_excel(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/intensity.xlsx',index_col=0)
    prot_int_norm_list = prot_int_norm_df.index.tolist()
    distance_df = pd.read_csv(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/distance_to_center_all_filter.tsv',sep='\t',index_col=0)
    distance_df = distance_df.fillna(0)
    new_df = pd.DataFrame(index=distance_df.index,columns=distance_df.columns)
    for row in distance_df.itertuples(index=True):
        prot = row[0]
        int_array = prot_int_norm_df.loc[prot,:].to_numpy() if prot in prot_int_norm_list else np.zeros(7)
        # int_array_fillna = np.where(np.isnan(int_array),0,int_array)
        normalized_int = (int_array-np.min(int_array))/(np.max(int_array)-np.min(int_array))
        multiply = np.multiply(np.array(row[1:]), normalized_int)
        # multiply = np.multiply(np.array(row[1:]), int_array)
        new_df.loc[prot,:] = multiply
    new_df.to_csv(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/distance_times_norm_intensity.tsv',sep='\t')


def normalize_distance():
    distance_df = pd.read_csv(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/distance_times_norm_intensity.tsv', sep='\t',index_col=0)
    normal_dist_df = pd.DataFrame(index=distance_df.index, columns=distance_df.columns[:-1])
    dist_2d = distance_df.to_numpy()[:,:-1]
    for prot, dist in zip(normal_dist_df.index, dist_2d):
        dist = np.where(np.isnan(dist), 0, dist)
        norm_dist = (dist-np.min(dist))/(np.max(dist)-np.min(dist))
        normal_dist_df.loc[prot,:] = norm_dist
    normal_dist_df = normal_dist_df.dropna()
    normal_dist_df.to_csv(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/distance_times_norm_intensity_norm01.tsv',sep='\t')


def output_intensity_totsv():
    distance_df = pd.read_excel('F:/native_digestion/01242023/analysis/distance_to_center_all.xlsx', index_col=0)
    protein_int_dict = pk.load(open('F:/native_digestion/01242023/time_points/prot_intensity_dict.p', 'rb'))
    int_df = pd.DataFrame(index=distance_df.index, columns=distance_df.columns.tolist()[:-2])
    for row in distance_df.itertuples(index=True):
        prot = row[0]
        int_array = protein_int_dict[prot][:-2]
        int_array_fillna = np.where(np.isnan(int_array), 0, int_array)
        normalized_int = (int_array_fillna - np.min(int_array_fillna)) / (
                    np.max(int_array_fillna) - np.min(int_array_fillna))
        int_df.at[prot,:] = normalized_int
    int_df = int_df.dropna()
    int_df.to_csv('F:/native_digestion/01242023/analysis/norm_intensity_removelast2.tsv',sep='\t')


def distance_intensity_():
    """
    plot distance on y axis, time on x axis, color as intensity
    :return:
    """

    protein_int_dict = pk.load(open('F:/native_digestion/01242023/time_points/prot_intensity_dict.p','rb'))
    distance_df = pd.read_csv('F:/native_digestion/01242023/analysis/distance_to_center_all_filter_nomorethan3zeros.tsv',sep='\t',index_col=0)
    new_df = pd.DataFrame(columns=['protein_id','average_distance','intensity','time'])
    protein_id_list = []
    for prot in distance_df.index:
        protein_id_list+=[prot]*len(distance_df.columns)
    ave_dist = distance_df.to_numpy().flatten()
    int_total_array = np.concatenate([protein_int_dict[prot] for prot in distance_df.index])
    time_list = list(range(1,len(distance_df.columns)+1))*distance_df.shape[0]
    for k,v in zip(['protein_id','average_distance','intensity','time'],
                   [protein_id_list,ave_dist,int_total_array,time_list]):
        new_df[k] = v
    new_df.to_csv('F:/native_digestion/01242023/analysis/prot_distance_intensity_1.tsv',sep='\t')


def distance_intensity_plot():
    """
    plot heat map, use data from distance_intensity_
    :return:
    """
    from functools import reduce
    from operator import add
    df = pd.read_csv('F:/native_digestion/01242023/analysis/prot_distance_intensity_2.tsv', sep='\t',index_col=0)
    df = df.copy().fillna(0)
    # df = df.copy().dropna()
    normalize_factor = 50
    time_range = list(range(1,14))
    # time_index_dict = {i:j for i,j in zip(time_range,chunk_(list(range(1,normalize_factor)),int(normalize_factor/13))[:-1])}
    # time_index_join = reduce(add,[v for v in time_index_dict.values()])
    # print (time_index_dict)
    plot_array = np.zeros((normalize_factor,13))
    prot_unique = df['protein_id'].unique()
    # df_plot = pd.DataFrame(columns=['norm_dist_index','int','time'])
    # norm_dist_index, int_color, time = [],[],[]
    for each in prot_unique:

        sub_df = df[df['protein_id']==each]
        dist_array = sub_df['deltaD'].to_numpy()
        # normalize cleavage to center distance

        normalized_dist = (dist_array - np.min(dist_array)) / (
                    np.max(dist_array) - np.min(dist_array)) * (normalize_factor-1)
        # dist_index = (normalize_factor-1)-normalized_dist  # used as the y-axis index on heatmap
        dist_index = normalized_dist
        # dist_index_list = [[j]*len(time_index_dict[i+1]) for i, j in enumerate(dist_index)]
        # dist_index_list = reduce(add,dist_index_list)

        # time as index on x-axis, intensity as color
        int_array,time_array = sub_df['intensity'], sub_df['time']
        # int_index_list = [[j]*len(time_index_dict[i+1]) for i, j in enumerate(int_array)]
        # int_index_list = reduce(add,int_index_list)
        # norm_dist_index+=dist_index_list
        # int_color+=int_index_list
        # time+=time_index_join

        for dist, intensity, time in zip(dist_index, int_array, time_array):
            try:

                plot_array[int(dist),time-1] += intensity
            except ValueError:
                print (each,dist_array)

    # for i, j in zip(['norm_dist_index','int','time'],[norm_dist_index,int_color,time]):
    #     df_plot[i]=j

    # print(df_plot.shape)
    # df_plot.to_csv('F:/native_digestion/01242023/analysis/scatter_plot_0225.tsv',sep='\t')
    # g = sns.scatterplot(data=df_plot,x='time',y='norm_dist_index',hue='int',s=5,palette='viridis')
    # plt.legend([], [], frameon=False)
    g = sns.heatmap(data=plot_array[:,:-2],cmap='viridis',robust=True)
    plt.show()


def delta_dist_cal():
    # calcuate D'= Dmax-D
    distance_dict = pk.load(open('F:/native_digestion/01242023/time_points/to_center_distance_dict.pkl','rb'))
    max_dict_dict = {prot:max([v for v in distance_dict[prot].values()]) for prot in distance_dict}
    prot_dist_int_df = pd.read_csv('F:/native_digestion/01242023/analysis/prot_distance_intensity_1.tsv',sep='\t',index_col=0)
    delta_dist = []
    for dist, prot in zip(prot_dist_int_df['average_distance'],prot_dist_int_df['protein_id']):
        delta = 0.95*max_dict_dict[prot]-dist
        if delta<0:
            delta_dist.append(0)
        else:
            delta_dist.append(delta)
    prot_dist_int_df['deltaD'] = delta_dist
    prot_dist_int_df.to_csv('F:/native_digestion/01242023/analysis/prot_distance_intensity_2.tsv',sep='\t')


def chunk_(list, n):
    geneartor = (list[i:n + i] for i in range(0, len(list), n))
    return [i for i in geneartor]


def filter_df():
    # filter out all zeros in df
    import time
    time.sleep(3)
    df = pd.read_excel('D:/data/native_protein_digestion/12072021/heat_shock_ionquant_MBR/distance_to_center_all.xlsx',index_col=0)
    df = df.copy().fillna(0)
    data = []
    index = []
    for row in df.itertuples():
        # filter rows with all 0s
        if np.count_nonzero([row[i]==0 for i in range(-7,0)])>=4:  # if there are more than 6 zeros
        # if all([row[i] == "0" for i in range(1, len(row))]):
            continue
        else:
            data.append([i for i in row][1:])
            index.append(row[0])
    new_df = pd.DataFrame(data, columns=df.columns,index=index)
    # after filtering, 2293 proteins retained (originally 3886)
    new_df.to_csv('D:/data/native_protein_digestion/12072021/heat_shock_ionquant_MBR/distance_to_center_all_filter.tsv',sep='\t')


def cluster_map():
    import matplotlib.pyplot as plt
    import seaborn as sns

    data_heatshock = pd.read_csv(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/distance_times_norm_intensity_norm01.tsv',sep='\t',index_col=0)
    protein_list = data_heatshock.index.tolist()
    data = pd.read_csv(r'F:\native_digestion\01242023\analysis/distance_to_center_times_normIntensity_filter_norm01_removelast2.tsv',sep='\t',index_col=0)
    data = data[data.index.isin(protein_list)].dropna()
    protein_list = data.index.tolist()
    protein_hex_list = disorder_to_hex(protein_list)  # add row colors on clustermap
    plddt_hex_list = plddt_to_hex(protein_list)
    row_colors = pd.DataFrame(dict(disorder=protein_hex_list,plddt=plddt_hex_list))
    prot_disorder_dict = pk.load(open('F:/native_digestion/01242023/analysis/prot_disorder_dict.p', 'rb'))

    # new_data = pd.DataFrame(index=data.index, columns=data.columns.to_list()[:-2])
    # normalize each protein between 0 and 1
    # for prot, array in zip(data.index, data.to_numpy()[:,:-2]):
    #     normalize_0_1 = (array-np.min(array))/(np.max(array)-np.min(array))
    #     new_data.loc[prot,:] = normalize_0_1
    # new_data.to_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity_filter_norm01_removelast2.tsv', sep='\t')

    fig, ax = plt.subplots(1, 1, figsize=(12,8))
    g = sns.clustermap(data=data.reset_index(drop=True),col_cluster=False,cmap="YlGnBu",yticklabels=False,row_colors=row_colors, figsize=(10,8),dendrogram_ratio=0.15,
                       cbar_kws={"orientation": "horizontal","shrink": 0.5})
    re_order_protein_index = g.dendrogram_row.reordered_ind
    re_order_protein_list = [protein_list[each] for each in re_order_protein_index]
    print ([(each, prot_disorder_dict[each]) if each in prot_disorder_dict else (each, 'None') for each in re_order_protein_list])
    # plt.tight_layout()
    # plt.show()
    # plt.savefig(r'D:\data\native_protein_digestion\12072021\heat_shock_ionquant_MBR/heatshock_cluster_disorder_plddt.png',dpi=300)


def disorder_to_hex(protein_list):
    """
    convert disorder ratio to hex color string, high disorder will show as darker red
    :return:
    """

    base_rgb = (255,255,255)
    prot_disorder_dict = pk.load(open('F:/native_digestion/01242023/analysis/prot_disorder_dict.p','rb'))
    protein_hexcolor_list = ['#%02x%02x%02x' % (255, int(210*(1-prot_disorder_dict[each])),int(210*(1-prot_disorder_dict[each])))
                             if each in prot_disorder_dict else '#b3b3b3' for each in protein_list]
    print (protein_hexcolor_list[:5])
    # df = pd.DataFrame(index=protein_list,columns=['Disorder'])
    # df['Disorder'] = protein_hexcolor_list
    return protein_hexcolor_list


def plddt_to_hex(protein_list):

    pal = sns.color_palette("magma",24)
    pal_list = pal.as_hex()
    print (pal_list)
    plddt_mean_dict = pk.load(open('D:/data/alphafold_pdb/pLDDT_human_mean.pkl','rb'))
    plddt_hex_list = [pal_list[int(plddt_mean_dict[each]/100*23)] if each in plddt_mean_dict else '#b3b3b3' for each in protein_list]
    # print (plddt_hex_list)
    # df = pd.DataFrame(index=protein_list,columns=['plddt'])
    # df['plddt'] = plddt_hex_list
    return plddt_hex_list


def umap_hdbscan():
    import umap
    import hdbscan
    data = pd.read_csv('F:/native_digestion/01242023/analysis/distance_to_center_times_normIntensity_filter_norm01.tsv',
                       sep='\t',index_col=0)
    clusterable_embedding = umap.UMAP(
        n_neighbors=20,
        min_dist=0.1,
        n_components=2,
        random_state=42,
    ).fit_transform(data)
    plt.scatter(clusterable_embedding[:, 0], clusterable_embedding[:, 1],
                 s=1.5,alpha=0.8)
    labels = hdbscan.HDBSCAN(
        min_samples=20,
        min_cluster_size=200, # 50
    ).fit_predict(clusterable_embedding)
    # print (labels)
    clustered = (labels >= 0)  # label = -1 for noise data
    # for i, j, k in zip(data.index,labels, clusterable_embedding[clustered,0]):
    #     if j == 0:
    #         print (i, data.loc[i,:],k)
    plt.scatter(clusterable_embedding[~clustered, 0],
                clusterable_embedding[~clustered, 1],
                color=(0.5, 0.5, 0.5),
                s=3,
                alpha=0.6)
    plt.scatter(clusterable_embedding[clustered, 0],
                clusterable_embedding[clustered, 1],
                c=labels[clustered],
                s=3,
                cmap='Spectral')

    plt.show()


def ion_quant_analysis():
    from collections import Counter
    df = pd.read_csv('F:/native_digestion/01202023/ionquant_search/combined_ion.tsv',sep='\t')
    all_have_prot, at_least6, at_least5, at_least4 = [],[],[],[]
    for row in df.itertuples():
        num_nonzero = np.count_nonzero(row[-13:])
        prot_id = row[12]
        print (prot_id)
        if num_nonzero == 13:
            all_have_prot.append(prot_id)
        if num_nonzero >= 12:
            at_least6.append(prot_id)
        if num_nonzero >= 11:
            at_least5.append(prot_id)
        if num_nonzero >= 10:
            at_least4.append(prot_id)
    print (len(set(all_have_prot)),len(set(at_least6)),len(set(at_least5)),len(set(at_least4)))


def cluster_analysis():
    """
    analyze results from cluster map
    :return:
    """
    from scipy import stats
    df = pd.read_excel('F:/native_digestion/01242023/analysis/cluster_heatshock_overlap_output.xlsx',index_col=0)
    plddt_mean_dict = pk.load(open('D:/data/alphafold_pdb/pLDDT_human_mean.pkl', 'rb'))

    clusters = ['cluster_1','cluster_2']
    cluster_1,cluster_2 = df.iloc[:150,:], df.iloc[150:,:]  # manually cut 1400
    cluster_1,cluster_2 = cluster_1[cluster_1['disorder_ratio']!='None'], cluster_2[cluster_2['disorder_ratio']!='None']
    clus1_disorder, clus2_disorder = cluster_1['disorder_ratio'].tolist(), cluster_2['disorder_ratio'].tolist()
    # clus1_disorder, clus2_disorder = [plddt_mean_dict[prot] for prot in cluster_1.index], [plddt_mean_dict[prot] for prot in cluster_2.index]
    t_test = stats.ttest_ind(clus1_disorder,clus2_disorder)
    df_plot = pd.DataFrame(columns=['cluster','disorder_ratio'])
    df_plot['cluster'] = ['cluster_1']*cluster_1.shape[0]+['cluster_2']*cluster_2.shape[0]
    df_plot['disorder_ratio'] = clus1_disorder+clus2_disorder
    print (f't test result {t_test}')
    fig,ax = plt.subplots(1, 1, figsize=(6, 5))
    box_plot = ax.boxplot([df_plot[df_plot.cluster==clu]['disorder_ratio'] for clu in clusters],
               positions=range(len(clusters)), showfliers=False,
               medianprops=dict(linestyle='-', linewidth=3,color='k'), widths=[0.4,0.4],
                          boxprops=dict(linestyle='-', linewidth=3, color='k'),
                          whiskerprops=dict(linestyle='-', linewidth=3),
                          capprops=dict(linestyle='-', linewidth=3))
    sns.swarmplot(x='cluster', y='disorder_ratio', data=df_plot,size=3)

    plt.show()


if __name__ == '__main__':
    import seaborn as sns
    # get_unique_peptides()
    # get_protein_set()
    # print (protein_intensity())
    # combine_distance_intensity()
    # filter_df()
    # cluster_map()
    # umap_hdbscan()
    # distance_intensity_()
    # distance_intensity_plot()
    # ion_quant_analysis()
    # distance_intensity_plot()
    # delta_dist_cal()
    # output_intensity_totsv()
    # normalize_distance()
    cluster_analysis()

