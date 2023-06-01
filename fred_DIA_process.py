"""
process ECM time lapsed DIA data, 12/1/2022, Blackjack at UIC
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from collections import defaultdict
import pickle as pk
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
from protein_coverage import fasta_reader
from statannot import add_stat_annotation
import json
from functools import reduce
import operator
from protein_coverage import fasta_reader_gene
import time
from tsv_reader import venn_diagram_gen, venn_diagram_gen2

ecm_class_color_dict = {"ECM Glycoproteins":"#133463","Collagens":"#0584B7","Proteoglycans":"#59D8E6",
                        "ECM-affiliated Proteins":"#F4511E","ECM Regulators":"#F9A287","Secreted Factors":"#FFE188"}

sort_category = ["ECM Glycoproteins","Collagens","Proteoglycans","ECM-affiliated Proteins","ECM Regulators",
                  "Secreted Factors"]

def qc_check():
    """
    check coefficient variation of gene-centric intensities between replicates, use output from DIANN report.gg_matrix
    :return:
    """
    gene_int_dict = defaultdict(list)
    base_path = 'F:/fred_time_lapse/'
    # intensity_files = glob(base_path+'search*/report.gg_matrix*.tsv')
    intensity_files = glob(base_path+'20230503/report.pg_matrix_gene.tsv')
    for each_f in intensity_files:
        print (each_f)
        df = pd.read_csv(each_f, sep='\t').fillna(0)
        columns = df.columns.tolist()[1:]
        for each_col in columns:
            for gene, intensity in zip(df['Genes'], df[each_col]):
                if gene != 0:
                    gene_list = [gene] if ';' not in gene else gene.split(';')
                    for g in gene_list:
                        gene_int_dict[g].append((each_col.split('\\')[-1].split('.mzML')[0],intensity))
    gene_f_int_dict = {}
    for gene in gene_int_dict:
        f_int_dict = defaultdict(float)
        for tp in gene_int_dict[gene]:
            f_int_dict[tp[0]]=tp[1]
        gene_f_int_dict[gene] = f_int_dict
    pk.dump(gene_f_int_dict,open(base_path+'20230503/analysis/gene_f_intesity_dict_of_dict_0524.p','wb'))

    # calculate cv (coefficient variance) of gene intensities among 5 replicates,normalize between 0 and 1 first
    """
    gene_list = [gene for gene in gene_f_int_dict.keys()]
    new_df = pd.DataFrame(index=gene_list)
    time_points = ['15', '30', '60', '120', '240']
    for each in gene_list:
        for time in time_points:
            for cell in ['144','145']:
                print (each,time,cell)
                intensity = np.array([gene_f_int_dict[each][cell+rep+'_'+time] for rep in ['A','B','C','D','E']])
                cv = np.std(intensity)/np.mean(intensity)*100
                new_df.loc[each,cell+'_'+time+'_CV'] = cv
        # process 1080 time points
        for cell in ['144','145']:
            intensity1080 = np.array([gene_f_int_dict[each][cell+rep+'_1080'] for rep in ['NA','NB','NC','ND','NE']])
            cv_1080 = np.std(intensity1080)/np.mean(intensity1080)*100
            new_df.loc[each,cell+'_1080_CV'] = cv_1080
    new_df.to_csv(base_path+'replicates_gene_intensity_cv.tsv',sep='\t')
    """

def cv_box_plot():
    """
    plot CV% in box or violin plot
    :return:
    """
    base_path = 'F:/fred_time_lapse/'
    df = pd.read_csv(base_path+'replicates_gene_intensity_cv.tsv',index_col=0,sep='\t')
    columns = df.columns.tolist()
    df_plot = pd.DataFrame(columns=['cv','sample'])
    df_plot['cv'] = df.to_numpy().flatten()

    df_plot['sample'] = columns*df.shape[0]

    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    sns.violinplot(data=df_plot, x='sample', y='cv',ax=ax)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_ylabel('CV%', fontsize=20)
    ax.set_xlabel('Sample',fontsize=20)
    plt.xticks(rotation=30)
    plt.show()


def pr_matrix_reader():
    """
    read precursor matrix report tsv from DIANN
    :return:
    """
    base_path = 'F:/fred_time_lapse/'
    prot_info_dict = defaultdict(list)  # {'uniprotid':[(file,psm),(file,psm)]}
    # psm_f = glob(base_path+'search*/report.pr_matrix*.tsv')
    psm_f = glob(base_path+'2023_05_25/report.pr_matrix_gene.tsv')
    for each_f in psm_f:
        print (each_f)
        df = pd.read_csv(each_f, sep='\t').fillna('null')
        columns = [c for c in df.columns if '.mzML' in c]
        # iterate each row, grab protein, file and psm
        for i in range(df.shape[0]):
            # protein_id_cell = df.loc[i,'Protein.Group']
            protein_id_cell = df.loc[i, 'Genes'] # gene-centric
            prot_id_list = [protein_id_cell] if ';' not in protein_id_cell else protein_id_cell.split(';')
            # read file and psm in tuple (file, psm) when intensity is not zero
            file_psm = [(col.split('\\')[-1].split('.mzML')[0], df.loc[i,'Stripped.Sequence'])
                        for col in columns if df.loc[i,col] != 'null']
            for prot in prot_id_list:
                prot_info_dict[prot]+=file_psm
    pk.dump(prot_info_dict,open(base_path+'2023_05_25/analysis/gene_file_psm_tuple_dict_0601.p','wb'))


def aggregate_psms():
    """
    aggregate PSMs between different time
    :return:
    """

    prot_f_psm_dict_of_dict = {}
    replicates_combined_dict = {}  # {'uniprotid':{'144_15':{combined_psm_set},'144_30':{combined_psm_set}}}
    aggregate_pep_dict = {}  # {'uniprotid':{'144_15':{aggregated_psm_set},'144_30':{aggregate_psm_set}}}

    # {'uniprotid':[(file,psm),(file,psm)]}
    prot_f_psm_tuple_dict = pk.load(open('F:/fred_time_lapse/2023_05_25/analysis/gene_file_psm_tuple_dict_0601.p','rb'))

    for prot in prot_f_psm_tuple_dict:
        f_psm_dict = defaultdict(list)
        for tup in prot_f_psm_tuple_dict[prot]:
            f_psm_dict[tup[0]].append(tup[1])  
        prot_f_psm_dict_of_dict[prot] = f_psm_dict
    # print (prot_f_psm_dict_of_dict['Q3TW96'])
    pk.dump(prot_f_psm_dict_of_dict,open('F:/fred_time_lapse/2023_05_25/analysis/gene_f_psm_dict_of_dict_0601.p','wb'))

    # combine all psms among replicates (A to E), union all psms in replicates for now for better coverage
    """
    prot_f_psm_dict_of_dict = pk.load(open('F:/fred_time_lapse/20230508/analysis/gene_f_psm_dict_of_dict_0509.p','rb'))
    time_points = ['15', '30', '60', '120', '240']
    for prot in prot_f_psm_dict_of_dict:
        replicates_combined_psm_dict = {}
        for time in time_points:
            for sample in ['144','145']:
                total_psm_set = set([psm for rep in ['A','B','C','D','E'] for psm in prot_f_psm_dict_of_dict[prot][sample+rep+'_'+time]
                                     if len(set(prot_f_psm_dict_of_dict[prot][sample+rep+'_'+time]))>1]) # num of peptides>1
                replicates_combined_psm_dict[sample+'_'+time] = total_psm_set
        for sample in ['144','145']:
            total_psm_set1080 = set([psm for rep in ['NA','NB','NC','ND','NE'] for psm in prot_f_psm_dict_of_dict[prot][sample+rep+'_1080']
                                     if len(set(prot_f_psm_dict_of_dict[prot][sample+rep+'_1080']))>1])
            replicates_combined_psm_dict[sample+'_1080'] = total_psm_set1080
        replicates_combined_dict[prot] = replicates_combined_psm_dict
    pk.dump(replicates_combined_dict, open('F:/fred_time_lapse/20230508/analysis/gene_f_rep_combined_peptide_dict_0509.p', 'wb'))

    # aggregate psms throughout time points
    for prot in replicates_combined_dict:
        aggregate_psm_dict = {}
        for ind,val in enumerate(time_points):
            time_points_slice = time_points[:ind+1]
            for sample in ['144','145']:
                # aggregate psms from 15min to 240min
                aggregate_psm_set = reduce(operator.or_,[replicates_combined_dict[prot][sample+'_'+t] for t in time_points_slice])
                aggregate_psm_dict[sample+'_'+val+'_aggregate'] = aggregate_psm_set
        aggregate_pep_dict[prot] = aggregate_psm_dict

    pk.dump(aggregate_pep_dict,open('F:/fred_time_lapse/20230508/analysis/gene_f_aggregate_peptide_dict_0509.p','wb'))
    """


def combine_psms_fraction():
    # for fractionation samples
    replicates_combined_dict = {}
    prot_f_psm_dict_of_dict = pk.load(open('F:/fred_time_lapse/2023_05_25/analysis/gene_f_psm_dict_of_dict_0601.p', 'rb'))
    # rep_list = ['145_nA1080_1_6', '145_nA1080_2_7', '145_nA1080_3_8', '145_nA1080_4_9', '145_nA1080_5_10']
    rep_list, bio_list = ['1_6','2_7','3_8','4_9','5_10'], ['B', 'C', 'D', 'E']
    for prot in prot_f_psm_dict_of_dict:
        replicates_combined_psm_dict = {}
        for bio in bio_list:
            total_psm_set1080 = set([psm for rep in rep_list for psm in
                                     prot_f_psm_dict_of_dict[prot][bio+rep]
                                     if len(set(prot_f_psm_dict_of_dict[prot][bio+rep])) > 1])
            replicates_combined_psm_dict[bio+'_frac_combined'] = total_psm_set1080
        replicates_combined_dict[prot] = replicates_combined_psm_dict
    pk.dump(replicates_combined_dict,
            open('F:/fred_time_lapse/2023_05_25/analysis/gene_f_rep_combined_peptide_dict_0601.p', 'wb'))


def coverage_calculation():
    """
    calculate aggregated coverage and 1080 coverage
    :return:
    """
    from protein_coverage import fasta_reader,fasta_reader_gene
    from tsv_reader import protein_info_from_fasta
    # load fasta information
    protein_seq_dict = fasta_reader_gene('F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')
    protein_info_dict = protein_info_from_fasta('F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')

    # load aggregated and combined 1080 peptides
    prot_f_aggregate_pep_dict = pk.load(open('F:/fred_time_lapse/20230508/analysis/gene_f_aggregate_peptide_dict_0509.p','rb')) # 25 files, 5*5
    # prot_f_aggregate_pep_dict = pk.load(open('F:/fred_time_lapse/analysis/ecm_gene_f_aggregate_peptide_dict_best_rep_0327.p', 'rb')) # 5 files total, median replicate for each time point
    prot_f_combined_pep_dict = pk.load(open('F:/fred_time_lapse/2023_05_25/analysis/gene_f_rep_combined_peptide_dict_0601.p','rb')) # key is '144_1080' and '145_1080'

    intensity_dict = pk.load(open('F:/fred_time_lapse/20230508/analysis/gene_f_intesity_dict_of_dict_0509.p', 'rb'))
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r'))  # ECM genes
    # calculate coverage
    time_points, samples = ['15', '30', '60', '120', '240'], ['144','145']
    # samples_1080 = ['144_1080','145_1080']
    samples_1080 = [bio+'_frac_combined' for bio in ['B','C','D','E']] # fractionation sample combined
    aggre_cov_df = pd.DataFrame()
    for prot in prot_f_combined_pep_dict:
        # if prot != 'null':
        if prot in gene_category_dict: # filter non-ECM genes out
            print (prot)
            seq = protein_seq_dict[prot]
            prot_length = len(seq)
            # aggre_cov_df.at[prot,'gene'] = protein_info_dict[prot][0]
            # aggre_cov_df.at[prot,'description'] = protein_info_dict[prot][1]
            aggre_cov_df.loc[prot, 'category'] = gene_category_dict[prot]["Category"]
            aggre_cov_df.loc[prot, 'Sub'] = gene_category_dict[prot]["Sub"]
            aggre_cov_df.at[prot, 'length'] = len(protein_seq_dict[prot])
            # load peptides set
            # for sample in samples:
            #     for time in time_points:
            #         np_arry = np.zeros(prot_length)
            #         peptides_set = prot_f_aggregate_pep_dict[prot][sample+'_'+time+'_aggregate']
            #         for pep in peptides_set:
            #             pep_loc = seq.find(pep)
            #             pep_end_loc = pep_loc+len(pep)
            #             np_arry[pep_loc:pep_end_loc]+=1
            #         agg_cov = np.count_nonzero(np_arry)/prot_length*100
            #         aggre_cov_df.at[prot,sample+'_'+time+'_aggre_cov'] = agg_cov
                    # aggre_cov_df.at[prot,sample+'_'+time+'_aggre_pep'] = str(peptides_set)
            for sample2 in samples_1080:
                np_arry = np.zeros(prot_length)
                peptides1080_set = prot_f_combined_pep_dict[prot][sample2]
                for pep in peptides1080_set:
                    pep_loc = seq.find(pep)
                    pep_end_loc = pep_loc + len(pep)
                    np_arry[pep_loc:pep_end_loc] += 1
                cov1080 = np.count_nonzero(np_arry)/prot_length*100
                aggre_cov_df.at[prot,sample2+'_cov'] = cov1080
                # aggre_cov_df.at[prot,sample2+'_pep'] = str(peptides1080_set)
    aggre_cov_df.to_csv('F:/fred_time_lapse/2023_05_25/analysis/ecm_gene_fraction_combined_cov_0601.tsv',sep='\t')


def abs_coverage_calculation():
    """
    calculate the absolute sequence coverage for time series data without aggregation
    :return:
    """
    protein_seq_dict = fasta_reader_gene(
        'F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r'))  # ECM genes
    gene_rep_combined_psm_dict = pk.load(open('F:/fred_time_lapse/20230508/analysis/gene_f_rep_combined_peptide_dict_0509.p', 'rb'))
    samples, times = ['144','145'], ['15', '30', '60', '120', '240']
    cov_df = pd.DataFrame()
    for prot in gene_rep_combined_psm_dict:
        if prot in gene_category_dict:
        # if prot != 'null':
            print(prot)
            seq = protein_seq_dict[prot]
            prot_length = len(seq)
            cov_df.loc[prot, 'category'] = gene_category_dict[prot]["Category"]
            cov_df.loc[prot, 'Sub'] = gene_category_dict[prot]["Sub"]
            cov_df.at[prot, 'length'] = len(protein_seq_dict[prot])
            for sample in samples:
                for time in times:
                    np_array = np.zeros(prot_length)
                    peptide_set = gene_rep_combined_psm_dict[prot][sample + '_' + time]
                    for pep in peptide_set:
                        pep_loc = seq.find(pep)
                        pep_end_loc = pep_loc + len(pep)
                        np_array[pep_loc:pep_end_loc] += 1
                    cov = np.count_nonzero(np_array) / prot_length * 100
                    cov_df.at[prot, sample + '_' + time + '_abs_cov'] = cov
    cov_df.to_csv('F:/fred_time_lapse/20230508/analysis/ECM_gene_time_series_absoluteCov_0521.tsv',sep='\t')


def abs_coverage_plot():
    # plot absolute ECM seq coverage
    import plotly.express as px

    df = pd.read_csv('F:/fred_time_lapse/20230508/analysis/ECM_aggre_cov_0515.tsv',sep='\t', index_col=0)
    sub_df = df.iloc[:,-7:-2]
    cols = sub_df.columns.tolist()
    cov = []
    for each in cols:
        cov += sub_df[each].tolist()
    time = []
    for i in ['15min', '30min', '1hour', '2hour', '4hour']:
        time += [i]*sub_df.shape[0]
    df_plot = pd.DataFrame(dict(coverage=cov, time=time))
    fig = px.violin(df_plot,y='coverage',color='time',points = 'all', box = True)
    fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)'})
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_traces(marker={'size': 4})
    fig.write_image('F:/fred_time_lapse/20230508/analysis/aggre_cov_145.png',scale=2)


def myplot(x, y, s, bins=1000):

    # scatter density plot
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    heatmap = gaussian_filter(heatmap, sigma=s)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap.T, extent


def coverage_plot():
    # scatter plot showing aggregated vs standard coverage
    ecm_protein_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')
    ecm_prot_list = [k for k in ecm_protein_dict.keys()]
    cov_df = pd.read_csv('F:/fred_time_lapse/analysis/gene_aggre_cov_1219.tsv',sep='\t',index_col=0)
    # ecm_cov_df = cov_df[cov_df.index.isin(ecm_prot_list)]
    # ecm_cov_df.to_csv('F:/fred_time_lapse/analysis/ecm_aggre_cov.tsv',sep='\t')
    # df_fillna = cov_df.fillna(0)  # all proteins
    df_fillna = cov_df.fillna(0) # only ECM proteins from mat.fasta

    x = np.log2(df_fillna['145_240_aggre_cov']+1)
    y = np.log2(df_fillna['145_1080_cov']+1)
    fig, ax = plt.subplots()
    img, extent = myplot(x, y, 16)
    ax.plot(x, y, 'k.', markersize=5)
    ax.imshow(img, extent=extent, origin='lower', cmap=cm.binary)
    ax.set_xlabel('log2 (4h aggre. cov%)', fontsize=8)
    ax.set_ylabel('log2 (18h std. cov%)', fontsize=8)
    # ax.set_xticks(range(7))
    # ax.set_xticklabels([str(i) for i in range(7)])
    # ax.set_yticks(range(7))
    # ax.set_yticklabels([str(i) for i in range(7)])
    # plt.xlim(0, 6.5)
    # plt.ylim(0, 6.5)
    ax.set_aspect('equal', adjustable='box')
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    plt.savefig('F:/fred_time_lapse/figures/ecm_cov_145_240agg_vs_1080_1222.png', dpi=300)

    plt.show()


def scatter_plot_cov():
    # plot scatter plot to compare aggre.cov between sned1 and control
    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/analysis/gene_aggre_cov_0107.tsv', sep='\t',
                             index_col=0).fillna(0)
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    sned1_agg, control_agg, category = ecm_cov_df['144_1080_cov'], ecm_cov_df['145_1080_cov'], ecm_cov_df['Sub']
    for x, y, cat in zip(control_agg, sned1_agg, category):
        ax.plot(x,y,'o',markersize=12,markeredgewidth=2,fillstyle='none',color=ecm_class_color_dict[cat])
    ax.axline((1, 1), slope=1, ls='--', c='.3')
    ax.tick_params(axis='both', labelsize=14)
    plt.savefig('F:/fred_time_lapse/figures/ecm_stdcov_sned1y_controlx_0329.png', dpi=300)

    # plt.show()


def scatter_plot_cov2():
    # plot scatter plot to show ECM sequence coverage increase ratio from time lapsed digestion
    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/20230508/analysis/ECM_aggre_cov_0515.tsv', sep='\t',
                             index_col=0)
    ecm_cov_df['ratio'] = (ecm_cov_df['145_240_aggre_cov']-ecm_cov_df['145_1080_cov'])/ecm_cov_df['145_1080_cov']*100
    ecm_cov_df = ecm_cov_df.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    # ecm_cov_df.to_csv('F:/fred_time_lapse/20230508/analysis/ECM_aggre_ratio_0521.tsv', sep='\t')
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    length, agg, ratio, category = ecm_cov_df['length'], ecm_cov_df['145_240_aggre_cov'], ecm_cov_df['ratio'], ecm_cov_df['Sub']
    for x, y, size, cat in zip(length,agg,ratio,category):
        ax.plot(x,y,'o',markersize=size/5,markeredgewidth=2,fillstyle='none',color=ecm_class_color_dict[cat])
    # plt.savefig('F:/fred_time_lapse/20230508/analysis/ECM_coverage_increase_ratio_0521.png', dpi=300)
    # plt.show()

def category_cov_plot():
    # plot sequence coverage for different categories

    # gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json','r'))
    # ecm_cov_df = pd.read_csv('F:/fred_time_lapse/analysis/ecm_aggre_cov.tsv',sep='\t',index_col=0)
    # gene_list = ecm_cov_df['gene'].tolist()
    # ecm_catgory = [gene_category_dict[gene]["Category"] for gene in gene_list]
    # ecm_sub_cat = [gene_category_dict[gene]["Sub"] for gene in gene_list]
    # ecm_cov_df['category'] = ecm_catgory
    # ecm_cov_df['sub_category'] = ecm_sub_cat
    # ecm_cov_df.to_csv('F:/fred_time_lapse/analysis/ecm_aggre_cov_1219.tsv',sep='\t')

    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/20230508/analysis/ECM_aggre_cov_0515.tsv',sep='\t').fillna(0)
    # violin plot for each category
    fig,axs = plt.subplots(2,3,figsize=(10,5))
    for each, ax in zip(sort_category, [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2]]):

        color = ecm_class_color_dict[each]
        sub_df = ecm_cov_df[ecm_cov_df['Sub'] == each]
        aggre_cov, std_cov = sub_df['144_240_aggre_cov'].tolist(), sub_df['144_1080_cov'].tolist()
        # hue condition
        condition = []
        for agg, std in zip(aggre_cov, std_cov):
            if agg>0 and std==0:
                condition.append('only in aggre.')
            elif std>0 and agg==0:
                condition.append('only in std.')
            else:
                condition.append('both')
        sub_df_plot = pd.DataFrame(
            dict(agg_or_standard=['144_4h_agg'] * sub_df.shape[0] + ['144_18h'] * sub_df.shape[0],
                 coverage=aggre_cov + std_cov,condition=condition*2))
        median_list = sub_df_plot.groupby(['agg_or_standard'])['coverage'].median().tolist()
        x_pos_list = range(len(median_list))
        n_list = ['n = %i' % (sub_df.shape[0]) for i in range(len(median_list))]
        g = sns.swarmplot(data=sub_df_plot, x='agg_or_standard', y='coverage', ax=axs[ax[0], ax[1]],
                          color=color,hue='condition',hue_order=['only in std.','only in aggre.','both'])
        g.legend_.remove()
        # label sample size
        for i in range(len(median_list)):
            axs[ax[0], ax[1]].text(x_pos_list[i] + 0.05, median_list[i], n_list[i])
        add_stat_annotation(ax=axs[ax[0],ax[1]], data=sub_df_plot, x='agg_or_standard', y='coverage',
                            box_pairs=[("144_4h_agg", "144_18h")],
                            test='Wilcoxon', text_format='star',loc='inside', verbose=2,comparisons_correction='bonferroni')
        axs[ax[0], ax[1]].set_xlabel('')
        axs[ax[0], ax[1]].set_ylabel('')
    plt.savefig('F:/fred_time_lapse/20230508/analysis/ecm_cov_aggre_vs_18h_0515.png',dpi=300)
    # plt.show()


def dot_plot_connecting():
    # plot dot plots and connecting the dots for sequence coverage
    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/20230508/analysis/ECM_aggre_cov_0515.tsv', sep='\t',index_col=0).fillna(0)
    fig, axs = plt.subplots(2, 3, figsize=(8, 5))
    for each, ax in zip(sort_category, [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2]]):

        color = ecm_class_color_dict[each]
        sub_df = ecm_cov_df[ecm_cov_df['Sub'] == each]
        aggre_cov_list, std_cov_list = sub_df['145_240_aggre_cov'].tolist(), sub_df['145_1080_cov'].tolist()
        sub_df_plot = pd.DataFrame(
            dict(agg_or_standard=['145_4h_agg'] * sub_df.shape[0] + ['145_18h'] * sub_df.shape[0],
                 coverage=aggre_cov_list + std_cov_list))
        median_list = sub_df_plot.groupby(['agg_or_standard'])['coverage'].median().tolist()
        x_pos_list = range(len(median_list))
        n_list = ['n = %i' % (sub_df.shape[0]) for i in range(len(median_list))]
        x_agg, x_std = np.random.normal(0,0.05,size=len(aggre_cov_list)), np.random.normal(1,0.05,size=len(std_cov_list))

        # boxplot
        axs[ax[0], ax[1]].boxplot([aggre_cov_list,std_cov_list],positions=[0,1], widths=(0.5,0.5), showfliers=False)

        # plot dots
        for x, y in zip([x_agg, x_std],[aggre_cov_list,std_cov_list]):
            axs[ax[0], ax[1]].plot(x,y,'o',markersize=5,alpha=0.6,color=color)

        # connecting dots
        for x1,y1,x2,y2 in zip(x_agg,aggre_cov_list,x_std,std_cov_list):
            axs[ax[0], ax[1]].plot([x1,x2],[y1,y2],color=color, linestyle='--', alpha=0.3)

        # label sample size
        # for i in range(len(median_list)):
        #         #     axs[ax[0], ax[1]].text(x_pos_list[i] + 0.05, median_list[i], n_list[i])
        axs[ax[0], ax[1]].text(0.01,0.99,n_list[0], ha='left',va='top',transform=axs[ax[0], ax[1]].transAxes)

        # label statistics
        add_stat_annotation(ax=axs[ax[0], ax[1]], data=sub_df_plot, x='agg_or_standard', y='coverage',
                            box_pairs=[("145_4h_agg", "145_18h")],
                            test='Wilcoxon', text_format='star', loc='inside', verbose=2,
                            comparisons_correction='bonferroni')
        axs[ax[0], ax[1]].set_xlim([-0.5, 1.5])
        axs[ax[0], ax[1]].set_xticks([0,1])
        axs[ax[0], ax[1]].set_xticklabels(['4h_agg','18h_std'], fontsize=10)
        axs[ax[0], ax[1]].set_xlabel('')
        axs[ax[0], ax[1]].set_ylabel('')
    plt.savefig('F:/fred_time_lapse/20230508/analysis/145_ecm_aggre_cov_vs_18h_0521.png', dpi=300)
    # plt.show()


def dots_plotting_all():
    # plot sequence coverage data for all ECMs
    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/analysis/gene_aggre_pep_5compare5_best_psm_0328.tsv', sep='\t',
                             index_col=0).fillna(0)
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    aggre_cov_list, std_cov_list = ecm_cov_df['144_240_aggre_cov'].tolist(), ecm_cov_df['144_1080_cov'].tolist()
    sub_df_plot = pd.DataFrame(
        dict(agg_or_standard=['144_4h_agg'] * ecm_cov_df.shape[0] + ['144_18h'] * ecm_cov_df.shape[0],
             coverage=aggre_cov_list + std_cov_list))
    x_agg, x_std = np.random.normal(0, 0.05, size=len(aggre_cov_list)), np.random.normal(1, 0.05,
                                                                                           size=len(std_cov_list))
    # boxplot
    ax.boxplot([aggre_cov_list, std_cov_list], positions=[0, 1], widths=(0.3, 0.3), showfliers=False)
    # plot dots
    for x, y in zip([x_agg, x_std], [aggre_cov_list, std_cov_list]):
        ax.plot(x, y, 'o', markersize=5, alpha=0.6, color='k')
    # connecting dots
    for x1, y1, x2, y2 in zip(x_agg, aggre_cov_list, x_std, std_cov_list):
        ax.plot([x1, x2], [y1, y2], color='k', linestyle='--', alpha=0.3)

    ax.text(0.01, 0.99, 'n = %i' % ecm_cov_df.shape[0], ha='left', va='top', transform=ax.transAxes)
    add_stat_annotation(ax=ax, data=sub_df_plot, x='agg_or_standard', y='coverage',
                        box_pairs=[("144_4h_agg", "144_18h")],
                        test='Wilcoxon', text_format='star', loc='inside', verbose=2,
                        comparisons_correction='bonferroni')
    ax.set_xlim([-0.5, 1.5])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['4h_agg', '18h_std'], fontsize=10)
    ax.set_xlabel('')
    ax.set_ylabel('')
    plt.savefig('F:/fred_time_lapse/figures/ecm_cov_144_240agg_vs_1080_all_0328_5vs5_best.png', dpi=300)
    # plt.show()


def qc_check_ecm_ratio():
    # check total intensity of ECM entries/all entries
    gene_int_files = glob('F:/fred_time_lapse/search*/report.gg_matrix*.tsv')
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json','r'))
    for each_f in gene_int_files:
        gg_df = pd.read_csv(each_f,sep='\t')
        gg_df_ecm = gg_df[gg_df['Genes'].isin(gene_category_dict)]
        columns = gg_df.columns[1:]
        for col in columns:
            totol_int = gg_df[col].sum()
            ecm_total_int = gg_df_ecm[col].sum()
            print (col, ecm_total_int/totol_int)


def table_output():
    """
    output tables with peptide IDs, intensity, coverage
    :return:
    """
    from collections import Counter
    base_path = 'F:/fred_time_lapse/'
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r')) # ECM genes
    ecm_genelist = [k for k in gene_category_dict.keys()]
    gene_seq_dict = fasta_reader_gene('F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')

    gene_f_psm_dict_of_dict = pk.load(open('F:/fred_time_lapse/20230503/analysis/gene_f_psm_dict_of_dict_0524.p', 'rb'))
    intensity_dict = pk.load(open(base_path+'20230503/analysis/gene_f_intesity_dict_of_dict_0524.p','rb'))
    df = pd.DataFrame()

    times, samples, replicates = ['15', '30', '60', '120', '240'], ['144','145'],['A','B','C','D','E']
    std_replicates = ['NA_1080','NB_1080','NC_1080','ND_1080','NE_1080']
    # samples, std_replicates = ['145_nA1080_'], ['1_6', '2_7', '3_8', '4_9', '5_10']

    # for gene in gene_category_dict:
    #     if gene in intensity_dict:
    for gene in gene_f_psm_dict_of_dict:

            print(gene)
            seq = gene_seq_dict[gene]
            # df.loc[gene, 'category'] = gene_category_dict[gene]["Category"]
            # df.loc[gene, 'Sub'] = gene_category_dict[gene]["Sub"]

            for sample in samples:
                for time in times:
                    for rep in replicates:
                #         # sequence coverage
                #         np_array = np.zeros(len(seq))
                        psm_list = gene_f_psm_dict_of_dict[gene][sample + rep + '_' + time]
                #         psm_count = Counter(psm_list)
                        if len(set(psm_list)) > 1:  # filter out genes with only one peptide
                #             # for each in psm_list:
                #             #     pep_loc = seq.find(each)
                #             #     pep_end_loc = pep_loc + len(each)
                #             #     np_array[pep_loc:pep_end_loc] += 1
                #             # df.loc[gene,sample+rep+'_'+time] = np.count_nonzero(np_array)/len(seq)*100
                            df.loc[gene, sample + rep + '_' + time] = intensity_dict[gene][sample + rep + '_' + time]
                #             df.loc[gene,sample+rep+'_'+time] = str([(pep, 'psm_'+str(psm_count[pep]),'start_'+str(seq.find(pep)+1),'end_'+str(seq.find(pep)+len(pep))) for pep in psm_count])
                #             # df.loc[gene, sample + rep + '_' + time] = len(psm_list)
                        else:
                            df.loc[gene, sample + rep + '_' + time] = 0
                for std in std_replicates:
                    np_array_std = np.zeros(len(seq))
                    psm_list = gene_f_psm_dict_of_dict[gene][sample + std]
                    # psm_count = Counter(psm_list)
                    if len(set(psm_list)) > 1:  # filter out genes with only one peptide
                        # for each in psm_list:
                        #     pep_loc = seq.find(each)
                        #     pep_end_loc = pep_loc + len(each)
                        #     np_array_std[pep_loc:pep_end_loc] += 1
                        # df.loc[gene,sample+std] = np.count_nonzero(np_array_std)/len(seq)*100
                        df.loc[gene, sample + std] = intensity_dict[gene][sample + std]
                        # df.loc[gene,sample+std] = str([(pep, 'psm_'+str(psm_count[pep]),'start_'+str(seq.find(pep)+1),'end_'+str(seq.find(pep)+len(pep))) for pep in psm_count])
                    # df.loc[gene,sample+std] = len(psm_list)
                    else:
                        df.loc[gene, sample + std] = 0

    df.to_csv(base_path+'20230503/analysis/All_gene_intensity_0524.tsv',sep='\t')


def compare5_to5():

    # compare 5 standard 18 hour files versus 5 time lapse file (each from 1 time point)
    from functools import reduce
    import operator
    base_path = 'F:/fred_time_lapse/'
    gene_f_psm_dict_of_dict = pk.load(open(base_path+'analysis/gene_f_psm_dict_of_dict_1219.p', 'rb'))
    intensity_dict = pk.load(open(base_path + 'analysis/gene_f_intesity_dict_of_dict_1219.p', 'rb'))
    times, samples, replicates = ['15', '30', '60', '120', '240'], ['144', '145'], ['A', 'B', 'C', 'D', 'E']

    file_psm_n_dict = defaultdict(int)  # calculate number total psm for each file
    for gene in gene_f_psm_dict_of_dict:
        if gene in intensity_dict: # only count ECM gene
            for file in gene_f_psm_dict_of_dict[gene]:
                psm_n = len(gene_f_psm_dict_of_dict[gene][file])
                file_psm_n_dict[file] += psm_n
    file_group_psm_n_dict = defaultdict(list) # list for each value contains psm number for each replicate
    for file in file_psm_n_dict:
        sample, time, rep = file.split('_')[0][:-1], file.split('_')[-1], file.split('_')[0][-1]
        file_group_psm_n_dict[sample+'_'+time].append((file_psm_n_dict[file], rep))
    psm_median_group = {}
    for each in file_group_psm_n_dict:
        # print (each, sorted(file_group_psm_n_dict[each]))
        psm_median_group[each]=each.split('_')[0]+sorted(file_group_psm_n_dict[each])[-1][-1]+'_'+each.split('_')[1]
    print (psm_median_group)

    # aggregate psms throughout time points
    time_points = ['15', '30', '60', '120', '240']
    aggregate_pep_dict = {}
    for prot in intensity_dict:
        aggregate_psm_dict = {}

        for ind, val in enumerate(time_points):
            time_points_slice = time_points[:ind + 1]
            for sample in ['144', '145']:
                # aggregate psms from 15min to 240min
                aggregate_list = []
                for t in time_points_slice:
                    if len(set(gene_f_psm_dict_of_dict[prot][psm_median_group[sample+'_'+t]]))>1:
                        aggregate_list += gene_f_psm_dict_of_dict[prot][psm_median_group[sample+'_'+t]]
                # aggregate_psm_set = reduce(operator.or_, [set(gene_f_psm_dict_of_dict[prot][psm_median_group[sample+'_'+t]]) for t in
                #                                               time_points_slice if len(set(gene_f_psm_dict_of_dict[prot][psm_median_group[sample+'_'+t]]))>1])
                aggregate_psm_dict[sample+'_'+val+'_aggregate'] = set(aggregate_list)
        aggregate_pep_dict[prot] = aggregate_psm_dict

    pk.dump(aggregate_pep_dict, open('F:/fred_time_lapse/analysis/ecm_gene_f_aggregate_peptide_dict_best_rep_0327.p', 'wb'))


def filter_df():
    # filter out all zeros in df
    import time
    time.sleep(3)
    base_path = 'F:/fred_time_lapse/'
    df = pd.read_csv(base_path+'20230508/analysis/ECM_gene_time_series_absoluteCov_0521.tsv',sep='\t',index_col=0)
    df = df.copy().fillna(0)
    data = []
    index = []
    for row in df.itertuples():
        # filter rows with all 0s
        # if all([row[i]==0 for i in range(-10,0)]):
        if all([row[i] == 0 for i in range(4, len(row))]):
            continue
        else:
            data.append([i for i in row][1:])
            index.append(row[0])
    new_df = pd.DataFrame(data, columns=df.columns,index=index)
    new_df.to_csv(base_path+'20230508/analysis/ECM_gene_time_series_absoluteCov_0521.tsv',sep='\t')


def nsaf_cal():
    gene_seq_dict = fasta_reader_gene('F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')
    df = pd.read_csv('F:/fred_time_lapse/analysis/All_gene_timelapsed_psm_count_0109_atleast2pep.tsv',sep='\t',index_col=0)
    columns = df.columns
    gene_list = df.index
    new_df = pd.DataFrame(index=df.index)
    # new_df['category'] = df['category']
    # new_df['Sub'] = df['Sub']

    for col in columns:
        nsaf_total = sum([psm_count/len(gene_seq_dict[gene])for gene, psm_count in zip(gene_list,df[col])])
        new_df[col] = [psm_count/len(gene_seq_dict[gene])/nsaf_total for gene, psm_count in zip(gene_list,df[col])]
    new_df.to_csv('F:/fred_time_lapse/analysis/All_gene_timelapsed_nsaf_0109_atleast2pep.tsv',sep='\t')


def cosine_sim_calculating(v1, v2):
    """
    calculate the cosine similarity beweeen two b/y ions binned vectors
    :param v1:
    :param v2:
    :return:
    """
    from scipy import spatial
    return 1-spatial.distance.cosine(v1,v2)


def dissim_index():
    """
    compare dissimlarity between SNED1+/-
    :return:
    """
    cov_df = pd.read_csv('F:/fred_time_lapse/analysis/gene_aggre_cov_0107.tsv',sep='\t',index_col=0)
    cov_df_copy = cov_df.copy()
    gene_list, cov_matrix = cov_df.index, cov_df.to_numpy()[:,3:13]
    for gene, cov_list in zip(gene_list,cov_matrix):
        # print (gene, cov_list[:5], cov_list[5:])
        euclidean_dist = np.linalg.norm(cov_list[:5]-cov_list[5:])  # compare aggregated coverage between 144 and 145
        # euclidean distance is positive if 144 is bigger
        euclidean_dist_value = euclidean_dist if cov_list[4]-cov_list[9]>=0 else -euclidean_dist
        cos_sim = cosine_sim_calculating(cov_list[:5],cov_list[5:])
        cov_df_copy.at[gene,'euclidean_dis+-/cosine_sim'] = euclidean_dist_value/cos_sim
    cov_df_copy.to_csv('F:/fred_time_lapse/analysis/dissim_index.tsv',sep='\t')


def dissim_index_plot():
    df_ecm_aggre = pd.read_csv('F:/fred_time_lapse/analysis/dissim_index.tsv',sep='\t',index_col=0)
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax = sns.boxplot(x='Sub', y='euclidean_dis+-/cosine_sim', hue='Sub',palette=ecm_class_color_dict, data=df_ecm_aggre,
                     linewidth=2.5,
                 order=["ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated Proteins",
                                "ECM Regulators","Secreted Factors"],
                 hue_order=["ECM Glycoproteins", "Collagens", "Proteoglycans", "ECM-affiliated Proteins",
                                "ECM Regulators","Secreted Factors"], dodge=False)
    ax.get_legend().remove()
    ax.set_ylabel('Dissimilarity index', fontsize=8)
    plt.xticks(rotation=30, fontsize=8)
    plt.tight_layout()
    plt.savefig('F:/fred_time_lapse/figures/dissim_index.png', dpi=300)
    # plt.show()


def upsetplot_data():
    """
    upset plot on peptides in time points
    :return:
    """
    f_peptide_dict = defaultdict(set)
    total_peptide_set = set()
    peptide_gene_category_dict = {} # {peptide:(gene,category)}
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r'))  # ECM genes
    gene_f_peptides = pk.load(open('F:/fred_time_lapse/analysis/gene_f_rep_combined_peptide_dict_0107.p','rb'))
    f_list = ['144_15','144_30','144_60','144_120','144_240']
    df = pd.DataFrame(columns=f_list+['gene','category'])

    for gene in gene_f_peptides:
        if gene in gene_category_dict:
            for f in f_list:
                peptide_set = gene_f_peptides[gene][f]
                f_peptide_dict[f].update(peptide_set)
                total_peptide_set.update(peptide_set)
                for each_pep in peptide_set:
                    peptide_gene_category_dict[each_pep] = (gene,gene_category_dict[gene]['Sub'])
    # output a df with peptide seq as index, boolean to show peptide existence in time points
    for peptide in total_peptide_set:
        for f in f_list:
            if peptide in f_peptide_dict[f]:
                df.at[peptide,f] = 1
            else:
                df.at[peptide,f] = 0
        df.at[peptide,'gene'] = peptide_gene_category_dict[peptide][0]
        df.at[peptide,'category'] = peptide_gene_category_dict[peptide][1]

    df.to_csv('F:/fred_time_lapse/analysis/upset_plot_peptide.tsv',sep='\t')


def upset_plot():
    # plot upset plot based on data from upsetplot_data()
    from upsetplot import UpSet, plot

    f_list = ['144_15', '144_30', '144_60', '144_120', '144_240']
    df = pd.read_csv('F:/fred_time_lapse/analysis/upset_plot_peptide.tsv',sep='\t')
    df = df.set_index(f_list)
    # print (df.head)
    # fig = plt.figure(figsize=(10, 8))
    upset = UpSet(df,
                  intersection_plot_elements=0,sort_by='cardinality')
    # upset = plot(df,fig=fig,element_size=None)
    upset.add_stacked_bars(by='category',colors=ecm_class_color_dict, elements=18)
    upset.plot()
    plt.legend('',frameon=False)
    plt.savefig('F:/fred_time_lapse/figures/144_peptide_upset.png',dpi=300)
    # plt.show()


def abundance_plot_sned1():
    # plot heatmap showing difference between SNED1+/-
    df_sned1 = pd.read_csv('F:/fred_time_lapse/analysis/ECM_gene_144std_intensity_0104_atleast2pep.tsv',sep='\t',index_col=0)
    df_no_sned1 = pd.read_csv('F:/fred_time_lapse/analysis/ECM_gene_145std_intensity_0104_atleast2pep.tsv',sep='\t',index_col=0)
    prot_144, prot_145 = df_sned1.index.tolist(), df_no_sned1.index.tolist()
    column_144, column_145 = df_sned1.columns.tolist(),df_no_sned1.columns.tolist()
    prot_list = set(prot_144+prot_145)
    df_sum = pd.DataFrame(index=list(prot_list),columns=['144_1080sum','145_1080sum'])
    for each in prot_list:
        # print (df_sned1.loc[each,column_144[2]:])
        sum_144 = df_sned1.loc[each, column_144[2]:].sum() if each in prot_144 else 0
        sum_145 = df_no_sned1.loc[each, column_145[2]:].sum() if each in prot_145 else 0
        df_sum.at[each,'144_1080sum'] = sum_144
        df_sum.at[each,'145_1080sum'] = sum_145
    normalized_df = (df_sum-df_sum.min())/(df_sum.max()-df_sum.min())
    plt_array = (normalized_df['144_1080sum']-normalized_df['145_1080sum'])/normalized_df['144_1080sum']*100
    df_plot = pd.DataFrame(index=normalized_df.index, columns=['difference%'])
    df_plot['difference%'] = plt_array.tolist()
    df_plot = df_plot.sort_values(by='difference%', ascending=False)
    print (df_plot.head)
    fig,ax = plt.subplots(1,1, figsize=(15,4))
    g = sns.lineplot(data=df_plot['difference%'],ax=ax,linewidth=2.5,color='k')
    plt.show()


def get_unique_time_point():
    # from the replicate combined data, get peptides unique for each time point,
    # and label each time point with peptides first IDed at this time point
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r')) # ECM genes
    gene_combined_pep = pk.load(open('F:/fred_time_lapse/analysis/gene_f_rep_combined_peptide_dict_0107.p', 'rb'))
    times, samples = ['15', '30', '60', '120', '240'], ['144','145']
    df = pd.DataFrame()
    for gene in gene_category_dict:
        if gene in gene_combined_pep:
            print (gene)
            df.loc[gene, 'category'] = gene_category_dict[gene]["Category"]
            df.loc[gene, 'Sub'] = gene_category_dict[gene]["Sub"]
            for sample in samples:
                total_peptides = set()
                for time in times:
                    # times_copy = ['15', '30', '60', '120', '240']
                    # times_copy.remove(time)
                    # peptide_set = gene_combined_pep[gene][sample+'_'+time]
                    # peptide_other_set = set([pep for time1 in times_copy for pep in gene_combined_pep[gene][sample+'_'+time1]])
                    # # compare the peptide set in one time point with all other time points to get unique ones
                    # unique_pep = peptide_set.difference(peptide_other_set)
                    # # print (time, peptide_set, peptide_other_set, unique_pep)
                    # df.loc[gene,sample+'_'+time+'_unique'] = str(unique_pep)

                    peptide_set = gene_combined_pep[gene][sample+'_'+time]
                    peptide_set_first_time = peptide_set.difference(total_peptides)
                    total_peptides.update(peptide_set_first_time)
                    df.loc[gene,sample+'_'+time] = str(peptide_set_first_time)
    df.to_csv('F:/fred_time_lapse/analysis/ecm_peptide_first_timepoint_0420.tsv',sep='\t')


def aggregate_cov_each_bio_rep():
    # aggregate coverage of all time points for each biological replicate
    base_path = 'F:/fred_time_lapse/'

    gene_f_psm_dict_of_dict = pk.load(open('F:/fred_time_lapse/20230508/analysis/gene_f_psm_dict_of_dict_0509.p', 'rb'))
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r'))  # ECM genes
    intensity_dict = pk.load(open(base_path + '20230508/analysis/gene_f_intesity_dict_of_dict_0509.p', 'rb'))
    gene_seq_dict = fasta_reader_gene('F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')

    df = pd.DataFrame()

    times, samples, replicates = ['15', '30', '60', '120', '240'], ['144', '145'], ['A', 'B', 'C', 'D', 'E']
    # std_replicates = ['NA_1080', 'NB_1080', 'NC_1080', 'ND_1080', 'NE_1080']
    for gene in gene_category_dict:
        if gene in intensity_dict:
    # for gene in intensity_dict:
            print(gene)
            seq = gene_seq_dict[gene]
            df.loc[gene, 'category'] = gene_category_dict[gene]["Category"]
            df.loc[gene, 'Sub'] = gene_category_dict[gene]["Sub"]
            for sample in samples:
                for rep in replicates:
                    np_array = np.zeros(len(seq))
                    # peptide_set = reduce(operator.or_,[set(gene_f_psm_dict_of_dict[gene][sample + rep + '_' + time]) for time in times
                    #                if len(set(gene_f_psm_dict_of_dict[gene][sample + rep + '_' + time]))>=1])
                    peptide_set = set()
                    for time in times:  # which time point to aggregate
                        peptide_time_set = set(gene_f_psm_dict_of_dict[gene][sample + rep + '_' + time])
                        if len(peptide_time_set)>1: # at least 2 peptide ids
                            peptide_set.update(peptide_time_set)
                    for each in peptide_set:
                        pep_loc = seq.find(each)
                        pep_end_loc = pep_loc + len(each)
                        np_array[pep_loc:pep_end_loc] += 1
                    df.loc[gene,sample+rep+'_aggregate_240'] = np.count_nonzero(np_array)/len(np_array)*100
    df.to_csv('F:/fred_time_lapse/20230508/analysis/ECM_gene_aggre240_cov_each_replicate_0517.tsv',sep='\t')


def check_ecm_id():
    from tsv_reader import protein_info_from_fasta, venn_diagram_gen
    import json
    fasta = 'F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta'
    protein_info_dict = protein_info_from_fasta(fasta)
    protein_df = pd.read_csv(r'F:\fred_time_lapse\2023_05_25/report.pr_matrix.tsv',sep='\t').copy()
    genes = []
    for prot in protein_df['Protein.Group']:
        if ';' not in prot:
            try:
                genes.append(protein_info_dict[prot][0])
            except:
                genes.append('')
        else:
            try:
                genes.append(protein_info_dict[prot.split(';')[0]][0])
            except:
                genes.append('')
    protein_df['Genes'] = genes
    protein_df.to_csv(r'F:\fred_time_lapse\2023_05_25/report.pr_matrix_gene.tsv',sep='\t',index=False)

    # df = pd.read_csv(r'F:\fred_time_lapse\20230508/report.pr_matrix_gene.tsv',sep='\t') # merged spectral lib
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r')) # ECM genes
    # merged_gene_peptide_dict = defaultdict(set)
    # for gene, pep, charge in zip(df.Genes,df['Modified.Sequence'],df['Precursor.Charge']):
    #     if gene in gene_category_dict:
    #         # counter += 1
    #         merged_gene_peptide_dict[gene].add(pep+str(charge))
    # pk.dump(merged_gene_peptide_dict,open('F:/fred_time_lapse/20230508/ecm_peptide_0527.p','wb'))
    #
    # file_list = [r'F:\fred_time_lapse/search/report.pr_matrix_1.tsv',r'F:\fred_time_lapse\search_11282022/report.pr_matrix_2.tsv',
    #              r'F:\fred_time_lapse\search_11292022/report.pr_matrix_3.tsv']
    # gene_peptide_dict = defaultdict(set)
    # for f in file_list:
    #     df = pd.read_csv(f,sep='\t')
    #     for gene, pep, charge in zip(df.Genes, df['Modified.Sequence'],df['Precursor.Charge']):
    #         if gene in gene_category_dict:
    #             gene_peptide_dict[gene].add(pep+str(charge))
    #
    # pk.dump(gene_peptide_dict, open('F:/fred_time_lapse/analysis/ecm_peptide_0527.p', 'wb'))

    ## get ecm peptides from raw searched data
    # inter_gene_peptide_dict = defaultdict(set)
    # df = pd.read_csv('F:/fred_time_lapse/20230504/report.pr_matrix_gene.tsv', sep='\t') # dda-based spectral library search
    # for gene, pep, charge in zip(df.Genes, df['Modified.Sequence'], df['Precursor.Charge']):
    #     if gene in gene_category_dict:
    #         inter_gene_peptide_dict[gene].add(pep + str(charge))
    # pk.dump(inter_gene_peptide_dict, open('F:/fred_time_lapse/20230504/ecm_peptide_0527.p', 'wb'))


    #         counter+=1
    # print (counter)
    # print ([gene for gene in ecm_gene_set if gene not in prosit_gene])
    # venn_dict = {'spec_library':ecm_gene_set,'prosit_predict_spec_library':prosit_gene}
    # venn_diagram_gen(venn_dict)
    ### find prosit-unique identifed peptides from prosit library and append to custom library
    # spec_library2 = pd.read_csv(r'F:\fred_time_lapse\spectral_library_MDB/library.tsv', sep='\t')
    # columns = spec_library2.columns
    # unique_pep = [each for each in prosit_total_pepset if each not in ecm_peptide_set]
    # new_spec_lib = []
    # spec_library = pd.read_csv(r'F:\fred_time_lapse\spectral_library/library_combine_20230508.tsv', sep='\t')
    # count = 0
    # for pep in unique_pep:
    #     pep,charge = pep[:-1],int(pep[-1])
    #     sub_df = spec_library[(spec_library['ModifiedPeptideSequence']==pep)&(spec_library['PrecursorCharge']==charge)]
    #     if not sub_df.empty:
    #         count += 1
    # print (count)
    #         for tp in sub_df.itertuples():
    #             new_spec_lib.append([tp.PrecursorMz, tp.ProductMz, tp.FragmentType+'^'+str(tp.FragmentCharge), tp.ProteinGroup, tp.Genes,
    #                    tp.PeptideSequence, tp.ModifiedPeptide, tp.PrecursorCharge, tp.LibraryIntensity*10000, tp.Tr_recalibrated,
    #                    '', tp.FragmentType, tp.FragmentCharge, tp.FragmentSeriesNumber,'',''])
    # new_spec_df = pd.DataFrame(new_spec_lib, columns=columns)
    # new_spec_df.to_csv(r'F:\fred_time_lapse\spectral_library/library_fromprosit_toadd.tsv', sep='\t')

    # total_pep = set(spec_library['PeptideSequence'].tolist())
    # count = 0
    # for pep in unique_pep:
        # print (pep)
        # sub_df = spec_library[spec_library['PeptideSequence']==pep]
        # if not sub_df.empty:
        #     count += 1
        #     print (sub_df)
        # if pep in total_pep:
        #     count += 1
    # print (count)


def combine_spec_library():
    # spec_library1 = pd.read_csv(r'F:\fred_time_lapse\spectral_library/library.tsv', sep='\t').copy()
    spec_library2 = pd.read_csv(r'F:\fred_time_lapse\20230504/library_combine.tsv', sep='\t').copy()
    # check_set = set([mod_seq+str(charge) for mod_seq, charge in
    #                  zip(spec_library2.ModifiedPeptideSequence,spec_library2.PrecursorCharge)])
    # print (spec_library1.columns)
    # matrix = []
    # for tp in spec_library1.itertuples():
    #     mod_seq, charge = tp.ModifiedPeptideSequence, tp.PrecursorCharge
    #     if mod_seq+str(charge) not in check_set: # add rows from spec library1 that are not in spec library 2
    #         matrix.append([each for each in tp][1:])
    #         print (tp[0])
    #
    # spec_library3 = pd.DataFrame(matrix, columns=spec_library1.columns)
    # spec_library3.to_csv(r'F:\fred_time_lapse\spectral_library/library_to_add.tsv',sep='\t')

    spec_library3 = pd.read_csv(r'F:\fred_time_lapse\spectral_library/library_fromprosit_toadd.tsv',sep='\t')
    concat = pd.concat([spec_library2,spec_library3],axis=0)
    print (spec_library2.shape, concat.shape)
    concat.to_csv(r'F:\fred_time_lapse\spectral_library/library_combine_20230508.tsv', sep='\t',index=False)


def compare_results_spectral_library():

    # file_list = [r'F:\fred_time_lapse/search/report.pr_matrix_1.tsv',r'F:\fred_time_lapse\search_11282022/report.pr_matrix_2.tsv',
    #              r'F:\fred_time_lapse\search_11292022/report.pr_matrix_3.tsv']
    # prosit_total_pepset = set()
    # for f in file_list:
    #     df = pd.read_csv(f,sep='\t')
    #     for pep, charge in zip(df['Modified.Sequence'],df['Precursor.Charge']):
    #         prosit_total_pepset.add(pep.replace('(','[').replace(')',']')+str(charge))
    spec_library = pd.read_csv(r'F:\fred_time_lapse\20230504/library_combine.tsv', sep='\t')
    pep_charge_spec = set([pep.replace('(','[').replace(')',']')+str(charge) for pep, charge in zip(spec_library['ModifiedPeptideSequence'], spec_library['PrecursorCharge'])])

    file1 = r'F:\fred_time_lapse\new_spec/spectraRT.tsv'
    df1 = pd.read_csv(file1, sep='\t')
    pep_charge_set = set([pep+str(charge) for pep, charge in zip(df1['peptide'], df1['charge'])])
    print (len([each for each in pep_charge_spec if each in pep_charge_set]))
    venn_dict = {'library combined': pep_charge_spec, 'new_spec': pep_charge_set}
    venn_diagram_gen(venn_dict)


def fraction_time_lapsed_cov():
    """
    from scipy import stats
    # compare fractionation and time lapsed cov, 5 raw versus 5 raw files
    # frac_cov_f = r'F:\fred_time_lapse\20230517\analysis/ecm_gene_fraction_combined_cov_0517.tsv'
    frac_cov_f_2 = r'F:\fred_time_lapse\2023_05_25\analysis/frac_cov_A_to_E.tsv'
    time_lapsed_f = r'F:\fred_time_lapse\20230508\analysis/ECM_gene_aggre240_cov_each_replicate_0517.tsv'
    df1 = pd.read_csv(frac_cov_f_2, sep='\t', index_col=0)
    df2 = pd.read_csv(time_lapsed_f, sep='\t', index_col=0)
    shared_gene = [each for each in df1.index if each in df2.index]
    print (len(shared_gene))
    df = pd.DataFrame(index=shared_gene, columns= ['category','Sub']+df2.columns.tolist()[-5:]+df1.columns.tolist()[-5:])
    # df = pd.DataFrame(index=shared_gene, columns= ['category','Sub']+['A_frac_combined_cov']+df2.columns.tolist()[-4:])

    for gene in shared_gene:
        df.loc[gene, 'category'] = df1.loc[gene,'category']
        df.loc[gene, 'Sub'] = df1.loc[gene,'Sub']
        df.loc[gene, '145A_aggregate_240':'145E_aggregate_240'] = df2.loc[gene, '145A_aggregate_240':'145E_aggregate_240']
        df.loc[gene, 'A_frac_combined_cov':'E_frac_combined_cov'] = df1.loc[gene, 'A_frac_combined_cov':'E_frac_combined_cov']
        # df.loc[gene, 'A_frac_combined_cov'] = df1.loc[gene, '145_frac_1080_cov']
    # df.to_csv('F:/fred_time_lapse/2023_05_25/analysis/aggre_and_frac_cov_A_to_E.tsv',sep='\t')

    # transform dataframe
    df = df.fillna(0)
    new_df = pd.DataFrame(columns=['cov', 'replicate','method'])
    cov, replicate, method = [], [], []
    for col, label in zip(df.columns.to_list()[-10:], [bio + '_aggre' for bio in ['A', 'B', 'C', 'D', 'E']]+
                                                      [bio + '_frac' for bio in ['A', 'B', 'C', 'D', 'E']]):
        cov += df[col].tolist()
        replicate += [label[0]]*df.shape[0]
        method += [label[2:]]*df.shape[0]
    new_df['cov'] = cov
    new_df['replicate'] = replicate
    new_df['method'] = method
    new_df.to_csv('F:/fred_time_lapse/2023_05_25/analysis/aggre_and_frac_cov_plot_df.tsv',sep='\t')
    """
    df = pd.read_csv('F:/fred_time_lapse/2023_05_25/analysis/aggre_and_frac_cov_plot_df.tsv', sep='\t', index_col=0)
    fig, ax = plt.subplots(figsize=(8,6))
    sns.boxplot(ax = ax,data = df, x='replicate', y='cov', hue='method')
    box_pairs = [((bio,'aggre'), (bio,'frac')) for bio in ['A','B','C','D','E']]
    add_stat_annotation(ax, data=df,x='replicate',y='cov',hue = 'method', box_pairs=box_pairs,
                    test='Wilcoxon', text_format='star',loc='outside', verbose=2,comparisons_correction='bonferroni')
    plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
    fig.tight_layout()
    plt.savefig('F:/fred_time_lapse/2023_05_25/analysis/ecm_aggre_frac_cov_compare_boxplot_0601.png',dpi=300)
    # plt.show()
    # df = pd.read_csv('F:/fred_time_lapse/20230517/analysis/time_lapsed_vs_frac_cov_0517.tsv',sep='\t',index_col=0)
    # df = df.fillna(0)
    # for col in df.columns.tolist()[2:-1]:
    #     print (col, stats.ttest_rel(df[col],df['145_frac_1080_cov']))
    # fig,ax = plt.subplots()
    # ax = df.boxplot(column = df.columns.tolist()[-6:],
    #              boxprops=dict(linestyle='-', linewidth=2, color='k'),
    #              flierprops=dict(linestyle='-', linewidth=2, color='k'),
    #              medianprops=dict(linestyle='-', linewidth=2, color='r'),
    #              whiskerprops=dict(linestyle='-', linewidth=2, color='k'),
    #              capprops=dict(linestyle='-', linewidth=2),
    #              showfliers=True, grid=False, rot=20)
    # fig.tight_layout()
    # plt.savefig('F:/fred_time_lapse/20230517/analysis/time_lapsed_vs_frac_cov_0517.png', dpi=300)


def compare_results_spectral_lib():
    # compare the results from prosit, dda-based, and merged spectral lib
    merged_gene_peptide_dict = pk.load(open('F:/fred_time_lapse/20230508/ecm_peptide_0527.p', 'rb'))
    gene_peptide_dict = pk.load(open('F:/fred_time_lapse/analysis/ecm_peptide_0527.p', 'rb'))
    inter_gene_peptide_dict = pk.load(open('F:/fred_time_lapse/20230503/analysis/ecm_peptide_0527.p', 'rb'))
    inter2_gene_peptide_dict = pk.load(open('F:/fred_time_lapse/20230504/ecm_peptide_0527.p','rb'))

    venn_dict = {}
    for each, name in zip([gene_peptide_dict,inter_gene_peptide_dict,inter2_gene_peptide_dict, merged_gene_peptide_dict],
                          ['prosit', 'DDA-based', 'DDA+MDB', 'Merged']):
        venn_dict[name] = [key for key in each.keys()]
    # print (venn_dict)
    for each in venn_dict:
        print (each, len(venn_dict[each]))
    venn_diagram_gen2(venn_dict, title='ECM proteins')


if __name__ == '__main__':
    from protein_coverage import fasta_reader_gene
    # qc_check()
    # cv_box_plot()
    # pr_matrix_reader()
    # aggregate_psms()
    # coverage_calculation()
    # coverage_plot()
    # category_cov_plot()
    # dot_plot_connecting()
    # dots_plotting_all()
    # scatter_plot_cov()
    # scatter_plot_cov2()
    # qc_check_ecm_ratio()
    # psm_dict = pk.load(open('F:/fred_time_lapse/analysis/prot_f_rep_combined_peptide_dict_1219.p', 'rb'))
    # print (psm_dict['Q8TER0']['145_1080'])
    # table_output()
    # filter_df()
    # nsaf_cal()
    # abs_coverage_calculation()
    # abs_coverage_plot()
    # dissim_index()
    # dissim_index_plot()
    # upsetplot_data()
    # upset_plot()
    # abundance_plot_sned1()
    # compare5_to5()
    # get_unique_time_point()
    # aggregate_cov_each_bio_rep()
    # check_ecm_id()
    # combine_spec_library()
    # compare_results_spectral_library()
    # compare_results_spectral_lib()
    # combine_psms_fraction()
    fraction_time_lapsed_cov()