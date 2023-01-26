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
from protein_coverage import fasta_reader_gene

ecm_class_color_dict = {"Collagens": '#0584B7', 'ECM-affiliated Proteins':'#F4511E',
                        'ECM Regulators':"#F9A287","Secreted Factors":"#FFE188",
                        "ECM Glycoproteins":"#133463", "Proteoglycans":"#59D8E6"}

sort_category = ["ECM Glycoproteins","Collagens","Proteoglycans","ECM-affiliated Proteins","ECM Regulators",
                  "Secreted Factors"]

def qc_check():
    """
    check coefficient variation of gene-centric intensities between replicates, use output from DIANN report.gg_matrix
    :return:
    """
    gene_int_dict = defaultdict(list)
    base_path = 'F:/fred_time_lapse/'
    intensity_files = glob(base_path+'search*/report.gg_matrix*.tsv')

    for each_f in intensity_files:
        print (each_f)
        df = pd.read_csv(each_f, sep='\t').fillna(0)
        columns = df.columns.tolist()[1:]
        for each_col in columns:
            for gene, intensity in zip(df['Genes'], df[each_col]):
                gene_list = [gene] if ';' not in gene else gene.split(';')
                for g in gene_list:
                    gene_int_dict[g].append((each_col.split('\\')[-1].split('.mzML')[0],intensity))
    gene_f_int_dict = {}
    for gene in gene_int_dict:
        f_int_dict = defaultdict(float)
        for tp in gene_int_dict[gene]:
            f_int_dict[tp[0]]=tp[1]
        gene_f_int_dict[gene] = f_int_dict
    pk.dump(gene_f_int_dict,open(base_path+'analysis/gene_f_intesity_dict_of_dict_1219.p','wb'))

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
    psm_f = glob(base_path+'search*/report.pr_matrix*.tsv')
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
    pk.dump(prot_info_dict,open(base_path+'analysis/gene_file_psm_tuple_dict_1219.p','wb'))


def aggregate_psms():
    """
    aggregate PSMs between different time
    :return:
    """
    from functools import reduce
    import operator
    prot_f_psm_dict_of_dict = {}
    replicates_combined_dict = {}  # {'uniprotid':{'144_15':{combined_psm_set},'144_30':{combined_psm_set}}}
    aggregate_pep_dict = {}  # {'uniprotid':{'144_15':{aggregated_psm_set},'144_30':{aggregate_psm_set}}}
    """
    # {'uniprotid':[(file,psm),(file,psm)]}
    prot_f_psm_tuple_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_file_psm_tuple_dict_1219.p','rb'))

    for prot in prot_f_psm_tuple_dict:
        f_psm_dict = defaultdict(list)
        for tup in prot_f_psm_tuple_dict[prot]:
            f_psm_dict[tup[0]].append(tup[1])  
        prot_f_psm_dict_of_dict[prot] = f_psm_dict
    # print (prot_f_psm_dict_of_dict['Q3TW96'])
    pk.dump(prot_f_psm_dict_of_dict,open('F:/fred_time_lapse/analysis/gene_f_psm_dict_of_dict_1219.p','wb'))
    """
    # combine all psms among replicates (A to E), union all psms in replicates for now for better coverage
    prot_f_psm_dict_of_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_f_psm_dict_of_dict_1219.p','rb'))
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
    pk.dump(replicates_combined_dict, open('F:/fred_time_lapse/analysis/gene_f_rep_combined_peptide_dict_0107.p', 'wb'))

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

    pk.dump(aggregate_pep_dict,open('F:/fred_time_lapse/analysis/gene_f_aggregate_peptide_dict_0107.p','wb'))


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
    prot_f_aggregate_pep_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_f_aggregate_peptide_dict_0107.p','rb'))
    prot_f_combined_pep_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_f_rep_combined_peptide_dict_0107.p','rb')) # key is '144_1080' and '145_1080'

    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r'))  # ECM genes
    # calculate coverage
    time_points, samples = ['15', '30', '60', '120', '240'], ['144','145']
    samples_1080 = ['144_1080','145_1080']
    aggre_cov_df = pd.DataFrame()
    for prot in prot_f_aggregate_pep_dict:
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
            for sample in samples:
                for time in time_points:
                    np_arry = np.zeros(prot_length)
                    peptides_set = prot_f_aggregate_pep_dict[prot][sample+'_'+time+'_aggregate']
                    # for pep in peptides_set:
                    #     pep_loc = seq.find(pep)
                    #     pep_end_loc = pep_loc+len(pep)
                    #     np_arry[pep_loc:pep_end_loc]+=1
                    # agg_cov = np.count_nonzero(np_arry)/prot_length*100
                    # aggre_cov_df.at[prot,sample+'_'+time+'_aggre_cov'] = agg_cov
                    aggre_cov_df.at[prot,sample+'_'+time+'_aggre_pep'] = str(peptides_set)
            for sample2 in samples_1080:
                np_arry = np.zeros(prot_length)
                peptides1080_set = prot_f_combined_pep_dict[prot][sample2]
                # for pep in peptides1080_set:
                #     pep_loc = seq.find(pep)
                #     pep_end_loc = pep_loc + len(pep)
                #     np_arry[pep_loc:pep_end_loc] += 1
                #     cov1080 = np.count_nonzero(np_arry)/prot_length*100
                #     aggre_cov_df.at[prot,sample2+'_cov'] = cov1080
                aggre_cov_df.at[prot,sample2+'_pep'] = str(peptides1080_set)
    aggre_cov_df.to_csv('F:/fred_time_lapse/analysis/gene_aggre_pep_0107.tsv',sep='\t')


def abs_coverage_calculation():
    """
    calculate the absolute sequence coverage for time series data without aggregation
    :return:
    """
    protein_seq_dict = fasta_reader_gene(
        'F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r'))  # ECM genes
    gene_rep_combined_psm_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_f_rep_combined_peptide_dict_0107.p', 'rb'))
    samples, times = ['144','145'], ['15', '30', '60', '120', '240']
    cov_df = pd.DataFrame()
    for prot in gene_rep_combined_psm_dict:
        # if prot in gene_category_dict:
        if prot != 'null':
            print(prot)
            seq = protein_seq_dict[prot]
            prot_length = len(seq)
            # cov_df.loc[prot, 'category'] = gene_category_dict[prot]["Category"]
            # cov_df.loc[prot, 'Sub'] = gene_category_dict[prot]["Sub"]
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
    cov_df.to_csv('F:/fred_time_lapse/analysis/All_gene_time_series_absoluteCov_0125.tsv',sep='\t')


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

    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/analysis/gene_aggre_cov_1219.tsv',sep='\t').fillna(0)
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
    plt.savefig('F:/fred_time_lapse/figures/ecm_cov_144_240agg_vs_1080_category_1223_swarm.png',dpi=300)
    # plt.show()


def dot_plot_connecting():
    # plot dot plots and connecting the dots for sequence coverage
    ecm_cov_df = pd.read_csv('F:/fred_time_lapse/analysis/gene_aggre_cov_0107.tsv', sep='\t',index_col=0).fillna(0)
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
    plt.savefig('F:/fred_time_lapse/figures/ecm_cov_145_240agg_vs_1080_category_0107_dotsconnect.png', dpi=300)
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

    base_path = 'F:/fred_time_lapse/'
    gene_category_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json', 'r')) # ECM genes
    ecm_genelist = [k for k in gene_category_dict.keys()]
    gene_seq_dict = fasta_reader_gene('F:/sned1_biotinalytion/uniprot-proteome_UP000000589_mouse_human_SNED1_BSA.fasta')

    gene_f_psm_dict_of_dict = pk.load(open('F:/fred_time_lapse/analysis/gene_f_psm_dict_of_dict_1219.p', 'rb'))
    intensity_dict = pk.load(open(base_path+'analysis/gene_f_intesity_dict_of_dict_1219.p','rb'))
    df = pd.DataFrame()

    times, samples, replicates = ['15', '30', '60', '120', '240'], ['144','145'],['A','B','C','D','E']
    std_replicates = ['NA_1080','NB_1080','NC_1080','ND_1080','NE_1080']
    # for gene in gene_category_dict:
    #     if gene in intensity_dict:
    for gene in intensity_dict:
        print(gene)
        seq = gene_seq_dict[gene]
        # df.loc[gene, 'category'] = gene_category_dict[gene]["Category"]
        # df.loc[gene, 'Sub'] = gene_category_dict[gene]["Sub"]

        for sample in samples:
            for time in times:
                for rep in replicates:
                    # sequence coverage
                    np_array = np.zeros(len(seq))
                    psm_list = gene_f_psm_dict_of_dict[gene][sample + rep + '_' + time]
                    if len(set(psm_list)) > 1:  # filter out genes with only one peptide
                        # for each in psm_list:
                        #     pep_loc = seq.find(each)
                        #     pep_end_loc = pep_loc + len(each)
                        #     np_array[pep_loc:pep_end_loc] += 1
                        # df.loc[gene,sample+rep+'_'+time] = np.count_nonzero(np_array)/len(seq)*100
                        # df.loc[gene, sample + rep + '_' + time] = intensity_dict[gene][sample + rep + '_' + time]
                        df.loc[gene,sample+rep+'_'+time] = str(set(psm_list))
                        # df.loc[gene, sample + rep + '_' + time] = len(psm_list)
                    else:
                        df.loc[gene, sample + rep + '_' + time] = 0
            # for std in std_replicates:
            #     np_array_std = np.zeros(len(seq))
            # psm_list = gene_f_psm_dict_of_dict[gene][sample + std]
            # if len(set(psm_list)) > 1:  # filter out genes with only one peptide
            # for each in psm_list:
            #     pep_loc = seq.find(each)
            #     pep_end_loc = pep_loc + len(each)
            #     np_array_std[pep_loc:pep_end_loc] += 1
            # df.loc[gene,sample+std] = np.count_nonzero(np_array_std)/len(seq)*100
            # df.loc[gene, sample + std] = intensity_dict[gene][sample + std]
            # df.loc[gene,sample+std] = str(set(psm_list))
            # df.loc[gene,sample+std] = len(psm_list)
            # else:
            #     df.loc[gene, sample + std] = 0

    df.to_csv(base_path+'analysis/All_gene_timelapsed_peptide_0109_atleast2pep.tsv',sep='\t')


def filter_df():
    # filter out all zeros in df
    import time
    time.sleep(3)
    base_path = 'F:/fred_time_lapse/'
    df = pd.read_csv(base_path+'analysis/All_gene_time_series_absoluteCov_0125.tsv',sep='\t',index_col=0)
    df = df.copy()
    data = []
    index = []
    for row in df.itertuples():
        # filter rows with all 0s
        if all([row[i]==0 for i in range(-10,0)]):
        # if all([row[i] == "0" for i in range(1, len(row))]):
            continue
        else:
            data.append([i for i in row][1:])
            index.append(row[0])
    new_df = pd.DataFrame(data, columns=df.columns,index=index)
    new_df.to_csv(base_path+'analysis/All_gene_time_series_absoluteCov_0125.tsv',sep='\t')


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
    # qc_check_ecm_ratio()
    # psm_dict = pk.load(open('F:/fred_time_lapse/analysis/prot_f_rep_combined_peptide_dict_1219.p', 'rb'))
    # print (psm_dict['Q8TER0']['145_1080'])
    # table_output()
    filter_df()
    # nsaf_cal()
    # abs_coverage_calculation()