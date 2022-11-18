"""
-----
data process and analysis of colon dataset from Cornell University
-----
Compare the MS2 spectra with Prosit predicted MS2 spectra:
1. get a list of PSM sequences of your interest together with charge, e.g. [PSM1, PSM2,...], [PSM1_charge, PSM2_charge,...]
2. output info from 1 into a csv for prosit prediction https://www.proteomicsdb.org/prosit/, go to spectral library tab
3. run prosit server to have predicted spectra (m/z, intensity) for each of the PSM, you will have a myPrositLib.msp file
4. convert raw files to .ms2 format, read the m/z, int pairs of your PSMs of interest into a data structure
   Now you have predicted m/z int pairs and real m/z int pairs of each of the PSMs
5. compare real pairs and predicted pairs for each of PSMs: a. bin m/z values based on theoretical mass,
   b. get sum of intensity of each m/z bin for both predicted/real spectra
   c. cos similarity compare both spectra
"""

from multiprocessing_naive_algorithym import *
from aho_corasick import automaton_trie,automaton_matching
from glob import glob
from protein_coverage import fasta_reader
import pandas as pd


def psmlist_todict(psm_list):
    regex_pat = '\w{1}\[\d+\.?\d+\]'
    psm_dict = defaultdict(list)

    for psm in psm_list:
        reg_sub = re.sub(regex_pat, my_replace, psm)
        psm_dict[reg_sub].append(psm)
    return psm_dict


def my_replace(match_obj):
    match_obj = match_obj.group()
    matched_aa = match_obj[0]
    if matched_aa != 'n':
        return matched_aa  # gives back the first element of matched object as string
    else:
        # if first match is n, then n acetylation, get rid of n
        return ''


def filter_pin(pin_file):
    """
    filter a pin file to find psms mapping to target membrane proteins (_Hagay)
    :param pin_file:
    :return:
    """
    psm_prot_spec = defaultdict(list)
    with open(pin_file,'r') as f_open:
        next(f_open)
        for line in f_open:
            if '_Hagay_' in line and 'rev_' not in line and 'sp|' not in line and 'tr|' not in line:
                line_split = line.split('\t')
                psm = line_split[19][2:-2]
                prot = line_split[20].rstrip('\n')
                spec = '.'.join(line_split[0].split('.')[:2])
                hyper_score = float(line_split[9])
                psm_prot_spec[psm].append((prot,spec,hyper_score))
    return psm_prot_spec


def filter_pin_target(pin_file,human_prot_dict):
    """
    filter a pin file from target database search (only includes membrane protein sequences, not human fasta) to find
    PSMs mapping to target proteins (_Hagay)
    :param pin_file:
    :param human_prot_dict: protein seq dictionary from human fasta file
    :return:
    """
    psm_dict = defaultdict(list)
    with open(pin_file,'r') as f_open:
        next(f_open)
        for line in f_open:
            if '_Hagay_' in line and 'rev_' not in line:
                line_split = line.split('\t')
                psm = line_split[19][2:-2]
                prot = line_split[20].rstrip('\n')
                spec = '.'.join(line_split[0].split('.')[:2])
                hyper_score = float(line_split[9])
                psm_dict[psm].append((prot,spec,hyper_score))
    clean_psm_dict = psmlist_todict([k for k in psm_dict.keys()]) # remove PTMs {'PEPTIDE':['PEP[100]TIDE']}

    # map all clean PSMs to human proteome to see if there's any match, if there is, delete the PSM
    seq_list = extract_UNID_and_seq(human_prot_dict)[1]
    seqline = creat_total_seq_line(seq_list, sep='|')
    aho_result = automaton_matching(automaton_trie([pep for pep in clean_psm_dict.keys()]), seqline)
    mapped_peptides = set([pos[2] for pos in aho_result])
    # print (f'{len(mapped_peptides)} mapped to human proteome')

    # get psm not mapped (unique mapping to membrane proteins)
    membrane_map = [each for each in clean_psm_dict.keys() if each not in mapped_peptides]

    return {psm:psm_dict[psm] for pep in membrane_map for psm in clean_psm_dict[pep]}


def compare_between2(combined_search_dict,target_search_dict,target_pin_files_dict):
    """
    compare results between two searches
    :param combined_search_dict: filtered search results from combined search
    :param target_search_dict: filtered seach results from target search
    :param target_pin_files_dict: target pin files info dictionary
    :return:
    """
    from tsv_reader import venn_diagram_gen
    combined_psms = set([psm+'.'+tp[1] for f in combined_search_dict for psm in combined_search_dict[f] for tp in combined_search_dict[f][psm]])
    target_psms = set([psm+'.'+tp[1] for f in target_search_dict for psm in target_search_dict[f] for tp in target_search_dict[f][psm]])
    # venn_diagram_gen({'PSMs of combined search':combined_psms,'PSMs of target search':target_psms},title='Colon membrane PSMs compare')

    # get overlapped psms
    overlaped = [each for each in target_psms if each in combined_psms]

    # write PSM into a csv
    # file, spec_no, psm, charge, prot, hyperscore
    df_array = [[each.split('.')[-2],each.split('.')[-1],'.'.join(each.split('.')[:-2]),
                 int(target_pin_files_dict[each].split('\t')[0].split('.')[-1].split('_')[0]),
                 target_pin_files_dict[each].split('\t')[-1].rstrip('\n'),
                 float(target_pin_files_dict[each].split('\t')[9])]
                for each in overlaped]
    df = pd.DataFrame(df_array,columns=['file_name','spec_no','PSM','charge','Protein','hyper_score'])
    df.to_csv('F:/Colon/overlapped_target_PSM_from_unfilter.csv',sep=',')
    # print (len(overlaped))
    # for psms in target result but not combined, compare the hyper score
    target_unique_psms = [each for each in target_psms if each not in combined_psms]


def read_pin_into_dict(pin_f_list:list,pickle_output):
    # read lines of pin files into dict with psm.file_name.spec_no as key, the entire line as value
    line_dict = {}
    for f in pin_f_list:
        with open(f,'r') as f_o:
            next(f_o)
            for line in f_o:
                if '_Hagay_' in line and 'rev_' not in line: # save only target hits (w/wo human proteome match)
                    line_split = line.split('\t')
                    psm = line_split[19][2:-2]
                    spec = '.'.join(line_split[0].split('.')[:2])
                    line_dict[psm+'.'+spec] = line
        print (f.split('\\')[-1])
    print (f'{len(line_dict)} target PSMs saved')
    pickle.dump(line_dict,open(pickle_output,'wb'))


def prosit_prepare():
    # output the msp file need for prosit prediction
    df = pd.read_csv('F:/Colon/overlapped_target_PSM_from_unfilter.csv')

    psm_list, charge_list = df['PSM'].tolist(), df['charge'].tolist()
    prosit_array = []
    for psm, charge in zip(psm_list, charge_list):
        clean_psm = re.sub('M\[15\.9949\]','M(ox)',psm)  # replace M oxidation with M(ox)
        clean_psm = re.sub('C\[57\.0215\]','C',clean_psm)  # replace C[] with C
        clean_psm = re.sub('n\[42\.0106\]','',clean_psm)  # delete n-term acetylation
        if 7<=len(clean_psm.replace('(ox)',''))<=30:
            prosit_array.append([clean_psm,30,charge])
        else:
            print (clean_psm)
    new_df = pd.DataFrame(prosit_array, columns=['modified_sequence', 'collision_energy', 'precursor_charge'])
    new_df = new_df.set_index('modified_sequence')
    new_df.to_csv('F:/Colon/overlapped_prosit.csv')


def prepare_data():
    # prepare a data structure for cos sim cal. {peptidecharge:{file:[spec1,spec2...]}}
    df = pd.read_csv('F:/Colon/overlapped_target_PSM_from_unfilter.csv')
    info_dict = defaultdict(list)
    for tp in df.itertuples():
        psm = tp.PSM
        clean_psm = re.sub('M\[15\.9949\]', 'M(ox)', psm)  # replace M oxidation with M(ox)
        clean_psm = re.sub('C\[57\.0215\]', 'C', clean_psm)  # replace C[] with C
        clean_psm = re.sub('n\[42\.0106\]', '', clean_psm)  # delete n-term acetylation
        charge = str(tp.charge)
        f_name, spec = tp.file_name, tp.spec_no
        info_dict[clean_psm+charge].append((f_name,spec))
    pep_file_spec_dict_of_dict = {}
    for pep in info_dict:
        f_spec_list_dict=defaultdict(list)
        for each_tp in info_dict[pep]:
            f_spec_list_dict[each_tp[0]].append(each_tp[1])
        pep_file_spec_dict_of_dict[pep]=f_spec_list_dict
    pickle.dump(pep_file_spec_dict_of_dict,open('F:/Colon/prosit/pep_file_spec_dict_of_dict.p','wb'))


def copy_files(target_f_list, output_folder):
    import os
    import shutil
    for f in target_f_list:
        f_name = f.split('\\')[-1]
        if not os.path.exists(os.path.join(output_folder,f_name)):
            print (f)
            shutil.copy(f,output_folder)


if __name__ == '__main__':
    import pickle
    # combined_dict = {}
    # pin_f_list = glob('F:/Colon/search*/*/*_edited.pin')
    # for f in pin_f_list:
    #     f_name = f.split('\\')[-1].split('_edited')[0]
    #     print (f_name)
    #     psm_dict = filter_pin(f)
    #     combined_dict[f_name] = psm_dict
    # pickle.dump(combined_dict, open('F:/Colon/combined_search_psm_filter_dict.p', 'wb'))

    # human_prot_dict = fasta_reader('F:/Colon/Homo_sapiens_0318_canonical_and_isoform.fasta')
    #
    # target_dict = {}
    # target_pin_f_list = glob('F:/Colon/target_search*/*/*_edited.pin')
    # for f in target_pin_f_list:
    #     f_name = f.split('\\')[-1].split('_edited')[0]
    #     print (f_name)
    #     psm_dict = filter_pin_target(f,human_prot_dict)
    #     target_dict[f_name] = psm_dict
    # pickle.dump(target_dict,open('F:/Colon/target_search_psm_filter_dict.p','wb'))

    # combined_psm_dict = pickle.load(open('F:/Colon/combined_search_psm_filter_dict.p','rb'))
    # target_psm_dict = pickle.load(open('F:/Colon/target_search_psm_filter_dict.p','rb'))
    # target_pin_dict = pickle.load(open('F:/Colon/target_search_pin_files.p','rb'))
    # compare_between2(combined_psm_dict,target_psm_dict,target_pin_dict)

    # read all pin files
    # pin_f_list = glob('F:/Colon/target_search*/*/*_edited.pin')
    # read_pin_into_dict(pin_f_list,pickle_output='F:/Colon/target_search_pin_files.p')

    # prepare_data()

    # copy pin files to a folder
    pin_f_list = glob('F:/Colon/target_search*/*/*_edited.pin')
    output_folder = 'F:/Colon/unfiltered_pin_files/target_search'
    copy_files(pin_f_list,output_folder)