import re
import numpy as np
from collections import defaultdict


def modified_peptide_from_psm(psm_path):

    psm_list = []
    with open(psm_path, 'r') as f_open:
        next(f_open)
        for line in f_open:
            line_split = line.split('\t')
            match = re.search('\w{1}\[\d+\.?\d+\]',line)
            if match:
                psm_list.append(line_split[3])
            else:
                psm_list.append(line_split[2])
    return psm_list


def my_replace(match_obj):
    match_obj = match_obj.group()
    return match_obj[0]  # gives back the first element of matched object as string


def freq_array_and_PTM_index_generator(peptide_list, protein_seq_string,regex_pat='\w{1}\[\d+\.?\d+\]'):
    """
    map single protein seq
    :param peptide_list:
    :param protein_seq_string:
    :param regex_pat:
    :return:
    """

    freq_array = np.zeros(len(protein_seq_string))
    PTM_sites_counting = defaultdict(int)
    PTM_loc_list = []

    # reformat the peptide with PTM numbers into characters only
    new_pep_list = [re.sub(regex_pat, my_replace, pep) for pep in peptide_list]
    PTM_list = [re.findall(regex_pat, pep) for pep in peptide_list]
    # print (PTM_list)
    # calculation

    for pep, new_pep, PTM in zip(peptide_list, new_pep_list,PTM_list):  # PTM_list is list of list
        if new_pep in protein_seq_string:
            # print(new_pep)
            start_pos = protein_seq_string.find(new_pep)
            end_pos = start_pos + len(new_pep) -1
            # print (start_pos,end_pos,new_pep)
            freq_array[start_pos:end_pos + 1] += 1
            if PTM:  # the peptide has ptm site
                for ele in PTM:

                    PTM_index = pep.find(ele)
                    PTM_sites_counting[ele] += 1
                    PTM_loc_list.append(start_pos+PTM_index)
    # print (PTM_sites_counting, PTM_loc_list)

    return freq_array, PTM_loc_list, PTM_sites_counting


def one_d_cov_bar(html_template,peptide_list,protein_seq,output_html_path=None):
    # generate a one d coverage bar html

    freq_array, ptm_loc = freq_array_and_PTM_index_generator(peptide_list,protein_seq)[:2]
    print (ptm_loc)
    gray_rgb, cov_rgb, ptm_rgb = 'rgba(191,191,191,1) ','rgba(255,0,0,1) ','rgba(0,255,255,1) '
    replace_str = ''

    for i in range(len(protein_seq)):
        percen = f'{i / float(len(protein_seq)) * 100:.1f}%'
        percen_next = f'{(i + 1) / float(len(protein_seq)) * 100:.1f}%'
        # print (i,i+1,len(protein_seq),percen,percen_next)

        if i in set(ptm_loc):
            replace_str += ptm_rgb
            replace_str += percen
            replace_str += ', '
            replace_str += ptm_rgb
            replace_str += percen_next
            replace_str += ', '

        if freq_array[i] == 0:
            replace_str += gray_rgb  # not covered part show as grey
            replace_str += percen
            replace_str += ', '
            replace_str += gray_rgb
            replace_str += percen_next
            replace_str += ', '
        elif freq_array[i] != 0:
            replace_str += cov_rgb
            replace_str += percen
            replace_str += ', '
            replace_str += cov_rgb
            replace_str += percen_next
            replace_str += ', '

    replace_str = replace_str[:-2]  # delete ', ' in the end

    if output_html_path:
        f_open = open(html_template, 'r')
        new_html = f_open.read().replace('#insert_replace_here', replace_str)
        f_write = open(output_html_path, 'w')
        f_write.write(new_html)
        f_write.close()

        f_open.close()
    return replace_str


if __name__=='__main__':
    from protein_coverage import fasta_reader
    import time
    fasta = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_only.fasta'
    html_template = 'D:/data/Naba_deep_matrisome/html_template.html'
    psm_list = modified_peptide_from_psm('D:/data/native_protein_digestion/11182021/search_result_XS/0060min_XS/psm.tsv')
    protein_dict = fasta_reader(fasta)
    protein_id = 'Q9GZV4'
    start = time.time()
    one_d_cov_bar(html_template,psm_list,protein_dict[protein_id],'D:/data/temp/test.html')
    print (time.time()-start)