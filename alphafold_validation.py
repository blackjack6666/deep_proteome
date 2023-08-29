### some functions and scripts used in Alphafold validation data processing (partially cleaned)

def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def get_unique_peptide(list_of_peptsv:list):
    """
    from pep tsv file only get unique peptides compared with previous ones, e.g. in 4 hour sample, filter out peptides
    in 1h,2h and only retain peptides uniquely identified in 4h
    :param list_of_peptide:
    :return:
    """

    unique_peptide_dict = {}
    peptide_list = []
    for idx, val in enumerate(list_of_peptsv):
        if '\\' in val:
            file_name = val.split('\\')[-2]
        else:
            file_name = val.split("/")[-2]
        print (file_name)
        unique_peptide_list = [each for each in peptide_counting(val) if each not in peptide_list]

        peptide_list += unique_peptide_list

        unique_peptide_dict[file_name] = unique_peptide_list

    return unique_peptide_dict


def peptide_counting(peptide_tsv_file):

    with open(peptide_tsv_file, 'r') as file_open:
        next(file_open)

        peptide_list = [line.split("\t")[0] for line in file_open]
    return peptide_list


def protein_tsv_reader(protein_tsv_file, protein_column=1):
    """

    :param protein_tsv_file:
    :param protein_column: the column number corresponding to UniprotID
    :return:
    """
    with open(protein_tsv_file, 'r') as file_open:
        next(file_open)
        return [line.split("\t")[protein_column] for line in file_open]


def mapping_KR_toarray(psm_list, protein_dict):
    """
    only map the cleavage sites which are start and end of peptide
    :param psm_list:
    :param protein_dict:
    :return:
    """
    id_KR_array_dict = {}
    id_KR_index_dict = {}

    # aho mapping
    id_list, seq_list = commons.extract_UNID_and_seq(protein_dict)
    seq_line = commons.creat_total_seq_line(seq_list, sep="|")
    zero_line = commons.zero_line_for_seq(seq_line)
    separtor_pos_array = commons.separator_pos(seq_line)

    aho_result = automaton_matching(automaton_trie([pep for pep in psm_list]), seq_line)
    for tp in aho_result:
        # matched_pep = tp[2]  # without ptm site
        # print (matched_pep,seq_line[tp[0]-1],seq_line[tp[1]])
        zero_line[tp[0]-1]+=1  # map the start and end of peptide to the array
        zero_line[tp[1]]+=1
    for i in range(len(separtor_pos_array)-1):
        zero_line_slice = zero_line[separtor_pos_array[i]+1:separtor_pos_array[i+1]]
        # percentage_cov = np.count_nonzero(zero_line_slice)/len(zero_line_slice)*100
        # if percentage_cov != 0:
        id_KR_array_dict[id_list[i]] = zero_line_slice
        id_KR_index_dict[id_list[i]] = np.nonzero(zero_line_slice)[0]

    return id_KR_array_dict, id_KR_index_dict


def plddt_retrieve(alphafold_pdb):
    """
    get pLDDT value from alphafold pdb file
    :param alphafold_pdb:
    :return:
    """
    with open(alphafold_pdb,'r') as f_o:
        file_split = f_o.read().split('\nATOM')[1:]

        return [float(re.search('\d+\.\d+(?=\s+[A-Z])',line).group()) for line in file_split]


def find_centroid(residue_atom_xyz, centroid_method='mean'):
    """
    find geometric center of protein given xyz coordinates
    :param residue_atom_xyz:
    :return: use median as centroid
    """
    atom_coordinates = np.array([coord for each in residue_atom_xyz for coord in residue_atom_xyz[each]])

    if centroid_method == 'mean':
        return np.mean(atom_coordinates,axis=0)
    elif centroid_method == 'median':
        return np.median(atom_coordinates,axis=0)


def residue_distance(residue_atom_xyz):
    """
    calculate the distance for each residue to center of 3d structure
    :param residue_atom_xyz:
    :return:
    """
    # zero_point = np.array([0,0,0])
    zero_point = find_centroid(residue_atom_xyz)
    residue_distance_dict = {}
    for each_pos in residue_atom_xyz:

        total_dist = sum([np.linalg.norm(np.array(each_atom)-zero_point)
                          for each_atom in residue_atom_xyz[each_pos]])

        average_dist = total_dist/len(residue_atom_xyz[each_pos])
        residue_distance_dict[each_pos] = average_dist
    return residue_distance_dict


def sasa_pdb(input_tuple, protease='trypsin'):
    """
    compute solvent accessible surface area (SASA) given a pdb file
    :param pdb_file:
    :return:
    """
    import pymol
    from commons import expasy_rules
    pdb_file, protein_seq = input_tuple
    cleavage_index = [m.end() for m in re.finditer(expasy_rules[protease], protein_seq)]

    # pymol.pymol_argv = ['pymol', '-qc']  # pymol launching: quiet (-q), without GUI (-c)
    # pymol.finish_launching()

    pymol.cmd.set('dot_solvent', 1)
    pymol.cmd.set('dot_density', 3)  # surface area
    pymol.cmd.set('solvent_radius', 15)  # same radius as trypsin

    pdb_name = os.path.split(pdb_file)[1]
    pymol.cmd.load(pdb_file, pdb_name)
    # residues = []
    # cmd.iterate('all', 'residues.append(resi)') # iterate residues and store into the list

    residue_sasa_dict = {}
    for i in cleavage_index:
        residue_sasa_dict[i] = pymol.cmd.get_area('resi %s' % i)
    pymol.cmd.delete(pdb_name)
    print(pdb_file + ' done')
    return {pdb_file.split('\\')[-1].split('-')[1]: residue_sasa_dict}


def inSphere2(ref, atoms, radius):
    """
    faster than inSphere, calculate multiple atoms at once
    :param ref: reference point
    :param atoms: atoms in the filtered cube space, 2D numpy array
    :param radius:
    :return: number of atoms in the sphere radius
    """
    ### filter surrounding atoms in a cube with equivalent radius
    # atoms = np.array(atoms) if not type(atoms) == np.ndarray else atoms

    cube_edge_upper, cubu_edge_lower = np.sum([ref, [radius, radius, radius]], axis=0), \
                                       np.sum([ref, [-radius, -radius, -radius]], axis=0)
    # get atoms inside the cube, lower edge<atom<upper edge

    filtered_atoms = [np.sum(each <= cube_edge_upper) == 3 and
                      np.sum(each >= cubu_edge_lower) == 3 for each in atoms]
    filtered_atoms_number = np.sum(filtered_atoms)
    # print(f'orginal number of atoms: {len(atoms)}, after filter: {filtered_atoms_number}')

    filtered_atoms_coords = atoms[np.array(filtered_atoms)]

    ### check if in the sphere, calculate euclidean distance
    distance_array = np.linalg.norm(filtered_atoms_coords - np.tile(ref, [filtered_atoms_number, 1]), axis=1)

    return len(np.where(distance_array <= radius)[0])

def residue_density_cal2(input_tuple, protease='trypsin', radius=29):
    """
    calculate number of atoms within certain range of a residue
    :param input_tuple:
    :param protease:
    :param radius_power2: radius power of protease
    :return:
    """

    # chymotrypsin radius in water =2.1 nm ,reference https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2567952/
    from commons import expasy_rules
    time_start = time.time()
    alphafold_pdb_file, protein_seq = input_tuple
    cleavage_density_dict = {}

    cleavage_index = [m.end() for m in re.finditer(expasy_rules[protease], protein_seq)]
    residue_atom_coord_dict = pdb_file_reader(alphafold_pdb_file)
    xyz_nparray = [each for v in residue_atom_coord_dict.values() for each in v]  # 2D

    xyz_2d_reshape = np.reshape(xyz_nparray, (-1, 3))  # reshape atom coords into 2d array

    for each in cleavage_index:
        ref = residue_atom_coord_dict[each][-1]
        # bool_array = [inSphere(i, ref, radius*radius) for i in xyz_nparray]
        # num_resi_inrange = np.count_nonzero(bool_array)
        num_resi_inrange = inSphere2(ref, xyz_2d_reshape, radius)
        cleavage_density_dict[each] = num_resi_inrange

    print(alphafold_pdb_file.split('\\')[-1].split('-')[1] + ' done.')
    # print (f'time used: {time.time()-time_start}')
    return {alphafold_pdb_file.split('\\')[-1].split('-')[1]: cleavage_density_dict}


def pdb_file_reader(pdb_file):
    """
    read pdb file and map xyz coordinates of each residue (alphafold pdbs)
    :param pdb_file:
    :return:
    """
    with open(pdb_file,'r') as f_o:
        file_split = f_o.read().split('\nATOM')[1:]

    residue_atom_xyz = defaultdict(list)

    # append xyz coords of each atom to residue position
    for line in file_split:
        # dont take hydrogen into account
        if line.split('           ')[1].split('  ')[0] != 'H':
            residue_atom_xyz[int(re.search('\d+(?=\s+[+-]?\d+\.)', line).group())].append(
                [float(i) for i in re.findall(r'[+-]?\d+\.\d{3}', line)])
        else:
            continue
    # print (residue_atom_xyz)
    return residue_atom_xyz


### protein dictionary and peptide lists
fasta_file = 'D:/data/pats/human_fasta/uniprot-proteome_UP000005640_sp_tr.fasta'
protein_dict = fasta_reader(fasta_file)

protein_tsv = 'D:/data/native_protein_digestion/12072021/control/combined_protein.tsv'
protein_list = protein_tsv_reader(protein_tsv, protein_column=3)
sub_protein_dict = {prot: protein_dict[prot] for prot in protein_list}

# base_path = 'F:/native_digestion/trypsin_lysc_5_25/search/'
# folders = [base_path + folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]
# time_points = [each.split('/')[-1] for each in folders]
# pep_path_list = [each + '/peptide.tsv' for each in folders]
# psm_path_list = [each + '/psm.tsv' for each in folders]
# unique_peptide_dict = get_unique_peptide(pep_path_list)

### calculate covered distance/average pLDDT and write to excel
### e.g. /Alphafold_validation/trypsin_1207/cov_dist_unique.xlsx
"""
df = pd.DataFrame(index=protein_list, columns=time_points)  # some protein entry does not have pdb

for pep_tsv in pep_path_list:
    print (pep_tsv)
    # peptide_list = peptide_counting(pep_tsv)
    peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]

    # if peptide_list:
    # freq_array_dict = freq_ptm_index_gen_batch_v2(peptide_list,protein_dict)[0]
    freq_array_dict = mapping_KR_toarray(peptide_list, sub_protein_dict)[0]
    for prot in protein_list:
        pdb_file_path = pdb_base + 'AF-' + prot + '-F1-model_v1.pdb'
        if os.path.exists(pdb_file_path):
            residue_dist_dict = residue_distance(pdb_file_reader(pdb_file_path))
            # plddt_dict = residue_plddt_retrieve(pdb_file_path)
            if len(residue_dist_dict) == len(protein_dict[prot]):  # filter out those really long proteins
                # if len(plddt_dict) == len(protein_dict[prot]):
                freq_array = freq_array_dict[prot]
                # print (np.count_nonzero(freq_array))
                cov_dist = cov_distance(freq_array, residue_dist_dict)
                # print (cov_dist)
                # ave_cov_plddt = cov_plddt(freq_array,plddt_dict)
                df.at[prot, pep_tsv.split('/')[-2]] = cov_dist
                # df.at[prot,pep_tsv.split('/')[-2]] = ave_cov_plddt
            else:
                print('%s protein len between pdb and fasta is not same' % prot)
        else:
            continue
    # else:
    #     for prot in protein_list:
    #         df.at[prot, pep_tsv.split('/')[-2]] = np.nan
df.to_excel('F:/native_digestion/trypsin_lysc_5_25/search/distance.xlsx')
"""

### calculate covered K/R density or solvent accessibility area and write to excel
### e.g. /Alphafold_validation/trypsin_1207/cov_KR_density_15A.xlsx
### or /Alphafold_validation/trypsin_1207/sasa.xlsx
"""
import pandas as pd

df = pd.DataFrame(index=protein_list, columns=time_points)  # some protein entry does not have pdb
# solven_acc_dict = pickle.load(open('D:/data/alphafold_pdb/1207control_protein_KR_sasa_dict.pkl','rb'))
KR_density_alpha_dict = pickle.load(open('D:/data/alphafold_pdb/human_file_KR_density_dict.pkl', 'rb'))
# chymo_cleav_density_dict = pickle.load(
#     open('D:/data/alphafold_pdb/688_prot_chymotry_cleave_density_dict.pkl', 'rb'))
for pep_tsv in pep_path_list:
    print(pep_tsv)
    # peptide_list = peptide_counting(pep_tsv)
    peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]
    freq_array_dict, freq_array_index_dict = mapping_KR_toarray(peptide_list, sub_protein_dict)
    for prot in protein_list:
        print (prot)
        pdb_file_path = pdb_base + 'AF-' + prot + '-F1-model_v1.pdb'
        if os.path.exists(pdb_file_path):
            # residue_dist_dict = residue_distance(pdb_file_reader(pdb_file_path))
            # plddt_dict = residue_plddt_retrieve(pdb_file_path)
            # solvent_access_dict = solven_acc_dict[prot]
            # if len(residue_dist_dict) == len(protein_dict[prot]):  # filter out those really long proteins
            if len(alphafold_protein_dict[pdb_file_path.split('/')[-1]]) == len(protein_dict[prot]):
                freq_array = freq_array_dict[prot]
                ave_KR_density = cov_KR_density(freq_array, KR_density_alpha_dict[prot])
                # ave_solvent_access = cov_KR_density(freq_array,solvent_access_dict)

                df.at[prot, pep_tsv.split('/')[-2]] = ave_KR_density
                # df.at[prot, pep_tsv.split('/')[-2]] = ave_solvent_access
                # df.at[prot, pep_tsv.split('/')[-2]] = freq_array_index_dict[prot]
            else:
                print('%s protein len between pdb and fasta is not same' % prot)
        else: 
            continue
df.to_excel('F:/native_digestion/trypsin_lysc_5_25/search/cov_KR_density_15A.xlsx')
"""

### spearman correlation analysis of cleaved K/R densities between control and heat shock
### e.g. /Alphafold_validation/trypsin_1207/spearman_corr_pval_nofill.xlsx
"""
df_control = pd.read_excel('F:/native_digestion/chymotrypsin_4_16/search/cov_chymo_density.xlsx',index_col=0)
# df_heatshock = pd.read_excel('D:/data/native_protein_digestion/12072021/heat_shock/cov_KR_density_heatshock.xlsx',index_col=0)
from scipy.stats import spearmanr

df_control_median = df_control.median().tolist()[:8]
# df_control_fill = df_control.fillna(df_control.median())

df_spearman = pd.DataFrame(index=df_control.index, columns=['spearman correlation', 'p value'])
for tp in df_control.itertuples():
    prot, kr_densities = tp[0], tp[1:9]
    try:
        corr, p_val = spearmanr(kr_densities,df_control_median,nan_policy='omit')
        df_spearman.at[prot,'spearman correlation'] = corr
        df_spearman.at[prot,'p value'] = p_val
    except: # incase spearman couldn't perform for low number of data points
        df_spearman.at[prot, 'spearman correlation'] = 0
        df_spearman.at[prot,'p value'] = 0
df_spearman.to_excel('F:/native_digestion/chymotrypsin_4_16/search/spearman_atom_dens_5_240min.xlsx')
"""

### linear regression analysis, calcuate the slope and R-square
### e.g.  /Alphafold_validation/trypsin_1207/KR_atoms_linear_reg.xlsx
"""
from scipy.stats import linregress
df_dist = pd.read_excel('F:/native_digestion/chymotrypsin_4_16/search/cov_chymo_density.xlsx',index_col=0)
df_out = pd.DataFrame(columns=['slope','r_square','p_val'])
for tp in df_dist.itertuples():
    prot, kr_densities = tp[0], np.array(tp[1:-2])
    cleaned_kr_densities = kr_densities[np.isfinite(kr_densities)]
    if len(cleaned_kr_densities) == 1 or len(cleaned_kr_densities) == 0:
        slope, r_squre, p_val = np.nan, np.nan, np.nan
    elif len(cleaned_kr_densities) == 2:
        slope, r_squre, p_val = np.nan, np.nan, np.nan
    else:
        x = range(0,len(cleaned_kr_densities))
        result = linregress(x,cleaned_kr_densities)
        slope, r_squre, p_val = result.slope, np.square(result.rvalue), result.pvalue
    df_out.at[prot,:] = slope, r_squre, p_val

df_out.to_excel('F:/native_digestion/chymotrypsin_4_16/search/atom_den_5_1200min_linear_reg.xlsx')
"""

### implement amino acid exposure data from structuremap (full exposure, value is number of neiboring amino acids)
### https://github.com/MannLabs/structuremap_analysis/blob/master/data_analysis_structuremap.ipynb
### export excel: /Alphafold_validation/trypsin_1207/aa_exposure_structuremap.xlsx
"""

# df_full_exp = pd.read_csv('D:/data/alphafold_pdb/full_sphere_exposure.csv') 

## split big csv into chunks and save individually for higher reading speed

# pd_dicts = {key: df_full_exp.loc[value] for key, value in df_full_exp.groupby("protein_id").groups.items()}
# for k in pd_dicts:
#     print (k)
#     pd_dicts[k].to_csv(os.path.join('D:/data/alphafold_pdb/full_sphere_expo_split/',k+'.csv'))
# print (df_full_exp.head)

protein_tsv = 'D:/data/native_protein_digestion/12072021/control/combined_protein.tsv'
protein_list = protein_tsv_reader(protein_tsv, protein_column=3)
sub_protein_dict = {prot:protein_dict[prot] for prot in protein_list}

base_path = 'D:/data/native_protein_digestion/12072021/control/'
folders = [base_path + folder for folder in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, folder))]
time_points = [each.split('/')[-1] for each in folders]
pep_path_list = [each + '/peptide.tsv' for each in folders]
psm_path_list = [each + '/psm.tsv' for each in folders]
unique_peptide_dict = get_unique_peptide(pep_path_list)

df_native_exposure = pd.DataFrame(index=protein_list,columns=time_points)

for pep_tsv in pep_path_list:
    print (pep_tsv)
    peptide_list = unique_peptide_dict[pep_tsv.split('/')[-2]]
    freq_array_dict, freq_array_index_dict = map_k_r(peptide_list,sub_protein_dict)

    for prot in protein_list:
        print (prot)
        try:
            df_expo = pd.read_csv('D:/data/alphafold_pdb/full_sphere_expo_split/'+prot+'.csv')

            freq_index = freq_array_index_dict[prot]
            if len(freq_index)!=0: # has at least one peptide matched
                total_expo = sum([df_expo.at[each,'nAA_24_180_pae'] for each in freq_index])
                ave_expo = total_expo/len(freq_index)
                df_native_exposure.at[prot,pep_tsv.split('/')[-2]] = ave_expo
            else:
                df_native_exposure.at[prot,pep_tsv.split('/')[-2]] = np.nan
        except FileNotFoundError:
            df_native_exposure.at[prot, pep_tsv.split('/')[-2]] = np.nan
df_native_exposure.to_excel('D:/data/native_protein_digestion/12072021/control/aa_exposure_structuremap.xlsx')
"""

### protein fragments length analysis in native digestion, get the longest fragment at each time point
### e.g. /Alphafold_validation/trypsin_1207/digestion_max_peptide_relative_length.xlsx
"""
df_cleav_index = pd.read_excel('D:/data/native_protein_digestion/12072021/control/cleavage_index_4_24.xlsx',index_col=0)
selected_cols = df_cleav_index.columns[1:-2]
df_new = pd.DataFrame(index=df_cleav_index.index,columns=selected_cols)

# print (selected_cols)
for ind in df_cleav_index.index:
    seq_len = len(df_cleav_index.at[ind,'sequence'])
    for each_col in selected_cols:
        cells = df_cleav_index.loc[ind, '0010min':each_col].tolist()  # multiple columns

        cleav_array = ''.join([each for each in cells]).replace('[',' ').replace(']',' ').replace('\n','')
        np_array = np.append(np.fromstring(cleav_array,dtype=int, sep=' '),seq_len)
        unique_sort_array = np.sort(np.unique(np_array))
        substract_array = np.insert(unique_sort_array, 0, 0)[:-1]
        # print (unique_sort_array,substract_array)

        frag_ratio = np.subtract(unique_sort_array,substract_array)/seq_len # normalized by length of protein
        max_ind = np.argmax(frag_ratio)
        max_idx_first, max_idx_latter = substract_array[max_ind],unique_sort_array[max_ind]
        print (max_idx_first,max_idx_latter)
        # df_new.at[ind,each_col] = np.max(frag_ratio)
        df_new.at[ind,each_col] = (max_idx_first,max_idx_latter)
# df_new.to_excel('D:/data/native_protein_digestion/12072021/control/digestion_max_peptide_relative_length.xlsx')
df_new.to_excel('D:/data/native_protein_digestion/12072021/control/digestion_max_peptide_index.xlsx')
"""

### DisProt DB analysis (protein disordered structure ratio)
### e.g. /Alphafold_validation/trypsin_1207/digestion_max_peptide_disorder.xlsx
"""
from collections import defaultdict

df_disorder = pd.read_csv('D:/data/native_protein_digestion/12072021/control/DisProt release_2022_03 with_ambiguous_evidences (1).tsv',
                          delimiter='\t', index_col=0)
df_disorder = df_disorder[df_disorder['organism'] =='Homo sapiens']
protein_disoder_dict = defaultdict(set)
for each_tp in df_disorder.itertuples():
    protein_disoder_dict[each_tp.Index].add((each_tp.start,each_tp.end))


protein_disoder_ratio_dict = {}

for prot in protein_disoder_dict:
    if prot in protein_dict:
        protein_len = len(protein_dict[prot])
        zeroline = np.zeros(protein_len)
        for tp in protein_disoder_dict[prot]:
            zeroline[tp[0]-1:tp[1]]+=1
        percent = np.count_nonzero(zeroline)/protein_len*100
        protein_disoder_ratio_dict[prot] = percent

df_max_pep = pd.read_excel('D:/data/native_protein_digestion/12072021/control/digestion_max_peptide_relative_length_slope.xlsx',index_col=0)
df_max_pep['disorder_ratio'] = [protein_disoder_ratio_dict[each] if each in protein_disoder_ratio_dict else 0 for each in df_max_pep.index]
df_max_pep.to_excel('D:/data/native_protein_digestion/12072021/control/digestion_max_peptide_disorder.xlsx')
"""
