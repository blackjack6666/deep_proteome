# fragger_params created by ygao at 6/25/2021

fragger_params_dict = {'database_name': '', 'num_threads': '30', 'precursor_mass_lower': '-50', 'precursor_mass_upper': '50', 'precursor_mass_units': '1', 'precursor_true_tolerance': '20',
                       'precursor_true_units': '1', 'fragment_mass_tolerance': '20', 'fragment_mass_units': '1', 'calibrate_mass': '0', 'decoy_prefix': 'Rev_', 'deisotope': '1', 'isotope_error': '0/1/2',
                       'mass_offsets': '0', 'precursor_mass_mode': 'selected', 'remove_precursor_peak': '0', 'remove_precursor_range': '-1.500000,1.500000', 'intensity_transform': '0', 'write_calibrated_mgf': '0',
                       'mass_diff_to_variable_mod': '0', 'localize_delta_mass': '0', 'delta_mass_exclude_ranges': '(-1.5,3.5)', 'fragment_ion_series': 'b,y', 'ion_series_definitions': '', 'labile_search_mode': 'off',
                       'restrict_deltamass_to': 'all', 'diagnostic_intensity_filter': '0', 'Y_type_masses': '', 'diagnostic_fragments': '', 'search_enzyme_name': 'stricttrypsin', 'search_enzyme_cutafter': 'KR',
                       'search_enzyme_butnotafter': '', 'num_enzyme_termini': '2', 'allowed_missed_cleavage': '2', 'clip_nTerm_M': '1', 'allow_multiple_variable_mods_on_residue': '0',
                       'max_variable_mods_per_peptide': '3', 'max_variable_mods_combinations': '5000', 'output_file_extension': 'pepXML', 'output_format': 'tsv_pepxml_pin', 'output_report_topN': '5',
                       'output_max_expect': '50', 'report_alternative_proteins': '1', 'precursor_charge': '1 4', 'override_charge': '0', 'digest_min_length': '7', 'digest_max_length': '50',
                       'digest_mass_range': '500.0 5000.0', 'max_fragment_charge': '2', 'track_zero_topN': '0', 'zero_bin_accept_expect': '0', 'zero_bin_mult_expect': '1', 'add_topN_complementary': '0',
                       'minimum_peaks': '15', 'use_topN_peaks': '150', 'min_fragments_modelling': '2', 'min_matched_fragments': '4', 'minimum_ratio': '0.01', 'clear_mz_range': '0.0 0.0',
                       'add_Cterm_peptide': '0.000000', 'add_Nterm_peptide': '0.000000', 'add_Cterm_protein': '0.000000', 'add_Nterm_protein': '0.000000', 'add_G_glycine': '0.000000', 'add_A_alanine': '0.000000',
                       'add_S_serine': '0.000000', 'add_P_proline': '0.000000', 'add_V_valine': '0.000000', 'add_T_threonine': '0.000000', 'add_C_cysteine': '57.021464', 'add_L_leucine': '0.000000',
                       'add_I_isoleucine': '0.000000', 'add_N_asparagine': '0.000000', 'add_D_aspartic_acid': '0.000000', 'add_Q_glutamine': '0.000000', 'add_K_lysine': '0.000000', 'add_E_glutamic_acid': '0.000000',
                       'add_M_methionine': '0.000000', 'add_H_histidine': '0.000000', 'add_F_phenylalanine': '0.000000', 'add_R_arginine': '0.000000', 'add_Y_tyrosine': '0.000000', 'add_W_tryptophan': '0.000000',
                       'add_B_user_amino_acid': '0.000000', 'add_J_user_amino_acid': '0.000000', 'add_O_user_amino_acid': '0.000000', 'add_U_user_amino_acid': '0.000000', 'add_X_user_amino_acid': '0.000000',
                       'add_Z_user_amino_acid': '0.000000'}

def gen_fragger_params(params_dict=fragger_params_dict):
    fragger_params_dict['database_name'] = '/data/database/test.fasta'  # change any of the above parameters by accessing the dictionary key
    # add variable modification
    diff_mod_dict = {"[^": "42.010565", "M": "15.9949"}  # modified residue : mass ; this is equal to the entry of " [^ 42.010565" and "M 15.9949"
    i = 1
    for each_diff in diff_mod_dict:
        var_mod_key = "variable_mod_%02d" % i
        var_mod_value = "%s %s %s" % (each_diff, "".join(diff_mod_dict[each_diff]), "3")
        fragger_params_dict[var_mod_key] = var_mod_value
        i += 1
    # to change anything else
    fragger_params_dict['num_enzyme_termini'] = 2
    fragger_full = "".join(['%s = %s\n' % (key, fragger_params_dict[key]) for key in fragger_params_dict])
    return fragger_full  # fragger_full is the complete text string for the fragger.params

if __name__ == '__main__':
    with open('test.params', 'w') as fout:
        params_str = gen_fragger_params()
        fout.write(params_str)