from GeneralUtilityMethods import GUM


class FoldX(object):

    def prepare_for_FoldX_BuildModel(self, mutant_name, pdb_name_of_repaired_pdb):
        path_to_runscript = './'
        repaired_pdbs = 'RepairPDB_' + pdb_name_of_repaired_pdb + '.pdb'
        show_sequence_detail = False
        action = '<BuildModel>#,individual_list.txt'
        print_networks = False
        calculate_stability = False
        GUM.build_runscript_for_pdbs(path_to_runscript, repaired_pdbs,
                                     show_sequence_detail, action, print_networks,
                                     calculate_stability)
        individual_list_for_this_mutant_only = open('individual_list.txt', 'w')
        individual_list_for_this_mutant_only.write(mutant_name + ';\n')
        individual_list_for_this_mutant_only.close()


    def prepare_for_FoldX_AnalyseComplex(self, repair_pdb_name):
        _0_1_2_pdbs = ['0.pdb,', '1.pdb,', '2.pdb,']
        repair_pdb_name_1_ = repair_pdb_name + '_1_'
        wt_repair_pdb_name_1_ = 'WT_' + repair_pdb_name_1_

        path_to_runscript = './'
        pdbs_to_analyse = repair_pdb_name_1_ + _0_1_2_pdbs[0] + \
                          repair_pdb_name_1_ + _0_1_2_pdbs[1] + \
                          repair_pdb_name_1_ + _0_1_2_pdbs[2] + \
                          wt_repair_pdb_name_1_ + _0_1_2_pdbs[0] + \
                          wt_repair_pdb_name_1_ + _0_1_2_pdbs[1] + \
                          wt_repair_pdb_name_1_ + _0_1_2_pdbs[2]
        show_sequence_detail = False
        action = '<AnalyseComplex>#'
        print_networks = False
        calculate_stability = False
        GUM.build_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse,
                                     show_sequence_detail, action, print_networks,
                                     calculate_stability)
