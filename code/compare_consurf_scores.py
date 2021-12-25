#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import numpy as np
from pathlib import Path
from text_tools import read_list_table, write_list_table
from stat_tools import sderror, t_test
from evolution_tools import read_consurf_score_file

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to orthology
    orthoDir = procDir / 'orthology'
    
    # directory of ConSurf scores for human proteins
    humanConsurfDir = orthoDir / 'human_protein_consurf'
    
    # directory of ConSurf scores for mouse proteins
    mouseConsurfDir = orthoDir / 'mouse_protein_consurf'
    
    # input data files
    natMutPPIFile = orthoDir / 'common_mutation_edgetic_PPIs.txt'
    disMutPPIFile = orthoDir / 'disease_mutation_edgetic_PPIs.txt'
    sequenceMapFile = orthoDir / 'mouse_human_map.txt'
    
    # output data files
    natMutPPIOutFile = orthoDir / 'common_mutation_edgetic_PPIs_interface_consurf.txt'
    disMutPPIOutFile = orthoDir / 'disease_mutation_edgetic_PPIs_interface_consurf.txt'
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    natPPIs = read_list_table (natMutPPIFile, ['Protein_interface'], [int])
    disPPIs = read_list_table (disMutPPIFile, ['Protein_interface'], [int])
    
    seqMapTable = read_list_table (sequenceMapFile, ["Qpos", "Spos"], [int, int])
    seqMap = {}
    for _, row in seqMapTable.iterrows():
        seqMap[(row.Subject, row.Query)] = (row.Spos, row.Qpos, row.Expect)
    
    ortholog_interfaces, eValues = [], []
    for _, row in natPPIs.iterrows():
        pSeq, oSeq, evalue = seqMap[(row.Protein, row.Protein_ortholog)]
        intf = [oSeq[pSeq.index(i)] for i in row.Protein_interface if i in pSeq]
        ortholog_interfaces.append(intf)
        eValues.append(evalue)
    natPPIs['Ortholog_interface'] = ortholog_interfaces
    natPPIs['E-value'] = eValues
    
    ortholog_interfaces, eValues = [], []
    for _, row in disPPIs.iterrows():
        pSeq, oSeq, evalue = seqMap[(row.Protein, row.Protein_ortholog)]
        intf = [oSeq[pSeq.index(i)] for i in row.Protein_interface if i in pSeq]
        ortholog_interfaces.append(intf)
        eValues.append(evalue)
    disPPIs['Ortholog_interface'] = ortholog_interfaces
    disPPIs['E-value'] = eValues
    
    human_avgScore, mouse_avgScore = [], []    
    for _, row in natPPIs.iterrows():
        consurfFile = humanConsurfDir / (row.Protein + '.grades.txt')
        if consurfFile.is_file():
            scores = read_consurf_score_file (consurfFile)
            human_intf_scores = [scores[pos]['SCORE'] for pos in row.Protein_interface]
            human_avgScore.append(np.mean(human_intf_scores))
        else:
            human_avgScore.append(np.nan)
        
        consurfFile = mouseConsurfDir / (row.Protein_ortholog + '.grades.txt')
        if consurfFile.is_file():
            scores = read_consurf_score_file (consurfFile)
            mouse_intf_scores = [scores[pos]['SCORE'] for pos in row.Ortholog_interface]
            mouse_avgScore.append(np.mean(mouse_intf_scores))
        else:
            mouse_avgScore.append(np.nan)
    
    natPPIs['Protein_intf_consurf_score'] = human_avgScore
    natPPIs['Ortholog_intf_consurf_score'] = mouse_avgScore
    
    human_avgScore, mouse_avgScore = [], []    
    for _, row in disPPIs.iterrows():
        consurfFile = humanConsurfDir / (row.Protein + '.grades.txt')
        if consurfFile.is_file():
            scores = read_consurf_score_file (consurfFile)
            human_intf_scores = [scores[pos]['SCORE'] for pos in row.Protein_interface]
            human_avgScore.append(np.mean(human_intf_scores))
        else:
            human_avgScore.append(np.nan)
        
        consurfFile = mouseConsurfDir / (row.Protein_ortholog + '.grades.txt')
        if consurfFile.is_file():
            scores = read_consurf_score_file (consurfFile)
            mouse_intf_scores = [scores[pos]['SCORE'] for pos in row.Ortholog_interface]
            mouse_avgScore.append(np.mean(mouse_intf_scores))
        else:
            mouse_avgScore.append(np.nan)
    
    disPPIs['Protein_intf_consurf_score'] = human_avgScore
    disPPIs['Ortholog_intf_consurf_score'] = mouse_avgScore
    
    natPPIs = natPPIs[(np.isnan(natPPIs['Protein_intf_consurf_score']) == False) &
                      (np.isnan(natPPIs['Ortholog_intf_consurf_score']) == False)].reset_index(drop=True)
    disPPIs = disPPIs[(np.isnan(disPPIs['Protein_intf_consurf_score']) == False) &
                      (np.isnan(disPPIs['Ortholog_intf_consurf_score']) == False)].reset_index(drop=True)
    
    write_list_table (natPPIs, ['Protein_interface', 'Ortholog_interface'], natMutPPIOutFile)
    write_list_table (disPPIs, ['Protein_interface', 'Ortholog_interface'], disMutPPIOutFile)
    
    print()
    print('Average score for PPI interfaces disrupted by common mutations:')
    print('Human proteins: %.3f (SE = %g, n = %d)' % (natPPIs['Protein_intf_consurf_score'].mean(),
                                                      sderror(natPPIs['Protein_intf_consurf_score'].values),
                                                      len(natPPIs)))
    print('Mouse orthologs: %.3f (SE = %g, n = %d)' % (natPPIs['Ortholog_intf_consurf_score'].mean(),
                                                       sderror(natPPIs['Ortholog_intf_consurf_score'].values),
                                                       len(natPPIs)))
    t_test (natPPIs['Protein_intf_consurf_score'].values, natPPIs['Ortholog_intf_consurf_score'].values)
    
    print()
    print('Average score for PPI interfaces disrupted by disease mutations:')
    print('Human proteins: %.3f (SE = %g, n = %d)' % (disPPIs['Protein_intf_consurf_score'].mean(),
                                                      sderror(disPPIs['Protein_intf_consurf_score'].values),
                                                      len(disPPIs)))
    print('Mouse orthologs: %.3f (SE = %g, n = %d)' % (disPPIs['Ortholog_intf_consurf_score'].mean(),
                                                       sderror(disPPIs['Ortholog_intf_consurf_score'].values),
                                                       len(disPPIs)))
    t_test (disPPIs['Protein_intf_consurf_score'].values, disPPIs['Ortholog_intf_consurf_score'].values)
    
if __name__ == "__main__":
    main()
