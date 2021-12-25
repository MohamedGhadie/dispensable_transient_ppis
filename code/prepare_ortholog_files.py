#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import io
import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table, write_list_table
from interactome_tools import read_single_interface_annotated_interactome

def main():
    
    # interactome names
    interactome_names = ['HuRI', 'IntAct']
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to orthology
    orthoDir = procDir / 'orthology'
    
    # input data files
    orthologFile = extDir / 'mart_export.txt'
    humanUniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    mouseUniprotIDmapFile = procDir / 'to_mouse_uniprotID_map.pkl'
    humanSequenceFile = procDir / 'human_reference_sequences.pkl'
    mouseSequenceFile = procDir / 'mouse_reference_sequences.pkl'
    
    # output data files
    natMutPPIFile = orthoDir / 'common_mutation_edgetic_PPIs.txt'
    disMutPPIFile = orthoDir / 'disease_mutation_edgetic_PPIs.txt'
    humanNatMutProteinListFile = orthoDir / 'human_edgetic_natMut_proteins.txt'
    humanDisMutProteinListFile = orthoDir / 'human_edgetic_disMut_proteins.txt'
    mouseNatMutOrthologListFile = orthoDir / 'mouse_edgetic_natMut_orthologs.txt'
    mouseDisMutOrthologListFile = orthoDir / 'mouse_edgetic_disMut_orthologs.txt'
    humanNatMutProteinSeqFile = orthoDir / 'human_edgetic_natMut_protein_sequences.txt'
    humanDisMutProteinSeqFile = orthoDir / 'human_edgetic_disMut_protein_sequences.txt'
    mouseNatMutOrthologSeqFile = orthoDir / 'mouse_edgetic_natMut_ortholog_sequences.txt'
    mouseDisMutOrthologSeqFile = orthoDir / 'mouse_edgetic_disMut_ortholog_sequences.txt'
    
    # create output directories if not existing
    if not orthoDir.exists():
        os.makedirs(orthoDir)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.DataFrame()
    diseaseMutations = pd.DataFrame()
    interfaces = {}
    
    for name in interactome_names:
        edgeticDir = procDir / name / 'physics' / (ddg_method + '_edgetics')
        naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
        diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
        
        natMut = read_list_table (naturalMutationsFile, ["partners", "perturbations"], [str, float])
        disMut = read_list_table (diseaseMutationsFile, ["partners", "perturbations"], [str, float])
        
        naturalMutations = naturalMutations.append(natMut, ignore_index=True)
        diseaseMutations = diseaseMutations.append(disMut, ignore_index=True)
        
        structuralInteractomeFile = procDir / name / 'structural_interactome.txt'
        interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
        for p1, p2, intf in interactome[['Protein_1', 'Protein_2', 'Interfaces']].values:
            if (p1, p2) not in interfaces:
                interfaces.update({(p1, p2):intf[0], (p2, p1):intf[1]})
    
    print()
    print('Total number of common mutations = %d ' % len(naturalMutations))
    print('Total number of disease mutations = %d ' % len(diseaseMutations))
    
    naturalMutations = naturalMutations.drop_duplicates (subset=["protein", "mut_position", "mut_res"])
    diseaseMutations = diseaseMutations.drop_duplicates (subset=["protein", "mut_position", "mut_res"])
    
    naturalMutations["perturbations"] = naturalMutations["perturbations"].apply(
                                            lambda x: [int(p) if not np.isnan(p) else p for p in x])
    diseaseMutations["perturbations"] = diseaseMutations["perturbations"].apply(
                                            lambda x: [int(p) if not np.isnan(p) else p for p in x])
    
    naturalMutations = naturalMutations [naturalMutations['edgotype'] == 'edgetic'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations['edgotype'] == 'edgetic'].reset_index(drop=True)
    
    print()
    print('Edgetic mutations after removing duplicates:')
    print('Common mutations = %d ' % len(naturalMutations))
    print('Disease mutations = %d ' % len(diseaseMutations))
    
    proteins, partners = [], []
    for _, row in naturalMutations.iterrows():
        p1 = row.protein
        pr = [p for p, pert in zip(row.partners, row.perturbations) if pert == 1]
        for p2 in pr:
            proteins.append(p1)
            partners.append(p2)
    natPPIs = pd.DataFrame(data={'Protein':proteins, 'Partner':partners})
    
    proteins, partners = [], []
    for _, row in diseaseMutations.iterrows():
        p1 = row.protein
        pr = [p for p, pert in zip(row.partners, row.perturbations) if pert == 1]
        for p2 in pr:
            proteins.append(p1)
            partners.append(p2)
    disPPIs = pd.DataFrame(data={'Protein':proteins, 'Partner':partners})
    
    disPPIs = disPPIs.drop_duplicates(keep='first')
    
    print()
    print('PPIs disrupted by common mutations = %d ' % len(natPPIs))
    print('PPIs disrupted by disease mutations = %d ' % len(disPPIs))
    
    allPPIs = disPPIs.append(natPPIs, ignore_index=True)
    allPPIs = allPPIs.drop_duplicates(keep='first')
    natPPIs = allPPIs[len(disPPIs):].reset_index(drop=True)
    
    print()
    print('Disrupted PPIs after removing duplicates among common and disease mutations:')
    print('Common mutations = %d ' % len(natPPIs))
    print('Disease mutations = %d ' % len(disPPIs))
    
    with open(humanUniprotIDmapFile, 'rb') as f:
        uniprotID_human = pickle.load(f)
    with open(mouseUniprotIDmapFile, 'rb') as f:
        uniprotID_mouse = pickle.load(f)
    
    with open(humanSequenceFile, 'rb') as f:
        humanSeq = pickle.load(f)
    with open(mouseSequenceFile, 'rb') as f:
        mouseSeq = pickle.load(f)
    
    ortho = pd.read_table(orthologFile, sep='\t')
    #ortho = ortho[ortho['Mouse homology type'] == 'ortholog_one2one']
    ortho = ortho.drop_duplicates(subset=['Gene stable ID'])
    ortho = ortho.drop_duplicates(subset=['Mouse gene stable ID'])
        
    ortholog = {}
    for _, row in ortho.iterrows():
        if row['Gene stable ID'] in uniprotID_human:
            humProt = uniprotID_human[row['Gene stable ID']]
            if humProt in humanSeq:
                if row['Mouse gene stable ID'] in uniprotID_mouse:
                    mouseProt = uniprotID_mouse[row['Mouse gene stable ID']]
                    if mouseProt in mouseSeq:
                        ortholog[humProt] = mouseProt
    
    print()
    print('Human proteins with mouse orthologs = %d' % len(ortholog))
    
    natPPIs = natPPIs[natPPIs['Protein'].apply(lambda x: x in ortholog)].reset_index(drop=True)
    disPPIs = disPPIs[disPPIs['Protein'].apply(lambda x: x in ortholog)].reset_index(drop=True)
    
    natPPIs['Protein_ortholog'] = [ortholog[p] for p in natPPIs['Protein'].values]
    #natPPIs['Partner_ortholog'] = [ortholog[p] for p in natPPIs['Partner'].values]
    disPPIs['Protein_ortholog'] = [ortholog[p] for p in disPPIs['Protein'].values]
    #disPPIs['Partner_ortholog'] = [ortholog[p] for p in disPPIs['Partner'].values]
    
    print()
    print('PPIs with disruptive mutations at the interface, and disrupted protein has ortholog:')
    print('Common mutations = %d ' % len(natPPIs))
    print('Disease mutations = %d ' % len(disPPIs))
    
    natPPIs['Protein_interface'] = natPPIs.apply(lambda x: interfaces[x['Protein'], x['Partner']], axis=1)
    disPPIs['Protein_interface'] = disPPIs.apply(lambda x: interfaces[x['Protein'], x['Partner']], axis=1)    
    
    write_list_table (natPPIs, ['Protein_interface'], natMutPPIFile)
    write_list_table (disPPIs, ['Protein_interface'], disMutPPIFile)
    
    humanNatMutProteins = sorted(set(natPPIs['Protein'].values))
    humanDisMutProteins = sorted(set(disPPIs['Protein'].values))
    mouseNatMutProteins = sorted(set(natPPIs['Protein_ortholog'].values))
    mouseDisMutProteins = sorted(set(disPPIs['Protein_ortholog'].values))    
    
    print()
    print('Unique human proteins with disruptive mutations at the interface:')
    print('Common mutations = %d ' % len(humanNatMutProteins))
    print('Disease mutations = %d ' % len(humanDisMutProteins))
    
    # write protein lists for human and mouse
    with io.open(humanNatMutProteinListFile, "w") as fout:
        for p in humanNatMutProteins:
            fout.write(p + '\n')
    
    with io.open(humanDisMutProteinListFile, "w") as fout:
        for p in humanDisMutProteins:
            fout.write(p + '\n')
    
    with io.open(mouseNatMutOrthologListFile, "w") as fout:
        for p in mouseNatMutProteins:
            fout.write(p + '\n')
    
    with io.open(mouseDisMutOrthologListFile, "w") as fout:
        for p in mouseDisMutProteins:
            fout.write(p + '\n')
    
    # write protein sequences for human and mouse
    with io.open(humanNatMutProteinSeqFile, "w") as fout:
        for p in humanNatMutProteins:
            fout.write('>' + p + '\n' + humanSeq[p] + '\n')
    
    with io.open(humanDisMutProteinSeqFile, "w") as fout:
        for p in humanDisMutProteins:
            fout.write('>' + p + '\n' + humanSeq[p] + '\n')
    
    with io.open(mouseNatMutOrthologSeqFile, "w") as fout:
        for p in mouseNatMutProteins:
            fout.write('>' + p + '\n' + mouseSeq[p] + '\n')
    
    with io.open(mouseDisMutOrthologSeqFile, "w") as fout:
        for p in mouseDisMutProteins:
            fout.write('>' + p + '\n' + mouseSeq[p] + '\n')
    
if __name__ == "__main__":
    main()
