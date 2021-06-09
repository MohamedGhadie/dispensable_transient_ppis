#----------------------------------------------------------------------------------------
# Re-write structural interactome
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of processed data files specific to interactome
    interactomeDir = dataDir / interactome_name
    
    # input data files
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    outFile = interactomeDir / 'structural_interactome_2.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    #------------------------------------------------------------------------------------
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile, sep='\t')
    
    interactome["Template_ID"] = interactome["Template_file_ID"].apply(lambda x: x.replace('!', ''))
        
    chain_1, chain_2 = [], []
    for _, row in interactome.iterrows():
        complexID, templateID = row.Alignment_file_ID.split('_')
        p1, p2 = complexID.split('=')
        pdbid, c1, c2 = templateID.split('-')
        c1 = c1.replace('!', '')
        c2 = c2.replace('!', '')
        mapping = {p1:c1, p2:c2}
        chain_1.append(mapping[row.Protein_1])
        chain_2.append(mapping[row.Protein_2])
    interactome["Protein_1_template_chain"] = chain_1
    interactome["Protein_2_template_chain"] = chain_2
    
    interactome["Protein_1_interface_residues"] = interactome["Interfaces"].apply(lambda x: x.split('+')[0])
    interactome["Protein_2_interface_residues"] = interactome["Interfaces"].apply(lambda x: x.split('+')[1])
    
    interactome.drop (columns = ['Complex_ID', 'Template_file_ID', 'Alignment_file_ID',
                                 'Interfaces', 'Chain_pairs'], inplace=True)
    
    interactome.to_csv (outFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
