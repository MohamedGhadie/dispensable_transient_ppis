#----------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from interactome_tools import read_chain_annotated_interactome

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    expDir = procDir / 'La_and_Mintseris'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    ppiTypeFile = expDir / 'structural_interactome.txt'
    chainAnnotatedInteractome = interactomeDir / 'chain_annotated_interactome.txt'
    
    # output data files
    outFile = expDir / (interactome_name + '_ppi_types.txt')
    
    # create output directories if not existing
    if not expDir.exists():
        os.makedirs(expDir)
    
    #------------------------------------------------------------------------------------
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    expPPIs = pd.read_table (ppiTypeFile, sep='\t')
    ppiTypes = {}
    for _, ppi in expPPIs.iterrows():
        k = '_'.join(sorted([ppi.Protein_1, ppi.Protein_2]))
        ppiTypes[k] = ppi.Type
    
    interactome = read_chain_annotated_interactome (chainAnnotatedInteractome)
    
    types, mapChains = [], []
    for _, ppi in interactome.iterrows():
        t, c = '-', '-'
        for c1, c2 in ppi.Mapping_chains:
            k = '_'.join(sorted([c1.replace('_', '-'), c2.replace('_', '-')]))
            if k in ppiTypes:
                t, c = ppiTypes[k], c1 + '+' + c2
                break
        types.append(t)
        mapChains.append(c)
        
    interactome["Type"] = types
    interactome["Mapping_chains"] = mapChains
    interactome = interactome[interactome["Type"] != '-']
    interactome.to_csv (outFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
