#----------------------------------------------------------------------------------------
# Produce jobs for PPI structures to be submitted to FoldX for interaction energy 
# calculations. Each job contains only one PPI. PPIs with known energy values in the 
# input file are skipped.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from energy_tools import read_unprocessed_energy_ppis, produce_foldx_and_beluga_jobs

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of ppi structural models
    #modelDir = interactomeDir / 'ppi_models'
    modelDir = Path('/Volumes/MG_Samsung/edgotype_fitness_effect_full_model/data/processed/' +
                    interactome_name + '/model_based/ppi_models')
    
    # directory of foldx output jobs
    outDir = interactomeDir / 'foldx'
    
    # input file containing mutations to submit to foldx
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output files
    energyFile = interactomeDir / 'ppi_energy_foldx.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    if not energyFile.is_file():
        ppis = pd.read_table (interactomeFile, sep='\t')
        ppis["Interaction_energy"] = '-'
        chain_1, chain_2 = zip(* ppis["Chain_pairs"].apply(lambda x: x.split('+')).values)
        ppis["Chain_1"] = [c.split('_')[1] for c in chain_1]
        ppis["Chain_2"] = [c.split('_')[1] for c in chain_2] 
        ppis[["Protein_1",
              "Protein_2",
              "Complex_ID",
              "Chain_1",
              "Chain_2",
              "Interaction_energy"]].to_csv (energyFile, index=False, sep='\t')
    
    ppiStructures = read_unprocessed_energy_ppis (energyFile)
    produce_foldx_and_beluga_jobs (ppiStructures,
                                   modelDir,
                                   outDir,
                                   'ppi_energy',
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem_per_cpu = '4G',
                                   outputfile = '/project/ctb-yxia/ghadie84/foldx/data/%x-%j.out',
                                   username = 'ghadie84',
                                   serverDataDir = '/project/ctb-yxia/ghadie84/foldx/data')

if __name__ == "__main__":
    main()
