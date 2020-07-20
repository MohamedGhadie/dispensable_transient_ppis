#----------------------------------------------------------------------------------------
# Produce jobs for PPI template structures to be submitted to FoldX for interaction energy 
# calculations. Each job contains only one PPI. PPIs with known energy values in the 
# input file are skipped.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import (produce_initial_ppi_template_energy_file,
                          read_unprocessed_energy_ppis,
                          produce_foldx_and_beluga_jobs)

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
    modelDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # directory of foldx output jobs
    outDir = interactomeDir / 'foldx'
    
    # input data files
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    energyFile = interactomeDir / 'ppi_template_energy_foldx.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    if not energyFile.is_file():
        produce_initial_ppi_template_energy_file (interactomeFile, energyFile)
    
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
