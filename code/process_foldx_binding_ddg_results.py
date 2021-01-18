#----------------------------------------------------------------------------------------
# Process binding ∆∆G results from FoldX. Since each result file 
# contains only one mutation, no second-round jobs should be produced.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import (read_foldx_results,
                          write_mutation_ddg_tofile,
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
    
    # directory of processed model-related data files specific to interactome
    modellingDir = Path('../../edgotype_fitness_effect/data/processed/') / interactome_name / 'model_based'
    
    # directory of PPI structural models
    modelDir = modellingDir / 'ppi_models'
    
    # directory of foldx results
    inDir = interactomeDir / 'foldx' / 'results'
    
    # directory of foldx output jobs
    outDir = interactomeDir / 'foldx'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_foldx_results (inDir)
    
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'nondis_mut_binding_ddg_foldx.txt',
                               interactomeDir / 'nondis_mut_binding_ddg_foldx_2.txt',
                               type = 'binding')
    write_mutation_ddg_tofile (processed,
                               interactomeDir / 'dis_mut_binding_ddg_foldx.txt',
                               interactomeDir / 'dis_mut_binding_ddg_foldx_2.txt',
                               type = 'binding')
    
    os.remove (interactomeDir / 'nondis_mut_binding_ddg_foldx.txt')
    os.remove (interactomeDir / 'dis_mut_binding_ddg_foldx.txt')
    os.rename (interactomeDir / 'nondis_mut_binding_ddg_foldx_2.txt',
               interactomeDir / 'nondis_mut_binding_ddg_foldx.txt')
    os.rename (interactomeDir / 'dis_mut_binding_ddg_foldx_2.txt',
               interactomeDir / 'dis_mut_binding_ddg_foldx.txt')
    
    produce_foldx_and_beluga_jobs (unprocessed,
                                   modelDir,
                                   outDir,
                                   'binding',
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
