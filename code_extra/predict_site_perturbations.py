#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import read_single_interface_annotated_interactome
from mutation_interface_edgotype import mutation_site_perturbations

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'human_site_annotated_interactome.txt'
    naturalPerturbsFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    diseasePerturbsFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
    
    # output data files
    sitePerturbsFile = interactomeDir / 'site_perturbs_geometry.pkl'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    naturalPerturbs = read_list_table (naturalPerturbsFile,
                                       ["partners", "perturbations"],
                                       [str, int])
    diseasePerturbs = read_list_table (diseasePerturbsFile,
                                       ["partners", "perturbations"],
                                       [str, int])
    
    sites, perturbations = mutation_site_perturbations (naturalPerturbs, interactome)
    naturalPerturbs["sites"], naturalPerturbs["site_perturbations"] = sites, perturbations
    
    sites, perturbations = mutation_site_perturbations (diseasePerturbs, interactome)
    diseasePerturbs["sites"], diseasePerturbs["site_perturbations"] = sites, perturbations
    
    with open(sitePerturbsFile, 'wb') as fOut:
        pickle.dump([naturalPerturbs, diseasePerturbs], fOut)
    
if __name__ == "__main__":
    main()
