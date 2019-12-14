#----------------------------------------------------------------------------------------
# Produce a structural interactome annotated with PPI binding sites merged from PPI
# interfaces.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from structural_annotation import produce_site_annotated_interactome

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # minimum fraction of interface overlap needed to merge interfaces into one binding site
    overlap_cutoff = 0.5
    
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
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    interactomeSiteFile = interactomeDir / 'human_site_annotated_interactome.txt'
    
    # create directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # produce site-annotated interactome from interface-annotated interactome
    #------------------------------------------------------------------------------------
    
    produce_site_annotated_interactome (structuralInteractomeFile,
                                        interactomeSiteFile,
                                        cutoff = overlap_cutoff)

if __name__ == "__main__":
    main()
