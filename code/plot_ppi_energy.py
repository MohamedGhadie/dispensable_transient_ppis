#----------------------------------------------------------------------------------------
# Plot distribution of PPI interaction energy.
#----------------------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
from plot_tools import multi_histogram_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
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
    
    # output files
    energyFile = interactomeDir / 'ppi_energy_foldx.txt'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    ppis = pd.read_table (energyFile, sep='\t')
    print('Average energy = %f' % ppis["Interaction_energy"].mean())
    multi_histogram_plot (ppis["Interaction_energy"].values,
                          xlabel = 'Interaction energy',
                          ylabel = 'Number of PPIs',
                          bins = 25,
                          fontsize = 24,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'ppi_energy_histogram')

if __name__ == "__main__":
    main()
