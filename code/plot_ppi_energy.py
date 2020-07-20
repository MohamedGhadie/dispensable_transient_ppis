#----------------------------------------------------------------------------------------
# Plot distribution of PPI interaction energy.
#----------------------------------------------------------------------------------------

import os
import numpy as np
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
    
    energy = {}
    for _, ppi in ppis.iterrows():
        energy[(ppi.Complex_ID,) + tuple(sorted([ppi.Chain_1, ppi.Chain_2]))] = ppi.Interaction_energy
    
    values = list(energy.values())
    print('Unique structures = %d' % len(energy))
    print('Average energy = %f' % np.mean(values))
    print(sum([i <= 0 for i in values]))
    multi_histogram_plot (values,
                          xlabel = 'Interaction energy',
                          ylabel = 'Number of PPIs',
                          bins = 25,
                          fontsize = 16,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'ppi_energy_histogram')

if __name__ == "__main__":
    main()
