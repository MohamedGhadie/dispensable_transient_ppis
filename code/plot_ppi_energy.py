#----------------------------------------------------------------------------------------
# Plot distribution of PPI interaction energy.
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from energy_tools import read_ppi_energy, read_ppi_chain_energy
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
    
    # input data files
    energyFile = interactomeDir / 'ppi_template_energy_foldx.txt'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    energy = read_ppi_energy (energyFile)
    
    values = list(energy.values())
    numNeg = sum([i <= 0 for i in values])
    print('Unique PPI structures = %d' % len(energy))
    print('Average energy = %f' % np.mean(values))
    print('Median energy = %f' % np.median(values))
    print('Structures with energy ≤ 0 = %d (%.1f%%)' % (numNeg, 100 * numNeg / len(values)))
    multi_histogram_plot (values,
                          xlabel = 'Interaction energy',
                          ylabel = 'Number of PPIs',
                          bins = 50,
                          fontsize = 16,
                          xlim = [-250, 250],
                          #xlim = [-299, 299],
                          #xlabels = [-200, 0, 200, 400, 600],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'ppi_template_energy_histogram')

if __name__ == "__main__":
    main()
