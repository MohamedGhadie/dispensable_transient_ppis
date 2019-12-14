#----------------------------------------------------------------------------------------
# Calculate the number of PPIs mediated per interfacial site.
#----------------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pathlib import Path
from interactome_tools import (read_single_interface_annotated_interactome,
                               num_partners,
                               num_sites)
from plot_tools import bar_plot, multi_histogram_plot

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
    referenceInteractomeFile = interactomeDir / 'human_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_site_annotated_interactome.txt'
    
    # create directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Calculate distributions of PPIs per interaction site
    #------------------------------------------------------------------------------------
    
    ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
    struc_interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    
    numPartners_ref = num_partners (ref_interactome)
    numPartners_struc = num_partners (struc_interactome)
    print()
    print('Average degree in reference interactome: %.1f' % np.mean(list(numPartners_ref.values())))
    print('Average degree in structural interactome: %.1f' % np.mean(list(numPartners_struc.values())))
    
    multi_histogram_plot (list(numPartners_ref.values()),
                          xlabel = 'Interaction degree',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 25,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'ref_interactome_degree')
    multi_histogram_plot (list(numPartners_struc.values()),
                          xlabel = 'Interaction degree',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 25,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'struc_interactome_degree')
    return
    numPPIs = {}
    for _, ppi in struc_interactome.iterrows():
        for site in ((ppi.Protein_1, ppi.Site_1), (ppi.Protein_2, ppi.Site_2)):
            numPPIs[site] = numPPIs[site] + 1 if site in numPPIs else 1
    
    maxVal = 10
    numPPIs_toplot = [min(n, maxVal) for n in numPPIs.values()]
    numPPIs_toplot = [numPPIs_toplot.count(i) for i in np.arange(1, maxVal + 1)]
    
    print()
    print('Fraction of interfacial sites mediating 1 PPI: %.1f%% (%d out of %d)' 
          % (100 * numPPIs_toplot[0] / sum(numPPIs_toplot),
             numPPIs_toplot[0],
             sum(numPPIs_toplot)))
    print('Fraction of interfacial sites mediating 2 PPIs: %.1f%% (%d out of %d)' 
          % (100 * numPPIs_toplot[1] / sum(numPPIs_toplot),
             numPPIs_toplot[1],
             sum(numPPIs_toplot)))
    print('Fraction of interfacial sites mediating ≥3 PPIs: %.1f%% (%d out of %d)' 
          % (100 * sum(numPPIs_toplot[2:]) / sum(numPPIs_toplot),
             sum(numPPIs_toplot[2:]),
             sum(numPPIs_toplot)))
    
    bar_plot (numPPIs_toplot,
              xlabels = list(map(str, np.arange(1, maxVal))) + [ '≥' + str(maxVal)],
              xlabel = 'PPIs per interfacial site',
              ylabel = 'Frequency',
              colors = 'blue',
              fontsize = 16,
              show = showFigs,
              figdir = figDir,
              figname = 'PPIs_per_site')
    
    maxVal = 10
    numSites = num_sites (struc_interactome)
    numSites_toplot = [min(n, maxVal) for n in numSites]
    numSites_toplot = [numSites_toplot.count(i) for i in np.arange(1, maxVal + 1)]
    
    print()
    print('Fraction of proteins with 1 PPI site: %.1f%% (%d out of %d)' 
          % (100 * numSites_toplot[0] / sum(numSites_toplot),
             numSites_toplot[0],
             sum(numSites_toplot)))
    print('Fraction of proteins with 2 PPI sites: %.1f%% (%d out of %d)' 
          % (100 * numSites_toplot[1] / sum(numSites_toplot),
             numSites_toplot[1],
             sum(numSites_toplot)))
    print('Fraction of proteins with ≥3 PPI sites: %.1f%% (%d out of %d)' 
          % (100 * sum(numSites_toplot[2:]) / sum(numSites_toplot),
             sum(numSites_toplot[2:]),
             sum(numSites_toplot)))
    
    bar_plot (numSites_toplot,
              xlabels = list(map(str, np.arange(1, maxVal))) + [ '≥' + str(maxVal)],
              xlabel = 'PPI sites per protein',
              ylabel = 'Frequency',
              colors = 'blue',
              fontsize = 16,
              show = showFigs,
              figdir = figDir,
              figname = 'sites_per_protein')
    
if __name__ == "__main__":
    main()
