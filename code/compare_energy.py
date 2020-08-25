#----------------------------------------------------------------------------------------
# Compare PPI interaction energy between templates and models.
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from energy_tools import read_ppi_energy
from stat_tools import pearson_corr
from plot_tools import curve_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
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
    modelEnergyFile = interactomeDir / 'ppi_energy_foldx.txt'
    templateEnergyFile = interactomeDir / 'ppi_template_energy_foldx.txt'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    modelEnergy = read_ppi_energy (modelEnergyFile)
    templateEnergy = read_ppi_energy (templateEnergyFile)
    
    modelValues, templateValues = [], []
    for p, v in modelEnergy.items():
        if p in templateEnergy:
            modelValues.append(v)
            templateValues.append(templateEnergy[p])
    
    print('Number of PPIs = %d' % len(modelValues))
    print('Average energy for models = %f' % np.mean(modelValues))
    print('Average energy for templates = %f' % np.mean(templateValues))
    print('Median energy for models = %f' % np.median(modelValues))
    print('Median energy for templates = %f' % np.median(templateValues))
    pearson_corr (modelValues, templateValues)
    
    modelValues = [y for x, y in zip(templateValues, modelValues) if x < 300]
    templateValues = [x for x in templateValues if x < 300]
    curve_plot (modelValues,
                xdata = templateValues,
                error = None,
                xlim = [-299, 299],
                ylim = [-299, 299],
                styles = '.b',
                fitstyles = '--r',
                capsize = 10,
                msize = 6,
                mwidth = 0,
                ewidth = 1,
                ecolors = 'k',   
                xlabel = 'Template energy',
                ylabel = 'Model energy',
                xticks = None,
                yticks = None,
                yMinorTicks = False,
                xticklabels = None,
                yticklabels = None,
                fontsize = 12,
                leg = None,
                compress = False,
                linefit = True,
                xstart = 0,
                binwidth = 0.1,
                perbin = 0,
                adjustBottom = False,
                shiftBottomAxis = None,
                xbounds = None,
                show = showFigs,
                figdir = figDir,
                figname = 'energy_scatter_plot')

if __name__ == "__main__":
    main()
