#----------------------------------------------------------------------------------------
# Plot the fraction of dispensable PPIs calculated from predictions and experiments.
#
# Run the following scripts before running this script:
# - dispensable_content_all_PPIs.py
# - dispensable_content_experiment.py
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from plot_tools import curve_plot, bar_plot

def main():
    
    # reference interactome names
    interactome_names = ['HuRI', 'IntAct', 'experiment']
    
    # structural interactome names for plot labels
    struc_interactome_names = ['Y2H-SI', 'Lit-SI', 'Experiment']
    
    # plot confidence interval for the fraction of dispensable PPIs
    plotConfidenceIntervals = True
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
        
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    dispensablePPIFile_name = 'fraction_disp_PPIs.pkl'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    pN_E_results, pN_E_bounds = [], []
    for interactome_name in interactome_names:
        dispensablePPIFile = procDir / interactome_name / dispensablePPIFile_name
        with open(dispensablePPIFile, 'rb') as f:
            all_results = pickle.load(f)
        pN_E_results.append( 100 * all_results['P(N|E)'] )
        if 'P(N|E)_CI' in all_results:
            lower, upper = all_results['P(N|E)_CI']
            pN_E_bounds.append( (100 * lower, 100 * upper) )
        else:
            pN_E_bounds.append( (0, 0) )
    
    numInteractomes = len(interactome_names)
    if plotConfidenceIntervals:
        upper = [p + upper for p, (lower, upper) in zip(pN_E_results, pN_E_bounds)]
        maxY = max(upper)
    else:
        maxY = max(pN_E_results)
    maxY = 10 * np.ceil(maxY / 10)
    maxY = max(maxY, 30)
    
    curve_plot (pN_E_results,
                error = pN_E_bounds if plotConfidenceIntervals else None,
                xlim = [0.8, numInteractomes + 0.1],
                ylim = [0, maxY],
                styles = '.k',
                capsize = 10 if plotConfidenceIntervals else 0,
                msize = 26,
                ewidth = 2,
                ecolors = 'k',
                ylabel = 'Fraction of dispensable PPIs (%)',
                yMinorTicks = 4,
                xticks = list(np.arange(1, numInteractomes + 1)),
                xticklabels = struc_interactome_names,
                yticklabels = list(np.arange(0, maxY + 5, 10)),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, numInteractomes),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs')
    
    bar_plot (pN_E_results,
                    error = pN_E_bounds if plotConfidenceIntervals else None,
                    xlabels = struc_interactome_names,
                    ylabels = list(np.arange(0, maxY + 5, 10)),
                    xlabel = None,
                    ylabel = 'Fraction of dispensable PPIs (%)',
                    colors = 'orange',
                    barwidth = 0.5,
                    capsize = 10 if plotConfidenceIntervals else 0,
                    fmt = '.k',
                    msize = 26,
                    ewidth = 2,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 20,
                    opacity = None,
                    #xlim = [0.8, numInteractomes + 0.1],
                    ylim = [0, maxY],
                    #xticks = list(np.arange(1, numInteractomes + 1)),
                    yMinorTicks = 4,
                    #adjustBottom = 0.2,
                    #shiftBottomAxis = -0.1,
                    #xbounds = (1, numInteractomes),
                    show = showFigs,
                    figdir = figDir,
                    figname = 'dispensable_content_allPPIs')

if __name__ == "__main__":
    main()
