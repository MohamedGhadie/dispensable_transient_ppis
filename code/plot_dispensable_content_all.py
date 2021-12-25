#----------------------------------------------------------------------------------------
# Plot dispensable content among multiple groups of PPIs.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from plot_tools import multi_bar_plot

def main():
    
    # reference interactome names
    interactome_names = ['HuRI', 'IntAct', 'Sahni']
    
    # structural interactome names for plot labels
    struc_interactome_names = ['Y2H-SI', 'Lit-SI', 'Experiment']
    
    # options: MutExcl, SingleInterface, Transient, Weak, Unbalanced
    tp = 'Weak'
    
    # gene expression database name
    # options: Illumina, Fantom5, GEO
    expr_db = 'Fantom5'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    groups = {'MutExcl':        ['Simultaneously possible PPIs', 'Moderately exclusive PPIs', 'Highly exclusive PPIs'],
              'SingleInterface':['Multi-interface proteins', 'Single-interface proteins'],
              'Transient':      ['Permanent PPIs', 'Transient PPIs'],
              'Weak':           ['All PPIs', 'Strong PPIs', 'Weak PPIs'],
              'Unbalanced':     ['Balanced PPIs', 'Unbalanced PPIs']}
    
    groupColors = {'MutExcl':           ['mediumseagreen', 'orange', 'red'],
                   'SingleInterface':   ['mediumseagreen', 'red'],
                   'Transient':         ['mediumseagreen', 'red'],
                   'Weak':              ['orange', 'mediumseagreen', 'red'],
                   'Unbalanced':        ['mediumseagreen', 'red']}
    
    groupHatches = {'MutExcl':          ['..', 'o', '/'],
                    'SingleInterface':  ['..', '/'],
                    'Transient':        ['..', '/'],
                    'Weak':             ['..', 'o', '/'],
                    'Unbalanced':       ['..', '/']}
    
    figName = {'MutExcl':           'dispensable_content_MutExcl',
               'SingleInterface':   'dispensable_content_SingleInterface',
               'Transient':         'dispensable_content_Transient_%s' % expr_db,
               'Weak':              'dispensable_content_Weak',
               'Unbalanced':        'dispensable_content_Unbalanced_%s' % expr_db}
    
    barWidth = {'MutExcl':          0.2,
                'SingleInterface':  0.3,
                'Transient':        0.3,
                'Weak':             0.2,
                'Unbalanced':       0.3}
    
    barGap = 0.05
    
    # plot confidence interval for the fraction of dispensable PPIs
    plotConfidenceIntervals = True
    
    # show figure legend
    showLeg = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
        
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # input data files
    inFile = {'MutExcl':        'dispensable_content_MutExcl.pkl',
              'SingleInterface':'dispensable_content_SingleInterface.pkl',
              'Transient':      'dispensable_content_Transient_%s.pkl' % expr_db,
              'Weak':           'dispensable_content_Weak.pkl',
              'Unbalanced':     'dispensable_content_Unbalanced_%s.pkl' % expr_db}
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    allresults = {}
    for name in interactome_names[:-1]:
        inPath = procDir / name / 'physics' / (ddg_method + '_edgetics') / inFile[tp]
        with open(inPath, 'rb') as f:
            allresults[name] = pickle.load(f)
    inPath = procDir / 'Sahni' / 'dispensable_content.pkl'
    with open(inPath, 'rb') as f:
        allresults['experiment'] = pickle.load(f)
    
    alldisp, allconf = [], []
    for g in groups[tp]:
        disp, conf = [], []
        for name in interactome_names:
            results = allresults[name]
            if 'DC' in results:
                disp.append(results['DC'][g] if g in results['DC'] else np.nan)
                conf.append(results['CI'][g] if g in results['CI'] else (np.nan, np.nan))
            else:
                disp.append(np.nan)
                conf.append((np.nan, np.nan))
        alldisp.append(disp)
        allconf.append(conf)
        
    numInteractomes = len(interactome_names)
    if plotConfidenceIntervals:
        flatDisp = [p for ls in alldisp for p in ls]
        flatConf = [c for ls in allconf for c in ls]
        upper = [p + upper for p, (lower, upper) in zip(flatDisp, flatConf)]
        maxY = max(upper)
    else:
        maxY = max(flatDisp)
    maxY = 10 * np.ceil(maxY / 10)
    maxY = max(maxY, 30)
    
    multi_bar_plot (alldisp,
                    errors = allconf if plotConfidenceIntervals else None,
                    xlabels = struc_interactome_names,
                    ylabels = list(np.arange(0, maxY + 5, 10)),
                    xlabel = None,
                    ylabel = 'Fraction of\ndispensable PPIs (%)',
                    colors = groupColors[tp],
                    hatches = groupHatches[tp],
                    barwidth = barWidth[tp],
                    linewidth = 2,
                    bargap = barGap,
                    capsize = 10 if plotConfidenceIntervals else 0,
                    fmt = '.k',
                    msize = 22,
                    ewidth = 2,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 24,
                    opacity = None,
                    #xlim = [0.8, numInteractomes + 0.1],
                    ylim = [0, maxY],
                    #xticks = list(np.arange(1, numInteractomes + 1)),
                    yMinorTicks = 4,
                    bottomPos = None,
                    #adjustBottom = 0.2,
                    #shiftBottomAxis = -0.1,
                    #xbounds = (1, numInteractomes),
                    overlap = False,
                    leg = groups[tp] if showLeg else None,
                    show = showFigs,
                    figdir = figDir,
                    figname = figName[tp])

if __name__ == "__main__":
    main()
