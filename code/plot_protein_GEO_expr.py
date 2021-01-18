#----------------------------------------------------------------------------------------
# Plot distribution of protein expression.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from protein_function import produce_geo_expr_dict
from plot_tools import multi_histogram_plot

def main():
    
    # gene expression database name
    # options: GEO
    expr_db = 'GEO'
    
    logScale = True
    
    zeroRange = np.inf
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures')
    
    # input data files
    geoDir = extDir / 'GEO' / 'datasets'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    
    # output data files
    proteinExprFile = procDir / 'protein_expr_GEO.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing protein tissue expression dictionary')
        produce_geo_expr_dict (gdsTypeFile,
                               uniprotIDmapFile,
                               geoDir,
                               proteinExprFile,
                               numPoints = 5,
                               avg = 'all')
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    datasets = list(expr.values())
    allexpr = []
    for d in datasets:
        for v in d.values():
            allexpr.extend(v)
    
    allexpr = [e for e in allexpr if not np.isnan(e)]
    if logScale:
        allexpr = [np.log10(e) for e in allexpr if e > 0]
    
    print()
    print('All expression data:')
    print('Proteins = %d' % len(expr))
    print('All values = %d' % len(allexpr))
    print('Mean = %g' % np.mean(allexpr))
    print('Median = %g' % np.median(allexpr))
    print('SD = %g' % np.std(allexpr))
    print('Range = (%g, %g)' % (min(allexpr), max(allexpr)))
    
    if zeroRange < np.inf:
        allexpr = [e for e in allexpr if -zeroRange < e < zeroRange]
        print()
        print('Values within a range of %g from zero:' % zeroRange)
        print('All values = %d' % len(allexpr))
        print('Mean = %g' % np.mean(allexpr))
        print('Median = %g' % np.median(allexpr))
        print('SD = %g' % np.std(allexpr))
        print('Range = (%g, %g)' % (min(allexpr), max(allexpr)))
        
    if logScale:
        figname = 'protein_expr_histogram_log10'
    else:
        figname = 'protein_expr_histogram'
    
    if zeroRange < np.inf:
        figname = figname + '_range'
    
    multi_histogram_plot (allexpr,
                          xlabel = 'Expression level',
                          ylabel = 'Frequency',
                          bins = 200,
                          fontsize = 16,
                          show = showFigs,
                          figdir = figDir,
                          figname = figname + '_' + expr_db)

if __name__ == "__main__":
    main()
