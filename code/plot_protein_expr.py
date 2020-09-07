#----------------------------------------------------------------------------------------
# Plot distribution of protein expression.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              produce_fantom5_expr_dict)
from plot_tools import multi_histogram_plot

def main():
    
    # gene expression database name
    # options: Illumina, GTEx, Fantom5
    expr_db = 'Fantom5'
    
    normalize = True
    
    logScale = False
    
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
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    
    # output data files
    if (expr_db is 'GTEx') or normalize:
        proteinExprFile = procDir / ('protein_expr_norm_%s.pkl' % expr_db)
    else:
        proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing protein tissue expression dictionary')
        if expr_db is 'Illumina':
            produce_illumina_expr_dict (illuminaExprFile,
                                        uniprotIDmapFile,
                                        proteinExprFile,
                                        normalize = normalize,
                                        headers = list(range(1, 18)))
        elif expr_db is 'GTEx':
            produce_gtex_expr_dict (gtexDir,
                                    uniprotIDmapFile,
                                    proteinExprFile,
                                    uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'Fantom5':
            produce_fantom5_expr_dict (fantomExprFile,
                                       uniprotIDmapFile,
                                       proteinExprFile,
                                       normalize = normalize,
                                       sampleTypes = 'tissues',
                                       sampleTypeFile = fantomSampleTypeFile,
                                       uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    values = list(expr.values())
    allexpr = []
    for v in values:
        allexpr.extend(v)
    
#     allexpr = []
#     for e in expr.values():
#         allexpr.append(e[3])
    
    #print(values[12000])
    #allexpr = [e for e in values[12000] if not np.isnan(e)]
#     allexpr = [e[0] for e in values if not np.isnan(e[0])]
    allexpr = [e for e in allexpr if not np.isnan(e)]
    if (expr_db is not 'GTEx') and (not normalize) and logScale:
        allexpr = [np.log10(e) for e in allexpr if e > 0]
    
    print()
    print('All expression data:')
    print('Tissues = %d' % len(values[0]))
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
        
    if (expr_db is 'GTEx') or normalize:
        figname = 'protein_expr_histogram_norm'
    elif logScale:
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
