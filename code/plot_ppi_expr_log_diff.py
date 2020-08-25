#----------------------------------------------------------------------------------------
# Plot distribution of PPI log difference in expression.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              produce_fantom5_expr_dict,
                              expr_log_diff)
from plot_tools import multi_histogram_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # gene expression database name
    # options: Illumina, GTEx, Fantom5
    expr_db = 'GTEx'
    
    # minimum number of expression point values required for protein pair 
    # expression ratios to be considered
    minPoints = 1
    
    logBase = 10
    
    zeroRange = 1e-15
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
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
                                       sampleTypes = 'tissues',
                                       sampleTypeFile = fantomSampleTypeFile,
                                       uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    if expr_db is 'GTEx':
        logBase = None
    interactome = pd.read_table (interactomeFile, sep='\t')
    interactome ["expr_log_diff"] = interactome.apply (lambda x: 
                                                       expr_log_diff (x["Protein_1"],
                                                                      x["Protein_2"],
                                                                      expr,
                                                                      minPts = minPoints,
                                                                      logBase = logBase,
                                                                      method = 'mean'), axis=1)
    
    interactome = interactome [np.isnan(interactome ["expr_log_diff"]) == False]
    
    values = interactome ["expr_log_diff"].values
    print()
    print('All values:')
    print('PPIs = %d' % len(values))
    print('Mean = %g' % np.mean(values))
    print('Median = %g' % np.median(values))
    print('SD = %g' % np.std(values))
    print('Range = (%g, %g)' % (min(values), max(values)))
    
    if zeroRange < np.inf:
        values = [v for v in values if -zeroRange < v < zeroRange]
        print()
        print('Values within a range of %g from zero:' % zeroRange)
        print('PPIs = %d' % len(values))
        print('Mean = %g' % np.mean(values))
        print('Median = %g' % np.median(values))
        print('SD = %g' % np.std(values))
        print('Range = (%g, %g)' % (min(values), max(values)))
    
    if expr_db is not 'GTEx':
        figname = 'ppi_expr_log_diff_histogram'
    else:
        figname = 'ppi_expr_diff_histogram'
    
    if zeroRange < np.inf:
        figname = figname + '_range'
    
    multi_histogram_plot (values,
                          xlabel = 'Mean difference in expression',
                          ylabel = 'Number of PPIs',
                          bins = 200,
                          fontsize = 16,
                          show = showFigs,
                          figdir = figDir,
                          figname = figname + '_' + expr_db)

if __name__ == "__main__":
    main()
