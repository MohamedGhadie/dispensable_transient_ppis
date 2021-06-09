#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import sys
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from simple_tools import sample_random_pairs
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              produce_hpa_expr_dict,
                              produce_fantom5_expr_dict,
                              produce_geo_expr_dict,
                              coexpr)
from plot_tools import multi_histogram_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # gene expression database name
    # options: Illumina, GTEx, HPA, Fantom5, GEO
    expr_db = 'GEO'
    
    # minimum number of expression point values required for protein pair tissue
    # co-expression to be considered
    minPoints = 5
    
    # number of random protein pairs
    numRandPairs = 100000
    
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
    
    # directory of GEO dataset files
    geoDir = extDir / 'GEO' / 'datasets'
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    hpaExprFile = extDir / 'normal_tissue.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    interactomeFile = interactomeDir / 'reference_interactome.txt'
    #interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
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
                                        proteinExprFile)
        elif expr_db is 'GTEx':
            produce_gtex_expr_dict (gtexDir,
                                    uniprotIDmapFile,
                                    proteinExprFile,
                                    uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'HPA':
            produce_hpa_expr_dict (hpaExprFile,
                                   uniprotIDmapFile,
                                   proteinExprFile)
        elif expr_db is 'Fantom5':
            produce_fantom5_expr_dict (fantomExprFile,
                                       uniprotIDmapFile,
                                       proteinExprFile,
                                       sampleTypes = 'tissues',
                                       sampleTypeFile = fantomSampleTypeFile,
                                       uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'GEO':
            produce_geo_expr_dict (gdsTypeFile,
                                   uniprotIDmapFile,
                                   geoDir,
                                   proteinExprFile,
                                   numPoints = 5,
                                   avg = 'all')
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    if expr_db is 'HPA':
        exprMap = {'Not detected':0, 'Low':1, 'Medium':2, 'High':3}
        for k, v in expr.items():
            expr[k] = np.array([(exprMap[e] if e in exprMap else np.nan) for e in v])
    
    #------------------------------------------------------------------------------------
    # Load interactome
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile, sep='\t')
    
    #-----------------------------------------------------------------------------------------
    # Generate random PPIs
    #-----------------------------------------------------------------------------------------
    
    print('Sampling random pairs...')
    proteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    randPairs = sample_random_pairs (proteins,
                                     numRandPairs,
                                     replacement = False,
                                     selfPair = False)
        
    #------------------------------------------------------------------------------------
    
    all = []
    n = len(randPairs)
    print('Calculating co-expression')
    for i, (p1, p2) in enumerate(randPairs):
        sys.stdout.write('  Pair %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        c = coexpr (p1,
                    p2,
                    expr,
                    minPts = minPoints,
                    singleExp = False if expr_db is 'GEO' else True)
        all.append(c)
    print()
    all = [c for c in all if not np.isnan(c)]
    
    print('Pairs = %d' % len(all))
    print('Mean = %g' % np.mean(all))
    print('Median = %g' % np.median(all))
    print('SD = %g' % np.std(all))
    print('Range = (%g, %g)' % (min(all), max(all)))
    
    multi_histogram_plot (all,
                          xlabel = 'Co-expression coefficient',
                          ylabel = 'Number of random pairs',
                          bins = 100,
                          fontsize = 18,
                          xlim = [-1.1, 1.1],
                          xlabels = [-1, -0.5, 0, 0.5, 1],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'rand_pair_coexpr_histogram_%s' % expr_db)
    
if __name__ == "__main__":
    main()
