#----------------------------------------------------------------------------------------
# Plot distribution of PPI expression ratio.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from protein_function import produce_hpa_expr_dict, expr_ratio_symbolic
from plot_tools import bar_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # gene expression database name
    # options: HPA
    expr_db = 'HPA'
    
    # minimum number of expression point values required for protein pair 
    # expression ratios to be considered
    minPoints = 1
        
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
    hpaExprFile = extDir / 'normal_tissue.tsv'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
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
        if expr_db is 'HPA':
            produce_hpa_expr_dict (hpaExprFile,
                                   uniprotIDmapFile,
                                   proteinExprFile)
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
        
    interactome = pd.read_table (interactomeFile, sep='\t')
    interactome ["expr_ratio"] = interactome.apply (lambda x:
                                                    expr_ratio_symbolic (x["Protein_1"],
                                                                         x["Protein_2"],
                                                                         expr,
                                                                         minPts = minPoints,
                                                                         values = ['Low', 'Medium', 'High'],
                                                                         high = [{'Low', 'High'},
                                                                                 {'Low', 'Medium'},
                                                                                 {'Medium', 'High'}],
                                                                         method = 'majority'), axis=1)
    
    interactome = interactome [interactome ["expr_ratio"] != '-']
    
    values = interactome ["expr_ratio"].values
    print('PPIs = %d' % len(values))
    data = [sum(values == 'low'), sum(values == 'high')]
    bar_plot (data,
              xlabels = ['Low', 'High'],
              ylabel = 'Number of PPIs',
              ewidth = 2,
              barwidth = 0.6,
              fontsize = 16,
              show = showFigs,
              figdir = figDir,
              figname = 'ppi_expr_ratio_histogram_%s' % expr_db)

if __name__ == "__main__":
    main()
