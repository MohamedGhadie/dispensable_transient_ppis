#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from protein_function import (produce_illumina_expr_dict,
                              produce_fantom5_expr_dict,
                              coexpr)

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # gene expression database name
    # options: Illumina, Fantom5
    expr_db = ['Illumina', 'Fantom5']
    
    # minimum number of expression point values required for protein pair tissue
    # co-expression to be considered
    minPoints = 5
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    coexprOutFile = interactomeDir / 'ppi_tissue_coexpression.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    for db in expr_db:
        proteinExprFile = procDir / ('protein_expr_%s.pkl' % db)
        if not proteinExprFile.is_file():
            print('\n' + 'producing %s protein tissue expression dictionary' % db)
            if db is 'Illumina':
                produce_illumina_expr_dict (illuminaExprFile,
                                            uniprotIDmapFile,
                                            proteinExprFile)
            elif db is 'Fantom5':
                produce_fantom5_expr_dict (fantomExprFile,
                                           uniprotIDmapFile,
                                           proteinExprFile,
                                           sampleTypes = 'tissues',
                                           sampleTypeFile = fantomSampleTypeFile,
                                           uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    #------------------------------------------------------------------------------------
    # Calculate tissue co-expression for all PPIs
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    
    for db in expr_db:
        with open(procDir / ('protein_expr_%s.pkl' % db), 'rb') as f:
            expr = pickle.load(f)
    
        interactome ['Coexpr_' + db] = interactome.apply (lambda x: 
                                                          coexpr (x["Protein_1"],
                                                                  x["Protein_2"],
                                                                  expr,
                                                                  minPts = minPoints), axis=1)
    
    interactome.to_csv (coexprOutFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
