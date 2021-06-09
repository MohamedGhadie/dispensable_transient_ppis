#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from protein_function import produce_illumina_expr_dict, produce_fantom5_expr_dict

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # gene expression database name
    # options: Illumina, Fantom5
    expr_db = 'Illumina'
    
    logBase = 10
    
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
    proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    exprDiffOutFile = interactomeDir / ('ppi_%s_expr_log_diff.txt' % expr_db)
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing %s protein tissue expression dictionary' % expr_db)
        if expr_db is 'Illumina':
            produce_illumina_expr_dict (illuminaExprFile,
                                        uniprotIDmapFile,
                                        proteinExprFile)
        elif expr_db is 'Fantom5':
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
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    headers = expr['Tissues']
    ppis = interactome[["Protein_1", "Protein_2"]].values
    
    for i, h in enumerate(headers):
        diff = []
        for p1, p2 in ppis:
            d = np.nan
            if (p1 in expr) and (p2 in expr):
                e1, e2 = expr[p1][i], expr[p2][i]
                if not (np.isnan(e1) or np.isnan(e2)):
                    if (e1 > 0) and (e2 > 0):
                        e1 = np.log10(e1) / np.log10(logBase)
                        e2 = np.log10(e2) / np.log10(logBase)
                        d = e1 - e2
            diff.append(d)
        interactome[h] = diff
    
    interactome.to_csv (exprDiffOutFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
