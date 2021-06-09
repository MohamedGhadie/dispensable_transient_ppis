#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from protein_function import produce_geo_expr_dict

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
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
    geoDir = extDir / 'GEO' / 'datasets'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    proteinExprFile = procDir / 'protein_expr_GEO.pkl'
    exprDiffOutFile = interactomeDir / 'ppi_GEO_expr_log_diff.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    #------------------------------------------------------------------------------------
    # Produce GEO time-course expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing GEO time-course expression dictionary')
        produce_geo_expr_dict (gdsTypeFile,
                               uniprotIDmapFile,
                               geoDir,
                               proteinExprFile,
                               numPoints = 5,
                               avg = 'all')
    
    #------------------------------------------------------------------------------------
    # calculate time-course co-expression for all PPIs
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    dsIDs = set()
    for p, v in expr.items():
        dsIDs.update(list(v.keys()))
    dsIDs = sorted(dsIDs)
    
    ppis = interactome[["Protein_1", "Protein_2"]].values
    for ds in dsIDs:
        diff = []
        for p1, p2 in ppis:
            d = ''
            if (p1 in expr) and (p2 in expr):
                if (ds in expr[p1]) and (ds in expr[p2]):
                    e1, e2 = expr[p1][ds], expr[p2][ds]
                    e1 = np.array([(e if e > 0 else np.nan) for e in e1])
                    e2 = np.array([(e if e > 0 else np.nan) for e in e2])
                    e1 = np.log10(e1) / np.log10(logBase)
                    e2 = np.log10(e2) / np.log10(logBase)
                    d = e1 - e2
            diff.append(', '.join(map(str, d)))
        interactome [ds] = diff
    
    interactome.to_csv (exprDiffOutFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
