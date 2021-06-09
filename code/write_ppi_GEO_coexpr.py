#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from protein_function import produce_geo_expr_dict, coexpr

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
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
    geoDir = extDir / 'GEO' / 'datasets'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    proteinExprFile = procDir / 'protein_expr_GEO.pkl'
    coexprOutFile = interactomeDir / 'ppi_GEO_coexpression.txt'
    
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
        dsCoexpr = []
        for p1, p2 in ppis:
            c = np.nan
            if (p1 in expr) and (p2 in expr):
                if (ds in expr[p1]) and (ds in expr[p2]):
                    c = coexpr (p1,
                                p2,
                                {p1:expr[p1][ds], p2:expr[p2][ds]},
                                minPts = minPoints)
            dsCoexpr.append(c)
        interactome [ds] = dsCoexpr
    
    interactome.to_csv (coexprOutFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
