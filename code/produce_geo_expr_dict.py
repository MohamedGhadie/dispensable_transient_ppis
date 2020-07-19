#----------------------------------------------------------------------------------------
# Produce a summary of GEO datasets.
#----------------------------------------------------------------------------------------

import pickle
import numpy as np
from pathlib import Path
from protein_function import produce_geo_expr_dict

def main():
        
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of all GEO files
    geoDir = extDir / 'GEO'
    
    # directory of GEO dataset files
    gdsDir = geoDir / 'datasets'
    
    # input data files
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    
    # output data files
    exprFile = procDir / 'protein_expr_GEO.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # 
    #------------------------------------------------------------------------------------
    
    produce_geo_expr_dict (gdsTypeFile,
                           uniprotIDmapFile,
                           gdsDir,
                           exprFile,
                           numPoints = 5,
                           avg = 'all')
    
    with open(exprFile, 'rb') as f:
        expr = pickle.load(f)
    
    num = {k:0 for k in np.arange(1, 64)}
    for p, v in expr.items():
        num[len(v.keys())] += 1
    
    for p, v in num.items():
        print('%d: %d' % (p, v))
    
    print(len(expr.keys()))
    
if __name__ == "__main__":
    main()
