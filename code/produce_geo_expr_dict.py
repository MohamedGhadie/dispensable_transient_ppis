#----------------------------------------------------------------------------------------
# Produce a summary of GEO datasets.
#----------------------------------------------------------------------------------------

import pickle
import numpy as np
from collections import Counter
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
    
    if not exprFile.is_file():
        produce_geo_expr_dict (gdsTypeFile,
                               uniprotIDmapFile,
                               gdsDir,
                               exprFile,
                               numPoints = 5,
                               avg = 'all')
    
    with open(exprFile, 'rb') as f:
        expr = pickle.load(f)
    
#    gdsSummary = pd.read_table(gdsTypeFile, sep='\t')
    
#    numDS = {k:0 for k in np.arange(1, len(gdsSummary))}
    gdsIDs = set()
    numDS = []
    for p, v in expr.items():
        ids = list(v.keys())
        numDS.append(len(ids))
        gdsIDs.update(ids)
    numDS = Counter(numDS)
    
    print('Proteins: %d' % len(expr.keys()))
    print('Datasets: %d' % len(gdsIDs))
    print()
    print('Number of datasets: number of proteins')
    for n in sorted(numDS.keys()):
        print('%d: %d' % (n, numDS[n]))

if __name__ == "__main__":
    main()
