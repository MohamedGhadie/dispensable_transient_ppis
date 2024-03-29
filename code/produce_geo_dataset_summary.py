#----------------------------------------------------------------------------------------
# Produce a summary of GEO datasets.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from geo_tools import produce_geo_dataset_summary

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
    gdsIDfile = geoDir / 'gds_accessions.txt'
    
    # output data files
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # 
    #------------------------------------------------------------------------------------
    
    produce_geo_dataset_summary (gdsIDfile, gdsDir, gdsTypeFile)

if __name__ == "__main__":
    main()
