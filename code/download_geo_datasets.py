#----------------------------------------------------------------------------------------
# Download GEO datasets given a list of accession IDs.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from geo_tools import download_gds_softfiles

def main():
        
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # directory of all GEO files
    geoDir = extDir / 'GEO'
    
    # directory of GEO dataset files
    gdsDir = geoDir / 'datasets'
    
    # input data files
    gdsIDfile = geoDir / 'gds_accessions.txt'
    
    # create output directories if not existing
    if not gdsDir.exists():
        os.makedirs(gdsDir)
    
    #------------------------------------------------------------------------------------
    # Download GEO dataset soft files
    #------------------------------------------------------------------------------------
    
    download_gds_softfiles (gdsIDfile, outDir = gdsDir, how = 'full')
    
if __name__ == "__main__":
    main()
