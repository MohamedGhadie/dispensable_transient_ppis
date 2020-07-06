#----------------------------------------------------------------------------------------
# Modules for processing GEO data.
#----------------------------------------------------------------------------------------

import os
import sys
import pickle
import warnings
from pathlib import Path
from subprocess import call
from GEOparse import GEOparse

def read_gds_softfile (gdsID, inDir = './'):
    
    if not gdsID.startswith('GDS'):
        gdsID = 'GDS' + gdsID
    filepath = inDir / (gdsID + '.soft')
    if filepath.is_file():
        return GEOparse.parse_GDS (str(filepath))
    else:
        return None

def download_gds_softfiles (inPath,
                            outDir = './GEO',
                            how = 'full',
                            silent = False):

    if not outDir.exists():
        os.makedirs(outDir)
    with open(inPath, 'r') as fin:
        gdsIDs = list(fin.read().split())
    failed = 0
    for i, id in enumerate(gdsIDs):
        filename = outDir / ('GDS' + id + '.soft')
        if not filename.is_file():
            try:
                filepath, objtype = download_geo_softfile ('GDS' + id,
                                                           outDir = outDir,
                                                           how = how,
                                                           silent = silent)
                if Path(filepath).is_file():
                    call(['gunzip' , filepath])
                else:
                    failed += 1
            except:
                failed += 1
    print('\n' + '%d files requested, %d file downloads failed' % (len(gdsIDs), failed))

def download_geo_softfile (geoID,
                           outDir = './',
                           how = 'full',
                           annotate_gpl = False,
                           include_data = False,
                           silent = False,
                           aspera = False):
    try:
        return GEOparse.get_GEO_file (geo = geoID,
                                      destdir = outDir,
                                      how = how,
                                      annotate_gpl = annotate_gpl,
                                      include_data = include_data,
                                      silent = silent,
                                      aspera = aspera)
    except:
        raise
             