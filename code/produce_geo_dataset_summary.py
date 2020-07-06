#----------------------------------------------------------------------------------------
# Produce a summary of GEO datasets.
#----------------------------------------------------------------------------------------

import os
import io
from pathlib import Path
from geo_tools import read_gds_softfile

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
    gdsTypeFile = geoDir / 'gds_subset_type.txt'
    
    # create output directories if not existing
    if not geoDir.exists():
        os.makedirs(geoDir)
    
    #------------------------------------------------------------------------------------
    # 
    #------------------------------------------------------------------------------------
    
    with open(gdsIDfile, 'r') as fin:
        gdsIDs = list(fin.read().split())
    
    numSubsets = {}
    types = {id:[] for id in gdsIDs}
    timepoints = {id:[] for id in gdsIDs}
    numFeatures, numSamples = {}, {}
    for id in gdsIDs:
        gds = read_gds_softfile (id, inDir = gdsDir)
        if gds:
            numFeatures[id] = gds.metadata['feature_count'][0]
            numSamples[id] = gds.metadata['sample_count'][0]
            numSubsets[id] = len(gds.subsets)
            for sbid, subset in gds.subsets.items():
                types[id].extend(subset.metadata['type'])
                if 'time' in subset.metadata['type']:
                    timepoints[id].extend(subset.metadata['description'])
    types = {t:set(v) for t, v in types.items()}
    
    with io.open(gdsTypeFile, "w") as fout:
        fout.write ('\t'.join(['GDS_ID',
                               'Feature_count',
                               'Sample_count',
                               'Subset_count',
                               'Time_point_count',
                               'Subset_types',
                               'Time_points']) + '\n')
        for id, tp in types.items():
            fout.write ('\t'.join([id,
                                   str(numFeatures[id]),
                                   str(numSamples[id]),
                                   str(numSubsets[id]),
                                   str(len(timepoints[id])),
                                   '; '.join(tp),
                                   '; '.join(timepoints[id])]) + '\n')
    
    allTypes = []
    for v in types.values():
        allTypes.extend(list(v))
    
    print('\n' + 'Number of datasets: %d' % len(types))
    print('\n' + 'Unique subset types (found in N datasets):')
    for t in set(allTypes):
        print('%s (%d datasets)' % (t, allTypes.count(t)))

if __name__ == "__main__":
    main()
