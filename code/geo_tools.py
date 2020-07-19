#----------------------------------------------------------------------------------------
# Modules for processing GEO data.
#----------------------------------------------------------------------------------------

import os
import io
from pathlib import Path
from subprocess import call
from GEOparse import GEOparse

def read_gds_softfile (gdsID, inDir = './'):
    
    gdsID = str(gdsID)
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

# def produce_geo_dataset_summary (inPath, gdsDir, outPath):
#     
#     with open(inPath, 'r') as fin:
#         gdsIDs = list(fin.read().split())
#     
#     types = {id:[] for id in gdsIDs}
#     timePoints = {id:[] for id in gdsIDs}
#     numFeatures, numSamples, numSubsets, timeMeasure = {}, {}, {}, {}
#     
#     for id in gdsIDs:
#         gds = read_gds_softfile (id, inDir = gdsDir)
#         if gds:
#             for subset in gds.subsets.values():
#                 types[id].extend(subset.metadata['type'])
#             if 'time' in types[id]:
#                 timeMeasure[id] = 'time'
#             elif 'development stage' in types[id]:
#                 timeMeasure[id] = 'development stage'
#             elif 'age' in types[id]:
#                 timeMeasure[id] = 'age'
#             else:
#                 timeMeasure[id] = '-'
#             
#             if timeMeasure[id] != '-':
#                 for subset in gds.subsets.values():
#                     if timeMeasure[id] in subset.metadata['type']:
#                         timePoints[id].extend(subset.metadata['description'])
#             
#             numFeatures[id] = gds.metadata['feature_count'][0]
#             numSamples[id] = gds.metadata['sample_count'][0]
#             numSubsets[id] = len(gds.subsets)
#     
#     types = {k:set(v) for k, v in types.items()}
#     
#     allTypes = []
#     for v in types.values():
#         allTypes.extend(list(v))
#     
#     print('\n' + 'Number of datasets: %d' % len(types))
#     print('\n' + 'Unique subset types (found in N datasets):')
#     for t in set(allTypes):
#         print('%s (%d datasets)' % (t, allTypes.count(t)))
#     
#     with io.open(outPath, "w") as fout:
#         fout.write ('\t'.join(['GDS_ID',
#                                'Feature_count',
#                                'Sample_count',
#                                'Subset_count',
#                                'Time_point_count',
#                                'Subset_types',
#                                'Time_measure',
#                                'Time_points']) + '\n')
#         for id, t in types.items():
#             fout.write ('\t'.join(['GDS' + id,
#                                    str(numFeatures[id]),
#                                    str(numSamples[id]),
#                                    str(numSubsets[id]),
#                                    str(len(timePoints[id])),
#                                    '; '.join(t),
#                                    timeMeasure[id],
#                                    '; '.join(timePoints[id])]) + '\n')

def produce_geo_dataset_summary (inPath, gdsDir, outPath):
    
    with open(inPath, 'r') as fin:
        gdsIDs = list(fin.read().split())
    
    types = {id:{} for id in gdsIDs}
    #timePoints = {id:[] for id in gdsIDs}
    numFeatures, numSamples, numSubsets, timeMeasure = {}, {}, {}, {}
    
    for id in gdsIDs:
        gds = read_gds_softfile (id, inDir = gdsDir)
        if gds:
            numFeatures[id] = gds.metadata['feature_count'][0]
            numSamples[id] = gds.metadata['sample_count'][0]
            numSubsets[id] = len(gds.subsets)
            
            for subset in gds.subsets.values():
                t = subset.metadata['type'][0]
                if t not in types[id]:
                    types[id][t] = subset.metadata['description']
                else:
                    types[id][t].extend(subset.metadata['description'])
                #types[id].extend(subset.metadata['type'])
            if 'time' in types[id]:
                timeMeasure[id] = 'time'
            elif 'development stage' in types[id]:
                timeMeasure[id] = 'development stage'
            elif 'age' in types[id]:
                timeMeasure[id] = 'age'
            else:
                timeMeasure[id] = '-'
            
#             if timeMeasure[id] != '-':
#                 for subset in gds.subsets.values():
#                     if timeMeasure[id] in subset.metadata['type']:
#                         timePoints[id].extend(subset.metadata['description'])
    
    #types = {k:set(v) for k, v in types.items()}
    
    allTypes = []
    for v in types.values():
        allTypes.extend([t for t in v.keys()])
    
    for d, tp in types.items():
        print('\n' + 'Dataset %s' % d)
        for t, v in tp.items():
            print(t + ' = ' + ', '.join(v))
    print('\n' + 'Number of datasets: %d' % len(types))
    print('Unique subset types (found in N datasets):')
    for t in set(allTypes):
        print('%s (%d datasets)' % (t, allTypes.count(t)))
    
    with io.open(outPath, "w") as fout:
        fout.write ('\t'.join(['GDS_ID',
                               'Feature_count',
                               'Sample_count',
                               'Subset_count',
                               'Time_point_count',
                               'Subset_types',
                               'Subset_values',
                               'Time_measure',
                               'Time_points']) + '\n')
        for id, tp in types.items():
            timePoints = types[id][timeMeasure[id]] if timeMeasure[id] in types[id] else []
            subsetValues = [t + ':' + ';'.join(v) for t, v in tp.items()]
            fout.write ('\t'.join(['GDS' + id,
                                   str(numFeatures[id]),
                                   str(numSamples[id]),
                                   str(numSubsets[id]),
                                   str(len(timePoints)),
                                   '; '.join(list(tp.keys())),
                                   '/'.join(subsetValues),
                                   timeMeasure[id],
                                   '; '.join(timePoints)]) + '\n')
