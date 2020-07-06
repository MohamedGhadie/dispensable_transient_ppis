#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pathlib import Path

def main():
        
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # input data files
    reGEOFile = extDir / 'upload_20181112 (original).csv'
    
    # output data files
    human_timecourse_file = procDir / 'human_timecourse.txt'
    human_timecourse_5pts_file = procDir / 'human_timecourse_5pts.txt'
    human_timecourse_10pts_file = procDir / 'human_timecourse_10pts.txt'
    human_timecourse_nondis_file = procDir / 'human_timecourse_nondis.txt'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # Load ReGEO data file
    #------------------------------------------------------------------------------------
    
    reGEO = pd.read_table (reGEOFile, sep=',')
    print('Total number of series = %d' % len(reGEO))
    
    human_timecourse = reGEO [(reGEO["timepoint_count"] > 0) & (reGEO["organism"] == 'Homo sapiens')]
    human_timecourse_5pts = human_timecourse [human_timecourse["timepoint_count"] >= 5]
    human_timecourse_10pts = human_timecourse [human_timecourse["timepoint_count"] >= 10]
    human_timecourse_nondis = human_timecourse [human_timecourse["doid_termids"].apply(type) == float]
    human_timecourse_nondis_5pts = human_timecourse_nondis [human_timecourse_nondis["timepoint_count"] >= 5]
    human_timecourse_nondis_10pts = human_timecourse_nondis [human_timecourse_nondis["timepoint_count"] >= 10]
    print('Number of series for human = %d' % sum(reGEO["organism"] == 'Homo sapiens'))
    print('Number of time-course series = %d' % sum(reGEO["timepoint_count"] > 0))
    print('Number of time-course series for human = %d' % len(human_timecourse))
    print('Number of time-course series (≥5 pts) for human = %d' % len(human_timecourse_5pts))
    print('Number of time-course series (≥10 pts) for human = %d' % len(human_timecourse_10pts))
    print('Number of time-course series for human with no disease association = %d' % len(human_timecourse_nondis))
    print('Number of time-course series (≥5 pts) for human with no disease association = %d' % len(human_timecourse_nondis_5pts))
    print('Number of time-course series (≥10 pts) for human with no disease association = %d' % len(human_timecourse_nondis_10pts))
    
    human_timecourse.to_csv (human_timecourse_file, index=False, sep='\t')
    human_timecourse_5pts.to_csv (human_timecourse_5pts_file, index=False, sep='\t')
    human_timecourse_10pts.to_csv (human_timecourse_10pts_file, index=False, sep='\t')
    human_timecourse_nondis.to_csv (human_timecourse_nondis_file, index=False, sep='\t')
    
if __name__ == "__main__":
    main()
