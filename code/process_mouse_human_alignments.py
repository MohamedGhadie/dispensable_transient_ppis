#----------------------------------------------------------------------------------------
# Process mouse-human alignments.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import parse_blast_file
from structural_annotation import locate_alignments, filter_chain_annotations

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'

    # directory of processed data files specific to orthology
    orthoDir = procDir / 'orthology'
    
    # input data files
    blastFile = orthoDir / 'mouse_human_alignment_e100'
    
    # output data files
    mapFile1 = orthoDir / 'mouse_human_alignment.txt'
    mapFile2 = orthoDir / 'mouse_human_alignment_coverage.txt'
    mapFile3 = orthoDir / 'mouse_human_map.txt'
    
    if not orthoDir.exists():
        os.makedirs(orthoDir)
    
    if not mapFile1.is_file():
        print('parsing BLAST alignment file')
        parse_blast_file (blastFile, mapFile1)
    
    if not mapFile2.is_file():
        print('filtering annotations')
        filter_chain_annotations (mapFile1, mapFile2, evalue = 100, prCov = 0, chCov = 0)
    
    if not mapFile3.is_file():
        print('locating aligned positions')
        locate_alignments (mapFile2, mapFile3)

if __name__ == "__main__":
    main()
