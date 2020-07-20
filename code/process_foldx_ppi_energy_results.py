#----------------------------------------------------------------------------------------
# Process interaction energy results from FoldX.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from energy_tools import read_foldx_results, write_ppi_energy_tofile

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of foldx results
    inDir = interactomeDir / 'foldx' / 'results'
    
    # input files
    energyFile = interactomeDir / 'ppi_template_energy_foldx.txt'
    
    # temporary output file
    tempFile = interactomeDir / 'ppi_energy_foldx_temp.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    processed = read_foldx_results (inDir, type = 'ppi_energy')
    write_ppi_energy_tofile (processed, energyFile, tempFile)
    energyFile.unlink()
    tempFile.rename (energyFile)

if __name__ == "__main__":
    main()
