#----------------------------------------------------------------------------------------
# This script transfers interaction energy values from one file to another.
# Values are transfered only for PPIs with no energy values in the destination file. 
#----------------------------------------------------------------------------------------

from pathlib import Path
from energy_tools import copy_ppi_energy

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    fromPath = procDir / 'HuRI' / 'ppi_template_energy_foldx.txt'
    toPath = procDir / 'IntAct' / 'ppi_template_energy_foldx.txt'
    outPath = procDir / 'IntAct' / 'ppi_template_energy_foldx_2.txt'
    copy_ppi_energy (fromPath, toPath, outPath)

if __name__ == "__main__":
    main()