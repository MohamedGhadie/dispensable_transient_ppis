#----------------------------------------------------------------------------------------
# This script transfers mutation ∆∆G values from one file to another.
# Values are transfered only for mutations with no ∆∆G values in the destination file. 
#
# Requirements:
# Files must be in format produced by script produce_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

from pathlib import Path
from energy_tools import copy_mutation_ddg

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    fromPath = procDir / 'HuRI' / 'dis_mut_binding_ddg_foldx.txt'
    toPath = procDir / 'HuRI' / 'dis_mut_binding_ddg_foldx_new.txt'
    outPath = procDir / 'HuRI' / 'dis_mut_binding_ddg_foldx_new_2.txt'
    
    copy_mutation_ddg (fromPath, toPath, outPath, 'binding')

if __name__ == "__main__":
    main()
