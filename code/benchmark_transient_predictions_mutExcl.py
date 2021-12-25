#----------------------------------------------------------------------------------------
# Benchmark transient PPI predictions based on mutually exclusive PPI calculations
#----------------------------------------------------------------------------------------

import numpy as np
from pathlib import Path
from stat_tools import fisher_test
from interactome_tools import (read_single_interface_annotated_interactome,
                               mutExcl_simult_partners,
                               max_num_mutExcl)

def main():
    
    # reference interactome name: La_and_Mintseris
    interactome_name = 'La_and_Mintseris'
    
    # max fraction of interface overlap allowed for simultaneous PPIs
    maxSimultOverlap = 0.1
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    
    #------------------------------------------------------------------------------------
    # Count the maximum number of mutually exclusive PPIs for each PPI
    #------------------------------------------------------------------------------------
    
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    
    mutExclusive, simultaneous = mutExcl_simult_partners (interactome, maxSimultOverlap = maxSimultOverlap)
    
    interactome["Mutually_exclusive_PPIs"] = interactome.apply (lambda x: 
                                                                max_num_mutExcl (x["Protein_1"],
                                                                                 x["Protein_2"],
                                                                                 mutExclusive), axis=1)
    
    ppis = interactome[np.isnan(interactome["Mutually_exclusive_PPIs"]) == False]
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Mutually_exclusive_PPIs"] > 0
    
    print()
    print('Interactome: %s' % interactome_name)
    print('Transient PPIs in experiments = %d' % sum(expTrans))
    print('Permanent PPIs in experiments = %d' % sum(expTrans == False))
    print()
    print('Predictions:')
    print('mutually exclusive: %d' % sum(prdTrans))
    print('simultaneously possible: %d' % sum(prdTrans == False))  
    
    TP = sum(prdTrans & expTrans)
    FP = sum(prdTrans & (expTrans == False))
    TN = sum((prdTrans == False) & (expTrans == False))
    FN = sum((prdTrans == False) & expTrans)
    
    Acc = (TP + TN) / (TP + TN + FP + FN)
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    
    print()
    print('Acc = %.2f' % Acc)
    print('FDR = %.2f' % FDR)
    print('TPR = %.2f' % TPR)
    print('FPR = %.2f' % FPR)
    
    fisher_test([TP, FN], [FP, TN])
    
if __name__ == "__main__":
    main()
