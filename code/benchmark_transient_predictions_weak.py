#----------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

import pandas as pd
from pathlib import Path
from stat_tools import fisher_test

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'La_and_Mintseris'
        
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    interactomeFile = interactomeDir / 'ppi_energy_foldx.txt'
    
    #------------------------------------------------------------------------------------
    # Label weak and strong PPIs
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile, sep='\t')
    ppis = interactome[interactome["Interaction_energy"] != 'X']
    
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Interaction_energy"] > -15
    
    print()
    print('Interactome: %s' % interactome_name)
    print('Transient PPIs in experiments = %d' % sum(expTrans))
    print('Permanent PPIs in experiments = %d' % sum(expTrans == False))
    print()
    print('Predictions:')
    print('Weak: %d' % sum(prdTrans))
    print('Strong: %d' % sum(prdTrans == False))  
    
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