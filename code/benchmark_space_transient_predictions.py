#----------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from stat_tools import fisher_test
from protein_function import (produce_illumina_expr_dict,
                              produce_fantom5_expr_dict,
                              coexpr,
                              expr_log_diff)

def main():
    
    # minimum number of expression point values required for protein pair tissue
    # co-expression to be considered
    minCoexprPoints = 5
    
    # minimum number of expression point values required for protein pair 
    # expression ratios to be considered
    minRatioPoints = 1
    
    # log base used for expression ratio
    logBase = 10
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # input data files
    ppiTypeFile = extDir / 'CRG.integrated.human.interactome.txt'
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    
    # output data files
    illuminaDictFile = procDir / 'protein_expr_Illumina.pkl'
    fantomDictFile = procDir / 'protein_expr_Fantom5.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # Produce expression dictionaries
    #------------------------------------------------------------------------------------
    
    if not illuminaDictFile.is_file():
        produce_illumina_expr_dict (illuminaExprFile,
                                    uniprotIDmapFile,
                                    illuminaDictFile)
    
    if not fantomDictFile.is_file():
        produce_fantom5_expr_dict (fantomExprFile,
                                   uniprotIDmapFile,
                                   fantomDictFile,
                                   sampleTypes = 'tissues',
                                   sampleTypeFile = fantomSampleTypeFile,
                                   uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    #------------------------------------------------------------------------------------
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (ppiTypeFile, sep='\t')
    
    total = len(interactome.columns) - 6
    t = 0.9 * total
    
    types = []
    for _,row in interactome.iterrows():
        expr = [e for e in row[6:] if not np.isnan(e)]
        if len(expr) == total:
            types.append('permanent' if sum(expr) > t else 'transient')
        else:
            types.append('-')
    interactome['Type'] = types
    interactome = interactome[interactome["Type"] != '-']
    
    fracTran = sum(interactome["Type"] == 'transient') / len(interactome)
    fracPerm = sum(interactome["Type"] == 'permanent') / len(interactome)
    
    print()
    print('PPIs before mapping to proteins = %d' % len(interactome))
    print('Transient PPIs = %.1f%%' % (100 * fracTran))
    print('Permanent PPIs = %.1f%%' % (100 * fracPerm))
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    interactome['Protein_1'] = interactome['Gene1'].apply(lambda x: uniprotID[x] if x in uniprotID else '-')
    interactome['Protein_2'] = interactome['Gene2'].apply(lambda x: uniprotID[x] if x in uniprotID else '-')
    interactome = interactome[(interactome["Protein_1"] != '-') & (interactome["Protein_2"] != '-')].reset_index(drop=True)
    
    fracTran = sum(interactome["Type"] == 'transient') / len(interactome)
    fracPerm = sum(interactome["Type"] == 'permanent') / len(interactome)
    
    print()
    print('PPIs after mapping to proteins = %d' % len(interactome))
    print('Transient PPIs = %.1f%%' % (100 * fracTran))
    print('Permanent PPIs = %.1f%%' % (100 * fracPerm))
        
    #------------------------------------------------------------------------------------
    # Label transient in space PPIs using Illumina data
    #------------------------------------------------------------------------------------
    
    with open(illuminaDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Coexpr_(Illumina)"] = interactome.apply (lambda x: 
                                                    coexpr (x["Protein_1"],
                                                            x["Protein_2"],
                                                            expr,
                                                            minPts = minCoexprPoints), axis=1)
    
    
    ppis = interactome[np.isnan(interactome["Coexpr_(Illumina)"]) == False]
    maxCoexpr = ppis["Coexpr_(Illumina)"].quantile(fracTran)
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Coexpr_(Illumina)"] < maxCoexpr
    
    print()
    print('********************************************************************')
    print('Transient and permanent in space (Illumina):')
    print('********************************************************************')
    print()
    print('Transient in space: %d' % sum(prdTrans))
    print('Permanent in space: %d' % sum(prdTrans == False))  
    
    TP = sum(prdTrans & expTrans)
    FP = sum(prdTrans & (expTrans == False))
    TN = sum((prdTrans == False) & (expTrans == False))
    FN = sum((prdTrans == False) & expTrans)
    
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    TNR = 1 - FPR
    Prec = 1 - FDR
    Acc = (TP + TN) / (TP + TN + FP + FN)
    bAcc = (TPR + TNR) / 2
    
    print()
    print('Acc = %.2f' % Acc)
    print('bAcc = %.2f' % bAcc)
    print('Prec = %.2f' % Prec)
    print('TPR = %.2f' % TPR)
    print('FPR = %.2f' % FPR)
    
    fisher_test([TP, FN], [FP, TN])
    
    #------------------------------------------------------------------------------------
    # Label transient in space PPIs using Fantom5 data
    #------------------------------------------------------------------------------------
    
    with open(fantomDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Coexpr_(Fantom5)"] = interactome.apply (lambda x: 
                                                    coexpr (x["Protein_1"],
                                                            x["Protein_2"],
                                                            expr,
                                                            minPts = minCoexprPoints), axis=1)
    
    ppis = interactome[np.isnan(interactome["Coexpr_(Fantom5)"]) == False]
    maxCoexpr = ppis["Coexpr_(Fantom5)"].quantile(fracTran)
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Coexpr_(Fantom5)"] < maxCoexpr
    
    print()
    print('********************************************************************')
    print('Transient and permanent in space (Fantom5):')
    print('********************************************************************')
    print()
    print('Transient in space: %d' % sum(prdTrans))
    print('Permanent in space: %d' % sum(prdTrans == False))  
    
    TP = sum(prdTrans & expTrans)
    FP = sum(prdTrans & (expTrans == False))
    TN = sum((prdTrans == False) & (expTrans == False))
    FN = sum((prdTrans == False) & expTrans)
    
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    TNR = 1 - FPR
    Prec = 1 - FDR
    Acc = (TP + TN) / (TP + TN + FP + FN)
    bAcc = (TPR + TNR) / 2
    
    print()
    print('Acc = %.2f' % Acc)
    print('bAcc = %.2f' % bAcc)
    print('Prec = %.2f' % Prec)
    print('TPR = %.2f' % TPR)
    print('FPR = %.2f' % FPR)
    
    fisher_test([TP, FN], [FP, TN])
    
    #------------------------------------------------------------------------------------
    # Label unbalanced over space PPIs using Illumina data
    #------------------------------------------------------------------------------------
    
    with open(illuminaDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["ExprDiff_(Illumina)"] = interactome.apply (lambda x: 
                                                            expr_log_diff (x["Protein_1"],
                                                                           x["Protein_2"],
                                                                           expr,
                                                                           minPts = minRatioPoints,
                                                                           logBase = logBase), axis=1)
    
    ppis = interactome[np.isnan(interactome["ExprDiff_(Illumina)"]) == False]
    minDiff = ppis["ExprDiff_(Illumina)"].quantile(1 - fracTran)
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["ExprDiff_(Illumina)"] >= minDiff
    
    print()
    print('********************************************************************')
    print('Unbalanced and balanced over space (Illumina):')
    print('********************************************************************')
    print()
    print('Unbalanced over space: %d' % sum(prdTrans))
    print('Balanced over space: %d' % sum(prdTrans == False))  
    
    TP = sum(prdTrans & expTrans)
    FP = sum(prdTrans & (expTrans == False))
    TN = sum((prdTrans == False) & (expTrans == False))
    FN = sum((prdTrans == False) & expTrans)
    
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    TNR = 1 - FPR
    Prec = 1 - FDR
    Acc = (TP + TN) / (TP + TN + FP + FN)
    bAcc = (TPR + TNR) / 2
    
    print()
    print('Acc = %.2f' % Acc)
    print('bAcc = %.2f' % bAcc)
    print('Prec = %.2f' % Prec)
    print('TPR = %.2f' % TPR)
    print('FPR = %.2f' % FPR)
    
    fisher_test([TP, FN], [FP, TN])
    
    #------------------------------------------------------------------------------------
    # Label unbalanced over space PPIs using Fantom5 data
    #------------------------------------------------------------------------------------
    
    with open(fantomDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["ExprDiff_(Fantom5)"] = interactome.apply (lambda x: 
                                                            expr_log_diff (x["Protein_1"],
                                                                           x["Protein_2"],
                                                                           expr,
                                                                           minPts = minRatioPoints,
                                                                           logBase = logBase), axis=1)
    
    ppis = interactome[np.isnan(interactome["ExprDiff_(Fantom5)"]) == False]
    minDiff = ppis["ExprDiff_(Fantom5)"].quantile(1 - fracTran)
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["ExprDiff_(Fantom5)"] >= minDiff
    
    print()
    print('********************************************************************')
    print('Unbalanced and balanced over space (Fantom5):')
    print('********************************************************************')
    print()
    print('Unbalanced over space: %d' % sum(prdTrans))
    print('Balanced over space: %d' % sum(prdTrans == False))  
    
    TP = sum(prdTrans & expTrans)
    FP = sum(prdTrans & (expTrans == False))
    TN = sum((prdTrans == False) & (expTrans == False))
    FN = sum((prdTrans == False) & expTrans)
    
    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    TNR = 1 - FPR
    Prec = 1 - FDR
    Acc = (TP + TN) / (TP + TN + FP + FN)
    bAcc = (TPR + TNR) / 2
    
    print()
    print('Acc = %.2f' % Acc)
    print('bAcc = %.2f' % bAcc)
    print('Prec = %.2f' % Prec)
    print('TPR = %.2f' % TPR)
    print('FPR = %.2f' % FPR)
    
    fisher_test([TP, FN], [FP, TN])
    
if __name__ == "__main__":
    main()
