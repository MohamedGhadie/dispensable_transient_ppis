#----------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from stat_tools import fisher_test
from energy_tools import read_ppi_energy
from interactome_tools import (read_single_interface_annotated_interactome,
                               mutExcl_simult_partners,
                               max_num_mutExcl)
from protein_function import (produce_illumina_expr_dict,
                              produce_fantom5_expr_dict,
                              produce_geo_expr_dict,
                              is_weak,
                              is_transient,
                              is_unbalanced)

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
    # type of data used for prediction: energy, coexpr_GEO, coexpr_Illumina, coexpr_Fantom5, 
    # stoich_GEO, stoich_Illumina, stoich_Fantom5, interface
    #predData = 'energy'
    
    # structure to use for PPI energy
    # options: template, model
    energyStructure = 'template'
    
    # use the median of PPI binding free energy as a cutoff instead of fixed cutoff
    medEnergy = False
    
    # median interaction free energy for all PPIs in structural interactome
    medianEnergy = {'HuRI':     {'template': -20, 'model': -14},
                    'IntAct':   {'template': -15, 'model': -9}}
    
    # minumum interaction free energy required for weak PPIs
    if medEnergy:
        minEnergy = medianEnergy[interactome_name][energyStructure]
    else:
        minEnergy = -8
    
    # minimum number of expression point values required for protein pair tissue
    # co-expression to be considered
    minCoexprPoints = 5
    
    # median of PPI co-expression distribution in structural interactome
    medianCoexpr = {'HuRI':    {'Illumina':0.39, 'GTEx':0.80, 'HPA':0.22, 'Fantom5':0.16, 'GEO':0.10},
                    'IntAct':  {'Illumina':0.45, 'GTEx':0.83, 'HPA':0.26, 'Fantom5':0.22, 'GEO':0.11}}
    
    # minimum number of expression point values required for protein pair 
    # expression ratios to be considered
    minRatioPoints = 1
    
    # log base used for expression ratio
    logBase = 10
    
    # median of log difference in expression for PPIs in structural interactome
    medianDiff = {'HuRI':    {'Illumina':0.63, 'Fantom5':0.68, 'GEO':0.38},
                  'IntAct':  {'Illumina':0.56, 'Fantom5':0.59, 'GEO':0.40}}
    
    # max fraction of interface overlap allowed for simultaneous PPIs
    maxSimultOverlap = 0.1
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed data files specific to La and Mintseris datasets
    expDir = procDir / 'La_and_Mintseris'
    
    # input data files
    ppiTypeFile = expDir / (interactome_name + '_ppi_types.txt')
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    geoDir = extDir / 'GEO' / 'datasets'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    strucInteractomeFile = interactomeDir / 'structural_interactome.txt'
    energyFile = interactomeDir / ('ppi_%s_energy_foldx.txt' % energyStructure)
    
    # output data files
    illuminaDictFile = procDir / 'protein_expr_Illumina.pkl'
    fantomDictFile = procDir / 'protein_expr_Fantom5.pkl'
    geoDictFile = procDir / 'protein_expr_GEO.pkl'
    
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
    
    if not geoDictFile.is_file():
        produce_geo_expr_dict (gdsTypeFile,
                               uniprotIDmapFile,
                               geoDir,
                               geoDictFile,
                               numPoints = 5,
                               avg = 'all')
    
    #------------------------------------------------------------------------------------
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (ppiTypeFile, sep='\t')
    
    print()
    print('Interactome: %s' % interactome_name)
    print('Transient PPIs in experiments = %d' % sum(interactome["Type"] == 'transient'))
    print('Permanent PPIs in experiments = %d' % sum(interactome["Type"] == 'permanent'))
    
    #------------------------------------------------------------------------------------
    # Label weak and strong PPIs
    #------------------------------------------------------------------------------------
    
    energy = read_ppi_energy (energyFile)
    interactome ["Strength"] = interactome.apply (lambda x: is_weak (x["Protein_1"],
                                                                     x["Protein_2"],
                                                                     energy,
                                                                     minEnergy = minEnergy), axis=1)
    
    ppis = interactome[interactome["Strength"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Strength"] == 'weak'
    
    print()
    print('********************************************************************')
    print('Weak and strong PPIs:')
    print('********************************************************************')
    print()
    print('Weak: %d' % sum(prdTrans))
    print('Strong: %d' % sum(prdTrans == False))  
    
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
    # Label transient in time PPIs using GEO data
    #------------------------------------------------------------------------------------
    
    with open(geoDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Transient_in_time"] = interactome.apply (lambda x: 
                                        is_transient (x["Protein_1"],
                                                      x["Protein_2"],
                                                      expr,
                                                      minPts = minCoexprPoints,
                                                      maxCoexpr = medianCoexpr[interactome_name]['GEO'],
                                                      singleExp = False), axis=1)
    
    ppis = interactome[interactome["Transient_in_time"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Transient_in_time"] == 'transient'
    
    print()
    print('********************************************************************')
    print('Transient and permanent in time:')
    print('********************************************************************')
    print()
    print('Transient in time: %d' % sum(prdTrans))
    print('Permanent in time: %d' % sum(prdTrans == False))  
    
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
    # Label transient in space PPIs using Illumina data
    #------------------------------------------------------------------------------------
    
    with open(illuminaDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Transient_in_space_(Illumina)"] = interactome.apply (lambda x: 
                                        is_transient (x["Protein_1"],
                                                      x["Protein_2"],
                                                      expr,
                                                      minPts = minCoexprPoints,
                                                      maxCoexpr = medianCoexpr[interactome_name]['Illumina'],
                                                      singleExp = True), axis=1)
    
    ppis = interactome[interactome["Transient_in_space_(Illumina)"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Transient_in_space_(Illumina)"] == 'transient'
    
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
    
    interactome ["Transient_in_space_(Fantom5)"] = interactome.apply (lambda x: 
                                        is_transient (x["Protein_1"],
                                                      x["Protein_2"],
                                                      expr,
                                                      minPts = minCoexprPoints,
                                                      maxCoexpr = medianCoexpr[interactome_name]['Fantom5'],
                                                      singleExp = True), axis=1)
    
    ppis = interactome[interactome["Transient_in_space_(Fantom5)"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Transient_in_space_(Fantom5)"] == 'transient'
    
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
    # Label unbalanced over time PPIs using GEO data
    #------------------------------------------------------------------------------------
    
    with open(geoDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Balance_over_time"] = interactome.apply (lambda x: 
                                        is_unbalanced (x["Protein_1"],
                                                       x["Protein_2"],
                                                       expr,
                                                       minPts = minRatioPoints,
                                                       logBase = logBase,
                                                       minDiff = medianDiff[interactome_name]['GEO'],
                                                       singleExp = False), axis=1)
    
    ppis = interactome[interactome["Balance_over_time"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Balance_over_time"] == 'unbalanced'
    
    print()
    print('********************************************************************')
    print('Unbalanced and balanced over time:')
    print('********************************************************************')
    print()
    print('Unbalanced over time: %d' % sum(prdTrans))
    print('Balanced over time: %d' % sum(prdTrans == False))  
    
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
    
    interactome ["Balance_over_space_(Illumina)"] = interactome.apply (lambda x: 
                                        is_unbalanced (x["Protein_1"],
                                                       x["Protein_2"],
                                                       expr,
                                                       minPts = minRatioPoints,
                                                       logBase = logBase,
                                                       minDiff = medianDiff[interactome_name]['Illumina'],
                                                       singleExp = True), axis=1)
    
    ppis = interactome[interactome["Balance_over_space_(Illumina)"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Balance_over_space_(Illumina)"] == 'unbalanced'
    
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
    
    interactome ["Balance_over_space_(Fantom5)"] = interactome.apply (lambda x: 
                                        is_unbalanced (x["Protein_1"],
                                                       x["Protein_2"],
                                                       expr,
                                                       minPts = minRatioPoints,
                                                       logBase = logBase,
                                                       minDiff = medianDiff[interactome_name]['Fantom5'],
                                                       singleExp = True), axis=1)
    
    ppis = interactome[interactome["Balance_over_space_(Fantom5)"] != '-']
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Balance_over_space_(Fantom5)"] == 'unbalanced'
    
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
    
    #------------------------------------------------------------------------------------
    # Count the maximum number of mutually exclusive PPIs for each PPI
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Mutually exclusive PPIs:')
    print('********************************************************************')
    print()
    
    strucInteractome = read_single_interface_annotated_interactome (strucInteractomeFile)
    mutExclusive, simultaneous = mutExcl_simult_partners (strucInteractome, maxSimultOverlap = maxSimultOverlap)
    
    interactome["Mutually_exclusive_PPIs"] = interactome.apply (lambda x: 
                                                                max_num_mutExcl (x["Protein_1"],
                                                                                 x["Protein_2"],
                                                                                 mutExclusive), axis=1)
    
    ppis = interactome[np.isnan(interactome["Mutually_exclusive_PPIs"]) == False]
    expTrans = ppis["Type"] == 'transient'
    prdTrans = ppis["Mutually_exclusive_PPIs"] > 0
    
    print()
    print('Mutually exclusive: %d' % sum(prdTrans))
    print('Simultaneously possible: %d' % sum(prdTrans == False))  
    
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
