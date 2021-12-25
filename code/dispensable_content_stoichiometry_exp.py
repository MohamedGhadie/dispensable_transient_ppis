#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from text_tools import read_list_table, write_list_table
from perturbation_tools import (unique_perturbation_mutations,
                                num_balanced_ppis_perturbed,
                                num_unbalanced_ppis_perturbed)
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              produce_fantom5_expr_dict,
                              produce_geo_expr_dict,
                              is_unbalanced)
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect
from plot_tools import pie_plot, curve_plot

def main():
    
    # reference interactome name
    interactome_name = 'Sahni'

    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # expression database name
    # options: Illumina, GTEx, Fantom5, GEO
    expr_db = 'GEO'
    
    # minimum number of expression point values required for protein pair 
    # expression ratios to be considered
    minPoints = 1
    
    logBase = 10
    
    if expr_db is 'GTEx':
        logBase = None
    
    # median of log difference (no log for GTEx) in expression for PPIs in structural interactome
    medianDiff = {'Illumina':0.63, 'GTEx':1.58e-17, 'Fantom5':0.68, 'GEO':0.38}
    
    # minimum difference in expression for unbalanced PPIs
    minDiff = medianDiff [expr_db]
    
    # calculate confidence interval for the fraction of dispensable PPIs
    computeConfidenceIntervals = True
    
    # % confidence interval
    CI = 95
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27
    
    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53
    
    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    pN_E_keys = ['Balanced PPIs', 'Unbalanced PPIs']
    
    pN_E, conf = {}, {}
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    geoDir = extDir / 'GEO' / 'datasets'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    
    # output data files
    if expr_db is 'GTEx':
        proteinExprFile = procDir / ('protein_expr_norm_%s.pkl' % expr_db)
    else:
        proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    natMutOutFile = interactomeDir / ('nondisease_mutation_unbalanced_perturbs_%s.txt' % expr_db)
    disMutOutFile = interactomeDir / ('disease_mutation_unbalanced_perturbs_%s.txt' % expr_db)
    dispensablePPIFile = interactomeDir / ('dispensable_content_Unbalanced_%s.pkl' % expr_db)
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Produce expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing protein tissue expression dictionary')
        if expr_db is 'Illumina':
            produce_illumina_expr_dict (illuminaExprFile,
                                        uniprotIDmapFile,
                                        proteinExprFile)
        elif expr_db is 'GTEx':
            produce_gtex_expr_dict (gtexDir,
                                    uniprotIDmapFile,
                                    proteinExprFile,
                                    uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'Fantom5':
            produce_fantom5_expr_dict (fantomExprFile,
                                       uniprotIDmapFile,
                                       proteinExprFile,
                                       sampleTypes = 'tissues',
                                       sampleTypeFile = fantomSampleTypeFile,
                                       uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'GEO':
            produce_geo_expr_dict (gdsTypeFile,
                                   uniprotIDmapFile,
                                   geoDir,
                                   proteinExprFile,
                                   numPoints = 5,
                                   avg = 'all')
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (naturalMutationsFile, ["partners", "perturbations"], [str, int])
    diseaseMutations = read_list_table (diseaseMutationsFile, ["partners", "perturbations"], [str, int])
    
    naturalMutations["Entrez_Gene_ID"] = naturalMutations["Entrez_Gene_ID"].apply(str)
    diseaseMutations["Entrez_Gene_ID"] = diseaseMutations["Entrez_Gene_ID"].apply(str)
    
    naturalMutations = naturalMutations [naturalMutations["Edgotype_class"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["Edgotype_class"] != '-'].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        naturalMutations = naturalMutations [unique_perturbation_mutations (naturalMutations)].reset_index(drop=True)
        diseaseMutations = diseaseMutations [unique_perturbation_mutations (diseaseMutations)].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Print summary of input parameters
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Summary of input parameters')
    print('********************************************************************')
    print()
    print('Interactome dataset = %s' % interactome_name)
    print('Gene expression database = %s' % expr_db)
    print('Gene expression log base = %s' % str(logBase))
    print('Minimum points required for expression difference calculation = %d' % minPoints)
    print('Expression difference cutoff for unbalanced PPIs = %f' % minDiff)
    print('Non-disease mutations = %d' % len(naturalMutations))
    print('Disease mutations = %d' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    # Count edgetic and non-edgetic mutations among all PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalMutations["Edgotype_class"] == 'Edgetic')
    numDiseaseMut_edgetic = sum(diseaseMutations["Edgotype_class"] == 'Edgetic')
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    #------------------------------------------------------------------------------------
    # Dispensable content among all PPIs
    #------------------------------------------------------------------------------------
        
    print('\n********************************************************************')
    print('Dispensable content among all PPIs:')
    print('********************************************************************\n')
    
    print('Fraction of predicted edgetic mutations:')
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic, numNaturalMut_considered),
               numNaturalMut_edgetic,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic, numDiseaseMut_considered),
               numDiseaseMut_edgetic,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    #------------------------------------------------------------------------------------
    # Label balanced and unbalanced PPIs for proteins carrying mutations
    #------------------------------------------------------------------------------------
    
    singleExp = False if expr_db is 'GEO' else True
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)

    naturalMutations ["expr_balance"] = naturalMutations.apply(
        lambda x: [is_unbalanced (uniprotID [x["Entrez_Gene_ID"]] if x["Entrez_Gene_ID"] in uniprotID else '-',
                                  uniprotID [p] if p in uniprotID else '-',
                                  expr,
                                  minPts = minPoints,
                                  logBase = logBase,
                                  minDiff = minDiff,
                                  singleExp = singleExp) for p in x["partners"]], axis=1)
    
    diseaseMutations ["expr_balance"] = diseaseMutations.apply(
        lambda x: [is_unbalanced (uniprotID [x["Entrez_Gene_ID"]] if x["Entrez_Gene_ID"] in uniprotID else '-',
                                  uniprotID [p] if p in uniprotID else '-',
                                  expr,
                                  minPts = minPoints,
                                  logBase = logBase,
                                  minDiff = minDiff,
                                  singleExp = singleExp) for p in x["partners"]], axis=1)
    
    write_list_table (naturalMutations, ["partners", "perturbations", "expr_balance"], natMutOutFile)
    write_list_table (diseaseMutations, ["partners", "perturbations", "expr_balance"], disMutOutFile)

    #------------------------------------------------------------------------------------
    # Count edgetic mutations that disrupt balanced and unbalanced PPIs
    #------------------------------------------------------------------------------------
    
    natMut_numPerturbs = naturalMutations["perturbations"].apply(np.nansum).apply(int)
    disMut_numPerturbs = diseaseMutations["perturbations"].apply(np.nansum).apply(int)
    
    natMut_numBalancedPerturbs = naturalMutations.apply(
        lambda x: num_balanced_ppis_perturbed (x["perturbations"], x["expr_balance"]), axis=1)
    disMut_numBalancedPerturbs = diseaseMutations.apply(
        lambda x: num_balanced_ppis_perturbed (x["perturbations"], x["expr_balance"]), axis=1)
    
    natMut_numUnbalancedPerturbs = naturalMutations.apply(
        lambda x: num_unbalanced_ppis_perturbed (x["perturbations"], x["expr_balance"]), axis=1)
    disMut_numUnbalancedPerturbs = diseaseMutations.apply(
        lambda x: num_unbalanced_ppis_perturbed (x["perturbations"], x["expr_balance"]), axis=1)
    
    numNaturalMut_edgetic_bal = sum((naturalMutations["Edgotype_class"] == 'Edgetic') &
                                    (natMut_numBalancedPerturbs > 0))
    numDiseaseMut_edgetic_bal = sum((diseaseMutations["Edgotype_class"] == 'Edgetic') &
                                    (disMut_numBalancedPerturbs > 0))
    
    numNaturalMut_edgetic_unbal = sum((naturalMutations["Edgotype_class"] == 'Edgetic') & 
                                      (natMut_numUnbalancedPerturbs == natMut_numPerturbs))
    numDiseaseMut_edgetic_unbal = sum((diseaseMutations["Edgotype_class"] == 'Edgetic') & 
                                      (disMut_numUnbalancedPerturbs == disMut_numPerturbs))
    
    numNaturalMut_considered = (numNaturalMut_nonedgetic + 
                                numNaturalMut_edgetic_bal + 
                                numNaturalMut_edgetic_unbal) 
                            
    numDiseaseMut_considered = (numDiseaseMut_nonedgetic + 
                                numDiseaseMut_edgetic_bal + 
                                numDiseaseMut_edgetic_unbal)
    
    #------------------------------------------------------------------------------------
    # Dispensable content among balanced PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among balanced PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted edgetic mutations:')
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_bal / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_bal, numNaturalMut_considered),
               numNaturalMut_edgetic_bal,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_bal / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_bal, numDiseaseMut_considered),
               numDiseaseMut_edgetic_bal,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_bal, numNaturalMut_considered - numNaturalMut_edgetic_bal],
                 [numDiseaseMut_edgetic_bal, numDiseaseMut_considered - numDiseaseMut_edgetic_bal])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_bal,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_bal,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Balanced PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Balanced PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among unbalanced PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among unbalanced PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted edgetic mutations:')
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_unbal / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_unbal, numNaturalMut_considered),
               numNaturalMut_edgetic_unbal,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_unbal / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_unbal, numDiseaseMut_considered),
               numDiseaseMut_edgetic_unbal,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_unbal, numNaturalMut_considered - numNaturalMut_edgetic_unbal],
                 [numDiseaseMut_edgetic_unbal, numDiseaseMut_considered - numDiseaseMut_edgetic_unbal])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_unbal,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_unbal,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Unbalanced PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Unbalanced PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Save dispensable content to file
    #------------------------------------------------------------------------------------
    
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump({'DC':pN_E, 'CI':conf}, fOut)
    
    #------------------------------------------------------------------------------------
    # Plot pie charts
    #------------------------------------------------------------------------------------
    
    pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic_unbal, numNaturalMut_edgetic_bal],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'nondisease_mutations_unbalanced_edgetic_%s' % expr_db)
    
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic_unbal, numDiseaseMut_edgetic_bal],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_unbalanced_edgetic_%s' % expr_db)
    
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    if computeConfidenceIntervals:
        maxY = max([pN_E[p] + conf[p][1] for p in pN_E.keys()])
    else:
        maxY = max(pN_E.values())
    maxY = 10 * np.ceil(maxY / 10)
    maxY = max(maxY, 30)
    
    curve_plot ([pN_E[p] for p in pN_E_keys if p in pN_E],
                error = [conf[p] for p in pN_E_keys if p in pN_E] if computeConfidenceIntervals else None,
                xlim = [0.8, 3.1],
                ylim = [0, maxY],
                styles = '.k',
                capsize = 10 if computeConfidenceIntervals else 0,
                msize = 26,
                ewidth = 2,
                ecolors = 'k',
                ylabel = 'Fraction of dispensable PPIs (%)',
                yMinorTicks = 4,
                xticks = [1, 2],
                xticklabels = [p.replace(' ', '\n') for p in pN_E_keys if p in pN_E],
                yticklabels = list(np.arange(0, maxY + 5, 10)),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 2),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_unbalanced_%s' % expr_db)

if __name__ == "__main__":
    main()
