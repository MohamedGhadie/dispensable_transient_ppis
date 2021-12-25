#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table, write_list_table
from perturbation_tools import (unique_perturbation_mutations,
                                num_permanent_ppis_perturbed,
                                num_transient_ppis_perturbed)
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect
from plot_tools import pie_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    #interactome_names = ['IntAct', 'HuRI']
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False

    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
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
    
    pN_E_keys = ['All PPIs', 'Permanent PPIs', 'Transient PPIs']
    
    pN_E, conf = {}, {}
    
    # show figures
    showFigs = False
    
    # write results to output files
    save_results = True
    
    # save figures
    save_figures = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to La and Mintseris datasets
    expDir = procDir / 'La_and_Mintseris'
    
    # figure directory
    figDir = Path('../figures') / 'La_and_Mintseris'
    
    # output data files
    natMutOutFile = expDir / (interactome_name + '_nondisease_mutation_transientPPI_perturbs.txt')
    disMutOutFile = expDir / (interactome_name + '_disease_mutation_transientPPI_perturbs.txt')
    dispensablePPIFile = expDir / (interactome_name + '_dispensable_content_Weak.pkl')
    
    # create output directories if not existing
    if not expDir.exists():
        os.makedirs(expDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Create edgotype label
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        eCol = 'mono-edgotype'
        edgotype = 'mono-edgetic'
    else:
        eCol = 'edgotype'
        edgotype = 'edgetic'
    
    #------------------------------------------------------------------------------------
    # Load PPI energy
    #------------------------------------------------------------------------------------
    
#     ppiTypes = {}
#     for name in interactome_names:
#         types = pd.read_table (expDir / (name + '_ppi_types.txt'))
#         for _, row in types.iterrows():
#             k = '-'.join(sorted([row.Protein_1, row.Protein_2]))
#             if k not in ppiTypes:
#                 ppiTypes[k] = row.Type
    
    ppiTypes = {}
    types = pd.read_table (expDir / (interactome_name + '_ppi_types.txt'))
    for _, row in types.iterrows():
        k = '-'.join(sorted([row.Protein_1, row.Protein_2]))
        if k not in ppiTypes:
            ppiTypes[k] = row.Type
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
#     naturalMutations = pd.DataFrame()
#     diseaseMutations = pd.DataFrame()
#     
#     for name in interactome_names:
#         edgeticDir = procDir / name / 'physics' / (ddg_method + '_edgetics')
#         
#         naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
#         diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
#         
#         natMut = read_list_table (naturalMutationsFile, ["partners", "perturbations"], [str, float])
#         disMut = read_list_table (diseaseMutationsFile, ["partners", "perturbations"], [str, float])
#         
#         naturalMutations = naturalMutations.append(natMut, ignore_index=True)
#         diseaseMutations = diseaseMutations.append(disMut, ignore_index=True)
#     
#     naturalMutations = naturalMutations.drop_duplicates(subset=['protein', 'mut_position', 'mut_res'])
#     diseaseMutations = diseaseMutations.drop_duplicates(subset=['protein', 'mut_position', 'mut_res'])
    
    edgeticDir = procDir / interactome_name / 'physics' / (ddg_method + '_edgetics')
    
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
    naturalMutations = read_list_table (naturalMutationsFile, ["partners", "perturbations"], [str, float])
    diseaseMutations = read_list_table (diseaseMutationsFile, ["partners", "perturbations"], [str, float])
    
    naturalMutations["perturbations"] = naturalMutations["perturbations"].apply(
                                                lambda x: [int(p) if not np.isnan(p) else p for p in x])
    diseaseMutations["perturbations"] = diseaseMutations["perturbations"].apply(
                                                lambda x: [int(p) if not np.isnan(p) else p for p in x])
    
    naturalMutations = naturalMutations [naturalMutations[eCol] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations[eCol] != '-'].reset_index(drop=True)
    
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
    print('Number of transient PPIs = %d' % sum([t == 'transient' for t in ppiTypes.values()]))
    print('Number of permanent PPIs = %d' % sum([t == 'permanent' for t in ppiTypes.values()]))
    print('Edgotype = %s' % edgotype)
    print('Non-disease mutations = %d' % len(naturalMutations))
    print('Disease mutations = %d' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    # Count edgetic and non-edgetic mutations among all PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalMutations[eCol] == edgotype)
    numDiseaseMut_edgetic = sum(diseaseMutations[eCol] == edgotype)
    
    numNaturalMut_nonedgetic = sum(naturalMutations[eCol] != edgotype)
    numDiseaseMut_nonedgetic = sum(diseaseMutations[eCol] != edgotype)
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
   #------------------------------------------------------------------------------------
    # Dispensable content among all PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among all PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
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
    
    if 'P(N|E)' in all_effects:
        pN_E['All PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['All PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Label transient and permanent PPIs for proteins carrying mutations
    #------------------------------------------------------------------------------------
    
    types = []
    for _, row in naturalMutations.iterrows():
        t = []
        for p in row.partners:
            k = '-'.join(sorted([row.protein, p]))
            t.append(ppiTypes[k] if k in ppiTypes else '-')
        types.append(t)
    naturalMutations ["transient_PPIs"] = types
    
    types = []
    for _, row in diseaseMutations.iterrows():
        t = []
        for p in row.partners:
            k = '-'.join(sorted([row.protein, p]))
            t.append(ppiTypes[k] if k in ppiTypes else '-')
        types.append(t)
    diseaseMutations ["transient_PPIs"] = types
    
    if save_results:
        write_list_table (naturalMutations, ["partners", "perturbations", "transient_PPIs"], natMutOutFile)
        write_list_table (diseaseMutations, ["partners", "perturbations", "transient_PPIs"], disMutOutFile)
    
    #------------------------------------------------------------------------------------
    # Count edgetic mutations that disrupt transient or permanent PPIs
    #------------------------------------------------------------------------------------
    
    natMut_numPerturbs = naturalMutations["perturbations"].apply(np.nansum).apply(int)
    disMut_numPerturbs = diseaseMutations["perturbations"].apply(np.nansum).apply(int)
    
    natMut_numPermPerturbs = naturalMutations.apply(
        lambda x: num_permanent_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    disMut_numPermPerturbs = diseaseMutations.apply(
        lambda x: num_permanent_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    
    natMut_numTransPerturbs = naturalMutations.apply(
        lambda x: num_transient_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    disMut_numTransPerturbs = diseaseMutations.apply(
        lambda x: num_transient_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    
    numNaturalMut_edgetic_perm = sum((naturalMutations[eCol] == edgotype) &
                                     (natMut_numPermPerturbs > 0))
    numDiseaseMut_edgetic_perm = sum((diseaseMutations[eCol] == edgotype) &
                                     (disMut_numPermPerturbs > 0))
    
    numNaturalMut_edgetic_tran = sum((naturalMutations[eCol] == edgotype) & 
                                     (natMut_numTransPerturbs == natMut_numPerturbs))
    numDiseaseMut_edgetic_tran = sum((diseaseMutations[eCol] == edgotype) & 
                                     (disMut_numTransPerturbs == disMut_numPerturbs))
    
    numNaturalMut_considered = (numNaturalMut_nonedgetic + 
                                numNaturalMut_edgetic_perm + 
                                numNaturalMut_edgetic_tran) 
                            
    numDiseaseMut_considered = (numDiseaseMut_nonedgetic + 
                                numDiseaseMut_edgetic_perm + 
                                numDiseaseMut_edgetic_tran)
    
    #------------------------------------------------------------------------------------
    # Dispensable content among permanent PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among permanent PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_perm / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_perm, numNaturalMut_considered),
               numNaturalMut_edgetic_perm,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_perm / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_perm, numDiseaseMut_considered),
               numDiseaseMut_edgetic_perm,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_perm, numNaturalMut_considered - numNaturalMut_edgetic_perm],
                 [numDiseaseMut_edgetic_perm, numDiseaseMut_considered - numDiseaseMut_edgetic_perm])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_perm,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_perm,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Permanent PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Permanent PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among transient PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among transient PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_tran / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_tran, numNaturalMut_considered),
               numNaturalMut_edgetic_tran,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_tran / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_tran, numDiseaseMut_considered),
               numDiseaseMut_edgetic_tran,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_tran, numNaturalMut_considered - numNaturalMut_edgetic_tran],
                 [numDiseaseMut_edgetic_tran, numDiseaseMut_considered - numDiseaseMut_edgetic_tran])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_tran,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_tran,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Transient PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Transient PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Save dispensable content to file
    #------------------------------------------------------------------------------------
    
    if save_results:
        with open(dispensablePPIFile, 'wb') as fOut:
            pickle.dump({'DC':pN_E, 'CI':conf}, fOut)
    
    #------------------------------------------------------------------------------------
    # Plot pie charts
    #------------------------------------------------------------------------------------
    
    if save_figures:
        pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic_tran, numNaturalMut_edgetic_perm],
                  angle = 90,
                  colors = ['lightsteelblue', 'red', 'mediumseagreen'],
                  edgewidth = 2,
                  show = showFigs,
                  figdir = figDir,
                  figname = '%s_nondisease_mutations_transient_%s' % (interactome_name, edgotype))
    
        pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic_tran, numDiseaseMut_edgetic_perm],
                  angle = 90,
                  colors = ['lightsteelblue', 'red', 'mediumseagreen'],
                  edgewidth = 2,
                  show = showFigs,
                  figdir = figDir,
                  figname = '%s_disease_mutations_transient_%s' % (interactome_name, edgotype))
    
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    if save_figures:
        if computeConfidenceIntervals:
            maxY = max([pN_E[p] + conf[p][1] for p in pN_E.keys()])
        else:
            maxY = max(pN_E.values())
        maxY = 10 * np.ceil(maxY / 10)
    
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
                    figname = 'Fraction_disp_PPIs_transient_%s' % interactome_name)

if __name__ == "__main__":
    main()
