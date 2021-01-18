#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import io
import pickle
import numpy as np
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import (read_single_interface_annotated_interactome,
                               is_single_interface)
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect
from plot_tools import pie_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False
    
    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # max fraction of interface overlap allowed for simultaneous PPIs
    maxInterfaceOverlap = 0.1
    
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
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
    
    pN_E_keys = ['Single-interface proteins', 'Multi-interface proteins']
    
    pN_E, conf = {}, {}
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of edgetic mutation calculation method
    edgeticDir = interactomeDir / 'physics' / (ddg_method + '_edgetics')
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    interfaceOutFile = edgeticDir / 'interactome_single_interface.txt'
    dispensablePPIFile = edgeticDir / 'dispensable_content_SingleInterface.pkl'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (naturalMutationsFile,
                                        ["partners", "perturbations"],
                                        [str, float])
    diseaseMutations = read_list_table (diseaseMutationsFile,
                                        ["partners", "perturbations"],
                                        [str, float])
    
    naturalMutations["perturbations"] = naturalMutations["perturbations"].apply(
                                            lambda x: [int(p) if not np.isnan(p) else p for p in x])
    diseaseMutations["perturbations"] = diseaseMutations["perturbations"].apply(
                                            lambda x: [int(p) if not np.isnan(p) else p for p in x])
    
    naturalMutations = naturalMutations [(naturalMutations["edgotype"] == 'edgetic') |
                                         (naturalMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    diseaseMutations = diseaseMutations [(diseaseMutations["edgotype"] == 'edgetic') |
                                         (diseaseMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        naturalMutations = naturalMutations [unique_perturbation_mutations (naturalMutations)].reset_index(drop=True)
        diseaseMutations = diseaseMutations [unique_perturbation_mutations (diseaseMutations)].reset_index(drop=True)
    
    print()
    print('%d non-disease mutations' % len(naturalMutations))
    print('%d disease mutations' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    
    interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
#     singleIntf = single_interface_proteins (interactome, partners, maxOverlap = maxInterfaceOverlap)
#     
#     with io.open(interfaceOutFile, "w") as fout:
#         fout.write('Protein' + '\t' + 'Single_interface' + '\n')
#         for p in sorted(singleIntf.keys()):
#             fout.write(p + '\t' + singleIntf[p]] + '\n')
    
    print('Calculating single- and multi-interface proteins')
    naturalMutations ["single_interface"] = naturalMutations.apply(
                                    lambda x: is_single_interface (x["protein"],
                                                                   x["partners"],
                                                                   interactome,
                                                                   maxOverlap = maxInterfaceOverlap),
                                                                   axis=1)
    diseaseMutations ["single_interface"] = diseaseMutations.apply(
                                    lambda x: is_single_interface (x["protein"],
                                                                   x["partners"],
                                                                   interactome,
                                                                   maxOverlap = maxInterfaceOverlap),
                                                                   axis=1)
                
    #------------------------------------------------------------------------------------
    # Dispensable content among all PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum(naturalMutations["mono-edgotype"] == 'mono-edgetic')
        numDiseaseMut_edgetic = sum(diseaseMutations["mono-edgotype"] == 'mono-edgetic')
    else:
        numNaturalMut_edgetic = sum(naturalMutations["edgotype"] == 'edgetic')
        numDiseaseMut_edgetic = sum(diseaseMutations["edgotype"] == 'edgetic')
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    
    print('\n********************************************************************')
    print('Dispensable content among all PPIs:')
    print('********************************************************************\n')
    
    print('Fraction of predicted %s mutations:' % label)
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
                                  pT_S = pE_S,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    #------------------------------------------------------------------------------------
    # Dispensable content among PPIs of single-interface proteins
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic_single = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                            (naturalMutations["single_interface"]))
        numDiseaseMut_edgetic_single = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                            (diseaseMutations["single_interface"]))
    else:
        numNaturalMut_edgetic_single = sum((naturalMutations["edgotype"] == 'edgetic') &
                                            (naturalMutations["single_interface"]))
        numDiseaseMut_edgetic_single = sum((diseaseMutations["edgotype"] == 'edgetic') &
                                            (diseaseMutations["single_interface"]))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic_single
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic_single

    numNaturalMut_considered = numNaturalMut_edgetic_single + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic_single + numDiseaseMut_nonedgetic
    
    print('\n********************************************************************')
    print('Dispensable content among PPIs of single-interface proteins:')
    print('********************************************************************\n')
    
    print('Fraction of predicted %s mutations:' % label)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_single / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_single, numNaturalMut_considered),
               numNaturalMut_edgetic_single,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_single / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_single, numDiseaseMut_considered),
               numDiseaseMut_edgetic_single,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_single, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic_single, numDiseaseMut_nonedgetic])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_single,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_single,
                                  numDiseaseMut_considered,
                                  pT_S = pE_S,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Single-interface proteins'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Single-interface proteins'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among PPIs of multi-interface proteins
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic_multi = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                            (naturalMutations["single_interface"] == False))
        numDiseaseMut_edgetic_multi = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                            (diseaseMutations["single_interface"] == False))
    else:
        numNaturalMut_edgetic_multi = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                            (naturalMutations["single_interface"] == False))
        numDiseaseMut_edgetic_multi = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                            (diseaseMutations["single_interface"] == False))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic_multi
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic_multi

    numNaturalMut_considered = numNaturalMut_edgetic_multi + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic_multi + numDiseaseMut_nonedgetic
    
    print('\n********************************************************************')
    print('Dispensable content among PPIs of multi-interface proteins:')
    print('********************************************************************\n')
    
    print('Fraction of predicted %s mutations:' % label)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_multi / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_multi, numNaturalMut_considered),
               numNaturalMut_edgetic_multi,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_multi / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_multi, numDiseaseMut_considered),
               numDiseaseMut_edgetic_multi,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_multi, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic_multi, numDiseaseMut_nonedgetic])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_multi,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_multi,
                                  numDiseaseMut_considered,
                                  pT_S = pE_S,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Multi-interface proteins'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Multi-interface proteins'] = 100 * lower, 100 * upper
    
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump({'DC':pN_E, 'CI':conf}, fOut)
    
    #------------------------------------------------------------------------------------
    # Plot pie charts
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = numNaturalMut_edgetic_multi + numNaturalMut_edgetic_single
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    
    numDiseaseMut_edgetic = numDiseaseMut_edgetic_multi + numDiseaseMut_edgetic_single
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
    
    pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic_single, numNaturalMut_edgetic_multi],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'nondisease_mutations_singleInterface_%s' % label)
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic_single, numDiseaseMut_edgetic_multi],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_singleInterface_%s' % label)
    
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
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
                fontsize = 18,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 2),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_single_interface')

if __name__ == "__main__":
    main()
