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
                               mutExcl_simult_partners,
                               max_num_mutExcl)
from perturbation_tools import num_perturbed_ppis_n_mutExcl
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
    maxSimultOverlap = 0.1
    
    # min number of mutually exclusive PPIs considered to be medium
    minMedMutExcl = 1
    
    # max number of mutually exclusive PPIs considered to be medium
    maxMedMutExcl = 4
    
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
    
#     pN_E_keys = ['<%d' % (minMedMutExcl + 1),
#                  '%d-%d' % (minMedMutExcl + 1, maxMedMutExcl + 1),
#                  '>%d' % (maxMedMutExcl + 1)]
    
    pN_E_keys = ['Simultaneously possible PPIs',
                 'Moderately exclusive PPIs',
                 'Highly exclusive PPIs']
    
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
    exclusivesOutFile = edgeticDir / 'interactome_mutExcl_simult.txt'
    dispensablePPIFile = edgeticDir / 'dispensable_content_MutExcl.pkl'

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
    mutExclusive, simultaneous = mutExcl_simult_partners (interactome, maxSimultOverlap = maxSimultOverlap)
    
    with io.open(exclusivesOutFile, "w") as fout:
        fout.write('\t'.join(['Protein_1',
                              'Protein_2',
                              'Protein_1_mutually_exclusive_partners',
                              'Protein_1_simultaneous_partners']) + '\n')
        for p in sorted(mutExclusive.keys()):
            for partner in sorted(mutExclusive[p].keys()):
                excl, simult = mutExclusive[p][partner], simultaneous[p][partner]
                excl = ','.join(sorted(excl)) if excl else '-'
                simult = ','.join(sorted(simult)) if simult else '-'
                fout.write('\t'.join([p, partner, excl, simult]) + '\n')
    
    mutExcl = []
    for p1, p2 in interactome[["Protein_1", "Protein_2"]].values:
        if (len(mutExclusive[p1][p2]) > 0) or (len(mutExclusive[p2][p1]) > 0):
            mutExcl.append(1)
        else:
            mutExcl.append(0)
    print('%d out of %d PPIs are mutually exclusive' % (sum(mutExcl), len(mutExcl)))
    
    
    naturalMutations["numMutExcl"] = naturalMutations.apply(
                                        lambda x: [max_num_mutExcl (x["protein"], p, mutExclusive) 
                                                   for p in x["partners"]], axis=1)
    diseaseMutations["numMutExcl"] = diseaseMutations.apply(
                                        lambda x: [max_num_mutExcl (x["protein"], p, mutExclusive)
                                                   for p in x["partners"]], axis=1)
    
    natMut_numLowMutExcl_disrupt = naturalMutations.apply(
                                    lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
                                                                            x["numMutExcl"],
                                                                            minN = 0,
                                                                            maxN = minMedMutExcl - 1),
                                                                            axis=1)
    disMut_numLowMutExcl_disrupt = diseaseMutations.apply(
                                    lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
                                                                            x["numMutExcl"],
                                                                            minN = 0,
                                                                            maxN = minMedMutExcl - 1),
                                                                            axis=1)
    
    natMut_numMedMutExcl_disrupt = naturalMutations.apply(
                                    lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
                                                                            x["numMutExcl"],
                                                                            minN = minMedMutExcl,
                                                                            maxN = maxMedMutExcl),
                                                                            axis=1)
    disMut_numMedMutExcl_disrupt = diseaseMutations.apply(
                                    lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
                                                                            x["numMutExcl"],
                                                                            minN = minMedMutExcl,
                                                                            maxN = maxMedMutExcl),
                                                                            axis=1)
    
    natMut_numHighMutExcl_disrupt = naturalMutations.apply(
                                    lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
                                                                            x["numMutExcl"],
                                                                            minN = maxMedMutExcl + 1),
                                                                            axis=1)
    disMut_numHighMutExcl_disrupt = diseaseMutations.apply(
                                    lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
                                                                            x["numMutExcl"],
                                                                            minN = maxMedMutExcl + 1),
                                                                            axis=1)
#     natMut_numHighMutExcl_disrupt = naturalMutations.apply(
#                                     lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
#                                                                             x["numMutExcl"],
#                                                                             minN = minMedMutExcl),
#                                                                             axis=1)
#     disMut_numHighMutExcl_disrupt = diseaseMutations.apply(
#                                     lambda x: num_perturbed_ppis_n_mutExcl (x["perturbations"],
#                                                                             x["numMutExcl"],
#                                                                             minN = minMedMutExcl),
#                                                                             axis=1)
#     naturalPerturbs ["perturbed_ppis_max_num_mutExcl"] = naturalPerturbs.apply(lambda x:
#                                                         perturbed_ppis_max_num_mutExcl (x["protein"],
#                                                                                         x["partners"],
#                                                                                         x["perturbations"],
#                                                                                         mutExclusive), axis=1)
#     diseasePerturbs ["perturbed_ppis_max_num_mutExcl"] = diseasePerturbs.apply(lambda x:
#                                                         perturbed_ppis_max_num_mutExcl (x["protein"],
#                                                                                         x["partners"],
#                                                                                         x["perturbations"],
#                                                                                         mutExclusive), axis=1)
        
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
    
    print('\n********************************************************************')
    print('Dispensable content among all PPIs:')
    print('********************************************************************\n')
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
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
    # Dispensable content among simultaneous PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic_lowME = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                          (natMut_numLowMutExcl_disrupt > 0) &
                                          (natMut_numMedMutExcl_disrupt == 0) &
                                          (natMut_numHighMutExcl_disrupt == 0))
        numDiseaseMut_edgetic_lowME = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                          (disMut_numLowMutExcl_disrupt > 0) &
                                          (disMut_numMedMutExcl_disrupt == 0) &
                                          (disMut_numHighMutExcl_disrupt == 0))
    else:
        numNaturalMut_edgetic_lowME = sum((naturalMutations["edgotype"] == 'edgetic') &
                                          (natMut_numLowMutExcl_disrupt > 0) &
                                          (natMut_numMedMutExcl_disrupt == 0) &
                                          (natMut_numHighMutExcl_disrupt == 0))
        numDiseaseMut_edgetic_lowME = sum((diseaseMutations["edgotype"] == 'edgetic') &
                                          (disMut_numLowMutExcl_disrupt > 0) &
                                          (disMut_numMedMutExcl_disrupt == 0) &
                                          (disMut_numHighMutExcl_disrupt == 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic_lowME
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic_lowME

    numNaturalMut_considered = numNaturalMut_edgetic_lowME + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic_lowME + numDiseaseMut_nonedgetic
    
    print('\n********************************************************************')
    print('Dispensable content among simultaneous PPIs:')
    print('********************************************************************\n')
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print('Fraction of predicted %s mutations:' % label)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_lowME / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_lowME, numNaturalMut_considered),
               numNaturalMut_edgetic_lowME,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_lowME / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_lowME, numDiseaseMut_considered),
               numDiseaseMut_edgetic_lowME,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_lowME, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic_lowME, numDiseaseMut_nonedgetic])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_lowME,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_lowME,
                                  numDiseaseMut_considered,
                                  pT_S = pE_S,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Simultaneously possible PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Simultaneously possible PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among moderately mutually exclusive PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic_medME = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                          (natMut_numLowMutExcl_disrupt == 0) &
                                          (natMut_numMedMutExcl_disrupt > 0) &
                                          (natMut_numHighMutExcl_disrupt == 0))
        numDiseaseMut_edgetic_medME = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                          (disMut_numLowMutExcl_disrupt == 0) &
                                          (disMut_numMedMutExcl_disrupt > 0) &
                                          (disMut_numHighMutExcl_disrupt == 0))
    else:
        numNaturalMut_edgetic_medME = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                          (natMut_numLowMutExcl_disrupt == 0) &
                                          (natMut_numMedMutExcl_disrupt > 0) &
                                          (natMut_numHighMutExcl_disrupt == 0))
        numDiseaseMut_edgetic_medME = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                          (disMut_numLowMutExcl_disrupt == 0) &
                                          (disMut_numMedMutExcl_disrupt > 0) &
                                          (disMut_numHighMutExcl_disrupt == 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic_medME
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic_medME

    numNaturalMut_considered = numNaturalMut_edgetic_medME + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic_medME + numDiseaseMut_nonedgetic
    
    print('\n*************************************************************************')
    print('Dispensable content among moderately mutually exclusive PPIs:')
    print('*************************************************************************\n')
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print('Fraction of predicted %s mutations:' % label)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_medME / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_medME, numNaturalMut_considered),
               numNaturalMut_edgetic_medME,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_medME / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_medME, numDiseaseMut_considered),
               numDiseaseMut_edgetic_medME,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_medME, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic_medME, numDiseaseMut_nonedgetic])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_medME,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_medME,
                                  numDiseaseMut_considered,
                                  pT_S = pE_S,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Moderately exclusive PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Moderately exclusive PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among highly mutually exclusive PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic_highME = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                           (natMut_numLowMutExcl_disrupt == 0) &
                                           (natMut_numMedMutExcl_disrupt == 0) &
                                           (natMut_numHighMutExcl_disrupt > 0))
        numDiseaseMut_edgetic_highME = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                           (disMut_numLowMutExcl_disrupt == 0) &
                                           (disMut_numMedMutExcl_disrupt == 0) &
                                           (disMut_numHighMutExcl_disrupt > 0))
    else:
        numNaturalMut_edgetic_highME = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                           (natMut_numLowMutExcl_disrupt == 0) &
                                           (natMut_numMedMutExcl_disrupt == 0) &
                                           (natMut_numHighMutExcl_disrupt > 0))
        numDiseaseMut_edgetic_highME = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                           (disMut_numLowMutExcl_disrupt == 0) &
                                           (disMut_numMedMutExcl_disrupt == 0) &
                                           (disMut_numHighMutExcl_disrupt > 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic_highME
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic_highME

    numNaturalMut_considered = numNaturalMut_edgetic_highME + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic_highME + numDiseaseMut_nonedgetic
    
    print('\n***********************************************************************')
    print('Dispensable content among highly mutually exclusive PPIs:')
    print('***********************************************************************\n')
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print('Fraction of predicted %s mutations:' % label)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_highME / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_highME, numNaturalMut_considered),
               numNaturalMut_edgetic_highME,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_highME / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_highME, numDiseaseMut_considered),
               numDiseaseMut_edgetic_highME,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_highME, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic_highME, numDiseaseMut_nonedgetic])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_highME,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_highME,
                                  numDiseaseMut_considered,
                                  pT_S = pE_S,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Highly exclusive PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Highly exclusive PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Save dispensable content to file
    #------------------------------------------------------------------------------------
    
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump({'DC':pN_E, 'CI':conf}, fOut)
    
    #------------------------------------------------------------------------------------
    # Plot pie charts
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = (numNaturalMut_edgetic_highME + 
                             numNaturalMut_edgetic_medME + 
                             numNaturalMut_edgetic_lowME)
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    
    numDiseaseMut_edgetic = (numDiseaseMut_edgetic_highME + 
                             numDiseaseMut_edgetic_medME + 
                             numDiseaseMut_edgetic_lowME)
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
    
    pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic_highME,
               numNaturalMut_edgetic_medME, numNaturalMut_edgetic_lowME],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'orange', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'nondisease_mutations_mutExcl_%s' % label)
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic_highME,
               numDiseaseMut_edgetic_medME, numDiseaseMut_edgetic_lowME],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'orange', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_mutExcl_%s' % label)
    
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
                xticks = [1, 2, 3],
                xticklabels = [p.replace(' ', '\n') for p in pN_E_keys if p in pN_E],
                yticklabels = list(np.arange(0, maxY + 5, 10)),
                fontsize = 18,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 3),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_mutExcl')

if __name__ == "__main__":
    main()
