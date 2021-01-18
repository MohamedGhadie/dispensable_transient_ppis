#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import io
import numpy as np
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import (read_single_interface_annotated_interactome,
                               mutExcl_simult_partners,
                               max_num_simult)
from perturbation_tools import num_perturbed_ppis_n_simult
from math_tools import fitness_effect
from plot_tools import curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False
    
    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # max fraction of interface overlap allowed for simultaneous PPIs
    maxSimultOverlap = 0.05
    
    # min number of mutually exclusive PPIs considered to be medium
    minMedSimult = 1
    
    # max number of mutually exclusive PPIs considered to be medium
    maxMedSimult = 1
    
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
    
    pN_E_keys = ['All PPIs',
                 '<%d' % minMedSimult,
                 '%d-%d' % (minMedSimult, maxMedSimult),
                 '>%d' % maxMedSimult]
    
    pN_E, conf = {}, {}
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'strict'
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    exclusivesOutFile = interactomeDir / 'interactome_mutExcl_simult.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
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
    
    naturalMutations["numSimult"] = naturalMutations.apply(
                                        lambda x: [max_num_simult (x["protein"], p, simultaneous) 
                                                   for p in x["partners"]], axis=1)
    diseaseMutations["numSimult"] = diseaseMutations.apply(
                                        lambda x: [max_num_simult (x["protein"], p, simultaneous)
                                                   for p in x["partners"]], axis=1)
    
    natMut_numLowSimult_disrupt = naturalMutations.apply(
                                    lambda x: num_perturbed_ppis_n_simult (x["perturbations"],
                                                                            x["numSimult"],
                                                                            minN = 0,
                                                                            maxN = minMedSimult - 1),
                                                                            axis=1)
    disMut_numLowSimult_disrupt = diseaseMutations.apply(
                                    lambda x: num_perturbed_ppis_n_simult (x["perturbations"],
                                                                            x["numSimult"],
                                                                            minN = 0,
                                                                            maxN = minMedSimult - 1),
                                                                            axis=1)
    
    natMut_numMedSimult_disrupt = naturalMutations.apply(
                                    lambda x: num_perturbed_ppis_n_simult (x["perturbations"],
                                                                            x["numSimult"],
                                                                            minN = minMedSimult,
                                                                            maxN = maxMedSimult),
                                                                            axis=1)
    disMut_numMedSimult_disrupt = diseaseMutations.apply(
                                    lambda x: num_perturbed_ppis_n_simult (x["perturbations"],
                                                                            x["numSimult"],
                                                                            minN = minMedSimult,
                                                                            maxN = maxMedSimult),
                                                                            axis=1)
    
    natMut_numHighSimult_disrupt = naturalMutations.apply(
                                    lambda x: num_perturbed_ppis_n_simult (x["perturbations"],
                                                                            x["numSimult"],
                                                                            minN = maxMedSimult + 1),
                                                                            axis=1)
    disMut_numHighSimult_disrupt = diseaseMutations.apply(
                                    lambda x: num_perturbed_ppis_n_simult (x["perturbations"],
                                                                            x["numSimult"],
                                                                            minN = maxMedSimult + 1),
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
    
    print(numNaturalMut_considered)
    print(numNaturalMut_edgetic)
    print(numDiseaseMut_considered)
    print(numDiseaseMut_edgetic)
    print()
    print('Fitness effect for disruption among %s:' % pN_E_keys[0])
    
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
    
    if 'P(N|E)' in all_effects:
        pN_E[pN_E_keys[0]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[0]] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among PPIs with low number of simultaneous PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (natMut_numLowSimult_disrupt > 0) &
                                    (natMut_numMedSimult_disrupt == 0) &
                                    (natMut_numHighSimult_disrupt == 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (disMut_numLowSimult_disrupt > 0) &
                                    (disMut_numMedSimult_disrupt == 0) &
                                    (disMut_numHighSimult_disrupt == 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') &
                                    (natMut_numLowSimult_disrupt > 0) &
                                    (natMut_numMedSimult_disrupt == 0) &
                                    (natMut_numHighSimult_disrupt == 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') &
                                    (disMut_numLowSimult_disrupt > 0) &
                                    (disMut_numMedSimult_disrupt == 0) &
                                    (disMut_numHighSimult_disrupt == 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print(numNaturalMut_considered)
    print(numNaturalMut_edgetic)
    print(numDiseaseMut_considered)
    print(numDiseaseMut_edgetic)
    print()
    print('Fitness effect for disruption of PPIs with %s simultaneous PPIs:' % pN_E_keys[1])
    
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
    
    if 'P(N|E)' in all_effects:
        pN_E[pN_E_keys[1]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[1]] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among PPIs with a medium number of simultaneous PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (natMut_numLowSimult_disrupt == 0) &
                                    (natMut_numMedSimult_disrupt > 0) &
                                    (natMut_numHighSimult_disrupt == 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (disMut_numLowSimult_disrupt == 0) &
                                    (disMut_numMedSimult_disrupt > 0) &
                                    (disMut_numHighSimult_disrupt == 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                    (natMut_numLowSimult_disrupt == 0) &
                                    (natMut_numMedSimult_disrupt > 0) &
                                    (natMut_numHighSimult_disrupt == 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                    (disMut_numLowSimult_disrupt == 0) &
                                    (disMut_numMedSimult_disrupt > 0) &
                                    (disMut_numHighSimult_disrupt == 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print(numNaturalMut_considered)
    print(numNaturalMut_edgetic)
    print(numDiseaseMut_considered)
    print(numDiseaseMut_edgetic)
    print()
    print('Fitness effect for disruption of PPIs with %s simultaneous PPIs:' % pN_E_keys[2])
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
    
    if 'P(N|E)' in all_effects:
        pN_E[pN_E_keys[2]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[2]] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among PPIs with a high number of simultaneous PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (natMut_numLowSimult_disrupt == 0) &
                                    (natMut_numMedSimult_disrupt == 0) &
                                    (natMut_numHighSimult_disrupt > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (disMut_numLowSimult_disrupt == 0) &
                                    (disMut_numMedSimult_disrupt == 0) &
                                    (disMut_numHighSimult_disrupt > 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                    (natMut_numLowSimult_disrupt == 0) &
                                    (natMut_numMedSimult_disrupt == 0) &
                                    (natMut_numHighSimult_disrupt > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                    (disMut_numLowSimult_disrupt == 0) &
                                    (disMut_numMedSimult_disrupt == 0) &
                                    (disMut_numHighSimult_disrupt > 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print(numNaturalMut_considered)
    print(numNaturalMut_edgetic)
    print(numDiseaseMut_considered)
    print(numDiseaseMut_edgetic)
    print()
    print('Fitness effect for disruption of PPIs with %s simultaneous PPIs:' % pN_E_keys[3])
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
    
    if 'P(N|E)' in all_effects:
        pN_E[pN_E_keys[3]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[3]] = 100 * lower, 100 * upper
        
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    numGroups = len(pN_E_keys)
    if computeConfidenceIntervals:
        maxY = max([pN_E[p] + conf[p][1] for p in pN_E.keys()])
    else:
        maxY = max(pN_E.values())
    maxY = 10 * np.ceil(maxY / 10)
    
    curve_plot ([(pN_E[p] if p in pN_E else np.nan) for p in pN_E_keys],
                error = [(conf[p] if p in pN_E else [np.nan, np.nan]) for p in pN_E_keys] if computeConfidenceIntervals else None,
                xlim = [0.8, numGroups + 0.1],
                ylim = [0, maxY],
                styles = '.k',
                capsize = 10 if computeConfidenceIntervals else 0,
                msize = 16,
                ewidth = 1.25,
                ecolors = 'k',
                xlabel = 'Number of simultaneous PPIs',
                ylabel = 'Fraction of dispensable PPIs (%)',
                yMinorTicks = 4,
                xticks = list(np.arange(1, numGroups + 1)),
                xticklabels = pN_E_keys,
                yticklabels = list(np.arange(0, maxY + 10, 10)),
                fontsize = 18,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, numGroups),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_Simult')

if __name__ == "__main__":
    main()
