#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from perturbation_tools import unique_perturbation_mutations
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect
from plot_tools import pie_plot

def main():
    
    # reference interactome name
    # options: HuRI, IntAct
    interactome_name = 'HuRI'
    
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
    
    # show figures
    showFigs = False
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27
    
    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53
    
    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
    
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
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    dispensablePPIFile = edgeticDir / ('fraction_disp_PPIs%s.pkl' % ('_monoedgetic' if mono_edgetic else ''))
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load mutation edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
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
    
    #------------------------------------------------------------------------------------
    # Print mutation numbers
    #------------------------------------------------------------------------------------
    
    print('\n' + 'Total number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print('\n' + 'Total number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    
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
    return
    fisher_test ([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic],
              angle = 90,
              colors = ['lightsteelblue', 'orange'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'nondisease_mutations_%s' % label)
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
              angle=90,
              colors = ['lightsteelblue', 'orange'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_%s' % label)
    
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
    
    dispContent = {}
    if 'P(N|E)' in all_effects:
        dispContent['P(N|E)'] = all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            pN_E_lower, pN_E_upper = all_effects['P(N|E)_CI']
            dispContent['P(N|E)_CI'] = [pN_E_lower, pN_E_upper]
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump(dispContent, fOut)

if __name__ == "__main__":
    main()
