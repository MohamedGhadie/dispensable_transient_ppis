#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from perturbation_tools import unique_perturbation_mutations
from stat_tools import sderror_on_fraction, fisher_test
from math_tools import fitness_effect
from plot_tools import pie_plot

def main():
    
    # reference interactome name
    interactome_name = 'Sahni'

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
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    
    # output data files
    dispensablePPIFile = interactomeDir / 'dispensable_content.pkl'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load mutation edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations = pd.read_table (naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table (diseaseMutationsFile, sep='\t')
    
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
    
    print()
    print('%d non-disease mutations' % len(naturalMutations))
    print('%d disease mutations' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    # dispensable content among all PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalMutations["Edgotype_class"] == 'Edgetic')
    numDiseaseMut_edgetic = sum(diseaseMutations["Edgotype_class"] == 'Edgetic')
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
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
    
    pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic],
              angle = 90,
              colors = ['lightsteelblue', 'orange'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'nondisease_mutations_edgetic')
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
              angle=90,
              colors = ['lightsteelblue', 'orange'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_edgetic')
    
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
        pN_E = {'All PPIs': 100 * all_effects['P(N|E)']}
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf = {'All PPIs': (100 * lower, 100 * upper)}
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump({'DC':pN_E, 'CI':conf}, fOut)
    
if __name__ == "__main__":
    main()