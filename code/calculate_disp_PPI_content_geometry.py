#----------------------------------------------------------------------------------------
# Calculate the fraction of dispensable protein-protein interactions (PPIs) in the structural
# interactome, i.e., the fraction of PPIs that are effectively neutral upon perturbation.
# The fraction of dispensable PPIs is calculated from geometry-based predictions of PPI 
# perturbations.
#
# Run the following script before running this script:
# - predict_perturbations_geometry.py
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from text_tools import write_list_table
from stat_tools import fisher_test, sderror_on_fraction
from plot_tools import pie_plot
from interactome_tools import num_partners
from mutation_interface_edgotype import (assign_edgotypes,
                                         unique_perturbation_mutations,
                                         perturbed_partner_max_degree)
from math_tools import fitness_effect

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False

    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # minimum interaction degree required for interaction partner perturbed by edgetic mutation
    minDegree = 1
    
    # calculate confidence interval for the fraction of dispensable PPIs
    computeConfidenceIntervals = True
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    referenceInteractomeFile = interactomeDir / 'human_interactome.txt'
    mutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
    
    # output data files
    natMutEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    disMutEdgotypeFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
    dispensablePPIFile = interactomeDir / ('fraction_disp_PPIs_geometry%s.pkl' % ('_monoedgetic' if mono_edgetic else ''))
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    with open(mutationPerturbsFile, 'rb') as f:
        naturalPerturbs, diseasePerturbs = pickle.load(f)
    
    ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
    numPartners = num_partners (ref_interactome)
    naturalPerturbs ["perturbed_partner_max_degree"] = naturalPerturbs.apply(lambda x:
                                                        perturbed_partner_max_degree (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      numPartners), axis=1)
    diseasePerturbs ["perturbed_partner_max_degree"] = diseasePerturbs.apply(lambda x:
                                                        perturbed_partner_max_degree (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      numPartners), axis=1)
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes:' )
    print( '%d non-disease mutations' % len(naturalPerturbs) )
    print( '%d disease mutations' % len(diseasePerturbs) )
    
    naturalPerturbs["edgotype"] = assign_edgotypes (naturalPerturbs["perturbations"].tolist(),
                                                    mono_edgetic = False)
    diseasePerturbs["edgotype"] = assign_edgotypes (diseasePerturbs["perturbations"].tolist(),
                                                    mono_edgetic = False)
    
    naturalPerturbs = naturalPerturbs [naturalPerturbs["edgotype"] != '-'].reset_index(drop=True)
    diseasePerturbs = diseasePerturbs [diseasePerturbs["edgotype"] != '-'].reset_index(drop=True)
    
    if mono_edgetic:
        print( '\n' + 'Labeling mono-edgetic mutations' )
        naturalPerturbs["mono-edgotype"] = assign_edgotypes (naturalPerturbs["perturbations"].tolist(),
                                                             mono_edgetic = True)
        diseasePerturbs["mono-edgotype"] = assign_edgotypes (diseasePerturbs["perturbations"].tolist(),
                                                             mono_edgetic = True)
        naturalPerturbs = naturalPerturbs [naturalPerturbs["mono-edgotype"] != '-'].reset_index(drop=True)
        diseasePerturbs = diseasePerturbs [diseasePerturbs["mono-edgotype"] != '-'].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        naturalPerturbs = naturalPerturbs [unique_perturbation_mutations (naturalPerturbs)].reset_index(drop=True)
        diseasePerturbs = diseasePerturbs [unique_perturbation_mutations (diseasePerturbs)].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Write selected predicted mutation edgotypes to tab-delimited file
    #------------------------------------------------------------------------------------
    
#     write_list_table (naturalPerturbs, ["partners", "perturbations"], natMutEdgotypeFile)
#     write_list_table (diseasePerturbs, ["partners", "perturbations"], disMutEdgotypeFile)
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------

    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalPerturbs["mono-edgotype"] == 'mono-edgetic') &
                                    (naturalPerturbs["perturbed_partner_max_degree"] >= minDegree))
        numDiseaseMut_edgetic = sum((diseasePerturbs["mono-edgotype"] == 'mono-edgetic') &
                                    (diseasePerturbs["perturbed_partner_max_degree"] >= minDegree))
    else:
        numNaturalMut_edgetic = sum((naturalPerturbs["edgotype"] == 'edgetic') & 
                                    (naturalPerturbs["perturbed_partner_max_degree"] >= minDegree))
        numDiseaseMut_edgetic = sum((diseasePerturbs["edgotype"] == 'edgetic') & 
                                    (diseasePerturbs["perturbed_partner_max_degree"] >= minDegree))
    
    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    fracNaturalMut_error = sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered)
    fracDiseaseMut_error = sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered)
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print( '\n' + 'Fraction of predicted %s mutations:' % label )
    print( 'Non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNaturalMut_edgetic,
                                                                    fracNaturalMut_error,
                                                                    numNaturalMut_edgetic,
                                                                    numNaturalMut_considered) )
    
    print( 'Disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDiseaseMut_edgetic,
                                                                fracDiseaseMut_error,
                                                                numDiseaseMut_edgetic,
                                                                numDiseaseMut_considered) )
    
    fisher_test([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    pie_plot([numNaturalMut_nonedgetic, numNaturalMut_edgetic],
             angle = 90,
             colors = ['mediumslateblue', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_%s_mutations_geometry' % label)
    
    pie_plot([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
             angle=90,
             colors = ['mediumslateblue', 'red'],
             edgewidth = 2,
             show = showFigs,
             figdir = figDir,
             figname = 'disease_%s_mutations_geometry' % label)
            
    #------------------------------------------------------------------------------------
    # apply Bayes' theorem to calculate the fraction of PPIs that are dispensable, i.e., 
    # effectively neutral under perturbation
    #------------------------------------------------------------------------------------
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27

    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53

    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Number of neutral mutations (N) that are edgetic (E)
    k_N = numNaturalMut_edgetic
    
    # Number of mildly deleterious mutations (M) that are edgetic (E)
    k_M = numDiseaseMut_edgetic
    
    # Total number of neutral mutations (N)
    n_N = numNaturalMut_considered
    
    # Total number of mildly deleterious mutations (M)
    n_M = numDiseaseMut_considered
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
    
    allresults = fitness_effect (pN,
                                 pM,
                                 pS,
                                 k_N,
                                 n_N,
                                 k_M,
                                 n_M,
                                 pT_S = pE_S,
                                 edgotype = 'edgetic',
                                 CI = 95,
                                 output = True)
    
    results = {k:allresults[k] for k in ('P(N|E)', 'P(N|E)_CI') if k in allresults}
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump(results, fOut)

if __name__ == "__main__":
    main()
