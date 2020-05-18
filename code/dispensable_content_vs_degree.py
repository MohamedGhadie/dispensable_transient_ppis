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
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import num_partners
from perturbation_tools import unique_perturbation_mutations, perturbed_partner_max_degree
from math_tools import fitness_effect
from plot_tools import bar_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HuRI'
    
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
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
    
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
    referenceInteractomeFile = interactomeDir / 'reference_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgetics.txt'
        
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutations
    #------------------------------------------------------------------------------------
    
    print('Reading mutations')
    naturalMutations = read_list_table (naturalMutationsFile,
                                        ["partners", "perturbations"],
                                        [str, float])
    diseaseMutations = read_list_table (diseaseMutationsFile,
                                        ["partners", "perturbations"],
                                        [str, float])
    
    naturalMutations = naturalMutations [(naturalMutations["edgotype"] == 'edgetic') |
                                         (naturalMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    diseaseMutations = diseaseMutations [(diseaseMutations["edgotype"] == 'edgetic') |
                                         (diseaseMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    
    print('Number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    print()
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print('Number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    print()
    
    #------------------------------------------------------------------------------------
    # Load interactome data
    #------------------------------------------------------------------------------------
    
    ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
    numPartners = num_partners (ref_interactome)
    
#     maxVal = 100
#     degree_dist = [min(n, maxVal) for n in numPartners.values()]
#     degree_dist = [degree_dist.count(i) for i in np.arange(1, maxVal + 1)]
#     
#     bar_plot (degree_dist,
#               xlabels = list(map(str, np.arange(1, maxVal))) + [ '≥' + str(maxVal)],
#               xlabel = 'Protein interaction degree',
#               ylabel = 'Frequency',
#               colors = 'blue',
#               fontsize = 16,
#               show = showFigs,
#               figdir = figDir,
#               figname = 'protein_degree_dist')

    #------------------------------------------------------------------------------------
    # Calculate the maximum interaction degree among proteins whose interactions are 
    # disrupted by each mutation
    #------------------------------------------------------------------------------------
        
    naturalMutations ["perturbed_partner_max_degree"] = naturalMutations.apply(lambda x:
                                                        perturbed_partner_max_degree (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      numPartners), axis=1)
    diseaseMutations ["perturbed_partner_max_degree"] = diseaseMutations.apply(lambda x:
                                                        perturbed_partner_max_degree (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      numPartners), axis=1)
        
    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        naturalMutations = naturalMutations [unique_perturbation_mutations (naturalMutations)].reset_index(drop=True)
        diseaseMutations = diseaseMutations [unique_perturbation_mutations (diseaseMutations)].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------
    
    degRange = [5] + list(np.arange(10, 101, 10))
    pN_E_high_degree, pN_E_high, conf_high = [], [], []
    for minDegree in degRange:
        if mono_edgetic:
            numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                        (naturalMutations["perturbed_partner_max_degree"] >= minDegree))
            numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                        (diseaseMutations["perturbed_partner_max_degree"] >= minDegree))
        else:
            numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                        (naturalMutations["perturbed_partner_max_degree"] >= minDegree))
            numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                        (diseaseMutations["perturbed_partner_max_degree"] >= minDegree))
    
        numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
        numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
    
        numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
        numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
        
        allresults = fitness_effect (pN,
                                     pM,
                                     pS,
                                     numNaturalMut_edgetic,
                                     numNaturalMut_considered,
                                     numDiseaseMut_edgetic,
                                     numDiseaseMut_considered,
                                     pT_S = pE_S,
                                     edgotype = 'edgetic',
                                     CI = 95 if computeConfidenceIntervals else None,
                                     output = False)
        
        if allresults:
            pN_E_high_degree.append(minDegree)
            pN_E_high.append(100 * allresults['P(N|E)'])
            if 'P(N|E)_CI' in allresults:
                lower, upper = allresults['P(N|E)_CI']
                conf_high.append( (100 * lower, 100 * upper) )
    
    pN_E_low_degree, pN_E_low, conf_low = [], [], []
    for maxDegree in degRange:
        if mono_edgetic:
            numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                        (naturalMutations["perturbed_partner_max_degree"] < maxDegree))
            numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                        (diseaseMutations["perturbed_partner_max_degree"] < maxDegree))
        else:
            numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                        (naturalMutations["perturbed_partner_max_degree"] < maxDegree))
            numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                        (diseaseMutations["perturbed_partner_max_degree"] < maxDegree))
    
        numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
        numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
    
        numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
        numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
        
        allresults = fitness_effect (pN,
                                     pM,
                                     pS,
                                     numNaturalMut_edgetic,
                                     numNaturalMut_considered,
                                     numDiseaseMut_edgetic,
                                     numDiseaseMut_considered,
                                     pT_S = pE_S,
                                     edgotype = 'edgetic',
                                     CI = 95 if computeConfidenceIntervals else None,
                                     output = False)
    
        if allresults:
            pN_E_low_degree.append(maxDegree)
            pN_E_low.append(100 * allresults['P(N|E)'])
            if 'P(N|E)_CI' in allresults:
                lower, upper = allresults['P(N|E)_CI']
                conf_low.append( (100 * lower, 100 * upper) )
    
    maxY_high = max([p + upper for p, (lower, upper) in zip(pN_E_high, conf_high)]) if conf_high else max(pN_E_high)
    maxY_low = max([p + upper for p, (lower, upper) in zip(pN_E_low, conf_low)]) if conf_low else max(pN_E_low)
    maxY = max([maxY_high, maxY_low])
    maxY = 5 * np.ceil(maxY / 5)
    # overwrite maxY
    maxY = 100
    curve_plot ([pN_E_low, pN_E_high],
                xdata = [pN_E_low_degree, pN_E_high_degree],
                error = [conf_low, conf_high],
                ylim = [0, maxY],
                styles = ['.k', '.r'],
                capsize = 10 if conf_low else 0,
                msize = 16,
                ewidth = 2,
                ecolors = ['k', 'r'],
                xlabel = 'Protein interaction degree cutoff (c)',
                ylabel = 'Fraction of PPIs neutral upon disruption (%)',
                leg = ['Degree < c', 'Degree ≥ c'],
                yMinorTicks = 4,
                yticklabels = list(np.arange(0, maxY + 5, 20)),
                fontsize = 16,
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_vs_deg%s%s' % ('_monoedgetic' if mono_edgetic else '',
                                                             '_unique_edgetic' if unique_edgetics else ''))

if __name__ == "__main__":
    main()
