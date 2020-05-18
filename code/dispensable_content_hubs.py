#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table, write_list_table
from interactome_tools import (num_partners,
                               is_hub_ppi,
                               read_single_interface_annotated_interactome)
from perturbation_tools import (unique_perturbation_mutations,
                                num_nonhub_ppis_perturbed,
                                num_hub_ppis_perturbed)
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect
from plot_tools import multi_histogram_plot, curve_plot

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
    
    # minimum degree required for a protein hub
    hubDegree = 20
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27
    
    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53
    
    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
    
    pN_E_keys = ['Non-hubs', 'Hubs']
    
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
    figDir = Path('../figures') / interactome_name
    
    # input data files
    referenceInteractomeFile = interactomeDir / 'reference_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    natMutOutFile = interactomeDir / 'nondisease_mutation_hub_perturbs.txt'
    disMutOutFile = interactomeDir / 'disease_mutation_hub_perturbs.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load reference and structural interactome, find hub and non-hub proteins
    #------------------------------------------------------------------------------------
        
    ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
    numPartners = num_partners (ref_interactome)
    
    proteins, degrees = zip(* [(k, v) for k, v in numPartners.items()])
    degrees = np.array(degrees)
    nonhubs = [p for p, v in zip(proteins, degrees) if v < hubDegree]
    hubs = [p for p, v in zip(proteins, degrees) if v >= hubDegree]
    
    print()
    print('Fraction of hub proteins in reference interactome: %.1f%% (%d out of %d)' 
          % (100 * len(hubs) / len(proteins),
             len(hubs),
             len(proteins)))
    print('Fraction of non-hub proteins in reference interactome: %.1f%% (%d out of %d)' 
          % (100 * len(nonhubs) / len(proteins),
             len(nonhubs),
             len(proteins)))
    
    multi_histogram_plot (degrees,
                          xlabel = 'Protein interaction degree',
                          ylabel = 'Frequency',
                          edgecolor = 'k',
                          fontsize = 16,
                          bins = 50,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'protein_degree_dist')
    
    struc_interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    struc_proteins = set(struc_interactome [["Protein_1", "Protein_2"]].values.flatten())
    
    print( '\n' + 'Structural interactome:' )
    print( '%d PPIs' % len(struc_interactome) )
    print( '%d proteins' % len(struc_proteins) )
    print()
    
    struc_nonhubs = struc_proteins & set(nonhubs)
    struc_hubs = struc_proteins & set(hubs)
    
    print()
    print('Fraction of hub proteins that are in structural interactome: %.1f%% (%d out of %d)' 
          % (100 * len(struc_hubs) / len(hubs),
             len(struc_hubs),
             len(hubs)))
    print('Fraction of non-hub proteins that are in structural interactome: %.1f%% (%d out of %d)' 
          % (100 * len(struc_nonhubs) / len(nonhubs),
             len(struc_nonhubs),
             len(nonhubs)))
    
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
    
    naturalMutations ["hub_interactions"] = naturalMutations.apply (lambda x: 
                                                                    [is_hub_ppi (x["protein"],
                                                                                 p,
                                                                                 numPartners,
                                                                                 hubDegree = hubDegree)
                                                                     for p in x["partners"]], axis=1)
    diseaseMutations ["hub_interactions"] = diseaseMutations.apply (lambda x: 
                                                                    [is_hub_ppi (x["protein"],
                                                                                 p,
                                                                                 numPartners,
                                                                                 hubDegree = hubDegree)
                                                                     for p in x["partners"]], axis=1)
    
    write_list_table (naturalMutations, ["partners", "perturbations", "hub_interactions"], natMutOutFile)
    write_list_table (diseaseMutations, ["partners", "perturbations", "hub_interactions"], disMutOutFile)
    
    natMut_numNonHubPerturbs = naturalMutations.apply (lambda x:
                                                       num_nonhub_ppis_perturbed (x["perturbations"],
                                                                                  x["hub_interactions"]), 
                                                                                  axis=1)
    disMut_numNonHubPerturbs = diseaseMutations.apply (lambda x:
                                                       num_nonhub_ppis_perturbed (x["perturbations"],
                                                                                  x["hub_interactions"]), 
                                                                                  axis=1)
    natMut_numHubPerturbs = naturalMutations.apply (lambda x:
                                                    num_hub_ppis_perturbed (x["perturbations"],
                                                                            x["hub_interactions"]), 
                                                                            axis=1)
    disMut_numHubPerturbs = diseaseMutations.apply (lambda x:
                                                    num_hub_ppis_perturbed (x["perturbations"],
                                                                            x["hub_interactions"]), 
                                                                            axis=1)        
    #------------------------------------------------------------------------------------
    # Fraction of dispensable PPIs among all proteins
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
    # Dispensable content among non-hub PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (natMut_numNonHubPerturbs > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (disMut_numNonHubPerturbs > 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                    (natMut_numNonHubPerturbs > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                    (disMut_numNonHubPerturbs > 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print('\n********************************************************************')
    print('Dispensable content among non-hub PPIs:')
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
    
    if 'P(N|E)' in all_effects:
        pN_E['Non-hubs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Non-hubs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among hub PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (natMut_numNonHubPerturbs == 0) &
                                    (natMut_numHubPerturbs > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (disMut_numNonHubPerturbs == 0) &
                                    (disMut_numHubPerturbs > 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                    (natMut_numNonHubPerturbs == 0) &
                                    (natMut_numHubPerturbs > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                    (disMut_numNonHubPerturbs == 0) &
                                    (disMut_numHubPerturbs > 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print('\n********************************************************************')
    print('Dispensable content among hub PPIs:')
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
    
    if 'P(N|E)' in all_effects:
        pN_E['Hubs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Hubs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    if computeConfidenceIntervals:
        maxY = max([pN_E[p] + conf[p][1] for p in pN_E.keys()])
    else:
        maxY = max(pN_E.values())
    maxY = 10 * np.ceil(maxY / 10)
    maxY = 40
    
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
                xticklabels = [p for p in pN_E_keys if p in pN_E],
                yticklabels = list(np.arange(0, maxY + 10, 10)),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 2),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_hubs')

if __name__ == "__main__":
    main()
