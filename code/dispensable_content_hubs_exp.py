#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table, write_list_table
from interactome_tools import num_partners, is_hub_ppi
from perturbation_tools import (unique_perturbation_mutations,
                                num_nonhub_ppis_perturbed,
                                num_hub_ppis_perturbed)
from math_tools import fitness_effect
from plot_tools import multi_histogram_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'experiment'

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
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    
    # output data files
    natMutOutFile = interactomeDir / 'nondisease_mutation_hub_perturbs.txt'
    disMutOutFile = interactomeDir / 'disease_mutation_hub_perturbs.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
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
    
    print()
    print('%d non-disease mutations' % len(naturalMutations))
    print('%d disease mutations' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    # Identify perturbation partner max degree
    #------------------------------------------------------------------------------------
    
    ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
    ref_interactome["EntrezID_1"] = ref_interactome["EntrezID_1"].apply(str)
    ref_interactome["EntrezID_2"] = ref_interactome["EntrezID_2"].apply(str)
    
    numPartners = num_partners (ref_interactome, colnames = ["EntrezID_1", "EntrezID_2"])
    
    partners = {p:set() for p in set(ref_interactome[["EntrezID_1", "EntrezID_2"]].values.flatten())}
    for p1, p2 in ref_interactome[["EntrezID_1", "EntrezID_2"]].values:
        partners[p1].add(p2)
        partners[p2].add(p1)
    
#     partners = {}
#     for _, mut in naturalMutations.iterrows():
#         if mut.Entrez_Gene_ID in partners:
#             partners[mut.Entrez_Gene_ID].update(mut.partners)
#         else:
#             partners[mut.Entrez_Gene_ID] = set(mut.partners)
#         for p in mut.partners:
#             if p in partners:
#                 partners[p].add(mut.Entrez_Gene_ID)
#             else:
#                 partners[p] = {mut.Entrez_Gene_ID}
#     for _, mut in diseaseMutations.iterrows():
#         if mut.Entrez_Gene_ID in partners:
#             partners[mut.Entrez_Gene_ID].update(mut.partners)
#         else:
#             partners[mut.Entrez_Gene_ID] = set(mut.partners)
#         for p in mut.partners:
#             if p in partners:
#                 partners[p].add(mut.Entrez_Gene_ID)
#             else:
#                 partners[p] = {mut.Entrez_Gene_ID}
#     numPartners = {p:len(pr) for p, pr in partners.items()}
    
    proteins, degrees = zip(* [(k, v) for k, v in numPartners.items()])
    degrees = np.array(degrees)
    hubs = [p for p, v in zip(proteins, degrees) if v >= hubDegree]
    nonhubs = [p for p, v in zip(proteins, degrees) if v < hubDegree]
    
    print()
    print('Fraction of hub proteins: %.1f%% (%d out of %d)' 
          % (100 * len(hubs) / len(proteins),
             len(hubs),
             len(proteins)))
    print('Fraction of non-hub proteins: %.1f%% (%d out of %d)' 
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
    
    naturalMutations ["hub_interactions"] = naturalMutations.apply (lambda x: 
                                                                    [is_hub_ppi (x["Entrez_Gene_ID"],
                                                                                 p,
                                                                                 numPartners,
                                                                                 hubDegree = hubDegree)
                                                                     for p in x["partners"]], axis=1)
    diseaseMutations ["hub_interactions"] = diseaseMutations.apply (lambda x: 
                                                                    [is_hub_ppi (x["Entrez_Gene_ID"],
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
    
#     naturalMutations ["perturbed_partner_max_degree"] = naturalMutations.apply(lambda x:
#                                                         perturbed_partner_max_degree (x["Entrez_Gene_ID"],
#                                                                                       x["partners"],
#                                                                                       x["perturbations"],
#                                                                                       numPartners), axis=1)
#     diseaseMutations ["perturbed_partner_max_degree"] = diseaseMutations.apply(lambda x:
#                                                         perturbed_partner_max_degree (x["Entrez_Gene_ID"],
#                                                                                       x["partners"],
#                                                                                       x["perturbations"],
#                                                                                       numPartners), axis=1)
    
    #------------------------------------------------------------------------------------
    # dispensable content among all PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalMutations["Edgotype_class"] == 'Edgetic')
    numDiseaseMut_edgetic = sum(diseaseMutations["Edgotype_class"] == 'Edgetic')
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
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
    
    pN_E_all = 100 * all_effects['P(N|E)']
    if 'P(N|E)_CI' in all_effects:
        lower, upper = all_effects['P(N|E)_CI']
        conf_all = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # dispensable content among non-hub PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum((naturalMutations["Edgotype_class"] == 'Edgetic') & 
                                (natMut_numNonHubPerturbs > 0))
    numDiseaseMut_edgetic = sum((diseaseMutations["Edgotype_class"] == 'Edgetic') & 
                                (disMut_numNonHubPerturbs > 0))
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
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
    
    pN_E_nonhub = 100 * all_effects['P(N|E)']
    if 'P(N|E)_CI' in all_effects:
        lower, upper = all_effects['P(N|E)_CI']
        conf_nohub = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # dispensable content among hub PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum((naturalMutations["Edgotype_class"] == 'Edgetic') & 
                                (natMut_numNonHubPerturbs == 0) &
                                (natMut_numHubPerturbs > 0))
    numDiseaseMut_edgetic = sum((diseaseMutations["Edgotype_class"] == 'Edgetic') & 
                                (disMut_numNonHubPerturbs == 0) &
                                (disMut_numHubPerturbs > 0))
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
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
    
    pN_E_hub = 100 * all_effects['P(N|E)']
    if 'P(N|E)_CI' in all_effects:
        lower, upper = all_effects['P(N|E)_CI']
        conf_hub = 100 * lower, 100 * upper
        
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    numGroups = 3
    if computeConfidenceIntervals:
        maxY = max([pN_E_all + conf_all[1], pN_E_hub + conf_hub[1], pN_E_nonhub + conf_nohub[1]])
    else:
        maxY = max([pN_E_all, pN_E_hub, pN_E_nonhub])
    maxY = 10 * np.ceil(maxY / 10)
    maxY = 40
    
    curve_plot ([pN_E_all, pN_E_nonhub, pN_E_hub],
                error = [conf_all, conf_nohub, conf_hub] if computeConfidenceIntervals else None,
                xlim = [0.8, numGroups + 0.1],
                ylim = [0, maxY],
                styles = '.k',
                capsize = 10 if computeConfidenceIntervals else 0,
                msize = 16,
                ewidth = 1.25,
                ecolors = 'k',
                ylabel = 'Fraction of dispensable PPIs (%)',
                yMinorTicks = 4,
                xticks = list(np.arange(1, numGroups + 1)),
                xticklabels = ["All PPIs", "Non-hub PPIs", "Hub PPIs"],
                yticklabels = list(np.arange(0, maxY + 10, 10)),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, numGroups),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_hubs')

if __name__ == "__main__":
    main()
