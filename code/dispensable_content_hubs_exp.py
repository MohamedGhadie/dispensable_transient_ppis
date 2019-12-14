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
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import num_partners
from mutation_interface_edgotype import unique_perturbation_mutations, perturbed_partner_max_degree
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
    hub_deg = 10
    
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
    referenceInteractomeFile = procDir / 'HI-II-14' / 'human_interactome.txt'
    natPerturbsFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    disPerturbsFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    naturalPerturbs = read_list_table (natPerturbsFile, ["partners", "perturbations"], [str, int])
    diseasePerturbs = read_list_table (disPerturbsFile, ["partners", "perturbations"], [str, int])
    
    naturalPerturbs["Entrez_Gene_ID"] = naturalPerturbs["Entrez_Gene_ID"].apply(str)
    diseasePerturbs["Entrez_Gene_ID"] = diseasePerturbs["Entrez_Gene_ID"].apply(str)
    
    naturalPerturbs = naturalPerturbs [naturalPerturbs["Edgotype_class"] != '-'].reset_index(drop=True)
    diseasePerturbs = diseasePerturbs [diseasePerturbs["Edgotype_class"] != '-'].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        naturalPerturbs = naturalPerturbs [unique_perturbation_mutations (naturalPerturbs)].reset_index(drop=True)
        diseasePerturbs = diseasePerturbs [unique_perturbation_mutations (diseasePerturbs)].reset_index(drop=True)
    
    print()
    print('%d non-disease mutations' % len(naturalPerturbs))
    print('%d disease mutations' % len(diseasePerturbs))
    
    #------------------------------------------------------------------------------------
    # Fraction of dispensable PPIs among all proteins
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalPerturbs["Edgotype_class"] == 'Edgetic')
    numDiseaseMut_edgetic = sum(diseasePerturbs["Edgotype_class"] == 'Edgetic')
    
    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic

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
                                  output = False)
    
    pN_E_all = 100 * all_effects['P(N|E)']
    if 'P(N|E)_CI' in all_effects:
        lower, upper = all_effects['P(N|E)_CI']
        conf_all = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Identify perturbation partner max degree
    #------------------------------------------------------------------------------------
    
#     ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
#     ref_interactome["EntrezID_1"] = ref_interactome["EntrezID_1"].apply(str)
#     ref_interactome["EntrezID_2"] = ref_interactome["EntrezID_2"].apply(str)
#     
#     numPartners = num_partners (ref_interactome, colnames = ["EntrezID_1", "EntrezID_2"])
#     
#     partners = {p:set() for p in set(ref_interactome[["EntrezID_1", "EntrezID_2"]].values.flatten())}
#     for p1, p2 in ref_interactome[["EntrezID_1", "EntrezID_2"]].values:
#         partners[p1].add(p2)
#         partners[p2].add(p1)
    
    partners = {}
    for _, mut in naturalPerturbs.iterrows():
        if mut.Entrez_Gene_ID in partners:
            partners[mut.Entrez_Gene_ID].update(mut.partners)
        else:
            partners[mut.Entrez_Gene_ID] = set(mut.partners)
        for p in mut.partners:
            if p in partners:
                partners[p].add(mut.Entrez_Gene_ID)
            else:
                partners[p] = {mut.Entrez_Gene_ID}
    for _, mut in diseasePerturbs.iterrows():
        if mut.Entrez_Gene_ID in partners:
            partners[mut.Entrez_Gene_ID].update(mut.partners)
        else:
            partners[mut.Entrez_Gene_ID] = set(mut.partners)
        for p in mut.partners:
            if p in partners:
                partners[p].add(mut.Entrez_Gene_ID)
            else:
                partners[p] = {mut.Entrez_Gene_ID}
    numPartners = {p:len(pr) for p, pr in partners.items()}
    
    proteins, degrees = zip(* [(k, v) for k, v in numPartners.items()])
    degrees = np.array(degrees)
    hubs = [p for p, v in zip(proteins, degrees) if v >= hub_deg]
    nonhubs = [p for p, v in zip(proteins, degrees) if v < hub_deg]
    
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
    
    naturalPerturbs ["perturbed_partner_max_degree"] = naturalPerturbs.apply(lambda x:
                                                        perturbed_partner_max_degree (x["Entrez_Gene_ID"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      numPartners), axis=1)
    diseasePerturbs ["perturbed_partner_max_degree"] = diseasePerturbs.apply(lambda x:
                                                        perturbed_partner_max_degree (x["Entrez_Gene_ID"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      numPartners), axis=1)
    
    #------------------------------------------------------------------------------------
    # Fraction of dispensable PPIs among hub proteins
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum((naturalPerturbs["Edgotype_class"] == 'Edgetic') & 
                                (naturalPerturbs["perturbed_partner_max_degree"] >= hub_deg))
    numDiseaseMut_edgetic = sum((diseasePerturbs["Edgotype_class"] == 'Edgetic') & 
                                (diseasePerturbs["perturbed_partner_max_degree"] >= hub_deg))
    
    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic

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
                                  output = False)
    
    pN_E_hub = 100 * all_effects['P(N|E)']
    if 'P(N|E)_CI' in all_effects:
        lower, upper = all_effects['P(N|E)_CI']
        conf_hub = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Fraction of dispensable PPIs among non-hub proteins
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum((naturalPerturbs["Edgotype_class"] == 'Edgetic') & 
                                (naturalPerturbs["perturbed_partner_max_degree"] < hub_deg))
    numDiseaseMut_edgetic = sum((diseasePerturbs["Edgotype_class"] == 'Edgetic') & 
                                (diseasePerturbs["perturbed_partner_max_degree"] < hub_deg))
    
    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic

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
                                  output = False)
    
    pN_E_nonhub = 100 * all_effects['P(N|E)']
    if 'P(N|E)_CI' in all_effects:
        lower, upper = all_effects['P(N|E)_CI']
        conf_nohub = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    numGroups = 3
    if computeConfidenceIntervals:
        maxY = max([pN_E_all + conf_all[1], pN_E_hub + conf_hub[1], pN_E_nonhub + conf_nohub[1]])
    else:
        maxY = max([pN_E_all, pN_E_hub, pN_E_nonhub])
    maxY = 10 * np.ceil(maxY / 10)
    
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
                xticklabels = ["All proteins", "Non-hubs", "Hubs"],
                yticklabels = list(np.arange(0, maxY + 10, 10)),
                fontsize = 16,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, numGroups),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_hubs')

if __name__ == "__main__":
    main()
