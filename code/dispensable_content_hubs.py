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
from interactome_tools import num_partners, read_single_interface_annotated_interactome
from mutation_interface_edgotype import (assign_edgotypes,
                                         unique_perturbation_mutations,
                                         perturbed_partner_max_degree)
from math_tools import fitness_effect
from plot_tools import multi_histogram_plot, curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
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
    referenceInteractomeFile = interactomeDir / 'human_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    mutationPerturbsFile = interactomeDir / 'unique_mutation_perturbs_geometry.pkl'
        
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
    
    struc_interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    struc_proteins = set(struc_interactome [["Protein_1", "Protein_2"]].values.flatten())
    struc_hubs = struc_proteins & set(hubs)
    struc_nonhubs = struc_proteins & set(nonhubs)
    
    print()
    print('Fraction of hub proteins that are in structural interactome: %.1f%% (%d out of %d)' 
          % (100 * len(struc_hubs) / len(hubs),
             len(struc_hubs),
             len(hubs)))
    print('Fraction of non-hub proteins that are in structural interactome: %.1f%% (%d out of %d)' 
          % (100 * len(struc_nonhubs) / len(nonhubs),
             len(struc_nonhubs),
             len(nonhubs)))
    
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
    # Fraction of dispensable PPIs among all proteins
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum(naturalPerturbs["mono-edgotype"] == 'mono-edgetic')
        numDiseaseMut_edgetic = sum(diseasePerturbs["mono-edgotype"] == 'mono-edgetic')
    else:
        numNaturalMut_edgetic = sum(naturalPerturbs["edgotype"] == 'edgetic')
        numDiseaseMut_edgetic = sum(diseasePerturbs["edgotype"] == 'edgetic')
    
    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('Fitness effect for PPI disruption among all proteins:')
    
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
    # Fraction of dispensable PPIs among hub proteins
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalPerturbs["mono-edgotype"] == 'mono-edgetic') &
                                    (naturalPerturbs["perturbed_partner_max_degree"] >= hub_deg))
        numDiseaseMut_edgetic = sum((diseasePerturbs["mono-edgotype"] == 'mono-edgetic') &
                                    (diseasePerturbs["perturbed_partner_max_degree"] >= hub_deg))
    else:
        numNaturalMut_edgetic = sum((naturalPerturbs["edgotype"] == 'edgetic') & 
                                    (naturalPerturbs["perturbed_partner_max_degree"] >= hub_deg))
        numDiseaseMut_edgetic = sum((diseasePerturbs["edgotype"] == 'edgetic') & 
                                    (diseasePerturbs["perturbed_partner_max_degree"] >= hub_deg))

    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('Fitness effect for PPI disruption among hub proteins:')
    
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
    # Fraction of dispensable PPIs among non-hub proteins
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalPerturbs["mono-edgotype"] == 'mono-edgetic') &
                                    (naturalPerturbs["perturbed_partner_max_degree"] < hub_deg))
        numDiseaseMut_edgetic = sum((diseasePerturbs["mono-edgotype"] == 'mono-edgetic') &
                                    (diseasePerturbs["perturbed_partner_max_degree"] < hub_deg))
    else:
        numNaturalMut_edgetic = sum((naturalPerturbs["edgotype"] == 'edgetic') & 
                                    (naturalPerturbs["perturbed_partner_max_degree"] < hub_deg))
        numDiseaseMut_edgetic = sum((diseasePerturbs["edgotype"] == 'edgetic') & 
                                    (diseasePerturbs["perturbed_partner_max_degree"] < hub_deg))

    numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('Fitness effect for PPI disruption among non-hub proteins:')
    
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
