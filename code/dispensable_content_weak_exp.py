#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from text_tools import read_list_table, write_list_table
from perturbation_tools import (unique_perturbation_mutations,
                                num_strong_ppis_perturbed,
                                num_weak_ppis_perturbed)
from energy_tools import read_ppi_energy
from protein_function import is_weak
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect
from plot_tools import pie_plot, curve_plot

def main():
    
    # reference interactome name
    interactome_name = 'Sahni'
    
    # structure to use for PPI energy
    # options: template, model
    structure = 'template'
    
    # use the median of PPI binding free energy as a cutoff instead of fixed cutoff
    medEnergy = False
    
    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # median interaction free energy for all PPIs in structural interactome
    medianEnergy = {'template': -20, 'model': -14}
    
    # minumum interaction free energy required for weak PPIs
    if medEnergy:
        minEnergy = medianEnergy[structure]
    else:
        minEnergy = -25
    
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
    
    pN_E_keys = ['All PPIs', 'Strong PPIs', 'Weak PPIs']
    
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
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    energyFile = interactomeDir / ('ppi_%s_energy_foldx.txt' % structure)
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    
    # output data files
    natMutOutFile = interactomeDir / 'nondisease_mutation_weakPPI_perturbs.txt'
    disMutOutFile = interactomeDir / 'disease_mutation_weakPPI_perturbs.txt'
    dispensablePPIFile = interactomeDir / 'dispensable_content_Weak.pkl'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Load PPI energy
    #------------------------------------------------------------------------------------
    
    energy = read_ppi_energy (energyFile)
    
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
    
    #------------------------------------------------------------------------------------
    # Print summary of input parameters
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Summary of input parameters')
    print('********************************************************************')
    print()
    print('Interactome dataset = %s' % interactome_name)
    print('Using median for energy cutoff = %s' % ('Yes' if medEnergy else 'No'))
    print('Minimum energy for transient PPIs = %f' % minEnergy)
    print('Non-disease mutations = %d' % len(naturalMutations))
    print('Disease mutations = %d' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    # Count edgetic and non-edgetic mutations among all PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalMutations["Edgotype_class"] == 'Edgetic')
    numDiseaseMut_edgetic = sum(diseaseMutations["Edgotype_class"] == 'Edgetic')
    
    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
   #------------------------------------------------------------------------------------
    # Dispensable content among all PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among all PPIs:')
    print('********************************************************************')
    print()
    
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
        pN_E['All PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['All PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Label weak and strong PPIs for proteins carrying mutations
    #------------------------------------------------------------------------------------
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    naturalMutations ["weak_PPIs"] = naturalMutations.apply(
        lambda x: [is_weak (uniprotID [x["Entrez_Gene_ID"]] if x["Entrez_Gene_ID"] in uniprotID else '-',
                            uniprotID [p] if p in uniprotID else '-',
                            energy,
                            minEnergy = minEnergy) for p in x["partners"]], axis=1)
    
    diseaseMutations ["weak_PPIs"] = diseaseMutations.apply(
        lambda x: [is_weak (uniprotID [x["Entrez_Gene_ID"]] if x["Entrez_Gene_ID"] in uniprotID else '-',
                            uniprotID [p] if p in uniprotID else '-',
                            energy,
                            minEnergy = minEnergy) for p in x["partners"]], axis=1)
    
    write_list_table (naturalMutations, ["partners", "perturbations", "weak_PPIs"], natMutOutFile)
    write_list_table (diseaseMutations, ["partners", "perturbations", "weak_PPIs"], disMutOutFile)
    
    #------------------------------------------------------------------------------------
    # Count edgetic mutations that disrupt weak or strong PPIs
    #------------------------------------------------------------------------------------
    
    natMut_numPerturbs = naturalMutations["perturbations"].apply(np.nansum).apply(int)
    disMut_numPerturbs = diseaseMutations["perturbations"].apply(np.nansum).apply(int)
    
    natMut_numStrongPerturbs = naturalMutations.apply(
        lambda x: num_strong_ppis_perturbed (x["perturbations"], x["weak_PPIs"]), axis=1)
    disMut_numStrongPerturbs = diseaseMutations.apply(
        lambda x: num_strong_ppis_perturbed (x["perturbations"], x["weak_PPIs"]), axis=1)
    
    natMut_numWeakPerturbs = naturalMutations.apply(
        lambda x: num_weak_ppis_perturbed (x["perturbations"], x["weak_PPIs"]), axis=1)
    disMut_numWeakPerturbs = diseaseMutations.apply(
        lambda x: num_weak_ppis_perturbed (x["perturbations"], x["weak_PPIs"]), axis=1)
    
    numNaturalMut_edgetic_strong = sum((naturalMutations["Edgotype_class"] == 'Edgetic') &
                                        (natMut_numStrongPerturbs > 0))
    numDiseaseMut_edgetic_strong = sum((diseaseMutations["Edgotype_class"] == 'Edgetic') &
                                        (disMut_numStrongPerturbs > 0))
    
    numNaturalMut_edgetic_weak = sum((naturalMutations["Edgotype_class"] == 'Edgetic') & 
                                     (natMut_numWeakPerturbs == natMut_numPerturbs))
    numDiseaseMut_edgetic_weak = sum((diseaseMutations["Edgotype_class"] == 'Edgetic') & 
                                     (disMut_numWeakPerturbs == disMut_numPerturbs))
    
    numNaturalMut_considered = (numNaturalMut_nonedgetic + 
                                numNaturalMut_edgetic_strong + 
                                numNaturalMut_edgetic_weak) 
                            
    numDiseaseMut_considered = (numDiseaseMut_nonedgetic + 
                                numDiseaseMut_edgetic_strong + 
                                numDiseaseMut_edgetic_weak)
    
#     print(sum((natMut_numPerturbs == 1) & (natMut_numWeakPerturbs == 1)))
#     print(sum((disMut_numPerturbs == 1) & (disMut_numWeakPerturbs == 1)))
    
    #------------------------------------------------------------------------------------
    # Dispensable content among strong PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among strong PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted edgetic mutations:')
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_strong / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_strong, numNaturalMut_considered),
               numNaturalMut_edgetic_strong,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_strong / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_strong, numDiseaseMut_considered),
               numDiseaseMut_edgetic_strong,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_strong, numNaturalMut_considered - numNaturalMut_edgetic_strong],
                 [numDiseaseMut_edgetic_strong, numDiseaseMut_considered - numDiseaseMut_edgetic_strong])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_strong,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_strong,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
        
    if 'P(N|E)' in all_effects:
        pN_E['Strong PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Strong PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among weak PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among weak PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted edgetic mutations:')
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_weak / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_weak, numNaturalMut_considered),
               numNaturalMut_edgetic_weak,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_weak / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_weak, numDiseaseMut_considered),
               numDiseaseMut_edgetic_weak,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_weak, numNaturalMut_considered - numNaturalMut_edgetic_weak],
                 [numDiseaseMut_edgetic_weak, numDiseaseMut_considered - numDiseaseMut_edgetic_weak])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_weak,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_weak,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    if 'P(N|E)' in all_effects:
        pN_E['Weak PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['Weak PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Save dispensable content to file
    #------------------------------------------------------------------------------------
    
    with open(dispensablePPIFile, 'wb') as fOut:
        pickle.dump({'DC':pN_E, 'CI':conf}, fOut)
    
    #------------------------------------------------------------------------------------
    # Plot pie charts
    #------------------------------------------------------------------------------------
    
    pie_plot ([numNaturalMut_nonedgetic, numNaturalMut_edgetic_weak, numNaturalMut_edgetic_strong],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'nondisease_mutations_weak_edgetic')
    
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic_weak, numDiseaseMut_edgetic_strong],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_weak_edgetic')
    
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
                xticks = [1, 2],
                xticklabels = [p.replace(' ', '\n') for p in pN_E_keys if p in pN_E],
                yticklabels = list(np.arange(0, maxY + 5, 10)),
                fontsize = 20,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, 2),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_weak')

if __name__ == "__main__":
    main()
