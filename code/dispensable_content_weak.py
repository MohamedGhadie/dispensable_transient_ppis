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
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # structure to use for PPI energy
    # options: template, model
    structure = 'template'
    
    # use the median of PPI binding free energy as a cutoff instead of fixed cutoff
    medEnergy = False
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False

    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # median interaction free energy for all PPIs in structural interactome
    medianEnergy = {'HuRI':     {'template': -20, 'model': -14},
                    'IntAct':   {'template': -15, 'model': -9}}
    
    # minumum interaction free energy required for weak PPIs
    if medEnergy:
        minEnergy = medianEnergy[interactome_name][structure]
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
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
    
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
    
    # directory of edgetic mutation calculation method
    edgeticDir = interactomeDir / 'physics' / (ddg_method + '_edgetics')
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    energyFile = interactomeDir / ('ppi_%s_energy_foldx.txt' % structure)
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    natMutOutFile = edgeticDir / 'nondisease_mutation_weakPPI_perturbs.txt'
    disMutOutFile = edgeticDir / 'disease_mutation_weakPPI_perturbs.txt'
    dispensablePPIFile = edgeticDir / 'dispensable_content_Weak.pkl'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Create edgotype label
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        eCol = 'mono-edgotype'
        edgotype = 'mono-edgetic'
    else:
        eCol = 'edgotype'
        edgotype = 'edgetic'
    
    #------------------------------------------------------------------------------------
    # Load PPI energy
    #------------------------------------------------------------------------------------
    
    energy = read_ppi_energy (energyFile)
    
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
    
#     naturalMutations = naturalMutations [(naturalMutations["edgotype"] == 'edgetic') |
#                                          (naturalMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
#     diseaseMutations = diseaseMutations [(diseaseMutations["edgotype"] == 'edgetic') |
#                                          (diseaseMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    
    naturalMutations = naturalMutations [naturalMutations[eCol] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations[eCol] != '-'].reset_index(drop=True)
    
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
    print('Edgotype = %s' % edgotype)
    print('Non-disease mutations = %d' % len(naturalMutations))
    print('Disease mutations = %d' % len(diseaseMutations))
    
#     if mono_edgetic:
#         numNaturalMut_edgetic = sum(naturalMutations["mono-edgotype"] == 'mono-edgetic')
#         numDiseaseMut_edgetic = sum(diseaseMutations["mono-edgotype"] == 'mono-edgetic')
#     else:
#         numNaturalMut_edgetic = sum(naturalMutations["edgotype"] == 'edgetic')
#         numDiseaseMut_edgetic = sum(diseaseMutations["edgotype"] == 'edgetic')
#     
#     numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
#     numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
    
    numNaturalMut_edgetic = sum(naturalMutations[eCol] == edgotype)
    numDiseaseMut_edgetic = sum(diseaseMutations[eCol] == edgotype)
    
    numNaturalMut_nonedgetic = sum(naturalMutations[eCol] != edgotype)
    numDiseaseMut_nonedgetic = sum(diseaseMutations[eCol] != edgotype)
    
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
    
    print('Fraction of predicted %s mutations:' % edgotype)
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
        pN_E['All PPIs'] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf['All PPIs'] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Label weak and strong PPIs for proteins carrying mutations
    #------------------------------------------------------------------------------------
    
    naturalMutations ["weak_PPIs"] = naturalMutations.apply (lambda x: 
                                                             [is_weak (x["protein"],
                                                                       p,
                                                                       energy,
                                                                       minEnergy = minEnergy)
                                                              for p in x["partners"]], axis=1)
    
    diseaseMutations ["weak_PPIs"] = diseaseMutations.apply (lambda x: 
                                                             [is_weak (x["protein"],
                                                                       p,
                                                                       energy,
                                                                       minEnergy = minEnergy)
                                                              for p in x["partners"]], axis=1)
    
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
    
    numNaturalMut_edgetic_strong = sum((naturalMutations[eCol] == edgotype) &
                                        (natMut_numStrongPerturbs > 0))
    numDiseaseMut_edgetic_strong = sum((diseaseMutations[eCol] == edgotype) &
                                        (disMut_numStrongPerturbs > 0))
    
    numNaturalMut_edgetic_weak = sum((naturalMutations[eCol] == edgotype) & 
                                     (natMut_numWeakPerturbs == natMut_numPerturbs))
    numDiseaseMut_edgetic_weak = sum((diseaseMutations[eCol] == edgotype) & 
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
    
#     if mono_edgetic:
#         numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
#                                     (natMut_numStrongPerturbs > 0))
#         numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
#                                     (disMut_numStrongPerturbs > 0))
#     else:
#         numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') &
#                                     (natMut_numStrongPerturbs > 0))
#         numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') &
#                                     (disMut_numStrongPerturbs > 0))
#     
#     numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
#     numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
# 
#     numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
#     numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('********************************************************************')
    print('Dispensable content among strong PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
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
                                  pT_S = pE_S,
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
    
#     if mono_edgetic:
#         numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
#                                     (natMut_numStrongPerturbs == 0) &
#                                     (natMut_numWeakPerturbs > 0))
#         numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
#                                     (disMut_numStrongPerturbs == 0) &
#                                     (disMut_numWeakPerturbs > 0))
#     else:
#         numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
#                                     (natMut_numStrongPerturbs == 0) &
#                                     (natMut_numWeakPerturbs > 0))
#         numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
#                                     (disMut_numStrongPerturbs == 0) &
#                                     (disMut_numWeakPerturbs > 0))
#     
#     numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
#     numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic
# 
#     numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
#     numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('********************************************************************')
    print('Dispensable content among weak PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
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
                                  pT_S = pE_S,
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
              figname = 'nondisease_mutations_weak_%s' % edgotype)
    pie_plot ([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic_weak, numDiseaseMut_edgetic_strong],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'disease_mutations_weak_%s' % edgotype)
    
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
