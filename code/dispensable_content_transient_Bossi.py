#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from text_tools import read_list_table
from perturbation_tools import (unique_perturbation_mutations,
                                num_permanent_ppis_perturbed,
                                num_transient_ppis_perturbed)
from stat_tools import fisher_test, sderror_on_fraction
from math_tools import fitness_effect

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'IntAct'
    
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
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27
    
    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53
    
    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20
    
    pN_E_keys = ['Permanent PPIs', 'Transient PPIs']
    
    pN_E, conf = {}, {}
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of edgetic mutation calculation method
    edgeticDir = interactomeDir / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    ppiTypeFile = extDir / 'CRG.integrated.human.interactome.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
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
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    ppiTypes = pd.read_table (ppiTypeFile, sep='\t')
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    total = len(ppiTypes.columns) - 6
    print('Total number of tissues = %d' % total)
    t = 0.9 * total
    
    types = {}
    for _,row in ppiTypes.iterrows():
        if (row.Gene1 in uniprotID) and (row.Gene2 in uniprotID):
            expr = [e for e in row[6:] if not np.isnan(e)]
            if len(expr) >= t:
                p1 = uniprotID[row.Gene1]
                p2 = uniprotID[row.Gene2]
                types[tuple(sorted([p1, p2]))] = 'permanent' if sum(expr) >= t else 'transient'
        
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
    print('Edgotype = %s' % edgotype)
    print('Non-disease mutations = %d' % len(naturalMutations))
    print('Disease mutations = %d' % len(diseaseMutations))
    
    #------------------------------------------------------------------------------------
    # Count edgetic and non-edgetic mutations among all PPIs
    #------------------------------------------------------------------------------------
    
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
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    #------------------------------------------------------------------------------------
    # Label transient and permanent PPIs for proteins carrying mutations
    #------------------------------------------------------------------------------------
    
    keys = naturalMutations.apply(lambda x: [tuple(sorted([x["protein"], p])) for p in x["partners"]], axis=1)
    naturalMutations ["transient_PPIs"] = keys.apply(lambda x: [(types[k] if k in types else '-') for k in x])
    
    keys = diseaseMutations.apply(lambda x: [tuple(sorted([x["protein"], p])) for p in x["partners"]], axis=1)
    diseaseMutations ["transient_PPIs"] = keys.apply(lambda x: [(types[k] if k in types else '-') for k in x])
    
    #------------------------------------------------------------------------------------
    # Count edgetic mutations that disrupt transient or permanent PPIs
    #------------------------------------------------------------------------------------
    
    natMut_numPerturbs = naturalMutations["perturbations"].apply(np.nansum).apply(int)
    disMut_numPerturbs = diseaseMutations["perturbations"].apply(np.nansum).apply(int)
    
    natMut_numPermPerturbs = naturalMutations.apply(
        lambda x: num_permanent_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    disMut_numPermPerturbs = diseaseMutations.apply(
        lambda x: num_permanent_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    
    natMut_numTransPerturbs = naturalMutations.apply(
        lambda x: num_transient_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    disMut_numTransPerturbs = diseaseMutations.apply(
        lambda x: num_transient_ppis_perturbed (x["perturbations"], x["transient_PPIs"]), axis=1)
    
    numNaturalMut_edgetic_perm = sum((naturalMutations[eCol] == edgotype) &
                                     (natMut_numPermPerturbs > 0))
    numDiseaseMut_edgetic_perm = sum((diseaseMutations[eCol] == edgotype) &
                                     (disMut_numPermPerturbs > 0))
    
    numNaturalMut_edgetic_tran = sum((naturalMutations[eCol] == edgotype) & 
                                     (natMut_numTransPerturbs == natMut_numPerturbs))
    numDiseaseMut_edgetic_tran = sum((diseaseMutations[eCol] == edgotype) & 
                                     (disMut_numTransPerturbs == disMut_numPerturbs))
    
    numNaturalMut_considered = (numNaturalMut_nonedgetic + 
                                numNaturalMut_edgetic_perm + 
                                numNaturalMut_edgetic_tran) 
                            
    numDiseaseMut_considered = (numDiseaseMut_nonedgetic + 
                                numDiseaseMut_edgetic_perm + 
                                numDiseaseMut_edgetic_tran)
    
    #------------------------------------------------------------------------------------
    # Dispensable content among permanent PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among permanent PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_perm / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_perm, numNaturalMut_considered),
               numNaturalMut_edgetic_perm,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_perm / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_perm, numDiseaseMut_considered),
               numDiseaseMut_edgetic_perm,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_perm, numNaturalMut_considered - numNaturalMut_edgetic_perm],
                 [numDiseaseMut_edgetic_perm, numDiseaseMut_considered - numDiseaseMut_edgetic_perm])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_perm,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_perm,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)
    
    #------------------------------------------------------------------------------------
    # Dispensable content among transient PPIs
    #------------------------------------------------------------------------------------
    
    print()
    print('********************************************************************')
    print('Dispensable content among transient PPIs:')
    print('********************************************************************')
    print()
    
    print('Fraction of predicted %s mutations:' % edgotype)
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic_tran / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic_tran, numNaturalMut_considered),
               numNaturalMut_edgetic_tran,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic_tran / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic_tran, numDiseaseMut_considered),
               numDiseaseMut_edgetic_tran,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic_tran, numNaturalMut_considered - numNaturalMut_edgetic_tran],
                 [numDiseaseMut_edgetic_tran, numDiseaseMut_considered - numDiseaseMut_edgetic_tran])
    
    print()
    all_effects = fitness_effect (pN,
                                  pM,
                                  pS,
                                  numNaturalMut_edgetic_tran,
                                  numNaturalMut_considered,
                                  numDiseaseMut_edgetic_tran,
                                  numDiseaseMut_considered,
                                  pT_S = 0,
                                  edgotype = 'edgetic',
                                  CI = 95 if computeConfidenceIntervals else None,
                                  output = True)

if __name__ == "__main__":
    main()
