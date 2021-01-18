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
from text_tools import write_list_table
from stat_tools import fisher_test, sderror_on_fraction, proportion_ratio_CI
from plot_tools import pie_plot
from interactome_tools import read_single_interface_annotated_interactome, num_sites, num_partners
from mutation_interface_edgotype import assign_sitotypes

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of edgetic mutations
    mono_edgetic = False
    
    unique_edgetics = False
    
    minPartners = 1
    
    minSites = 1
    
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
    structuralInteractomeFile = interactomeDir / 'human_site_annotated_interactome.txt'
    mutationPerturbsFile = interactomeDir / 'site_perturbs_geometry.pkl'
    
    # output data files
    natMutSitotypeFile = interactomeDir / 'nondisease_mutation_site_edgotype_geometry.txt'
    disMutSitotypeFile = interactomeDir / 'disease_mutation_site_edgotype_geometry.txt'
    dispensableSiteFile = interactomeDir / ('fraction_disp_sites_geometry%s.pkl' % ('_monoedgetic' if mono_edgetic else ''))
    
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
    
    #------------------------------------------------------------------------------------
    # Assign mutation sitotypes
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes:' )
    print( '%d non-disease mutations' % len(naturalPerturbs) )
    print( '%d disease mutations' % len(diseasePerturbs) )
    
    naturalPerturbs["sitotype"] = assign_sitotypes (naturalPerturbs["site_perturbations"].tolist(),
                                                    mono_edgetic = False)
    diseasePerturbs["sitotype"] = assign_sitotypes (diseasePerturbs["site_perturbations"].tolist(),
                                                    mono_edgetic = False)
    
    nat_mono_edgotype = assign_sitotypes (naturalPerturbs["site_perturbations"].tolist(), mono_edgetic = True)
    dis_mono_edgotype = assign_sitotypes (diseasePerturbs["site_perturbations"].tolist(), mono_edgetic = True)
    
    if mono_edgetic:
        print( '\n' + 'Labeling mono-edgetic mutations' )
        naturalPerturbs["mono-sitotype"] = nat_mono_edgotype
        diseasePerturbs["mono-sitotype"] = dis_mono_edgotype
    
    # write predicted mutation edgotypes to tab-delimited file
    write_list_table (naturalPerturbs, ["partners", "perturbations", "sites", "site_perturbations"], natMutSitotypeFile)
    write_list_table (diseasePerturbs, ["partners", "perturbations", "sites", "site_perturbations"], disMutSitotypeFile)
    
    if mono_edgetic:
        naturalPerturbs = naturalPerturbs [(naturalPerturbs["mono-sitotype"] == 'mono-sitic') | 
                                           (naturalPerturbs["mono-sitotype"] == 'non-sitic')].reset_index(drop=True) 
        diseasePerturbs = diseasePerturbs [(diseasePerturbs["mono-sitotype"] == 'mono-sitic') | 
                                           (diseasePerturbs["mono-sitotype"] == 'non-sitic')].reset_index(drop=True) 
    
    ref_interactome = pd.read_table (referenceInteractomeFile, sep='\t')
    numPartners = num_partners (ref_interactome)
    naturalPerturbs = naturalPerturbs [naturalPerturbs["protein"].apply(
                                        lambda x: numPartners[x] >= minPartners)].reset_index(drop=True)
    diseasePerturbs = diseasePerturbs [diseasePerturbs["protein"].apply(
                                        lambda x: numPartners[x] >= minPartners)].reset_index(drop=True)
        
    struc_interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    numSites = num_sites (struc_interactome)
    naturalPerturbs = naturalPerturbs [naturalPerturbs["protein"].apply(
                                        lambda x: numSites[x] >= minSites)].reset_index(drop=True)
    diseasePerturbs = diseasePerturbs [diseasePerturbs["protein"].apply(
                                        lambda x: numSites[x] >= minSites)].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        for mutations in (naturalPerturbs, diseasePerturbs):
            perturbs = {p:{} for p in set(mutations["protein"].values)}
            for _, mut in mutations.iterrows():
                mutPerturbs = set([p for p, pert in zip(mut.sites, mut.site_perturbations) if pert > 0])
                perturbs[mut.protein][(mut.mut_position, mut.mut_res)] = mutPerturbs
        
            mutations["keep"] = True
            for i, mut in mutations.iterrows():
                otherPerturbs = set()
                for k, val in perturbs[mut.protein].items():
                    if k != (mut.mut_position, mut.mut_res):
                        otherPerturbs.update(val)
                mutPerturbs = perturbs[mut.protein][(mut.mut_position, mut.mut_res)]
                if mutPerturbs:
                    if len(mutPerturbs - otherPerturbs) == 0:
                        mutations.loc[i, "keep"] = False
                        del perturbs[mut.protein][(mut.mut_position, mut.mut_res)]
        naturalPerturbs = naturalPerturbs [naturalPerturbs["keep"]].reset_index(drop=True)
        diseasePerturbs = diseasePerturbs [diseasePerturbs["keep"]].reset_index(drop=True)
        del naturalPerturbs["keep"], diseasePerturbs["keep"]
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum(naturalPerturbs["mono-sitotype"] == 'mono-sitic')
        numNaturalMut_nonedgetic = sum(naturalPerturbs["mono-sitotype"].apply(lambda x: 
                                                        x in ('non-sitic', 'sitic')))
        numDiseaseMut_edgetic = sum(diseasePerturbs["mono-sitotype"] == 'mono-sitic')
        numDiseaseMut_nonedgetic = sum(diseasePerturbs["mono-sitotype"].apply(lambda x: 
                                                        x in ('non-sitic', 'sitic')))
    else:
        numNaturalMut_edgetic = sum(naturalPerturbs["sitotype"] == 'sitic')
        numDiseaseMut_edgetic = sum(diseasePerturbs["sitotype"] == 'sitic')
        numNaturalMut_nonedgetic = sum(naturalPerturbs["sitotype"] == 'non-sitic')
        numDiseaseMut_nonedgetic = sum(diseasePerturbs["sitotype"] == 'non-sitic')
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    fracNaturalMut_error = sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered)
    fracDiseaseMut_error = sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered)
    
    label = 'mono-sitic' if mono_edgetic else 'sitic'
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
    
    if not mono_edgetic:
        print( '\n' + 'Fraction of mono-sitic mutations among non-disease sitic mutations:' ) 
        print( '%.3f (%d out of %d)' % (nat_mono_edgotype.count('mono-sitic') / numNaturalMut_edgetic,
                                        nat_mono_edgotype.count('mono-sitic'),
                                        numNaturalMut_edgetic) )
        print( 'Fraction of mono-sitic mutations among disease sitic mutations:' )
        print( '%.3f (%d out of %d)' % (dis_mono_edgotype.count('mono-sitic') / numDiseaseMut_edgetic,
                                        dis_mono_edgotype.count('mono-sitic'),
                                        numDiseaseMut_edgetic) )
        print( 'Fraction of mono-sitic mutations among all sitic mutations:' )
        print( '%.3f (%d out of %d)' % ((nat_mono_edgotype.count('mono-sitic') + dis_mono_edgotype.count('mono-edgetic')) 
                                        / (numNaturalMut_edgetic + numDiseaseMut_edgetic),
                                        nat_mono_edgotype.count('mono-sitic') + dis_mono_edgotype.count('mono-edgetic'),
                                        numNaturalMut_edgetic + numDiseaseMut_edgetic) )
    
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
    
    # Probability for neutral mutations (N) to be edgetic (E)
    pE_N = fracNaturalMut_edgetic
    
    # Probability for mildly deleterious mutations (M) to be edgetic (E)
    pE_M = fracDiseaseMut_edgetic
    
    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0
        
    # Probability for a new missense mutation to be edgetic
    pE = (pE_N * pN) + (pE_M * pM) + (pE_S * pS)
    
    # Probability for edgetic mutations to be effectively neutral
    pN_E = pE_N * pN / pE
    
    allresults = {}
    allresults['P(N|E)'] = pN_E
    
    print('')
    print( 'P(N) = %.1f %%' % (100 * pN) )
    print( 'P(M) = %.1f %%' % (100 * pM) )
    print( 'P(S) = %.1f %%' % (100 * pS) )
    print( 'P(E|N) = %.1f %%' % (100 * pE_N) )
    print( 'P(E|M) = %.1f %%' % (100 * pE_M) )
    print( 'P(E|S) = %.1f %%' % (100 * pE_S) )
    print( 'P(E) = P(E|N)*P(N) + P(E|M)*P(M) + P(E|S)*P(S) = %.1f %%' % (100 * pE) )
    print( 'Fraction of dispensable PPIs P(N|E) = P(E|N)*P(N)/P(E) = %.1f %%' % (100 * pN_E) )
    
    # calculate 95% confidence interval
    if computeConfidenceIntervals:
        n_N, n_M = numNaturalMut_considered, numDiseaseMut_considered
        k_obs_N, k_obs_M = numNaturalMut_edgetic, numDiseaseMut_edgetic
        pE_M_pE_N_lower, pE_M_pE_N_upper = proportion_ratio_CI (k_obs_M,
                                                                n_M,
                                                                k_obs_N,
                                                                n_N,
                                                                conf = CI)
        pN_E_lower = 1 / ( pE_M_pE_N_upper * (pM / pN) + 1 )
        pN_E_upper = 1 / ( pE_M_pE_N_lower * (pM / pN) + 1 )
        print( '%.1f%% confidence interval for P(N|E) = (%f, %f)' 
                % (CI, 100 * pN_E_lower, 100 * pN_E_upper) )
        if not (np.isnan(pN_E_lower) or np.isnan(pN_E_upper)):
            allresults['P(N|E)_CI'] = [pN_E - pN_E_lower, pN_E_upper - pN_E]
    with open(dispensableSiteFile, 'wb') as fOut:
        pickle.dump(allresults, fOut)

if __name__ == "__main__":
    main()
