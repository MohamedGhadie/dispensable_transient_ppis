#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import (read_single_interface_annotated_interactome,
                               mutExcl_simult_partners)
from perturbation_analysis import perturbed_ppis_max_num_mutExcl
from math_tools import fitness_effect
from plot_tools import curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False
    
    # max fraction of interface overlap for simultaneous PPIs
    overlap_cutoff = 0.5
    
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
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    natPerturbsFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    disPerturbsFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
        
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
    
    interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    mutExclusive, simultaneous = mutExcl_simult_partners (interactome, cutoff = overlap_cutoff)
    
    naturalPerturbs ["perturbed_ppis_max_num_mutExcl"] = naturalPerturbs.apply(lambda x:
                                                        perturbed_ppis_max_num_mutExcl (x["protein"],
                                                                                        x["partners"],
                                                                                        x["perturbations"],
                                                                                        mutExclusive), axis=1)
    diseasePerturbs ["perturbed_ppis_max_num_mutExcl"] = diseasePerturbs.apply(lambda x:
                                                        perturbed_ppis_max_num_mutExcl (x["protein"],
                                                                                        x["partners"],
                                                                                        x["perturbations"],
                                                                                        mutExclusive), axis=1)
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------
    
    degRange = [1] + list(np.arange(5, 21, 5))
    pN_E_high_degree, pN_E_high, conf_high = [], [], []
    for minDegree in degRange:
        if mono_edgetic:
            numNaturalMut_edgetic = sum((naturalPerturbs["mono-edgotype"] == 'mono-edgetic') &
                                        (naturalPerturbs["perturbed_ppis_max_num_mutExcl"] >= minDegree))
            numDiseaseMut_edgetic = sum((diseasePerturbs["mono-edgotype"] == 'mono-edgetic') &
                                        (diseasePerturbs["perturbed_ppis_max_num_mutExcl"] >= minDegree))
        else:
            numNaturalMut_edgetic = sum((naturalPerturbs["edgotype"] == 'edgetic') & 
                                        (naturalPerturbs["perturbed_ppis_max_num_mutExcl"] >= minDegree))
            numDiseaseMut_edgetic = sum((diseasePerturbs["edgotype"] == 'edgetic') & 
                                        (diseasePerturbs["perturbed_ppis_max_num_mutExcl"] >= minDegree))
    
        numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
        numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic
    
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
            numNaturalMut_edgetic = sum((naturalPerturbs["mono-edgotype"] == 'mono-edgetic') &
                                        (naturalPerturbs["perturbed_ppis_max_num_mutExcl"] < maxDegree))
            numDiseaseMut_edgetic = sum((diseasePerturbs["mono-edgotype"] == 'mono-edgetic') &
                                        (diseasePerturbs["perturbed_ppis_max_num_mutExcl"] < maxDegree))
        else:
            numNaturalMut_edgetic = sum((naturalPerturbs["edgotype"] == 'edgetic') & 
                                        (naturalPerturbs["perturbed_ppis_max_num_mutExcl"] < maxDegree))
            numDiseaseMut_edgetic = sum((diseasePerturbs["edgotype"] == 'edgetic') & 
                                        (diseasePerturbs["perturbed_ppis_max_num_mutExcl"] < maxDegree))
    
        numNaturalMut_nonedgetic = len(naturalPerturbs) - numNaturalMut_edgetic
        numDiseaseMut_nonedgetic = len(diseasePerturbs) - numDiseaseMut_edgetic
    
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
                xlabel = 'Disrupted PPI number of mutually exclusives (c)',
                ylabel = 'Fraction of PPIs neutral upon disruption (%)',
                leg = ['Number of mutually exclusives < c', 'Number of mutually exclusives ≥ c'],
                xticklabels = degRange,
                yMinorTicks = 4,
                yticklabels = list(np.arange(0, maxY + 5, 20)),
                fontsize = 16,
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_vs_numMutExcl%s_%.2f' % ('_monoedgetic' if mono_edgetic else '',
                                                                       overlap_cutoff))

if __name__ == "__main__":
    main()
