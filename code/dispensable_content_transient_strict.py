#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from text_tools import read_list_table
from perturbation_tools import (unique_perturbation_mutations,
                                num_transient_perturbed_ppis,
                                num_permanent_perturbed_ppis)
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              produce_hpa_expr_dict,
                              produce_fantom5_expr_dict)
from math_tools import fitness_effect
from plot_tools import curve_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # set to True to calculate dispensable PPI content using fraction of mono-edgetic mutations 
    # instead of all edgetic mutations
    mono_edgetic = False

    # set to True to remove mutations that have no unique PPI perturbation
    unique_edgetics = False
    
    # tissue expression database name
    # options: Illumina, GTEx, HPA, Fantom5
    expr_db = 'Fantom5'
    
    # minimum number of tissue expression values required for protein pair tissue
    # co-expression to be considered
    minTissues = 5
    
    # maximum co-expression level (not inclusive) for transient PPIs
    maxCoexpr = 0.1
    
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
    
    pN_E_keys = ['All PPIs', 'Permanent PPIs', 'Transient PPIs']
    
    pN_E, conf = {}, {}
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'strict'
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513.tsv.txt'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    hpaExprFile = extDir / 'normal_tissue.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    naturalMutationsFile = interactomeDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = interactomeDir / 'disease_mutation_edgetics.txt'
    
    # output data files
    proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing protein tissue expression dictionary')
        if expr_db is 'Illumina':
            produce_illumina_expr_dict (illuminaExprFile,
                                        uniprotIDmapFile,
                                        proteinExprFile,
                                        headers = list(range(1, 18)))
        elif expr_db is 'GTEx':
            produce_gtex_expr_dict (gtexDir,
                                    uniprotIDmapFile,
                                    proteinExprFile,
                                    uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'HPA':
            produce_hpa_expr_dict (hpaExprFile,
                                   uniprotIDmapFile,
                                   proteinExprFile)
        elif expr_db is 'Fantom5':
            produce_fantom5_expr_dict (fantomExprFile,
                                       uniprotIDmapFile,
                                       proteinExprFile,
                                       sampleTypes = 'tissues',
                                       sampleTypeFile = fantomSampleTypeFile,
                                       uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    if expr_db is 'HPA':
        exprMap = {'Not detected':0, 'Low':1, 'Medium':2, 'High':3}
        for k, v in expr.items():
            expr[k] = np.array([(exprMap[e] if e in exprMap else np.nan) for e in v])
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
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
    
    naturalMutations ["num_transient_perturbed_ppis"] = naturalMutations.apply(lambda x:
                                                        num_transient_perturbed_ppis (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      expr,
                                                                                      minTissues = minTissues,
                                                                                      maxCoexpr = maxCoexpr), axis=1)
    diseaseMutations ["num_transient_perturbed_ppis"] = diseaseMutations.apply(lambda x:
                                                        num_transient_perturbed_ppis (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      expr,
                                                                                      minTissues = minTissues,
                                                                                      maxCoexpr = maxCoexpr), axis=1)
    
    naturalMutations ["num_permanent_perturbed_ppis"] = naturalMutations.apply(lambda x:
                                                        num_permanent_perturbed_ppis (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      expr,
                                                                                      minTissues = minTissues,
                                                                                      minCoexpr = maxCoexpr), axis=1)
    diseaseMutations ["num_permanent_perturbed_ppis"] = diseaseMutations.apply(lambda x:
                                                        num_permanent_perturbed_ppis (x["protein"],
                                                                                      x["partners"],
                                                                                      x["perturbations"],
                                                                                      expr,
                                                                                      minTissues = minTissues,
                                                                                      minCoexpr = maxCoexpr), axis=1)

    #------------------------------------------------------------------------------------
    # Remove mutations with no unique PPI perturbation
    #------------------------------------------------------------------------------------
    
    if unique_edgetics:
        naturalMutations = naturalMutations [unique_perturbation_mutations (naturalMutations)].reset_index(drop=True)
        diseaseMutations = diseaseMutations [unique_perturbation_mutations (diseaseMutations)].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Dispensable content among all PPIs
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
    
    print()
    print('Fitness effect for disruption among all PPIs:')
    
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
        pN_E[pN_E_keys[0]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[0]] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among permanent PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (naturalMutations["num_permanent_perturbed_ppis"] > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (diseaseMutations["num_permanent_perturbed_ppis"] > 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') &
                                    (naturalMutations["num_permanent_perturbed_ppis"] > 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') &
                                    (diseaseMutations["num_permanent_perturbed_ppis"] > 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('Fitness effect for disruption of permanent PPIs:')
    
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
        pN_E[pN_E_keys[1]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[1]] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Dispensable content among transient PPIs
    #------------------------------------------------------------------------------------
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum((naturalMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (naturalMutations["num_transient_perturbed_ppis"] > 0) &
                                    (naturalMutations["num_permanent_perturbed_ppis"] == 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["mono-edgotype"] == 'mono-edgetic') &
                                    (diseaseMutations["num_transient_perturbed_ppis"] > 0) &
                                    (diseaseMutations["num_permanent_perturbed_ppis"] == 0))
    else:
        numNaturalMut_edgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                    (naturalMutations["num_transient_perturbed_ppis"] > 0) &
                                    (naturalMutations["num_permanent_perturbed_ppis"] == 0))
        numDiseaseMut_edgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                    (diseaseMutations["num_transient_perturbed_ppis"] > 0) &
                                    (diseaseMutations["num_permanent_perturbed_ppis"] == 0))

    numNaturalMut_nonedgetic = len(naturalMutations) - numNaturalMut_edgetic
    numDiseaseMut_nonedgetic = len(diseaseMutations) - numDiseaseMut_edgetic

    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    print()
    print('Fitness effect for disruption of transient PPIs:')
    print(numNaturalMut_edgetic)
    print(numNaturalMut_considered)
    print(numDiseaseMut_edgetic)
    print(numDiseaseMut_considered)
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
        pN_E[pN_E_keys[2]] = 100 * all_effects['P(N|E)']
        if 'P(N|E)_CI' in all_effects:
            lower, upper = all_effects['P(N|E)_CI']
            conf[pN_E_keys[2]] = 100 * lower, 100 * upper
    
    #------------------------------------------------------------------------------------
    # Plot dispensable PPI content
    #------------------------------------------------------------------------------------
    
    numGroups = len(pN_E.keys())
    numGroups = 3
    if computeConfidenceIntervals:
        maxY = max([pN_E[p] + conf[p][1] for p in pN_E.keys()])
    else:
        maxY = max(pN_E.values())
    maxY = 10 * np.ceil(maxY / 10)
    
    curve_plot ([(pN_E[p] if p in pN_E else np.nan) for p in pN_E_keys],
                error = [(conf[p] if p in pN_E else [np.nan, np.nan]) for p in pN_E_keys] if computeConfidenceIntervals else None,
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
                xticklabels = pN_E_keys,
                yticklabels = list(np.arange(0, maxY + 10, 10)),
                fontsize = 18,
                adjustBottom = 0.2,
                shiftBottomAxis = -0.1,
                xbounds = (1, numGroups),
                show = showFigs,
                figdir = figDir,
                figname = 'Fraction_disp_PPIs_transient_%s' % expr_db)

if __name__ == "__main__":
    main()
