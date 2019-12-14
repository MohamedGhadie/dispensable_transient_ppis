#----------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import read_single_interface_annotated_interactome, mutExcl_simult_partners
from stat_tools import fisher_test, sderror_on_fraction
from plot_tools import bar_plot

# reference interactome name
# options: HI-II-14, IntAct
interactome_name = 'IntAct'

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

def main():
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    naturalPerturbsFile = interactomeDir / 'nondisease_mutation_edgotype_geometry.txt'
    diseasePerturbsFile = interactomeDir / 'disease_mutation_edgotype_geometry.txt'
    
    # create directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    naturalPerturbs = read_list_table (naturalPerturbsFile, ["partners", "perturbations"], [str, int])
    diseasePerturbs = read_list_table (diseasePerturbsFile, ["partners", "perturbations"], [str, int])
    
#     naturalPerturbs = naturalPerturbs [naturalPerturbs["edgotype"] == 'edgetic'].reset_index(drop=True)
#     diseasePerturbs = diseasePerturbs [diseasePerturbs["edgotype"] == 'edgetic'].reset_index(drop=True)
    
    interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    mutExclusive, simultaneous = mutExcl_simult_partners (interactome)
    
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a PPI having another simultaneous
    #------------------------------------------------------------------------------------
    
    numNatural_simult_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_simult_perturbs.append(num)
      
    numDisease_simult_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_simult_perturbs.append(num)
    
    print_results ([p > 0 for p in numNatural_simult_perturbs],
                   [p > 0 for p in numDisease_simult_perturbs],
                   title = ['Fraction that disrupt a PPI', 'having another simultaneous'],
                   fig_name = 'frac_mut_disrupt_1>ppi_has_simult',
                   yticks = [0, 0.05, 0.10, 0.15, 0.20, 0.25])
        
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a PPI having mutually exclusives 
    #------------------------------------------------------------------------------------
        
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(mutExclusive[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_mutExcl_perturbs.append(num)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(mutExclusive[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_mutExcl_perturbs.append(num)
    
    print_results ([p > 0 for p in numNatural_mutExcl_perturbs],
                   [p > 0 for p in numDisease_mutExcl_perturbs],
                   title = ['Fraction that disrupt a PPI', 'having mutually exclusives'],
                   fig_name = 'frac_mut_disrupt_1>ppi_has_mutExcl',
                   yticks = [0, 0.2, 0.4, 0.6, 0.8])
    
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a proper subset of mutually exclusive PPIs
    #------------------------------------------------------------------------------------
    
    nat_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                for p2 in mutExclusive[mut.protein][p]:
                    if (p2, 0) in perturbations:
                        found = True
        nat_mutExcl_perturbs.append(found)
    
    dis_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                for p2 in mutExclusive[mut.protein][p]:
                    if (p2, 0) in perturbations:
                        found = True
        dis_mutExcl_perturbs.append(found)
    
    print_results (nat_mutExcl_perturbs,
                   dis_mutExcl_perturbs,
                   title = ['Fraction that disrupt a proper subset', 'of mutually exclusive PPIs'],
                   fig_name = 'frac_mut_disrupt_proper_subset_of_mutExcl',
                   yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5])
    
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a PPI and all its mutually exclusives (if any)
    #------------------------------------------------------------------------------------
    
    nat_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [(p2, 1) in perturbations for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = True
        nat_mutExcl_perturbs.append(found)
    
    dis_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [(p2, 1) in perturbations for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = True
        dis_mutExcl_perturbs.append(found)
    
    print_results (nat_mutExcl_perturbs,
                   dis_mutExcl_perturbs,
                   title = ['Fraction that disrupt a PPI', 'and all mutually exclusives'],
                   fig_name = 'frac_mut_disrupt_ppi_and_all_mutExcl',
                   yticks = [0, 0.2, 0.4, 0.6, 0.8])
    
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a PPI having mutually exclusives but no simultaneous
    #------------------------------------------------------------------------------------
        
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and 
                   (len(mutExclusive[mut.protein][partner]) > 0) and 
                   (len(simultaneous[mut.protein][partner]) == 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_mutExcl_perturbs.append(num)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and 
                   (len(mutExclusive[mut.protein][partner]) > 0) and 
                   (len(simultaneous[mut.protein][partner]) == 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_mutExcl_perturbs.append(num)
    
    print_results ([p > 0 for p in numNatural_mutExcl_perturbs],
                   [p > 0 for p in numDisease_mutExcl_perturbs],
                   title = ['Fraction that disrupt a PPI', 'with mutually exclusives and no simultaneous'],
                   fig_name = 'frac_mut_disrupt_1>ppi_has_mutExcl_no_simult',
                   yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
            
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a PPI and its mutually exclusives (if any),
    # with none having simultaneous PPIs
    #------------------------------------------------------------------------------------
    
    nat_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if (pert > 0) and (len(simultaneous[mut.protein][p]) == 0):
                otherPerturbs = [((p2, 1) in perturbations) and
                                 (len(simultaneous[mut.protein][p2]) == 0)
                                 for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = True
        nat_mutExcl_perturbs.append(found)
    
    dis_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [((p2, 1) in perturbations) and 
                                 (len(simultaneous[mut.protein][p2]) == 0) 
                                 for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = True
        dis_mutExcl_perturbs.append(found)
    
    print_results (nat_mutExcl_perturbs,
                   dis_mutExcl_perturbs,
                   title = ['Fraction that disrupt a PPI', 'and mutually exclusives', 'with none having simultaneous'],
                   fig_name = 'frac_mut_disrupt_1>ppi_and_all_mutExcl_no_simult',
                   yticks = [0, 0.2, 0.4, 0.6, 0.8])
    
    #------------------------------------------------------------------------------------
    # Edgetic mutations that disrupt a PPI and its mutually exclusives (if any),
    # with at least one having a simultaneous PPI
    #------------------------------------------------------------------------------------
    
    nat_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if (pert > 0):
                otherPerturbs = [(p2, 1) in perturbations for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    for p2 in (mutExclusive[mut.protein][p] | {p}):
                        if len(simultaneous[mut.protein][p2]) > 0:
                            found = True
        nat_mutExcl_perturbs.append(found)
    
    dis_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = False
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [(p2, 1) in perturbations for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    for p2 in (mutExclusive[mut.protein][p] | {p}):
                        if len(simultaneous[mut.protein][p2]) > 0:
                            found = True
        dis_mutExcl_perturbs.append(found)
    
    print_results (nat_mutExcl_perturbs,
                   dis_mutExcl_perturbs,
                   title = ['Fraction that disrupt a PPI', 'and mutually exclusives', 'with at least one having simultaneous'],
                   fig_name = 'frac_mut_disrupt_1>ppi_and_all_mutExcl_1>_simult',
                   yticks = [0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12])
    
def print_results (nat_perturbs,
                   dis_perturbs,
                   title = 'Fraction',
                   fig_name = None,
                   yticks = None):
    
    if isinstance(title, str):
        title = [title]
    
    numNatural_perturbs, total_nat = sum(nat_perturbs), len(nat_perturbs)
    numDisease_perturbs, total_dis = sum(dis_perturbs), len(dis_perturbs)
    
    frac_natural_perturbs = numNatural_perturbs / total_nat
    frac_disease_perturbs = numDisease_perturbs / total_dis
    
    sderror_nat = sderror_on_fraction (numNatural_perturbs, total_nat)
    sderror_dis = sderror_on_fraction (numDisease_perturbs, total_dis)
    
    print()
    print(' '.join(title))
    print('Non-disease mutations: %.1f%% (%d out of %d), SE = %g' % (100 * frac_natural_perturbs,
                                                                     numNatural_perturbs,
                                                                     total_nat,
                                                                     sderror_nat))
    print('Disease mutations: %.1f%% (%d out of %d), SE = %g' % (100 * frac_disease_perturbs,
                                                                 numDisease_perturbs,
                                                                 total_dis,
                                                                 sderror_dis))
    
    fisher_test ([numNatural_perturbs, total_nat - numNatural_perturbs],
                 [numDisease_perturbs, total_dis - numDisease_perturbs])
    
    if fig_name:
        bar_plot ([frac_natural_perturbs, frac_disease_perturbs],
                  error = [sderror_nat, sderror_dis],
                  xlabels = ('\n'.join(['Non-disease', 'mutations']),
                             '\n'.join(['Disease', 'mutations'])),
                  ylabel = '\n'.join(title),
                  ylabels = yticks,
                  ylim = [0, max(yticks)] if yticks else None,
                  colors = ('turquoise', 'red'),
                  edgecolor = 'black',
                  ewidth = 2.5,
                  barwidth = 0.6,
                  fontsize = 18,
                  capsize = 10,
                  msize = 26,
                  show = showFigs,
                  figdir = figDir,
                  figname = fig_name)

if __name__ == "__main__":
    main()
