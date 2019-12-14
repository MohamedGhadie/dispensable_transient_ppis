#----------------------------------------------------------------------------------------
# 
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from text_tools import read_list_table
from interactome_tools import read_single_interface_annotated_interactome, mutExcl_simult_partners
from stat_tools import fisher_test

def main():
    
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
    
    naturalPerturbs = naturalPerturbs [naturalPerturbs["edgotype"] == 'edgetic'].reset_index(drop=True)
    diseasePerturbs = diseasePerturbs [diseasePerturbs["edgotype"] == 'edgetic'].reset_index(drop=True)
    
    interactome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    mutExclusive, simultaneous = mutExcl_simult_partners (interactome)
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt at least one PPI that has another 
    # simultaneous PPI
    #------------------------------------------------------------------------------------
    
    numNatural_simult_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_simult_perturbs.append(num)
    numNatural_simult_perturbs = np.array(numNatural_simult_perturbs)
      
    numDisease_simult_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_simult_perturbs.append(num)
    numDisease_simult_perturbs = np.array(numDisease_simult_perturbs)
    
    print()
    print('Fraction of edgetic mutations disrupting at least 1 simultaneous PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_simult_perturbs > 0) / len(numNatural_simult_perturbs),
             sum(numNatural_simult_perturbs > 0),
             len(numNatural_simult_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_simult_perturbs > 0) / len(numDisease_simult_perturbs),
             sum(numDisease_simult_perturbs > 0),
             len(numDisease_simult_perturbs)))
    fisher_test ([sum(numNatural_simult_perturbs > 0), sum(numNatural_simult_perturbs == 0)],
                 [sum(numDisease_simult_perturbs > 0), sum(numDisease_simult_perturbs == 0)])
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt at least 1 PPI that has no other 
    # simultaneous PPI
    #------------------------------------------------------------------------------------
    
    numNatural_perturbs_nosimult = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) == 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_perturbs_nosimult.append(num)
    numNatural_perturbs_nosimult = np.array(numNatural_perturbs_nosimult)
    
    numDisease_perturbs_nosimult = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) == 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_perturbs_nosimult.append(num)
    numDisease_perturbs_nosimult = np.array(numDisease_perturbs_nosimult)
    
    print()
    print('Fraction of edgetic mutations disrupting at least 1 PPI that has no other simultaneous PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_perturbs_nosimult > 0) / len(numNatural_perturbs_nosimult),
             sum(numNatural_perturbs_nosimult > 0),
             len(numNatural_perturbs_nosimult)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_perturbs_nosimult > 0) / len(numDisease_perturbs_nosimult),
             sum(numDisease_perturbs_nosimult > 0),
             len(numDisease_perturbs_nosimult)))
    fisher_test ([sum(numNatural_perturbs_nosimult > 0), sum(numNatural_perturbs_nosimult == 0)],
                 [sum(numDisease_perturbs_nosimult > 0), sum(numDisease_perturbs_nosimult == 0)])
        
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt at least one PPI among multiple 
    # mutually exclusive PPIs 
    #------------------------------------------------------------------------------------
        
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(mutExclusive[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_mutExcl_perturbs.append(num)
    numNatural_mutExcl_perturbs = np.array(numNatural_mutExcl_perturbs)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(mutExclusive[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_mutExcl_perturbs.append(num)
    numDisease_mutExcl_perturbs = np.array(numDisease_mutExcl_perturbs)
    
    print()
    print('Fraction of edgetic mutations disrupting at least 1 mutually exclusive PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_mutExcl_perturbs > 0) / len(numNatural_mutExcl_perturbs),
             sum(numNatural_mutExcl_perturbs > 0),
             len(numNatural_mutExcl_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_mutExcl_perturbs > 0) / len(numDisease_mutExcl_perturbs),
             sum(numDisease_mutExcl_perturbs > 0),
             len(numDisease_mutExcl_perturbs)))
    fisher_test ([sum(numNatural_mutExcl_perturbs > 0), sum(numNatural_mutExcl_perturbs == 0)],
                 [sum(numDisease_mutExcl_perturbs > 0), sum(numDisease_mutExcl_perturbs == 0)])
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt at least one PPI among multiple 
    # mutually exclusive PPIs that has no other simultaneous PPI
    #------------------------------------------------------------------------------------
        
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and 
                   (len(mutExclusive[mut.protein][partner]) > 0) and 
                   (len(simultaneous[mut.protein][partner]) == 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_mutExcl_perturbs.append(num)
    numNatural_mutExcl_perturbs = np.array(numNatural_mutExcl_perturbs)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and 
                   (len(mutExclusive[mut.protein][partner]) > 0) and 
                   (len(simultaneous[mut.protein][partner]) == 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_mutExcl_perturbs.append(num)
    numDisease_mutExcl_perturbs = np.array(numDisease_mutExcl_perturbs)
    
    print()
    print('Fraction of edgetic mutations disrupting at least 1 mutually exclusive PPI' + 
          ' that has no other simultaneous PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_mutExcl_perturbs > 0) / len(numNatural_mutExcl_perturbs),
             sum(numNatural_mutExcl_perturbs > 0),
             len(numNatural_mutExcl_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_mutExcl_perturbs > 0) / len(numDisease_mutExcl_perturbs),
             sum(numDisease_mutExcl_perturbs > 0),
             len(numDisease_mutExcl_perturbs)))
    fisher_test ([sum(numNatural_mutExcl_perturbs > 0), sum(numNatural_mutExcl_perturbs == 0)],
                 [sum(numDisease_mutExcl_perturbs > 0), sum(numDisease_mutExcl_perturbs == 0)])
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt at least one but not all mutually 
    # exclusive PPIs 
    #------------------------------------------------------------------------------------
    
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = 0
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                for p2 in mutExclusive[mut.protein][p]:
                    if (p2, 0) in perturbations:
                        found = 1
        numNatural_mutExcl_perturbs.append(found)
    numNatural_mutExcl_perturbs = np.array(numNatural_mutExcl_perturbs)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = 0
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                for p2 in mutExclusive[mut.protein][p]:
                    if (p2, 0) in perturbations:
                        found = 1
        numDisease_mutExcl_perturbs.append(found)
    numDisease_mutExcl_perturbs = np.array(numDisease_mutExcl_perturbs)
    
    print()
    print('Fraction of edgetic mutations disrupting at least 1 but not all mutually exclusive PPIs:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_mutExcl_perturbs > 0) / len(numNatural_mutExcl_perturbs),
             sum(numNatural_mutExcl_perturbs > 0),
             len(numNatural_mutExcl_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_mutExcl_perturbs > 0) / len(numDisease_mutExcl_perturbs),
             sum(numDisease_mutExcl_perturbs > 0),
             len(numDisease_mutExcl_perturbs)))
    fisher_test ([sum(numNatural_mutExcl_perturbs > 0), sum(numNatural_mutExcl_perturbs == 0)],
                 [sum(numDisease_mutExcl_perturbs > 0), sum(numDisease_mutExcl_perturbs == 0)])
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt all mutually exclusive PPIs
    #------------------------------------------------------------------------------------
    
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = 0
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [(p2, 1) in perturbations for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = 1
        numNatural_mutExcl_perturbs.append(found)
    numNatural_mutExcl_perturbs = np.array(numNatural_mutExcl_perturbs)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = 0
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [(p2, 1) in perturbations for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = 1
        numDisease_mutExcl_perturbs.append(found)
    numDisease_mutExcl_perturbs = np.array(numDisease_mutExcl_perturbs)
    
    print()
    print('Fraction of edgetic mutations disrupting all mutually exclusive PPIs:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_mutExcl_perturbs > 0) / len(numNatural_mutExcl_perturbs),
             sum(numNatural_mutExcl_perturbs > 0),
             len(numNatural_mutExcl_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_mutExcl_perturbs > 0) / len(numDisease_mutExcl_perturbs),
             sum(numDisease_mutExcl_perturbs > 0),
             len(numDisease_mutExcl_perturbs)))
    fisher_test ([sum(numNatural_mutExcl_perturbs > 0), sum(numNatural_mutExcl_perturbs == 0)],
                 [sum(numDisease_mutExcl_perturbs > 0), sum(numDisease_mutExcl_perturbs == 0)])
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations that disrupt all mutually exclusive PPIs and have no
    # other simultaneous PPI
    #------------------------------------------------------------------------------------
    
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        found = 0
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if (pert > 0) and (len(simultaneous[mut.protein][p]) == 0):
                otherPerturbs = [((p2, 1) in perturbations) and
                                 (len(simultaneous[mut.protein][p2]) == 0)
                                 for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = 1
        numNatural_mutExcl_perturbs.append(found)
    numNatural_mutExcl_perturbs = np.array(numNatural_mutExcl_perturbs)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        found = 0
        perturbations = list(zip(mut.partners, mut.perturbations))
        for p, pert in perturbations:
            if pert > 0:
                otherPerturbs = [((p2, 1) in perturbations) and 
                                 (len(simultaneous[mut.protein][p2]) == 0) 
                                 for p2 in mutExclusive[mut.protein][p]]
                if sum(otherPerturbs) == len(otherPerturbs):
                    found = 1
        numDisease_mutExcl_perturbs.append(found)
    numDisease_mutExcl_perturbs = np.array(numDisease_mutExcl_perturbs)
    
    print()
    print('Fraction of edgetic mutations disrupting all mutually exclusive PPIs' + 
          ' and have no other simultaneous PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_mutExcl_perturbs > 0) / len(numNatural_mutExcl_perturbs),
             sum(numNatural_mutExcl_perturbs > 0),
             len(numNatural_mutExcl_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_mutExcl_perturbs > 0) / len(numDisease_mutExcl_perturbs),
             sum(numDisease_mutExcl_perturbs > 0),
             len(numDisease_mutExcl_perturbs)))
    fisher_test ([sum(numNatural_mutExcl_perturbs > 0), sum(numNatural_mutExcl_perturbs == 0)],
                 [sum(numDisease_mutExcl_perturbs > 0), sum(numDisease_mutExcl_perturbs == 0)])
    
    
    
    return
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations on proteins with multiple partners that disrupt at 
    # least one PPI among multiple mutually exclusive PPIs 
    #------------------------------------------------------------------------------------
    
    naturalPerturbs = naturalPerturbs [naturalPerturbs["partners"].apply(len) > 1].reset_index(drop=True)
    diseasePerturbs = diseasePerturbs [diseasePerturbs["partners"].apply(len) > 1].reset_index(drop=True)
    
    numNatural_mutExcl_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(mutExclusive[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_mutExcl_perturbs.append(num)
    numNatural_mutExcl_perturbs = np.array(numNatural_mutExcl_perturbs)
    
    numDisease_mutExcl_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(mutExclusive[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_mutExcl_perturbs.append(num)
    numDisease_mutExcl_perturbs = np.array(numDisease_mutExcl_perturbs)
    
    print()
    print('Fraction of edgetic mutations on proteins with multiple PPIs disrupting at least 1 mutually exclusive PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_mutExcl_perturbs > 0) / len(numNatural_mutExcl_perturbs),
             sum(numNatural_mutExcl_perturbs > 0),
             len(numNatural_mutExcl_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_mutExcl_perturbs > 0) / len(numDisease_mutExcl_perturbs),
             sum(numDisease_mutExcl_perturbs > 0),
             len(numDisease_mutExcl_perturbs)))
    fisher_test ([sum(numNatural_mutExcl_perturbs > 0), sum(numNatural_mutExcl_perturbs == 0)],
                 [sum(numDisease_mutExcl_perturbs > 0), sum(numDisease_mutExcl_perturbs == 0)])
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations on proteins with multiple partners that disrupt at 
    # least one PPI among multiple simultaneous PPIs 
    #------------------------------------------------------------------------------------
    
    numNatural_simult_perturbs = []
    for _, mut in naturalPerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numNatural_simult_perturbs.append(num)
    numNatural_simult_perturbs = np.array(numNatural_simult_perturbs)
    
    numDisease_simult_perturbs = []
    for _, mut in diseasePerturbs.iterrows():
        num = sum([(pert > 0) and (len(simultaneous[mut.protein][partner]) > 0)
                   for partner, pert in zip(mut.partners, mut.perturbations)])
        numDisease_simult_perturbs.append(num)
    numDisease_simult_perturbs = np.array(numDisease_simult_perturbs)
    
    print()
    print('Fraction of edgetic mutations on proteins with multiple PPIs disrupting at least 1 simultaneous PPI:')
    print('Non-disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numNatural_simult_perturbs > 0) / len(numNatural_simult_perturbs),
             sum(numNatural_simult_perturbs > 0),
             len(numNatural_simult_perturbs)))
    print('Disease mutations: %.1f%% (%d out of %d)' 
          % (100 * sum(numDisease_simult_perturbs > 0) / len(numDisease_simult_perturbs),
             sum(numDisease_simult_perturbs > 0),
             len(numDisease_simult_perturbs)))
    fisher_test ([sum(numNatural_simult_perturbs > 0), sum(numNatural_simult_perturbs == 0)],
                 [sum(numDisease_simult_perturbs > 0), sum(numDisease_simult_perturbs == 0)])

if __name__ == "__main__":
    main()
