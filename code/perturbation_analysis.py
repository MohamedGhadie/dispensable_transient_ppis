import numpy as np

def perturbed_ppis_max_num_mutExcl (protein, partners, perturbations, simultaneous):
    
    num_mutExcl = perturbed_ppis_num_mutExcl (protein, partners, perturbations, simultaneous)
    if num_mutExcl:
        return max([n for p, n in num_mutExcl])
    else:
        return np.nan

def perturbed_ppis_num_mutExcl (protein, partners, perturbations, simultaneous):
    
    num = num_mutExcl (protein, partners, simultaneous)
    return [(p, n) for n, p, pert in zip(num, partners, perturbations) if pert > 0]

def num_mutExcl (protein, partners, simultaneous):
    
    num_mutexcl = []
    for p in partners:
        if p in simultaneous[protein]:
            num_mutexcl.append(len(simultaneous[protein][p]))
        else:
            num_mutexcl.append(np.nan)
    return num_mutexcl

def perturbed_ppis_max_num_simult (protein, partners, perturbations, simultaneous):
    
    num_simult = perturbed_ppis_num_simult (protein, partners, perturbations, simultaneous)
    if num_simult:
        return max([n for p, n in num_simult])
    else:
        return np.nan

def perturbed_ppis_num_simult (protein, partners, perturbations, simultaneous):
    
    num = num_simult (protein, partners, simultaneous)
    return [(p, n) for n, p, pert in zip(num, partners, perturbations) if pert > 0]

def num_simult (protein, partners, simultaneous):
    
    num_simult = []
    for p in partners:
        if p in simultaneous[protein]:
            num_simult.append(len(simultaneous[protein][p]))
        else:
            num_simult.append(np.nan)
    return num_simult
