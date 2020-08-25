import numpy as np
from protein_function import is_transient

def unique_perturbation_mutations (mutations):
    
    perturbs = {p:{} for p in set(mutations["protein"].values)}
    for _, mut in mutations.iterrows():
        mutPerturbs = set([p for p, pert in zip(mut.partners, mut.perturbations) if pert > 0])
        perturbs[mut.protein][(mut.mut_position, mut.mut_res)] = mutPerturbs

    uniques = []
    for _, mut in mutations.iterrows():
        keep = True
        otherPerturbs = set()
        for k, val in perturbs[mut.protein].items():
            if k != (mut.mut_position, mut.mut_res):
                otherPerturbs.update(val)
        mutPerturbs = perturbs[mut.protein][(mut.mut_position, mut.mut_res)]
        if mutPerturbs:
            if len(mutPerturbs - otherPerturbs) == 0:
                keep = False
                del perturbs[mut.protein][(mut.mut_position, mut.mut_res)]
        uniques.append(keep)
    return uniques

def num_hub_ppis_perturbed (perturbs, hubPPIs):
    
    return sum([hub == 'hub PPI' for pert, hub in zip(perturbs, hubPPIs) if pert > 0])

def num_nonhub_ppis_perturbed (perturbs, hubPPIs):
    
    return sum([hub == 'non-hub PPI' for pert, hub in zip(perturbs, hubPPIs) if pert > 0])

def perturbed_partner_max_degree (protein, partners, perturbs, degree):
    
    d = perturbed_partner_degrees (protein, partners, perturbs, degree)
    if d:
        return max(d.values())
    else:
        return np.nan

def perturbed_partner_degrees (protein, partners, perturbs, degree):
    
    if num_perturbed_ppis (perturbs) > 0:
        proteins = {protein} | {p for p, pert in zip(partners, perturbs) if pert > 0}
        return {p:degree[p] for p in proteins if p in degree}
    else:
        return None

def num_transient_ppis_perturbed (perturbs, transients):
    
    return sum([trans == 'transient' for pert, trans in zip(perturbs, transients) if pert > 0])

def num_permanent_ppis_perturbed (perturbs, transients):
    
    return sum([trans == 'permanent' for pert, trans in zip(perturbs, transients) if pert > 0])

# def num_transient_ppis_perturbed (protein,
#                                   partners,
#                                   perturbs,
#                                   expr,
#                                   minTissues = 3,
#                                   maxCoexpr = 0.05):
#     
#     num = 0
#     for p, pert in zip(partners, perturbs):
#         if pert > 0:
#             if is_transient (protein, p, expr, minTissues = minTissues, maxCoexpr = maxCoexpr):
#                 num += 1
#     return num
# 
# def num_permanent_perturbed_ppis (protein,
#                                   partners,
#                                   perturbs,
#                                   expr,
#                                   minTissues = 3,
#                                   minCoexpr = 0.05):
#     
#     num = 0
#     for p, pert in zip(partners, perturbs):
#         if pert > 0:
#             if is_permanent (protein, p, expr, minTissues = minTissues, minCoexpr = minCoexpr):
#                 num += 1
#     return num

def num_perturbed_ppis (perturbations):
    
    return sum([p > 0 for p in perturbations if not np.isnan(p)])

def perturbed_ppis_max_num_mutExcl (protein, partners, perturbations, simultaneous):
    
    num_mutExcl = perturbed_ppis_num_mutExcl (protein, partners, perturbations, simultaneous)
    if num_mutExcl:
        return max([n for p, n in num_mutExcl])
    else:
        return np.nan

def num_perturbed_ppis_n_mutExcl (perturbs, numMutExcl, minN = 1, maxN = np.inf):
    
    return sum([minN <= n <= maxN for p, n in zip(perturbs, numMutExcl) if p > 0])

def num_perturbed_ppis_n_simult (perturbs, numSimult, minN = 1, maxN = np.inf):
    
    return sum([minN <= n <= maxN for p, n in zip(perturbs, numSimult) if p > 0])

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

def num_weak_ppis_perturbed (perturbs, strength):
    
    return sum([s == 'weak' for p, s in zip(perturbs, strength) if p > 0])

def num_strong_ppis_perturbed (perturbs, strength):
    
    return sum([s == 'strong' for p, s in zip(perturbs, strength) if p > 0])

def num_unbalanced_ppis_perturbed (perturbs, balance):
    
    return sum([b == 'unbalanced' for p, b in zip(perturbs, balance) if p > 0])

def num_balanced_ppis_perturbed (perturbs, balance):
    
    return sum([b == 'balanced' for p, b in zip(perturbs, balance) if p > 0])
