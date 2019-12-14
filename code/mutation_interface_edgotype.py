#----------------------------------------------------------------------------------------
# Modules for mutation edgotyping.
#----------------------------------------------------------------------------------------

import io
import numpy as np
import pandas as pd
from random import seed, random
from simple_tools import reverseTuples
from interactome_tools import partners_per_site

def mutation_PPI_interface_perturbations (mutations,
                                          interactome,
                                          maxInterfaces = np.inf,
                                          dist = 0):
    """Predict PPI perturbations by mutations located at PPI interfaces.

    Args:
        mutations (dataframe): mutation table.
        interactome (dataframe): interface-annotated interactome.
        maxInterfaces (int): maximum number of interfaces allowed per interaction partner.
        dist (int): max number of positions in protein sequence a PPI-perturbing 
                    mutation may be from interaction interface. Set to 0 to predict only 
                    mutations located strictly at interface as PPI-perturbing.
    
    Returns:
        list: PPI perturbations per mutation.
    
    """
    perturbations = []
    for protein, mut_pos in mutations[["protein", "mut_position"]].values:
        perturbs = single_mutation_PPI_perturbs(protein,
                                                mut_pos,
                                                interactome,
                                                maxInterfaces = maxInterfaces,
                                                dist = dist)
        perturbations.append(perturbs)
    return perturbations

def single_mutation_PPI_perturbs (protein,
                                  pos,
                                  interactome,
                                  maxInterfaces = np.inf,
                                  dist = 0):
    """Predict PPI perturbations by single mutation.

    Args:
        protein (str): protein ID.
        pos (int): mutation position on protein sequence.
        interactome (dataframe): interface-annotated interactome.
        maxInterfaces (int): maximum number of interfaces allowed per interaction partner.
        dist (int): max number of positions in protein sequence a PPI-perturbing 
                    mutation may be from interaction interface. Set to 0 to predict only 
                    mutations located strictly at interface as PPI-perturbing.
    
    Returns:
        list, array: interaction partners, number of interfaces perturbed per partner.
    
    """
    PPIs = interactome[(interactome[["Protein_1","Protein_2"]]==protein).any(1)
                       & interactome["Interfaces"].apply(lambda x:
                                                         0 < len(x) <= maxInterfaces)
                      ].reset_index(drop=True)
    if not PPIs.empty:
        inCol2 = (PPIs["Protein_2"] == protein)
        if sum(inCol2) > 0:
            PPIs.loc[inCol2, "Protein_2"] = PPIs.loc[inCol2, "Protein_1"].values
            PPIs.loc[inCol2, "Protein_1"] = protein
            PPIs.loc[inCol2, "Interfaces"] = PPIs.loc[inCol2, "Interfaces"].apply(reverseTuples)
        
        PPIs_Perturbed = PPIs.apply(lambda x:
                                    single_mutation_PPI_perturb(x["Protein_1"],
                                                                pos,
                                                                x["Interfaces"],
                                                                dist = dist),
                                    axis=1)
        return PPIs["Protein_2"].tolist(), PPIs_Perturbed.values
    else:
        return [], np.array([])

def single_mutation_PPI_perturb (protein,
                                 pos,
                                 interfaces,
                                 dist = 0):
    """Predict single PPI perturbation by single mutation.

    Args:
        protein (str): protein ID.
        pos (int): mutation position on protein sequence.
        interfaces (list): interaction interfaces.
        dist (int): max number of positions in protein sequence a PPI-perturbing 
                    mutation may be from interaction interface. Set to 0 to predict only 
                    mutations located strictly at interface as PPI-perturbing.
    
    Returns:
        int: number of interfaces perturbed.
    
    """
    leftside = [i for i, j in interfaces]
    mutationNeighbor = range(pos - dist, pos + dist + 1)
    onInterface = []
    for interface in leftside:
        onInterface.append(any([i in mutationNeighbor for i in interface]))
    return sum(onInterface)

def energy_based_perturbation (perturbs, ddg, cutoff, probabilistic = False):
    """Re-predict PPI perturbations for predicted PPI-perturbing mutations using ∆∆G.

    Args:
        perturbs (dataframe): mutation perturbations.
        ddg (dict): mutation ∆∆G values.
        cutoff (numeric): minimum ∆∆G required for PPI perturbation.
        probabilistic (bool): if True, predict perturbations for PPIs with unknown ∆∆Gs 
                                using the probability of perturbation calculated from 
                                known ∆∆Gs.
    
    Returns:
        list, int, int: perturbations,
                        number of PPIs with known ∆∆G,
                        number of PPIs with unknown ∆∆G.
    
    """
    knownDDG = unknownDDG = 0
    pertProb = sum([d > cutoff for _, _, _, _, d in ddg.values()]) / len(ddg)
    seed()
    allPred = []
    for protein, pos, mut_res, partners, perturbations in perturbs[["protein",
                                                                    "mut_position",
                                                                    "mut_res",
                                                                    "partners",
                                                                    "perturbations"]].values:
        pred = []
        for p, pert in zip(partners, perturbations):
            if pert == 0:
                pred.append(0)
            else:
                k = protein, p, pos, mut_res
                if k in ddg:
                    knownDDG += 1
                    pdb_id, chainID, ch_partner, chain_mut, ddgVal = ddg[k]
                    pred.append(1 if ddgVal > cutoff else 0)
                else:
                    unknownDDG += 1
                    if probabilistic:
                        pred.append(1 if random() < pertProb else 0)
                    else:
                        pred.append(np.nan)
        allPred.append(pred)
    return allPred, knownDDG, unknownDDG

def create_perturbed_network (interactome,
                              perturbations,
                              network_outPath,
                              nodeColor_outPath):
    """Create network from perturbed interactome for plotting by Cytoscape software.

    Args:
        interactome (dataframe): interactome with all edges.
        perturbations (dataframe): interactome perturbations.
        network_outPath (Path): file path to save edges for Cytoscape plotting.
        nodeColor_outPath (Path): file path to save node colors for Cytoscape plotting.

    Returns:
        list, list, list, list: nodes, edges, node colors, edge colors.
    
    """
    network = pd.DataFrame()
    edges = [tuple(e) for e in interactome[["Protein_1", "Protein_2"]].values]
    edgeColors = {e:'black' for e in edges}
    
    network["Node_1"], network["Node_2"] = zip(*edges)
    nodes = list(set(network.values.flatten()))
    numMut = {p:0 for p in nodes}
    network["Edge"] = network["Node_1"] + ' (pp) ' + network["Node_2"]
    
    # assign colors to perturbed edges
    for _, row in perturbations.iterrows():
        if row.protein in numMut:
            numMut[row.protein] += 1 
        for partner, perturb in zip(row.partners, row.perturbations):
            if perturb > 0:
                if (row.protein, partner) in edgeColors:
                    edgeColors[(row.protein, partner)] = 'red'
                elif (partner, row.protein) in edgeColors:
                    edgeColors[(partner, row.protein)] = 'red'
    network["Edge_color"] = [edgeColors[e] for e in edges]
    network.to_csv(network_outPath, index=False, sep='\t')

    # assign colors to nodes
    with io.open(nodeColor_outPath, "w") as fout:
        fout.write('Node' + '\t' + 'Color' + '\n')
        nodeColors = []
        for n in nodes:
            if numMut[n] > 9:
                nodeColors.append('blue')
            elif numMut[n] > 2:
                nodeColors.append('mediumpurple')
            elif numMut[n] > 0:
                nodeColors.append('magenta')
            else:
                nodeColors.append('grey')
            fout.write(n + '\t' + nodeColors[-1] + '\n')
    
    return nodes, edges, nodeColors, network["Edge_color"].tolist()

def assign_edgotypes (perturbs, mono_edgetic = False):
    """Assign edgotypes to mutation with PPI perturbation predictions.

    Args:
        perturbs (list): mutation perturbations.
        mono_edgetic (bool): if True, distinguish mono-edgetic mutations from other 
                                edgetic mutations.
    
    Returns:
        list: mutation edgotypes.
    
    """
    edgotypes = []
    for pert in perturbs:
        etype = assign_edgotype (pert, mono_edgetic = mono_edgetic)
        edgotypes.append(etype)
    return edgotypes

def assign_edgotype (perturbs, mono_edgetic = False):
    """Assign edgotype to single mutation with PPI perturbation predictions.

    Args:
        perturbs (list): mutation perturbations.
        mono_edgetic (bool): if True, distinguish mono-edgetic mutations from other 
                                edgetic mutations.
    
    Returns:
        str: mutation edgotype.
    
    """
    pert = []
    for p in perturbs:
        pert.append(p > 0 if not np.isnan(p) else np.nan)
    numPerturbs = pert.count(True)
    if (numPerturbs >= 2) or ((numPerturbs == 1) and not mono_edgetic):
        return 'edgetic'
    elif pert.count(np.nan) > 0:
        return '-'
    elif numPerturbs == 1:
        return 'mono-edgetic'
    else:
        return 'non-edgetic'

def mutation_site_perturbations (ppi_perturbations, interactome):
    """Predict PPI perturbations by mutations located at PPI interfaces.

    Args:
        mutations (dataframe): mutation table.
        interactome (dataframe): interface-annotated interactome.
        maxInterfaces (int): maximum number of interfaces allowed per interaction partner.
        dist (int): max number of positions in protein sequence a PPI-perturbing 
                    mutation may be from interaction interface. Set to 0 to predict only 
                    mutations located strictly at interface as PPI-perturbing.
    
    Returns:
        list: PPI perturbations per mutation.
    
    """
    allsites, allperturbs = [], []
    partnersPerSite = partners_per_site (interactome)
    for protein, partners, perturbs in ppi_perturbations[["protein",
                                                          "partners",
                                                          "perturbations"]].values:
        sites, sitePerturbs = single_mutation_site_perturbations (partnersPerSite[protein],
                                                                  protein,
                                                                  partners,
                                                                  perturbs)
        allsites.append(sites)
        allperturbs.append(sitePerturbs)
    return allsites, allperturbs
    
def single_mutation_site_perturbations (partnersPerSite, protein, partners, perturbs):
    
    perturbed_partners = set([p for p, pert in zip(partners, perturbs) if pert > 0])
    sites, sitePerturbs = [], []
    for site, site_partners in partnersPerSite.items():
        sites.append(site)
        if len(site_partners - perturbed_partners) == 0:
            sitePerturbs.append(1)
        else:
            sitePerturbs.append(0)
    return sites, sitePerturbs

def assign_sitotypes (perturbs, mono_edgetic = False):
    """Assign edgotypes to mutation with PPI perturbation predictions.

    Args:
        perturbs (list): mutation perturbations.
        mono_edgetic (bool): if True, distinguish mono-edgetic mutations from other 
                                edgetic mutations.
    
    Returns:
        list: mutation edgotypes.
    
    """
    edgotypes = []
    for pert in perturbs:
        etype = assign_sitotype (pert, mono_edgetic = mono_edgetic)
        edgotypes.append(etype)
    return edgotypes

def assign_sitotype (perturbs, mono_edgetic = False):
    """Assign edgotype to single mutation with PPI perturbation predictions.

    Args:
        perturbs (list): mutation perturbations.
        mono_edgetic (bool): if True, distinguish mono-edgetic mutations from other 
                                edgetic mutations.
    
    Returns:
        str: mutation edgotype.
    
    """
    pert = []
    for p in perturbs:
        pert.append(p > 0 if not np.isnan(p) else np.nan)
    numPerturbs = pert.count(True)
    if (numPerturbs >= 2) or ((numPerturbs == 1) and not mono_edgetic):
        return 'sitic'
    elif pert.count(np.nan) > 0:
        return '-'
    elif numPerturbs == 1:
        return 'mono-sitic'
    else:
        return 'non-sitic'

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

#     maxdegree = []
#     for _, mut in mutations.iterrows():
#         if num_perturbed_ppis (mut.perturbations) > 0:
#             proteins = {mut.protein} | {p for p, pert in zip(mut.partners, mut.perturbations) if pert > 0}
#             m = max([degree[p] for p in proteins])
#             maxdegree.append(m)
#         else:
#             maxdegree.append(0)        
#     return maxdegree

def num_perturbed_ppis (perturbations):
    
    return sum([p > 0 for p in perturbations if not np.isnan(p)])
