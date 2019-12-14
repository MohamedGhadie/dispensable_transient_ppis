#----------------------------------------------------------------------------------------
# Modules for processing interactomes.
#----------------------------------------------------------------------------------------

import pickle
import pandas as pd
import numpy as np
import networkx as nx
from text_tools import write_fasta_file
from simple_tools import str_to_tuples, str_to_list, list_to_str, sequence_overlap

def remove_interactions_reported_once (interactions):
    """Remove protein-protein interactions occuring only one time.

    Args:
        interactions (DataFrame):  protein-protein interactions.
    
    Returns:
        DataFrame
        
    """
    lessthan = interactions["Protein_2"] < interactions["Protein_1"]
    list1 = list(interactions.loc[lessthan, "Protein_1"].values)
    list2 = list(interactions.loc[lessthan, "Protein_2"].values)
    interactions.loc[lessthan, "Protein_1"] = list2
    interactions.loc[lessthan, "Protein_2"] = list1
    duplicates = interactions.duplicated(subset=["Protein_1", "Protein_2"], keep=False)
    interactions.loc[lessthan, "Protein_1"] = list1
    interactions.loc[lessthan, "Protein_2"] = list2
    duplicateInteractions = interactions[duplicates].reset_index(drop=True)
    return duplicateInteractions
    
def duplicated_PPIs (PPIs, minCount = 2):
    """Return protein-protein interactions occuring at least <minCount> times.

    Args:
        PPIs (DataFrame): protein-protein interactions.
        minCount (int): minimum number of occurrences for retained PPIs.
    
    Returns:
        DataFrame
        
    """
    PPIstrings = PPIs.apply(lambda x: '-'.join(sorted([ x["Protein_1"], x["Protein_2"]])), axis=1)
    counts = PPIstrings.apply(lambda x: sum(PPIstrings == x))
    return PPIs[counts >= minCount].reset_index(drop=True)

def remove_duplicate_PPIs (PPIs):
    """Remove duplicate protein-protein interactions (PPIs), keeping first occurence.

    Args:
        PPIs (DataFrame): protein-protein interactions.
    
    Returns:
        DataFrame
        
    """
    PPIstrings = PPIs.apply(lambda x: '-'.join(sorted([x["Protein_1"], x["Protein_2"]])), axis=1)
    duplicate = PPIstrings.duplicated()
    return PPIs[duplicate == False].reset_index(drop=True)

def write_interactome_sequences (inPath, sequenceFile, outPath):
    """Write interactome protein sequences to fasta file.

    Args:
        inPath (Path): path to file containing interactome.
        sequenceFile (Path): path to file containing protein sequences.
        outPath(Path): path to write protein sequences to.

    """
    interactome = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    proteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    interactomeSequences = pd.DataFrame(index=range(len(proteins)),
                                        columns=['ID', 'Sequence'])
    interactomeSequences["ID"] = proteins
    interactomeSequences["Sequence"] = interactomeSequences["ID"].apply(lambda x:
                                                                        sequences.loc[sequences["ID"]==x,
                                                                                      "Sequence"].item())
    write_fasta_file (interactomeSequences, "ID", "Sequence", outPath)

def write_chain_annotated_interactome_to_excel (interactome,
                                                outPath,
                                                sheet_name = 'Sheet1'):
    """Write interactome annotated with chain pairs to excel spreadsheet.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.
        sheet_name (str): name of spreadsheet where interactome will be written.

    """
    write_chain_annotated_interactome (interactome,
                                       outPath,
                                       toExcel = True,
                                       sheet_name = sheet_name )

def write_chain_annotated_interactome (interactome,
                                       outPath,
                                       toExcel = False,
                                       sheet_name = 'Sheet1' ):
    """Write interactome annotated with chain pairs to file.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.
        toExcel (bool): write to excel spreadsheet, otherwise to tab-delimited text file.
        sheet_name (str): name of spreadsheet where interactome will be written. Only
                            relevent if toExcel is True.

    """
    interactomeCp = interactome.copy()
    interactomeCp["Mapping_chains"] = interactomeCp["Mapping_chains"].apply(lambda x:
                                                                            ', '.join(['-'.join(pair)
                                                                                       for pair in x]))
    if toExcel:
        writer = pd.ExcelWriter(str(outPath))
        interactomeCp.to_excel(writer, sheet_name)
        writer.save()
    else:
        interactomeCp.to_csv (outPath, index=False, sep='\t')

def read_chain_annotated_interactome (inPath):
    """Read interactome annotated with chain pairs from a file.

    Args:
        inPath (Path): path to file containing chain-annotated interactome.
    
    Returns:
        DataFrame

    """
    interactome = pd.read_table(inPath, sep='\t')
    interactome["Mapping_chains"] = interactome["Mapping_chains"].apply( str_to_tuples )
    return interactome

def write_unmerged_interface_annotated_interactome (interactome, outPath):
    """Write interactome annotated with multiple interfaces per interaction to file.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.

    """
    interactomeCp = interactome.copy()
    if "Mapping_frac" in interactomeCp.columns:
        interactomeCp["Mapping_frac"] = interactomeCp["Mapping_frac"].apply(lambda x:
                                                                            list_to_str(x, ['|', '+', '/']))
    
    write_interface_annotated_interactome (interactomeCp,
                                           outPath,
                                           delm = ['|', '+', '/', ','],
                                           chain_delm = ['|', '+', '/'])

def read_unmerged_interface_annotated_interactome (inPath):
    """Read interactome annotated with multiple interfaces per interaction from file.

    Args:
        inPath (Path): file directory containing interface-annotated interactome.
    
    Returns:
        DataFrame

    """
    delm = ['|', '+', '/', ',']
    interactome = read_interface_annotated_interactome ( inPath, delm = delm)
    if "Mapping_frac" in interactome.columns:
        interactome["Mapping_frac"] = interactome["Mapping_frac"].apply(lambda x: 
                                                                        str_to_list(x, delm[ : -1 ], float))
    return interactome

def write_single_interface_annotated_interactome_to_excel (interactome,
                                                           outPath,
                                                           sheet_name = 'Sheet1'):
    """Write interactome annotated with a single interface per interaction to excel spreadsheet.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.
        sheet_name (str): name of spreadsheet where interactome will be written.

    """
    write_interface_annotated_interactome_to_excel (interactome,
                                                    outPath,
                                                    delm = ['+', ','],
                                                    sheet_name = sheet_name )

def write_single_interface_annotated_interactome (interactome, outPath):
    """Write interactome annotated with a single interface per interaction to file.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.

    """
    write_interface_annotated_interactome (interactome,
                                           outPath,
                                           delm = ['+', ','] )

def read_single_interface_annotated_interactome (inPath):
    """Read interactome annotated with a single interface per interaction from file.

    Args:
        inPath (Path): path to file containing interface-annotated interactome.
    
    Returns:
        DataFrame

    """
    interactome = read_interface_annotated_interactome ( inPath )
    interactome["Interfaces"] = interactome["Interfaces"].apply( lambda x: x[ 0 ] )
    return interactome

def write_interface_annotated_interactome_to_excel (interactome,
                                                    outPath,
                                                    delm = None,
                                                    sheet_name = 'Sheet1'):
    """Write interactome annotated with interfaces to excel spreadsheet.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.
        delm (list): delimiters used to separate interfaces before writing.
        sheet_name (str): name of spreadsheet where interactome will be written.

    """
    write_interface_annotated_interactome (interactome,
                                           outPath,
                                           delm = delm,
                                           toExcel = True,
                                           sheet_name = sheet_name)

def write_interface_annotated_interactome (interactome,
                                           outPath,
                                           delm = None,
                                           chain_delm = None,
                                           toExcel = False,
                                           sheet_name = 'Sheet1'):
    """Write interactome annotated with interfaces to file.

    Args:
        interactome (DataFrame): interactome to be written to file.
        outPath (Path): output file path to write interactome to.
        delm (list): delimiters used to separate interface annotations before writing.
        chain_delm (list): delimiters used to separate chain-pair annotations before writing.
        toExcel (bool): write to excel spreadsheet, otherwise to tab-delimited text file.
        sheet_name (str): name of spreadsheet where interactome will be written. Only
                            relevent if toExcel is True.

    """
    if not delm:
        delm = ['|', '+', ',']
    if not chain_delm:
        chain_delm = ['|', '+']
    interactomeCp = interactome.copy()
    interactomeCp["Interfaces"] = interactomeCp["Interfaces"].apply(lambda x: list_to_str(x, delm))
    if "Chain_pairs" in interactomeCp.columns:
        interactomeCp["Chain_pairs"] = interactomeCp["Chain_pairs"].apply(lambda x:
                                                                          list_to_str(x, chain_delm)
                                                                          if type(x) == list
                                                                          else x )
    if toExcel:
        writer = pd.ExcelWriter(str(outPath))
        interactomeCp.to_excel(writer, sheet_name)
        writer.save()
    else:
        interactomeCp.to_csv(outPath, index=False, sep='\t')

def read_interface_annotated_interactome (inPath, delm = None):
    """Read interactome annotated with interfaces from a file using a list of specified delimiters.

    Args:
        inPath (Path): file directory containing interface-annotated interactome.
        delm (list): delimiters used to split interface annotations after reading.
    
    Returns:
        DataFrame

    """
    if not delm:
        delm = ['|', '+', ',']
    interactome = pd.read_table(inPath, sep='\t')
    interactome["Interfaces"] = interactome["Interfaces"].apply(lambda x: str_to_list(x, delm, int))
    if "Chain_pairs" in interactome.columns:
        interactome["Chain_pairs"] = interactome["Chain_pairs"].apply(lambda x: str_to_list(x, delm[:-1], str))
    return interactome

def get_interface (interactome, p1, p2):
    
    ppi = interactome.loc[((interactome["Protein_1"] == p1) & (interactome["Protein_2"] == p2)) |
                          ((interactome["Protein_1"] == p2) & (interactome["Protein_2"] == p1))]
    if not ppi.empty:
        protein_1, protein_2, interface = ppi[["Protein_1", "Protein_2", "Interfaces"]].values[0]
    if protein_1 == p1:
        return interface
    else:
        return tuple(reversed(interface))

def create_ppi_graph (inPath):
    """Create a networkx graph from an interactome.

    Args:
        inPath (Path): path to dictionary of protein interaction partners.

    Returns:
        Networkx Graph

    """
    with open(inPath, 'rb') as f:
        partners = pickle.load(f)
    g = nx.Graph()
    proteins = list(partners.keys())
    g.add_nodes_from(proteins)
    for p in proteins:
        neighbors = partners[p]
        for nb in neighbors:
            g.add_edge(p, nb)
    return g

def find_all_shortest_paths (g, outPath):
    """Find shortest path between each pair of nodes in a graph.

    Args:
        g (Graph): Networkx graph.
        outPath (Path): path to save dictionary of shortest paths.

    """
    paths = dict(nx.all_pairs_shortest_path(g))
    with open(outPath, 'wb') as f:
        pickle.dump(paths, f)

def between (paths, n):
    """Calculate the fraction of paths in a list that node n falls on.

    Args:
        paths (list): paths to check. Each path is a list of node IDs.
        n: node ID.

    """
    subPaths = [p for p in paths if (p[0] != n) and (p[-1] != n)]
    return sum([n in p for p in subPaths]) / len(subPaths)
    
def shortest_path (p1, p2, partners, path = None):
    """Find the shortest path between two nodes.

    Args:
        p1: ID of first node.
        p2: ID of second node.
        partners (dict): node neighbors.
        path (list): starting path (used for recursive calls). Keep as default value.

    """
    if not path:
        path = []
    if p1 == p2:
        return [p1]
    else:
        minPath = []
        minLen = np.inf
        neighbors = partners[p1]
        for nb in neighbors:
            if nb not in path:
                pathExtension = [p1] + shortest_path(nb, p2, partners, path = path + [p1])
                if 1 < len(pathExtension) < minLen:
                    minPath = pathExtension
                    minLen = len(minPath)
        return minPath

def produce_protein_interaction_dict (inPath, outPath):
    """Make a dictionary of protein interaction partners.

    Args:
        inPath (Path): path to tab-delimited file containing protein-protein interactions.
        outPath (Path): path to save output pickle file to.

    """ 
    interactome = pd.read_table(inPath, sep="\t")
    partners = interaction_partners (interactome)
    with open(outPath, 'wb') as fOut:
        pickle.dump(partners, fOut)

def interaction_partners (interactome):
    
    proteins =  set(interactome[["Protein_1", "Protein_2"]].values.flatten())
    proteinPartners = {}
    for protein in proteins:
        partners = set(interactome.loc[(interactome[["Protein_1", "Protein_2"]]==protein).any(1),
                                ["Protein_1", "Protein_2"]].values.flatten()) - {protein}
        if sum((interactome[["Protein_1", "Protein_2"]]==protein).all(1)) > 0:
            partners.add(protein)
        proteinPartners[protein] = partners
    return proteinPartners

def partners_per_site (interactome):
    
    partners = {p:{} for p in set(interactome[["Protein_1", "Protein_2"]].values.flatten())}
    for protein_1, protein_2, site_1, site_2 in interactome[["Protein_1",
                                                             "Protein_2",
                                                             "Site_1",
                                                             "Site_2"]].values:
        if site_1 in partners[protein_1]:
            partners[protein_1][site_1].add(protein_2)
        else:
            partners[protein_1][site_1] = {protein_2}
        
        if site_2 in partners[protein_2]:
            partners[protein_2][site_2].add(protein_1)
        else:
            partners[protein_2][site_2] = {protein_1}
    return partners

def num_partners (interactome, colnames = None):
    
    if not colnames:
        colnames = ["Protein_1", "Protein_2"]
    partners = {p:set() for p in set(interactome[colnames].values.flatten())}
    for p1, p2 in interactome[colnames].values:
        partners[p1].add(p2)
        partners[p2].add(p1)
    return {p:len(pr) for p, pr in partners.items()}

def num_sites (interactome):
    
    sites = {p:set() for p in set(interactome[["Protein_1", "Protein_2"]].values.flatten())}
    for _, ppi in interactome.iterrows():
        for protein, site in ((ppi.Protein_1, ppi.Site_1), (ppi.Protein_2, ppi.Site_2)):
            if site not in sites[protein]:
                sites[protein].add(site)
    return {p:len(s) for p, s in sites.items()}

def mutExcl_simult_partners (interactome, cutoff = 0):
    
    partners = interaction_partners (interactome)
    mutExcl = {p:{} for p in partners.keys()}
    simult = {p:{} for p in partners.keys()}
    for p in partners.keys():
        mutExcl[p] = {pr:set() for pr in partners[p]}
        simult[p] = {pr:set() for pr in partners[p]}
        interface = {}
        for pr in partners[p]:
            interface[pr], _ = get_interface (interactome, p, pr)
        for pr in partners[p]:
            for pr2 in partners[p] - {pr}:
                overlap, _ = sequence_overlap (interface[pr], interface[pr2])
                if overlap > cutoff:
                    mutExcl[p][pr].add(pr2)
                else:
                    simult[p][pr].add(pr2)
    return mutExcl, simult
