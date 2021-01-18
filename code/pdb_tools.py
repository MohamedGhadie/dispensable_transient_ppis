#----------------------------------------------------------------------------------------
# Modules for processing PDB structure data.
#----------------------------------------------------------------------------------------

import os
import sys
import re
import pickle
import warnings
import numpy as np
import pandas as pd
from Bio import Seq, SeqIO
from Bio.PDB import *
from Bio.SeqUtils import seq1
from urllib import error
from itertools import compress
from text_tools import write_fasta_file

#-----------------------------------------
# Global variables modified by modules
#-----------------------------------------

# allow PDB structure downloading
downloadStructures = True

# suppress PDB warnings
suppressWarnings = True

# maximum number of currently loaded structures
MAX_STRUC = 5

# IDs of currently loaded PDB structures
strucList = []

# dictionary of loaded structures
structures = {}

# dictionary of PDB chain sequences
chainSeq = {}

# dictionary of labels for chain sequence positions associated with 3D coordinates
chainStrucResLabel = {}

# mapping of chain sequence positions to residue IDs 
resPosToID = {}

#-----------------------------------------
# Modules
#-----------------------------------------

def allow_pdb_downloads (download):
    """Set global variable to allow PDB downloads.

    Args:
        download (bool): if true, allow PDB downloads.

    """
    global downloadStructures
    downloadStructures = download

def suppress_pdb_warnings (suppress):
    """Set global variable to suppress PDB warnings.

    Args:
        suppress (bool): if true, suppress PDB warnings.

    """
    global suppressWarnings
    suppressWarnings = suppress

def pdbfile_id (strucid):
    """Return structure file ID from structure ID by appending ! to lowercase chain 
        letters, if any.

    Args:
        strucid (str): structure ID.

    Returns:
        str: structure file ID.

    """
    m = re.search(r'-|_', strucid)
    if m:
        ind = m.start(0)
        return strucid[:ind] + re.sub(r'((?<!!)[a-z])', r'!\1', strucid[ind:])
    else:
        return strucid

def pdbfile_name (strucid):
    """Return stucture file name from structure ID by appending ! to lowercase chain 
        letters, if any.

    Args:
        strucid (str): structure ID.

    Returns:
        str: structure file name.

    """
    return 'pdb' + pdbfile_id (strucid) + '.ent'

def solve_pdbfile_id (id):
    """Resolve structure ID from structure file ID by removing ! from lowercase chain 
        letters, if any.

    Args:
        id (str): structure file ID.

    Returns:
        str: structure ID.

    """
    return id.replace('!', '')

def load_pdbtools_chain_sequences (inPath):
    """Load PDB chain sequences into global variable.

    Args:
        inPath (Path): path to file containing dictionary of chain sequences.

    """
    global chainSeq
    with open( inPath, 'rb' ) as f:
        chainSeq = pickle.load(f)

def load_pdbtools_chain_strucRes_labels (inPath):
    """Load labels for chain sequence positions with 3D coordinates into global variable.

    Args:
        inPath (Path): path to file containing dictionary of labels.

    """
    global chainStrucResLabel
    with open( inPath, 'rb' ) as f:
        chainStrucResLabel = pickle.load(f)

def clear_structures ():
    """Clear dictionary of PDB structures.
    
    """
    global structures, strucList
    structures.clear()
    strucList = []

def load_structure (pdbid, structureFile):
    """Load PDB structure from file.

    Args:
        pdbid (str): structure ID.
        structureFile (Path): path to file containing structure.

    Returns:
        Structure

    """
    if structureFile.is_file():
        try:
            return PDBParser( QUIET=suppressWarnings ).get_structure( pdbid, str(structureFile) )
        except:
            warnings.warn('Structure file not loaded')
            return None
    else:
        warnings.warn('Structure file not found')
        return None

def add_structure (pdbid, struc):
    """Add structure to structure dictionary.

    Args:
        pdbid (str): structure ID.
        struc (structure): structure to be added.

    """
    global structures, strucList
    if len(strucList) == MAX_STRUC:
        del structures[strucList.pop(0)]
    structures[pdbid] = struc
    strucList.append(pdbid)

def retrieve_structure (pdbid, pdbDir):
    """Retrieve PDB structure from local directory if available, otherwise donwload 
        from PDB database.

    Args:
        pdbid (str): structure ID.
        pdbDir (Path): path to local file directory where PDB structures are saved.

    """
    if pdbid not in structures:
        structureFile = pdbDir / pdbfile_name (pdbid)
        if downloadStructures and not structureFile.is_file():
            try:
                download_structure (pdbid, pdbDir)
            except error.URLError:
                pass
        if structureFile.is_file():
            struc = load_structure (pdbid, structureFile)
            if struc:
                add_structure (pdbid, struc)

def return_structure (pdbid, pdbDir):
    """Return PDB structure.

    Args:
        pdbid (str): structure ID.
        pdbDir (Path): path to local file directory where PDB structures are saved.

    """
    retrieve_structure (pdbid, pdbDir)
    if pdbid in structures:
        return structures[pdbid]
    else:
        return None

def return_chain (pdbID, chainID, pdbDir):
    """Return chain object for given PDB chain ID.

    Args:
        pdbID (str): structure ID.
        chainID (str): chain letter.
        pdbDir (Path): path to local file directory where PDB structures are saved.

    """
    struc = return_structure (pdbID, pdbDir)
    if struc:
        model = struc[0]
        if model.has_id(chainID):
            return model[chainID]
    return None

def download_structures (inPath, outDir):
    """Download PDB structures associated with a list of pdb IDs.

    Args:
        inPath (Path): path to file containing list of PDB IDs.
        outDir (Path): path to directory where PDB files are saved.

    """
    with open(inPath, 'r') as fin:
        pdbIDs = list( fin.read().split() )
    n = len(pdbIDs)
    failed = 0
    for i, id in enumerate(pdbIDs):
        filename = outDir / ('pdb' + id + '.ent')
        if not filename.is_file():
            try:
                print()
                download_structure (id, outDir)
            except error.URLError:
                failed += 1
        sys.stdout.write('  %d out of %d IDs processed (%f%%), %d downloads failed \r' 
                            % (i+1, n, 100*(i+1)/n, failed))
        sys.stdout.flush()
    print()

def download_structure (pdbid, outDir, fileFormat = 'pdb'):
    """Download PDB structure file.

    Args:
        pdbid (str): structure ID.
        outDir (Path): path to directory where PDB file is saved.
        fileFormat (str): format of PDB file.

    """
    try:
        PDBList().retrieve_pdb_file (pdbid,
                                     pdir = outDir,
                                     file_format = fileFormat)
        filename = outDir / ('pdb' + pdbid + '.pdb')
        if (fileFormat is 'pdb') and filename.is_file():
            base = os.path.splitext(filename)[0]
            os.rename(filename, base + '.ent')
    except error.URLError:
        raise

def return_chain_sequence (chainID):
    """Return PDB chain sequence.

    Args:
        chainID (str): chain ID.

    """
    if len(chainSeq):
        if chainID in chainSeq:
            return chainSeq[chainID]
    else:
        warnings.warn("Chain sequence dictionary is empty")
    return None

def return_chain_strucRes_label (chainID):
    """Return labels for chain sequence positons having 3D coordinates.

    Args:
        chainID (str): chain ID.

    """
    if len(chainStrucResLabel):
        if chainID in chainStrucResLabel:
            return chainStrucResLabel[chainID]
    else:
        warnings.warn("Chain structured residue label dictionary is empty")
    return None

def residueID (residue):
    """Return residue ID in string format.

    Args:
        residue (Bio.PDB.Residue): residue whose ID is to be returned.
    
    Returns:
        str: residue ID.
    
    """
    het, pos, inCode = residue.get_id()
    return (het + '_' + str(pos) + '_' + inCode).strip()

def all_atoms (chain, residues):
    """Return all atoms combined for all residues in a chain.

    Args:
        chain (Bio.PDB.Chain): chain whose atoms are to be returned.
        residues (dict): residues per chain.
    
    Returns:
        list: list of atoms.
    
    """
    chainID = chain.get_id()
    allatoms = []
    atoms = atoms_per_residue (residues[chainID])
    for residue in residues[chainID]:
        resID = residueID (residue)
        allatoms.extend(atoms[resID])
    return allatoms

def atoms_per_residue (residues):
    """Return all atoms per residue for a list of PDB residues.

    Args:
        residues (list): list of PDB residues.
    
    Returns:
        dict: list of atoms by residue.
    
    """
    allatoms = {}
    for residue in residues:
        resID = residueID (residue)
        atoms = residue.get_unpacked_list()
        if len(atoms) > 0:
            allatoms[resID] = atoms
    return allatoms
    
def residues_per_chain (model):
    """Return all residues per chain for all chains in a PDB model.

    Args:
        model (Bio.PDB.Model): PDB model.
    
    Returns:
        dict: list of all residues for each chain, including hetero-residues.
    
    """
    allresidues = {}
    for chain in model.get_chains():
        residues = [res for res in chain.get_residues()]
        if len(residues) > 0:
            allresidues[chain.get_id()] = residues
    return allresidues

def structured_residues_per_chain (pdbid, pdbDir):
    """Return all SEQRES residues per chain for all chains in a PDB model.

    Args:
        pdbid (str): PDB ID for model parent structure.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        dict: list of all SEQRES residues for each chain, including hetero-residues.
    
    """
    resPerChain = {}
    struc = return_structure (pdbid, pdbDir)
    if struc:
        model = struc[0]
        for chain in model.get_chains():
            residues = structured_chain_residues (pdbid, chain.get_id(), pdbDir)
            if residues:
                resPerChain[chain.get_id()] = residues
    return resPerChain

def nonhet_residues_per_chain (model):
    """Return all non-hetero residues per chain for all chains in a PDB model. 

    Args:
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all non-hetero residues for each chain.
    
    """
    allresidues = {}
    for chain in model.get_chains():
        residues = []
        for res in chain.get_residues():
            het, resnum, incode = res.get_id()
            if het == ' ':
                residues.append(res)
        if len(residues) > 0:
            allresidues[chain.get_id()] = residues
    return allresidues

def structured_chain_residues (pdbid, chainID, pdbDir):
    """Return all SEQRES residues of a chain. 

    Args:
        pdbid (str): PDB ID for model parent structure.
        chainID (str): chain ID.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        dict: list of all SEQRES residues, including hetero-residues.
    
    """
    struc = return_structure (pdbid, pdbDir)
    if struc:
        model = struc[0]
        if model.has_id(chainID):
            residues = chain_residues (model, chainID)
            pos = return_chain_res_IDsToPos (pdbid, chainID, pdbDir)
            return [res for res in residues if res.get_id() in pos]
    return []

def chain_residues (model, chainID):
    """Return all residues of a chain that have 3D coordinates. 

    Args:
        model (Bio.PDB.Model): PDB model.
        chainID (str): chain ID.
    
    Returns:
        list: all chain residues, including hetero-residues.
    
    """
    chain = model[chainID]
    return [res for res in chain.get_residues()]

def get_chain_IDs (pdbid, pdbDir):
    """Return all chain ID letters for a given structure. 

    Args:
        pdbid (str): structure ID.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        list
    
    """
    struc = return_structure (pdbid, pdbDir)
    if struc:
        model = struc[0]
        return [chain.get_id() for chain in model.get_chains()]
    return []

def structured_residue_sequence (pdbid, chainID, pdbDir):
    """Return sequence of chain residues that have 3D cooridnates. 

    Args:
        pdbid (str): structure ID.
        chainID (str): chain letter.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        list
    
    """
    residues = structured_chain_residues (pdbid, chainID, pdbDir)
    if residues:
        return ''.join([seq1(res.get_resname(), undef_code='X') for res in residues])
    else:
        return ''

def structured_residue_IDs (pdbid, chainID, pdbDir):
    """Return residues IDs that have 3D cooridnates. 

    Args:
        pdbid (str): structure ID.
        chainID (str): chain letter.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        list
    
    """
    residues = structured_chain_residues (pdbid, chainID, pdbDir)
    return [res.get_id() for res in residues]

def structured_residue_ID (resOrder, pdbid, chainID, pdbDir):
    """Return residue ID in specific sequence position among residues with 3D cooridnates. 

    Args:
        resOrder (int, 'first', 'last'): position of residue.
        pdbid (str): structure ID.
        chainID (str): chain letter.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        list
    
    """
    resIDs = structured_residue_IDs (pdbid, chainID, pdbDir)
    toNum = {'first':1, 'last':len(resIDs)}
    if resOrder in toNum:
        resOrder = toNum[resOrder]
    try:
        return resIDs[resOrder - 1]
    except:
        warnings.warn("Residue order requested out of range")
        return ()
    
def is_contact (res, other_residues, maxDist):
    """Check if residue is in contact with a group of residues. 

    Args:
        res (Bio.PDB.Residue): residue to be checked.
        other_residues (list): other residues to check against.
        maxDist (numeric): maximum distance in Angstroms allowed for contact.
    
    Returns:
        boolean
    
    """
    for r in other_residues:
        if get_distance(res, r) < maxDist:
            return True
    return False

def get_distance (res1, res2):
    """Calculate Euclidean distance between two residues. 

    Args:
        res1 (Bio.PDB.Residue): first residue.
        res2 (Bio.PDB.Residue): second residue.
    
    Returns:
        float
    
    """
    dist = [a1 - a2 for a1 in res1.get_unpacked_list() for a2 in res2.get_unpacked_list()]
    return min(dist)

def count_res_neighbors (res, otherRes, nbDist):
    """Count neighbors of a residue within a group of residues.

    Args:
        res (Residue): residue to count neighbors for.
        otherRes (list): other residues in the group.
        nbDist (numeric): cutoff distance in Angstroms for neighbor residues.
    
    Returns:
        int
    
    """
    nb = [get_distance(res, res2) < nbDist for res2 in otherRes]
    return sum(nb)

def get_interface_by_chainIDs (pdbDir,
                               chain_id1,
                               chain_id2,
                               maxDist = 5):
    """Calculate binding interface between two chains.

    Args:
        pdbDir (Path): directory containing PDB structure files.
        chain_id1 (str): ID of first chain.
        chain_id2 (str): ID of second chain.
        maxDist (numeric): maximum distance in Angstroms allowed between binding residues.
    
    Returns:
        list, list: positions of interface residues on chain sequences.

    """
    pdbID1, id1 = chain_id1.split('_')
    pdbID2, id2 = chain_id2.split('_')
    if pdbID1 == pdbID2:
        pdbID = pdbID1
        struc = return_structure (pdbID, pdbDir)
        if struc:
            residues1 = structured_chain_residues (pdbID, id1, pdbDir)
            residues2 = structured_chain_residues (pdbID, id2, pdbDir)
            interfaceIndex1, interfaceIndex2 = get_interface_indices (residues1,
                                                                      residues2,
                                                                      maxDist = maxDist)
            interfaceResID1 = [residues1[i].get_id() for i in interfaceIndex1]
            interfaceResID2 = [residues2[i].get_id() for i in interfaceIndex2]
            chain1ResIDsToPos = return_chain_res_IDsToPos (pdbID, id1, pdbDir)
            chain2ResIDsToPos = return_chain_res_IDsToPos (pdbID, id2, pdbDir)
            interface1 = [chain1ResIDsToPos[id] for id in interfaceResID1 if id in chain1ResIDsToPos]
            interface2 = [chain2ResIDsToPos[id] for id in interfaceResID2 if id in chain2ResIDsToPos]
            return sorted(interface1), sorted(interface2)
        else:
            warnings.warn('Structure file for chains %s and %s not found' % (chain_id1, chain_id2))
            return [], []
    else:
        return [], []

def get_interface_indices (residues1, residues2, maxDist = 5):
    """Return the indices of residues at the interface between two groups of residues.

    Args:
        residues1 (list): first list of residues.
        residues2 (list): second list of residues.
        maxDist (numeric): maximum distance in Angstroms allowed between binding residues. 
    
    Returns:
        list, list: indices of interface residues in each group.

    """
    numRes1 = len(residues1)
    numRes2 = len(residues2)
    resPairs = [(r1, r2) for r1 in range(numRes1) for r2 in range(numRes2)]
    distances = [get_distance(res1, res2) for res1 in residues1 for res2 in residues2]
    interfaceIndices = [pair for pair, d in zip(resPairs, distances) if d < maxDist]
    if len(interfaceIndices) > 0:
        interfaceIndex1, interfaceIndex2 = zip(*interfaceIndices)
        return sorted(set(interfaceIndex1)), sorted(set(interfaceIndex2))
    else:
        return [], []

def count_neighbors_by_chainIDs (pdbid,
                                 resPos,
                                 hostchainID,
                                 chainIDs,
                                 pdbDir,
                                 nbDist = 5):
    """Count neighbors of a residue within a group of PDB chains.

    Args:
        pdbid (str): structure ID.
        resPos (int): residue position on chain sequence.
        hostchainID (str): ID of chain carrying the residue.
        chainIDs (list): IDs of other chains to consider in the structure.
        pdbDir (Path): directory path containing PDB structures.
        nbDist (numeric): cutoff distance in Angstroms for neighbor residues.
    
    Returns:
        numeric

    """
    _ , hostchainID = hostchainID.split('_')
    resID = return_chain_res_posToID (pdbid, hostchainID, resPos, pdbDir)
    if resID:
        struc = return_structure (pdbid, pdbDir)
        if struc:
            residues = structured_residues_per_chain (pdbid, pdbDir)
            res = struc[0][hostchainID][resID]
            otherRes = [r for r in residues[hostchainID] if r.get_id() != resID]
            for id in {id.split('_')[1] for id in chainIDs}:
                if id != hostchainID:
                    otherRes.extend(residues[id])
            return count_res_neighbors (res, otherRes, nbDist)
    return np.nan

def return_multichain_res_posToIDs (pdbid, chainResPos, pdbDir):
    """Return residue IDs in multiple chains using sequence positions. 

    Args:
        pdbid (str): structure ID.
        chainResPos (list): tuples of chain IDs with residue positions.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        list: residue IDs in the form (chainID, pos, resID).
    
    """
    IDs = []
    for chainID, pos in chainResPos:
        resID = return_chain_res_posToID (pdbid, chainID, pos, pdbDir)
        IDs.append( (chainID, pos, resID) )
    return IDs

def return_chain_res_posToID (pdbid, chainID, resPos, pdbDir):
    """Return residue ID using sequence position. 

    Args:
        pdbid (str): structure ID.
        chainID (str): chain ID.
        resPos (int): residue sequence position.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        tuple: residue ID.
    
    """
    IDmaps = return_chain_res_posToIDs (pdbid, chainID, pdbDir)
    if resPos in IDmaps:
        return IDmaps[resPos]
    else:
        return None

def return_chain_res_IDToPos (pdbid, chainID, resID, pdbDir):
    """Return sequence position using residue ID. 

    Args:
        pdbid (str): structure ID.
        chainID (str): chain ID.
        resID (tuple): residue ID.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        int: residue position.
    
    """
    posmaps = return_chain_res_IDsToPos (pdbid, chainID, pdbDir)
    if resID in posmaps:
        return posmaps[resID]
    else:
        return None
    
def return_chain_res_IDsToPos (pdbid, chainID, pdbDir):
    """Return mappings of all chain residue IDs to sequence positions. 

    Args:
        pdbid (str): structure ID.
        chainID (str): chain ID.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        dict
    
    """
    idmaps = return_chain_res_posToIDs (pdbid, chainID, pdbDir)
    posmaps = {}
    for pos, id in idmaps.items():
        posmaps[id] = pos
    return posmaps

def return_chain_res_posToIDs (pdbid, chainID, pdbDir):
    """Return mappings of all chain sequence positions to residue IDs.

    Args:
        pdbid (str): structure ID.
        chainID (str): chain ID.
        pdbDir (Path): path to file directory containing PDB structures.
    
    Returns:
        dict
    
    """
    get_chain_res_posToIDs (pdbid, chainID, pdbDir)
    fullID = pdbid + '_' + chainID
    if fullID in resPosToID:
        return resPosToID[fullID]
    else:
        return {}

def get_chain_res_posToIDs (pdbid, chainID, pdbDir):
    """Map chain sequence positions to residue IDs, and save to global dictionary.

    Args:
        pdbid (str): structure ID.
        chainID (str): chain ID.
        pdbDir (Path): path to file directory containing PDB structures.

    """
    global resPosToID
    fullID = pdbid + '_' + chainID
    if fullID not in resPosToID:
        struc = return_structure (pdbid, pdbDir)
        if struc:
            maps = allchain_res_pos_to_IDs (pdbid, struc)
            for id, posmaps in maps.items():
                resPosToID[id] = posmaps

def allchain_res_pos_to_IDs (pdbid, struc):
    """Return mappings of sequence positions to residue IDs for all chains in a structure.
        First model in structure is used.

    Args:
        pdbid (str): structure ID.
        struc (Bio.PDB.Structure): PDB structure.

    Returns:
        dict

    """
    allmaps = {}
    model = struc[0]
    for chain in model:
        idmaps = chain_res_pos_to_IDs (pdbid, chain)
        if idmaps:
            fullID = pdbid + '_' + chain.get_id()
            allmaps[fullID] = idmaps
    return allmaps

def chain_res_pos_to_IDs (pdbid, chain, minFracValidRes = 0.8):
    """Return mapping of all chain sequence positions to residue IDs.

    Args:
        pdbid (str): structure ID.
        chain (Bio.PDB.Chain): PDB chain.
        minFracValidRes (numeric): minimum required fraction of chain residue names 
                                    successfully read by Biopython.
    
    Returns:
        dict
    
    """
    IDs = {}
    fullID = pdbid + '_' + chain.get_id()
    sequence = return_chain_sequence (fullID)
    strucResLabel = return_chain_strucRes_label (fullID)
    mStart, resIDs = validated_strucRes_IDs (chain,
                                             sequence,
                                             strucResLabel,
                                             minFracValidRes = minFracValidRes)
    if resIDs:
        for seqPos, label in enumerate(strucResLabel):
            if (seqPos >= mStart) and (label == '-'):
                IDs[seqPos + 1] = resIDs.pop(0)
    return IDs

def valid_strucRes (pdbid, chainID, pdbDir, minFracValidRes = 0.8):
    """Check whether chain sequence read by Biopython matches a continuous subset of 
        SEQRES chain sequence.

    Args:
        pdbid (str): structure ID.
        chainID (str): chain letter.
        pdbDir (Path): path to file directory containing PDB structures.
        minFracValidRes (numeric): minimum required fraction of chain residue names 
                                    successfully read by Biopython.
    
    Returns:
        bool
    
    """
    fullID = pdbid + '_' + chainID
    sequence = return_chain_sequence (fullID)
    strucResLabel = return_chain_strucRes_label (fullID)
    struc = return_structure (pdbid, pdbDir)
    if struc:
        chain = struc[0][chainID]
        mStart, resIDs = validated_strucRes_IDs (chain,
                                                 sequence,
                                                 strucResLabel,
                                                 minFracValidRes = minFracValidRes)
        if resIDs:
            return True
    return False

def validated_strucRes_IDs (chain, seq, strucResLabel, minFracValidRes = 0.8):
    """Return residue IDs of chain sequence read by Biopython matching a continuous subset of 
        SEQRES chain sequence.

    Args:
        chain (Bio.PDB.Chain): chain object.
        seq (str): full chain sequence in SEQRES records.
        strucResLabel (str): labels for chain sequence residues with 3D coordinates.
        minFracValidRes (numeric): minimum required fraction of chain residue names 
                                    successfully read by Biopython.
    
    Returns:
        numeric: SEQRES start position where Biopython sequence matches, list: matching residue IDs
    
    """
    if seq and strucResLabel:
        if len(seq) == len(strucResLabel):
            # construct chain sequence of structured residues based on SEQRES record 
            seqres = ''.join(compress(seq, [c == '-' for c in strucResLabel]))
            seqresLen = len(seqres)
        
            # construct chain sequence of structured residues as read by BioPython
            coordresIDs = [ res.get_id() for res in chain.get_residues() ]
            coordres = [ seq1(res.get_resname(), undef_code='.') for res in chain.get_residues() ]
            coordresIDs = coordresIDs[ : seqresLen ]
            coordres = ''.join( coordres[ : seqresLen ] )
            # check if fraction of valid structured residue names read by Biopython 
            # meets cutoff
            if coordres:
                if ( 1 - coordres.count('.') / len(coordres) ) >= minFracValidRes:
                    # check if the sequence of structured residues read by Biopython 
                    # is a subset of SEQRES
                    match = re.search(coordres, seqres)
                    if match:
                        return match.start(), coordresIDs
    return np.nan, []

def produce_chain_list (inPath, outPath):
    """Produce list of PDB chain IDs from a fasta file of chain sequences.

    Args:
        inPath (Path): path to fasta file containing chain sequences.
        outPath (Path): file path to save chain IDs to.

    """
    chainList = []
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    chainList = map(str.strip, [row.id for row in s])
    with open(outPath, 'w') as fout:
        for id in sorted(set(chainList)):
            fout.write(id + '\n')

def produce_pdb_ids (inPath, outPath):
    """Produce list of PDB IDs from a list of chain IDs.

    Args:
        inPath (Path): path to file containing list of chain IDs.
        outPath (Path): file path to save PDB IDs to.

    """
    with open(inPath, 'r') as fin:
        chainIDs = list(fin.read().split())
    pdbIDs = [id[ : id.find('_')] for id in chainIDs if '_' in id]
    with open(outPath, 'w') as fout:
        for id in sorted(set(pdbIDs)):
            fout.write(id + '\n')

def produce_chain_dict (inPath, outPath):
    """Produce a dictionary of chain IDs associated with each PDB ID. 

    Args:
        inPath (Path): path to file containing list of chain IDs.
        outPath (Path): path to save pickle file to.

    """
    with open(inPath, 'r') as fin:
        chainIDs = list(fin.read().split())
    chains = {}
    for chainid in chainIDs:
        pdbid, id = chainid.split('_')
        if pdbid in chains:
            chains[pdbid].add(id)
        else:
            chains[pdbid] = {id}
    with open(outPath, 'wb') as fOut:
        pickle.dump(chains, fOut)

def produce_chain_struc_sequences (chainIDs, pdbDir, outPath):
    """Produce sequences of chain structure residues that are also part chain SEQRES 
        sequences in fasta format.

    Args:
        chainIDs (list): list of tuples with template ID and chain ID.
        pdbDir (Path): path to file directory containing PDB structures.
        outPath (Path): file path to save sequences to.

    """
    IDs, seqs, n = [], [], len(chainIDs)
    for i, (templateID, chainID) in enumerate(chainIDs):
        sys.stdout.write('  chain %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        pdbid = templateID.split('-')[0]
        seq = structured_residue_sequence (pdbid, chainID, pdbDir)
        if seq:
            IDs.append('_'.join([templateID, chainID]))
            seqs.append(seq)
    df = pd.DataFrame(data={"ChainID":IDs, "Sequence":seqs})
    write_fasta_file (df, "ChainID", "Sequence", outPath)

def write_partial_structure (pdbid, chainIDs, pdbDir, outPath, resIDs = None):
    """Produce a partial PDB structure file for a subset of chains.

    Args:
        pdbid (str): ID of full structure.
        chainIDs (list): IDs for chains to be included in partial structure file.
        pdbDir (Path): path to file directory containing PDB structures.
        outPath (Path): file path to save partial structure to.
        resIDs (dict): residue IDs to include from each chain.

    """
    struc = return_structure (pdbid, pdbDir)
    if struc:
        pdbio = PDBIO()
        pdbio.set_structure(struc)
        pdbio.save(str(outPath), select = ChainSelect(chainIDs, resIDs = resIDs))
    else:
        warnings.warn("Structure file for PDB ID %s not found in order to write parital structure" % pdbid)

class ChainSelect (Select):
    """Class for chain collection used for creating partial PDB structures.

    Args:
        chain_letters (list): chain letter IDs.
        resIDs (dict): residue IDs to include from each chain.

    """
    def __init__ (self, chain_letters, resIDs = None):
        self.chain_letters = chain_letters
        self.resIDs = resIDs
    
    def accept_chain (self, chain):
        return (chain.get_id() in self.chain_letters)
    
    def accept_residue (self, residue):
        chain_letter = residue.parent.get_id()
        if not self.resIDs:
            return True
        elif chain_letter in self.resIDs:
            return residue.get_id() in self.resIDs[chain_letter]
        else:
            return False
