#----------------------------------------------------------------------------------------
# Modules for simple processing of residue relative solvent accessibility (RSA) data.
#----------------------------------------------------------------------------------------

import os
import io
import re
import pickle
import numpy as np
import pandas as pd
from Bio.PDB import DSSP
from pdb_tools import pdbfile_id, return_structure, write_partial_structure

#-----------------------------------------
# Global variables modified by modules
#-----------------------------------------

# Maximum residue solvent accessibility
maxResAcc = {'Wilke': {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, 'C': 167.0,
                       'E': 223.0, 'Q': 225.0, 'G': 104.0, 'H': 224.0, 'I': 197.0,
                       'L': 201.0, 'K': 236.0, 'M': 224.0, 'F': 240.0, 'P': 159.0,
                       'S': 155.0, 'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0},
             'Empirical': {}}

#-----------------------------------------
# Modules
#-----------------------------------------

def load_empirical_maxAcc (inPath):
    """Load empirically calculated residue maximum solvent accessibility (MSA).

    Args:
        inPath (Path): path to file containing dictionary of MSA.

    """
    global maxResAcc
    with open( inPath, 'rb' ) as f:
        maxResAcc['Empirical'] = pickle.load(f)

def calculate_structure_rsa (accDir, pdbDir, pdbid, selChains = None):
    """Return a dictionary of RSA values for all residues in a structure.

    Args:
        accDir (Path): file directory where residue solvent accessibility values are saved.
        pdbDir (Path): path to local file directory where PDB structures are saved.
        pdbid (str): structure ID.
        selChains (list): letter IDs for chains to be included in calculations.
                            If None, all chais are included.

    Returns:
        dict

    """
    if selChains:
        pdbFileID = pdbfile_id (pdbid + '_' + '_'.join(sorted(selChains)))
    else:
        pdbFileID = pdbfile_id (pdbid)
    accFile = accDir / (pdbFileID + '.txt')
    if not accFile.is_file():
        make_acc_file (accDir, pdbDir, pdbid, selChains = selChains)
    if accFile.is_file():
        acc = read_acc_file (accFile)
        return acc_to_rsa (acc)
    else:
        return {}

def make_acc_file (accDir, pdbDir, pdbid, tempDir = None, selChains = None):
    """Write solvent accessibility data for a structure to file.

    Args:
        accDir (Path): file directory where residue solvent accessibility values are saved.
        pdbDir (Path): path to local file directory where PDB structures are saved.
        pdbid (str): structure ID.
        tempDir (Path): temporary file directory where partial structure files are saved.
                        If None, pdbDir is used. 
        selChains (list): letter IDs for chains to be included in calculations.
                            If None, all chais are included.

    """
    if selChains:
        if not tempDir:
            tempDir = pdbDir
        pdbfileID = pdbfile_id (pdbid + '_' + '_'.join(sorted(selChains)))
        strucFile = tempDir / ('pdb' + pdbfileID + '.ent')
        accFile = accDir / (pdbfileID + '.txt')
        write_partial_structure (pdbid, selChains, pdbDir, strucFile)
        struc = return_structure (pdbfileID, tempDir)
    else:
        pdbfileID = pdbfile_id (pdbid)
        strucFile = pdbDir / ('pdb' + pdbfileID + '.ent')
        accFile = accDir / (pdbfileID + '.txt')
        struc = return_structure (pdbfileID, pdbDir)
    if struc:
        try:
            rsa = DSSP (struc[0], strucFile, acc_array = 'Wilke')
            acc = rsa_to_acc (rsa, maxAccTable = 'Wilke')
            write_acc_file (acc, accFile)
        except:
            pass
    if selChains and strucFile.is_file():
        os.remove(strucFile)

def acc_to_rsa (acc, maxAccTable = 'Empirical'):
    """Return dictionary of residue RSA calculated from residue absolute solvent accessibility (ASA).

    Args:
        acc (dict): dictionary of residue ASA values.
        maxAccTable (str): name of maximum ASA table to use, either 'Wilke' or 'Empirical'.

    Returns:
        dict

    """
    rsa = {}
    maxAcc = maxResAcc[maxAccTable]
    for k in acc.keys():
        val = list(acc[k])
        if val[1] in maxAcc:
            val[3] /= maxAcc[val[1]]
        else:
            val[3] = np.nan
        rsa[k] = val
    return rsa

def rsa_to_acc (rsa, maxAccTable = 'Empirical'):
    """Return dictionary of residue absolute solvent accessibility (ASA) calculated from RSA.

    Args:
        rsa (dict): dictionary of residue RSA values.
        maxAccTable (str): name of maximum ASA table to use, either 'Wilke' or 'Empirical'.

    Returns:
        dict

    """
    acc = {}
    maxAcc = maxResAcc[maxAccTable]
    for k in rsa.keys():
        val = list(rsa[k])
        if val[1] in maxAcc:
            val[3] *= maxAcc[val[1]]
        else:
            val[3] = np.nan
        acc[k] = val
    return acc

def write_acc_file (dsspDict, outPath):
    """Write residue absolute solvent accessibility (ASA) data to file in DSSP format.

    Args:
        dsspDict (dict): dictionary of residue ASA values.
        outPath (Path): file path to save ASA data.

    """
    with io.open(outPath, "w") as fout:
        fout.write ('\t'.join(['dssp_index', 'chain', 'hetero', 'residue_number', 'insCode',
                               'residue', 'sec_structure', 'Acc', 'Phi', 'Psi',
                               'NH_to_O_1_relidx', 'NH_to_O_1_energy', 'O_to_NH_1_relidx',
                               'O_to_NH_1_energy', 'NH_to_O_2_relidx', 'NH_to_O_2_energy',
                               'O_to_NH_2_relidx', 'O_to_NH_2_energy']) + '\n')
        for k, dssp in dsspDict.items():
            fout.write ('\t'.join(map(str, [dssp[0], k[0]] + 
                                           list(k[1]) + 
                                           dssp[1:])) + '\n')

def read_acc_file (inPath):
    """Read residue absolute solvent accessibility (ASA) data from file in DSSP format.

    Args:
        inPath (Path): path to file containing ASA data.

    Returns:
        dict

    """
    accDict = {}
    acc = pd.read_table (inPath, dtype={'chain':str}, sep='\t')
    for _, row in acc.iterrows():
        k = (row.chain, (row.hetero, row.residue_number, row.insCode))
        accDict[k] = (row.dssp_index, row.residue, row.sec_structure, row.Acc, row.Phi,
                      row.Psi, row.NH_to_O_1_relidx, row.NH_to_O_1_energy,
                      row.O_to_NH_1_relidx, row.O_to_NH_1_energy, row.NH_to_O_2_relidx,
                      row.NH_to_O_2_energy, row.O_to_NH_2_relidx, row.O_to_NH_2_energy )
    return accDict

def read_dssp_file (inPath):
    """Read residue absolute solvent accessibility (ASA) data from file provided by DSSP.

    Args:
        inPath (Path): path to file containing ASA data.

    Returns:
        dict

    """
    resAcc = {}
    with io.open(inPath, 'r', encoding='utf-8') as f:
        for headerLines, line in enumerate(f):
            if line.strip().startswith('#'):
                break
    with io.open(inPath, "r", encoding="utf-8") as f:
        for _ in range(headerLines + 1):
            next(f)
        for line in f:
            if '!' not in line:
                ind, res, chain, rem = line[:5], line[5:11], line[11], line[13:]
                if chain == '>':
                    chain = line[159:]
                res, insCode = res[:-1].strip(), res[-1]
                if res.isdigit():
                    key = (chain.strip(), (' ', int( res ), insCode))
                else:
                    m = re.match(r"([\D]+)([0-9]+)", res)
                    if m:
                        het, num = m.groups()
                        key = (chain.strip(), (het.strip(), int( num ), insCode))
                    else:
                        print('invalid residue number')
                        continue
                
                aa, structure, bp1, bp2, acc, rem = (rem[ 0 ],
                                                     rem[ 3 : 12 ],
                                                     rem[ 12 : 16 ],
                                                     rem[ 16: 21 ],
                                                     rem[ 21: 25 ],
                                                     rem[ 25 : ])                                
                resAcc[key] = [int(ind), aa, structure, float(acc), bp1, bp2, rem]
    return resAcc
