#----------------------------------------------------------------------------------------
# Modules for calculating and processing structural energy.
#----------------------------------------------------------------------------------------

import os
import io
import re
import sys
import pandas as pd
from pathlib import Path
from text_tools import write_beluga_job
from pdb_tools import pdbfile_id, solve_pdbfile_id, clear_structures, write_partial_structure

def produce_initial_ppi_energy_file (inPath, outPath):
    """Write initial PPI structure interaction energy file with no energy values.

    Args:
        inPath (Path): path to file containing structural interactome.
        outPath (Path): path to output energy file.

    """
    ppis = pd.read_table (inPath, sep='\t')
    ppis["Interaction_energy"] = '-'
    chain_1, chain_2 = zip(* ppis["Chain_pairs"].apply(lambda x: x.split('+')).values)
    ppis["Chain_1"] = [c.split('_')[1] for c in chain_1]
    ppis["Chain_2"] = [c.split('_')[1] for c in chain_2] 
    ppis[["Protein_1", "Protein_2", "Complex_ID", "Chain_1", "Chain_2", "Interaction_energy"]
        ].to_csv (outPath, index=False, sep='\t')

def produce_initial_ppi_template_energy_file (inPath, outPath):
    """Write initial PPI template interaction energy file with no energy values.

    Args:
        inPath (Path): path to file containing structural interactome.
        outPath (Path): path to output energy file.

    """
    ppis = pd.read_table (inPath, sep='\t')
    ppis["Interaction_energy"], ppis["Chain_1"], ppis["Chain_2"] = '-', '-', '-'
    for i, row in ppis.iterrows():
        complexID, templateID = solve_pdbfile_id(row.Alignment_file_ID).split('_')
        p1, p2 = complexID.split('=')
        pdbid, c1, c2 = templateID.split('-')
        mapping = {p1:c1, p2:c2}
        ppis.loc[i, "Complex_ID"] = pdbid
        ppis.loc[i, "Chain_1"] = mapping[row.Protein_1]
        ppis.loc[i, "Chain_2"] = mapping[row.Protein_2] 
    ppis[["Protein_1", "Protein_2", "Complex_ID", "Chain_1", "Chain_2", "Interaction_energy"]
        ].to_csv (outPath, index=False, sep='\t')

def read_unprocessed_ddg_mutations (inPath, type = 'binding'):
    """Read PDB chain mutations with missing ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    mutations = {}
    done = set()
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) >= 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, mut, ch_partner = strsplit[:8]
                if type is 'binding':
                    mutation = '-'.join( [protein, partner, pr_pos, mut[-1]] )
                elif type is 'folding':
                    mutation = '-'.join( [protein, pr_pos, mut[-1]] )
                if mutation not in done:
                    if len(strsplit) == 8:
                        done.add(mutation)
                        if type is 'binding':
                            struc = (pdbid,) + tuple(sorted([chainID, ch_partner]))
                        elif type is 'folding':
                            struc = pdbid, chainID
                        if struc in mutations:
                            mutations[struc].add(mut)
                        else:
                            mutations[struc] = {mut}
                    elif strsplit[8] is not 'X':
                        done.add(mutation)
    return {k:list(v) for k, v in mutations.items()}

def read_unprocessed_energy_ppis (inPath):
    """Read PPI structures with missing interaction energy values from file.

    Args:
        inPath (Path): path to file containing PPI structure energy values.

    Returns:
        list

    """
    ppis = pd.read_table (inPath, sep='\t')
    ppis = ppis[ppis["Interaction_energy"] == '-']
    structures = {(strucID,) + tuple(sorted([c1, c2]))
                    for strucID, c1, c2 in ppis[["Complex_ID", "Chain_1", "Chain_2"]].values}
    return list(structures)

def produce_foldx_and_beluga_jobs (data,
                                   pdbDir,
                                   outDir,
                                   type,
                                   foldxParam = None,
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem = None,
                                   mem_per_cpu = '4G',
                                   outputfile = '%x-%j.out',
                                   errorfile = None,
                                   username = '',
                                   extraCommands = None,
                                   serverDataDir = '../data'):
    """Produce jobs for the FoldX method for ∆∆G calculations using several commands.

    Args:
        data (dict or list): dict of mutations associated with each structure, or list of
                             structure-chain triplets.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX data files and Beluga jobs to.
        type (str): type of FoldX calculation; 'binding' for interface ∆∆G, 
                    'folding' for protein folding ∆∆G, 'ppi_energy' for interaction energy.
        foldxParam (dict): parameters used for each FoldX job, otherwise default.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    if type is 'folding':
        produce_foldx_buildmodel_jobs (data,
                                       pdbDir,
                                       outDir / 'data',
                                       parameters = foldxParam)
        produce_beluga_foldx_jobs (data,
                                   type,
                                   outDir / 'jobs',
                                   account = account,
                                   walltime = walltime,
                                   ntasks = ntasks,
                                   nodes = nodes,
                                   ntasks_per_node = ntasks_per_node,
                                   cpus_per_task = cpus_per_task,
                                   mem = mem,
                                   mem_per_cpu = mem_per_cpu,
                                   outputfile = outputfile,
                                   errorfile = errorfile,
                                   username = username,
                                   extraCommands = extraCommands,
                                   serverDataDir = serverDataDir)
    elif type is 'binding':
        produce_foldx_pssm_jobs (data,
                                 pdbDir,
                                 outDir / 'data',
                                 parameters = foldxParam)
        produce_beluga_foldx_jobs (data,
                                   type,
                                   outDir / 'jobs',
                                   account = account,
                                   walltime = walltime,
                                   ntasks = ntasks,
                                   nodes = nodes,
                                   ntasks_per_node = ntasks_per_node,
                                   cpus_per_task = cpus_per_task,
                                   mem = mem,
                                   mem_per_cpu = mem_per_cpu,
                                   outputfile = outputfile,
                                   errorfile = errorfile,
                                   username = username,
                                   extraCommands = extraCommands,
                                   serverDataDir = serverDataDir)
    elif type is 'ppi_energy':
        produce_foldx_analyseComplex_jobs (data,
                                           pdbDir,
                                           outDir / 'data',
                                           parameters = foldxParam)
        produce_beluga_foldx_jobs (data,
                                   type,
                                   outDir / 'jobs',
                                   account = account,
                                   walltime = walltime,
                                   ntasks = ntasks,
                                   nodes = nodes,
                                   ntasks_per_node = ntasks_per_node,
                                   cpus_per_task = cpus_per_task,
                                   mem = mem,
                                   mem_per_cpu = mem_per_cpu,
                                   outputfile = outputfile,
                                   errorfile = errorfile,
                                   username = username,
                                   extraCommands = extraCommands,
                                   serverDataDir = serverDataDir)

def produce_foldx_buildmodel_jobs (mutations, pdbDir, outDir, parameters = None):
    """Produce jobs for the FoldX method for ∆∆G calculation using 'BuildModel' command.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX jobs to.
        parameters (dict): parameters used for each FoldX job, otherwise default.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    default_param = {'temp':298,
                     'ph':7,
                     'ionStrength':0.05,
                     'water':'-IGNORE',
                     'vdwDesign':2}
    if not parameters:
        parameters = {}
    mutations = mutations.items()
    n = len(mutations)
    print('Writing FoldX BuildModel files:')
    for i, (struc, mutList) in enumerate(mutations):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        strucFileID = strucfile_id (struc)
        strucDir = outDir / strucFileID
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'individual_list.txt'
        mutList = ['%s;' % mutList.pop(0)] + ['\n%s;' % mut for mut in mutList]
        with io.open(mutListFile, "w") as fout:
            for mut in mutList:
                fout.write(mut)
        pdbid, chainID = struc[:2]
        write_partial_structure (pdbid,
                                 [chainID],
                                 pdbDir,
                                 strucDir / (strucFileID + '.pdb'))
        write_foldx_config (strucDir / 'config_repairPDB.cfg',
                            'RepairPDB',
                            pdb_dir = '../data/%s' % strucFileID,
                            output_dir = '../data/%s' % strucFileID,
                            pdb_file = '%s.pdb' % strucFileID)
        mutParam = {k:v for k, v in default_param.items()}
        if struc in parameters:
            for k, v in parameters[struc].items():
                mutParam[k] = v
        write_foldx_config (strucDir / 'config_buildModel.cfg',
                            'BuildModel',
                            pdb_dir = '../data/%s' % strucFileID,
                            output_dir = '../data/%s' % strucFileID,
                            pdb_file = '%s_Repair.pdb' % strucFileID,
                            mutant_file = '../data/%s/individual_list.txt' % strucFileID,
                            temp = mutParam['temp'],
                            ph = mutParam['ph'],
                            ionStrength = mutParam['ionStrength'],
                            water = mutParam['water'],
                            vdwDesign = mutParam['vdwDesign'])
    print()

def produce_foldx_pssm_jobs (mutations, pdbDir, outDir, parameters = None):
    """Produce jobs for the FoldX method for ∆∆G calculation using PSSM command.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX jobs to.
        parameters (dict): parameters used for each FoldX job, otherwise default.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    default_param = {'temp':298,
                     'ph':7,
                     'ionStrength':0.05,
                     'water':'-IGNORE',
                     'vdwDesign':2}
    if not parameters:
        parameters = {}
    
    mutations = mutations.items()
    n = len(mutations)
    print('Writing FoldX PSSM files:')
    for i, (struc, mutList) in enumerate(mutations):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        strucFileID = strucfile_id (struc)
        pdbid, chainID1, chainID2 = struc[:3]
        for mut in mutList:
            mutID = '_'.join([strucFileID, mut])
            mutDir = outDir / mutID
            if not mutDir.exists():
                os.makedirs(mutDir)
            write_partial_structure (pdbid,
                                     [chainID1, chainID2],
                                     pdbDir,
                                     mutDir / (strucFileID + '.pdb'))
            write_foldx_config (mutDir / 'config_repairPDB.cfg',
                                'RepairPDB',
                                pdb_dir = '../data/%s' % mutID,
                                output_dir = '../data/%s' % mutID,
                                pdb_file = '%s.pdb' % strucFileID)
            mutParam = {k:v for k, v in default_param.items()}
            if (struc, mut) in parameters:
                for k, v in parameters[(struc, mut)].items():
                    mutParam[k] = v
            write_foldx_config (mutDir / 'config_pssm.cfg',
                                'Pssm',
                                other_cmd = ['analyseComplexChains=%s,%s' % (chainID1, chainID2),
                                             'aminoacids=%s' % mut[-1],
                                             'positions=%sa' % mut[:-1]],
                                pdb_dir = '../data/%s' % mutID,
                                output_dir = '../data/%s' % mutID,
                                pdb_file = '%s_Repair.pdb' % strucFileID,
                                temp = mutParam['temp'],
                                ph = mutParam['ph'],
                                ionStrength = mutParam['ionStrength'],
                                water = mutParam['water'],
                                vdwDesign = mutParam['vdwDesign'])
    print()

def produce_foldx_analyseComplex_jobs (structures, pdbDir, outDir, parameters = None):
    """Produce jobs for the FoldX interaction energy calculation using AnalyseComplex command.

    Args:
        structures (list): tuples of structure ID and chain letters.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX jobs to.
        parameters (dict): parameters used for each FoldX job, otherwise default.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    default_param = {'temp':298,
                     'ph':7,
                     'ionStrength':0.05,
                     'water':'-IGNORE',
                     'vdwDesign':2}
    if not parameters:
        parameters = {}
    
    n = len(structures)
    print('Writing FoldX AnalyseComplex files:')
    for i, (pdbid, chainID1, chainID2) in enumerate(structures):
        sys.stdout.write('  Structure %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
        sys.stdout.flush()
        struc = pdbid, chainID1, chainID2
        strucFileID = strucfile_id (struc)
        strucDir = outDir / strucFileID
        if not strucDir.exists():
            os.makedirs(strucDir)
        
        write_partial_structure (pdbid,
                                 [chainID1, chainID2],
                                 pdbDir,
                                 strucDir / (strucFileID + '.pdb'))
        
        write_foldx_config (strucDir / 'config_repairPDB.cfg',
                            'RepairPDB',
                            pdb_dir = '../data/%s' % strucFileID,
                            output_dir = '../data/%s' % strucFileID,
                            pdb_file = '%s.pdb' % strucFileID)
        
        strucParam = {k:v for k, v in default_param.items()}
        if struc in parameters:
            for k, v in parameters[struc].items():
                strucParam[k] = v
        
        write_foldx_config (strucDir / 'config_analyseComplex.cfg',
                            'AnalyseComplex',
                            other_cmd = ['analyseComplexChains=%s,%s' % (chainID1, chainID2)],
                            pdb_dir = '../data/%s' % strucFileID,
                            output_dir = '../data/%s' % strucFileID,
                            pdb_file = '%s_Repair.pdb' % strucFileID,
                            temp = strucParam['temp'],
                            ph = strucParam['ph'],
                            ionStrength = strucParam['ionStrength'],
                            water = strucParam['water'],
                            vdwDesign = strucParam['vdwDesign'])
    print()

def write_foldx_config (outPath,
                        command,
                        other_cmd = None,
                        pdb_dir = None,
                        output_dir = None,
                        pdb_file = None,
                        mutant_file = None,
                        temp = 298,
                        ph = 7,
                        ionStrength = 0.05,
                        water = '-IGNORE',
                        vdwDesign = 2):
    """Produce FoldX configuration file.

    Args:
        outPath (Path): file path to save FoldX configurations.
        command (str): main command to be run by FoldX.
        other_cmd (list): additional commands in string format to be run by FoldX.
        pdb_dir (Path): file directory containing PDB structures.
        output_dir (Path): file directory where FoldX results are saved.
        pdb_file (Path): path to file containing PDB structure to be processed.
        mutant_file (Path): path to file containing list of mutations if applicable.
        temp (numeric): temperature in Kelvins.
        ph (numeric): PH value.
        ionStrength (numeric): ionic strength of the solution in Moles.
        water (str): how to handle water molecules, '-CRYSTAL', '-PREDICT', '-IGNORE' or '-COMPARE'.
        vdwDesign (numeric): VDW design of the experiment, 0 very soft, 1 medium soft, 2 strong.

    """
    with io.open(outPath, "w") as fout:
        fout.write('command=%s' % command)
        if other_cmd:
            fout.write('\n' + '\n'.join(other_cmd))
        if pdb_dir:
            fout.write('\n' + 'pdb-dir=%s' % pdb_dir)
        if output_dir:
            fout.write('\n' + 'output-dir=%s' % output_dir)
        if pdb_file:
            fout.write('\n' + 'pdb=%s' % pdb_file)
        if mutant_file:
            fout.write('\n' + 'mutant-file=%s' % mutant_file)
        fout.write('\n' + 'temperature=%.1f' % temp)
        fout.write('\n' + 'pH=%.1f' % ph)
        fout.write('\n' + 'ionStrength=%f' % ionStrength)
        fout.write('\n' + 'water=%s' % water)
        fout.write('\n' + 'vdwDesign=%d' % vdwDesign)

def produce_beluga_foldx_jobs (data,
                               type,
                               outDir,
                               account = 'ctb-yxia',
                               walltime = '1-00',
                               ntasks = 1,
                               nodes = 1,
                               ntasks_per_node = 1,
                               cpus_per_task = 1,
                               mem = None,
                               mem_per_cpu = '4G',
                               outputfile = '%x-%j.out',
                               errorfile = None,
                               username = '',
                               extraCommands = None,
                               serverDataDir = '../data'):
    """Produce Beluga server job files specific to FoldX ∆∆G calculations.
        See https://docs.computecanada.ca/wiki/Béluga/en

    Args:
        data (dict or list): dict of mutations associated with each structure, or list of
                             structure-chain triplets.
        type (str): type of FoldX calculation; 'binding' for interface ∆∆G, 
                    'folding' for protein folding ∆∆G, 'ppi_energy' for interaction energy.
        outDir (Path): file directory to save jobs to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    if not outDir.exists():
        os.makedirs(outDir)
    
    print('Writing Beluga job files')
    if type is 'folding':
        for struc, _ in data.items():
            strucFileID = strucfile_id (struc)
            commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, strucFileID),
                        '../foldx -f %s/%s/config_buildModel.cfg' % (serverDataDir, strucFileID)]
            if extraCommands:
                commands = extraCommands + commands
            write_beluga_job (outDir / (strucFileID + '_job.sh'),
                              account = account,
                              walltime = walltime,
                              ntasks = ntasks,
                              nodes = nodes,
                              ntasks_per_node = ntasks_per_node,
                              cpus_per_task = cpus_per_task,
                              mem = mem,
                              mem_per_cpu = mem_per_cpu,
                              outputfile = outputfile,
                              errorfile = errorfile,
                              jobname = '%s_foldx_buildmodel_%s' % (username, strucFileID),
                              commands = commands)
    elif type is 'binding':
        for struc, mutList in data.items():
            strucFileID = strucfile_id (struc)
            for mut in mutList:
                mutID = '_'.join([strucFileID, mut])
                commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, mutID),
                            '../foldx -f %s/%s/config_pssm.cfg' % (serverDataDir, mutID)]
                if extraCommands:
                    commands = extraCommands + commands
                write_beluga_job (outDir / (mutID + '_job.sh'),
                                  account = account,
                                  walltime = walltime,
                                  ntasks = ntasks,
                                  nodes = nodes,
                                  ntasks_per_node = ntasks_per_node,
                                  cpus_per_task = cpus_per_task,
                                  mem = mem,
                                  mem_per_cpu = mem_per_cpu,
                                  outputfile = outputfile,
                                  errorfile = errorfile,
                                  jobname = '%s_foldx_pssm_%s' % (username, mutID),
                                  commands = commands)
    elif type is 'ppi_energy':
        for struc in data:
            strucFileID = strucfile_id (struc)
            commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, strucFileID),
                        '../foldx -f %s/%s/config_analyseComplex.cfg' % (serverDataDir, strucFileID)]
            if extraCommands:
                commands = extraCommands + commands
            write_beluga_job (outDir / (strucFileID + '_job.sh'),
                              account = account,
                              walltime = walltime,
                              ntasks = ntasks,
                              nodes = nodes,
                              ntasks_per_node = ntasks_per_node,
                              cpus_per_task = cpus_per_task,
                              mem = mem,
                              mem_per_cpu = mem_per_cpu,
                              outputfile = outputfile,
                              errorfile = errorfile,
                              jobname = '%s_foldx_%s' % (username, strucFileID),
                              commands = commands)

def read_foldx_results (inDir, type = 'binding'):
    """Read mutation ∆∆G results produced by FoldX commands.

    Args:
        inDir (Path): file directory containing FoldX results.
        type (str): type of FoldX calculation; 'binding' for interface ∆∆G, 
                    'folding' for protein folding ∆∆G, 'ppi_energy' for interaction energy.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    if type is 'folding':
        return read_foldx_buildmodel_results (inDir)
    elif type is 'binding':
        return read_foldx_pssm_results (inDir)
    elif type is 'ppi_energy':
        return read_foldx_analyseComplex_results (inDir)

def read_foldx_buildmodel_results (inDir):
    """Read mutation ∆∆G results produced by FoldX BuildModel command.

    Args:
        inDir (Path): file directory containing FoldX results.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    processed, unprocessed = {}, {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucFileID in strucDir:
        struc = tuple((solve_pdbfile_id(strucFileID)).split('_'))
        if re.match(r'\D\S\d+\D', struc[-1]):
            struc = struc[:-1]
        mutListFile = inDir / strucFileID / 'individual_list.txt'
        with io.open(mutListFile, "r") as f:
            mutList = list( map(str.strip, f.read().split(';')) )
        mutList.remove('')
        
        resultFile = None
        for filename in os.listdir(inDir / strucFileID):
            if filename.startswith('Average'):
                resultFile = filename
        
        results = []
        if resultFile:
            resultFile = inDir / strucFileID / resultFile
            headers = ["Pdb", "SD", "total energy", "Backbone Hbond", "Sidechain Hbond",
                       "Van der Waals", "Electrostatics", "Solvation Polar", "Solvation Hydrophobic", 
                       "Van der Waals clashes", "entropy sidechain", "entropy mainchain", 
                       "sloop_entropy", "mloop_entropy", "cis_bond", "torsional clash", "backbone clash", 
                       "helix dipole", "water bridge", "disulfide", "electrostatic kon", 
                       "partial covalent bonds", "energy Ionisation", "Entropy Complex"]
            with io.open(resultFile, "r") as f:
                for line in f:
                    linesplit = list(map(str.strip, line.split('\t')))
                    if linesplit == headers:
                        break
                for line in f:
                    linesplit = list(map(str.strip, line.split('\t')))
                    if len(linesplit) > 3:
                        results.append(linesplit[2])
                if len(results) == len(mutList):
                    for mut, ddg in zip(mutList, results):
                        processed[struc + (mut,)] = float(ddg)
        
        if not results:
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                for mut in mutList:
                    unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def read_foldx_pssm_results (inDir):
    """Read mutation ∆∆G results produced by FoldX PSSM command.

    Args:
        inDir (Path): file directory containing FoldX results.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    processed, unprocessed = {}, {}
    strucDirs = os.listdir(inDir)
    strucDirs = [dir for dir in strucDirs if os.path.isdir(inDir / dir)]
    for strucFileID in strucDirs:
        strucDir = inDir / strucFileID
        struc = tuple((solve_pdbfile_id(strucFileID)).split('_'))
        if len(struc) > 3:
            struc = struc[:-1]
        
        mutListFile = resultFile = None
        for filename in os.listdir(strucDir):
            if filename.startswith('individual_list'):
                mutListFile = strucDir / filename
                with io.open(mutListFile, "r") as f:
                    mutList = list(map(str.strip, f.read().split(';')))
                    mutList.remove('')
            elif filename.startswith('_'.join(('PSSM',) + struc)):
                resultFile = strucDir / filename
        
        if resultFile:
            with io.open(resultFile, "r") as f:
                mutResidues = list(map(str.strip, f.readline().strip().split('\t')))
                for line in f:
                    ddgList = list(map(str.strip, line.split('\t')))
                    if len(ddgList) > 1:
                        wt = ddgList[0]
                        for mt, ddg in zip(mutResidues, ddgList[1:]):
                            processed[struc + (wt + mt,)] = float(ddg) if len(ddg) else 'X'
                            
            for mut in mutList:
                if struc + (mut,) not in processed:
                    processed[struc + (mut,)] = 'X'
        else:
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                for mut in mutList:
                    unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def read_foldx_analyseComplex_results (inDir):
    """Read interaction energy results produced by FoldX AnalyseComplex command.

    Args:
        inDir (Path): file directory containing FoldX results.

    Returns:
        dict: processed structure interaction energies.

    """
    processed = {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucFileID in strucDir:
        struc = tuple((solve_pdbfile_id(strucFileID)).split('_'))
        
        resultFile = None
        for filename in os.listdir(inDir / strucFileID):
            if filename.startswith('Interaction_' + strucFileID) and filename.endswith('.fxout'):
                resultFile = filename
        
        if resultFile:
            resultFile = inDir / strucFileID / resultFile
            firstHeaders = ["Pdb", "Group1", "Group2", "IntraclashesGroup1", 
                               "IntraclashesGroup2", "Interaction Energy"]
            with io.open(resultFile, "r") as f:
                for line in f:
                    if line.startswith('\t'.join(firstHeaders)):
                        line = next(f, '')
                        if line:
                            linesplit = list(map(str.strip, line.split('\t')))
                            if len(linesplit) > 5:
                                processed[struc] = float(linesplit[5])
        
        if struc not in processed:
            processed[struc] = 'X'
    
    return processed

def read_protein_mutation_ddg (inPath, type = 'binding'):
    """Read protein mutation ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations with ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                if linesplit[8] is not 'X':
                    protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = linesplit
                    if type is 'binding':
                        k = protein, partner, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, ch_partner, ch_mut, float(ddg)
                    elif type is 'folding':
                        k = protein, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, ch_mut, float(ddg) 
                    if k not in ddgDict:
                        ddgDict[k] = val
    return ddgDict

def read_chain_mutation_ddg (inPath, type = 'binding'):
    """Read PDB chain mutation ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations with ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = linesplit[:9]
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                ddgDict[k] = ddg
    return ddgDict

def copy_mutation_ddg (inPath1, inPath2, outPath, type = 'binding'):
    """Transfer mutation ∆∆G values from one file to another file based on similar mutation 
        structure mappings.

    Args:
        inPath1 (Path): path to file containing mutation ∆∆G values.
        inPath2 (Path): path to file containing mutations with missing ∆∆G values.
        outPath (Path): file path to save a copy of file in inPath2 with updated ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    """
    ddg = read_chain_mutation_ddg (inPath1, type)
    write_mutation_ddg_tofile (ddg, inPath2, outPath, type)

def write_mutation_ddg_tofile (ddg, inPath, outPath, type = 'binding'):
    """Update file with mutation ∆∆G values.

    Args:
        ddg (dict): mutation ∆∆G values.
        inPath (Path): path to file whose mutations will be updated with ∆∆G values.
        outPath (Path): file path to save a copy of file in inPath with updated ∆∆G values.
        type (str): type of ∆∆G values; 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    """
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        fout.write(f.readline().strip() + '\n')
        for line in f:
            strsplit = list(map(str.strip, line.split('\t')))
            if len(strsplit) == 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner = strsplit
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                if k in ddg:
                    strsplit.append(str(ddg[k]))
            fout.write('\t'.join(map(str, strsplit)) + '\n')

def write_ppi_energy_tofile (energy, inPath, outPath):
    """Update PPI energy values in file.

    Args:
        energy (dict): PPI structure energy values.
        inPath (Path): path to file whose PPIs will be updated with energy values.
        outPath (Path): file path to save updated energy values.

    """
    ppis = pd.read_table (inPath, sep='\t')
    for i, row in ppis.iterrows():
        if row.Interaction_energy == '-':
            k = tuple([row.Complex_ID] + sorted([row.Chain_1, row.Chain_2]))
            if k in energy:
                ppis.loc[i, "Interaction_energy"] = energy[k]
    ppis.to_csv (outPath, index=False, sep='\t')

def append_mutation_ddg_files (inPath1, inPath2, outPath):
    """Append two ∆∆G files together and save to anothor file.

    Args:
        inPath1 (Path): path to first input file.
        inPath2 (Path): path to second input file.
        outPath (Path): file path to save appended output to.

    """
    with io.open(outPath, "w") as fout:
        with io.open(inPath1, "r", encoding="utf-8") as f1:
            fout.write(f1.readline().strip() + '\n')
            for line in f1:
                fout.write(line)
        with io.open(inPath2, "r", encoding="utf-8") as f2:
            next(f2)
            for line in f2:
                fout.write(line)

def strucfile_id (struc):
    """Return structure file ID from structure, chains and mutation ID tuple.

    Args:
        struc (tuple): structure, chains and mutation ID tuple.

    Returns:
        str: structure file ID.

    """
    if re.match(r'\D\S\d+\D', struc[-1]):
        return pdbfile_id ('_'.join(struc[:-1])) + '_' + struc[-1]
    else:
        return pdbfile_id ('_'.join(struc))
