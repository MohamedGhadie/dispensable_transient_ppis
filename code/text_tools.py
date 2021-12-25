#----------------------------------------------------------------------------------------
# Modules for text processing.
#----------------------------------------------------------------------------------------

import io
import time
import pandas as pd
import copy as cp
from pathlib import Path
from Bio import Seq, SeqIO

def parse_fasta_file (inPath, outPath):
    """Read sequences from fasta file.

    Args:
        inPath (Path): path to FASTA file containing sequences.
        outPath (Path): path to save tab-delimited sequence file to.

    """
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(['ID', 'Length', 'Sequence'])  + '\n')
        for i, row in enumerate(s):
            fout.write('\t'.join([row.id, str(len(row.seq)), str(row.seq)])  + '\n')
    print('%d sequences extracted from file ' % (i + 1) + str(inPath) +
                ' and written to file ' + str(outPath))

def reduce_fasta_headers (inPath, delimiter, first, last, outPath):
    """Reduce headers in FASTA file to specific subset of elements.

    Args:
        inPath (Path): path to FASTA file.
        delimiter (str): delimiter used to split FASTA headers.
        first (numeric): position of the first header element to include.
        last (numeric): position of the last header element to include.
        outPath (Path): path to save FASTA file with reduced headers to.
        
    """
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    with io.open(outPath, "w") as fout:
        for _, row in enumerate(s):
            idsplit = list(map(str.strip, row.id.split(delimiter)))
            id = '|'.join(idsplit[ first - 1 : last ])
            fout.write('>' + id + '\n')
            fout.write(str(row.seq) + '\n')

def write_fasta_file (data, idCol, seqCol, outPath):
    """Produce FASTA file from DataFrame columns.

    Args:
        data (DataFrame): table to be written to FASTA file.
        idCol (str): name of column containing IDs to be written as FASTA headers.
        seqCol (str): name of column containing sequences to be written to FASTA file.
        outPath (Path): path to save FASTA file to.
        
    """
    with io.open(outPath, "w") as fout:
        for _, row in data.iterrows():
            fout.write('>' + row[idCol] + '\n')
            fout.write(row[seqCol] + '\n')

def produce_item_list (inPath, cols, outPath):
    """Write unique items from specific table columns to file.

    Args:
        inPath (Path): path to file containing table to be written.
        cols (list, str): names of columns to be saved to file.
        outPath (Path): path to save list of items to.

    """
    df = pd.read_table(inPath, sep='\t')
    items = sorted(set(df[cols].values.flatten()))
    with open(outPath, 'w') as fout:
        for item in items:
            fout.write('%s\n' % item)

def read_list_table (inPath, cols, dtyp, delm = '\t'):
    """Read a table from a file into a DataFrame with specified column entries converted 
        to lists.

    Args:
        inPath (Path): path to file containing table to be read.
        cols (list, str): names of columns whose entries are to be converted to lists.
        dtyp (list, datatype): data types to convert each column list elements to.
        delm (str): delimiter used to read table.

    """
    df = pd.read_table(inPath, sep = delm)
    if not isinstance(cols, (list, tuple)):
        cols = [cols]
        dtyp = [dtyp]
    for col, typ in zip(cols, dtyp):
        strcol = df[col].apply(str)
        df[col] = strcol.apply(lambda x: list(map(typ, map(str.strip, x.split(',')))))
    return df

def write_list_table (df, cols, outPath, delm = '\t'):
    """Write a DataFrame to text file with specified columns with list entries converted 
        to comma-separated strings.

    Args:
        df (DataFrame): table to be writen to file.
        cols (list, str): names of columns whose list entries are to be converted to strings.
        outPath (Path): path to save table to.
        delm (str): delimiter used to write table to file.

    """
    table = cp.deepcopy(df)
    if isinstance(cols, (list, tuple)):
        for col in cols:
            table[col] = table[col].apply(lambda x: ','.join(map(str, x)))
    else:
        table[cols] = table[cols].apply(lambda x: ','.join(map(str, x)))
    table.to_csv(outPath, index=False, sep=delm)

def write_guillimin_job (outPath,
                         nodes = 1,
                         ppn = 1,
                         pmem = 7700,
                         walltime = '1:00:00:00',
                         outputfile = 'outputfile',
                         errorfile = 'errorfile',
                         rapid = None,
                         jobid = None,
                         commands = None):
    """Create a job file to be run by McGill HPC server.
        See http://www.hpc.mcgill.ca

    Args:
        outPath (Path): path to save job file to.
        nodes (numeric): number of server nodes to be allocated.
        ppn (numeric): total number of CPU cores to be allocated.
        pmem (numeric): default random access memory (RAM) in MB to be reserved per core.
        walltime (str: 'days:hr:min:sec'): maximum time allowed for job to run.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        rapid (str): resource allocation project identifier (RAPid).
        jobid ('str'): job name.
        commands (list): additional commands in string format to be written to job file.

    """
    with io.open(outPath, "w") as fout:
        fout.write('#!/bin/bash')
        fout.write('\n' + '#PBS -l nodes=%d:ppn=%d,pmem=%dm,walltime=%s' % (nodes, ppn, pmem, walltime))
        if rapid:
            fout.write('\n' + '#PBS -A %s' % rapid)
        fout.write('\n' + '#PBS -o %s' % outputfile)
        fout.write('\n' + '#PBS -e %s' % errorfile)
        if jobid:
            fout.write('\n' + '#PBS -N %s' % jobid)
        fout.write('\n\n' + 'cd $PBS_O_WORKDIR')
        if commands:
            for cmd in commands:
                fout.write('\n' + cmd)

def write_beluga_job (outPath,
                      account = 'ctb-yxia',
                      walltime = '1-00',
                      ntasks = 1,
                      nodes = 1,
                      ntasks_per_node = 1,
                      cpus_per_task = 1,
                      mem = None,
                      mem_per_cpu = '4000M',
                      outputfile = '%x-%j.out',
                      errorfile = None,
                      jobname = None,
                      commands = None):
    """Create a job file to be run on Compute Canada Beluga cluster.
        See https://docs.computecanada.ca/wiki/Beluga/en

    Args:
        outPath (Path): path to save job file to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run in the format days:hr:min:sec.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (numeric): memory per node.
        mem_per_cpu (numeric): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        jobname (str): job name.
        commands (list): commands in string format to be executed by job.

    """
    with io.open(outPath, "w") as fout:
        fout.write('#!/bin/bash')
        fout.write('\n' + '#SBATCH --account=%s' % account)
        fout.write('\n' + '#SBATCH --time=%s' % walltime)
        fout.write('\n' + '#SBATCH --ntasks=%d' % ntasks)
        fout.write('\n' + '#SBATCH --nodes=%d' % nodes)
        fout.write('\n' + '#SBATCH --ntasks-per-node=%d' % ntasks_per_node)
        fout.write('\n' + '#SBATCH --cpus-per-task=%d' %cpus_per_task)
        if mem:
            fout.write('\n' + '#SBATCH --mem=%s' % mem)
        fout.write('\n' + '#SBATCH --mem-per-cpu=%s' % mem_per_cpu)
        fout.write('\n' + '#SBATCH --output=%s' % outputfile)
        if errorfile:
            fout.write('\n' + '#SBATCH --error=%s' % errorfile)
        if jobname:
            fout.write('\n' + '#SBATCH --job-name=%s' % jobname)
        if commands:
            fout.write('\n')
            for cmd in commands:
                fout.write('\n' + cmd)

def parse_blast_file (inPath,
                      outPath,
                      encoding = 'us-ascii',
                      pausetime = 0):
    """Read alignments from BLAST output file.

    Args:
        inPath (Path): path to BLAST output file.
        outPath (Path): path to save tab-delimited alignment table to.
        encoding (str): encoding for reading BLAST output file.
        pausetime (numeric): time in seconds to pause after processing 50 million lines.

    """
    allignKeys = ['Query',
                  'Qlen',
                  'Subject',
                  'Slen',
                  'Score',
                  'Expect',
                  'Identities',
                  'Positives',
                  'Gaps',
                  'Qstart',
                  'Qend',
                  'Sstart',
                  'Send',
                  'Qseq',
                  'Sseq',
                  'Match']
    allignInfo = {k:k for k in allignKeys}
    currentInfo = {k:'-' for k in allignKeys[:4]}
    previous = '-'
    
    with io.open(inPath, "r", encoding=encoding) as f:
        for numLines, _ in enumerate(f):
            pass
            
    with io.open(inPath, "r", encoding=encoding) as fin, io.open(outPath, "w") as fout:
        for i, line in enumerate(fin):
            if not (i % 1e6):
                print('\t\t%d million out of %g million lines parsed' % (i // 1e6, numLines / 1e6))
            if not ((i+1) % 50e6):
                print('\t\tpausing for %d seconds' % pausetime)
                time.sleep(pausetime)
            line = line.replace('\n','')
            stripline = line.strip()
            if stripline.startswith('Query='):
                linesplit = list(map(str.strip, stripline.split('=')))
                currentInfo["Query"] = linesplit[1]
                previous = 'query'
            elif stripline.startswith('>'):
                currentInfo["Subject"] = stripline[1:].strip()
                previous = 'subject'
            elif stripline.startswith('Length'):
                linesplit = list(map(str.strip, stripline.split('=')))
                if previous=='query':
                    currentInfo["Qlen"] = linesplit[1]
                elif previous=='subject':
                    currentInfo["Slen"] = linesplit[1]
            elif any(x in stripline for x in ['Score =', 'Except =', 'Identities =', 'Positives =', 'Gaps =']):
                linesplit = list(map(str.strip, stripline.split(',')))
                entrysplit = [list(map(str.strip, x.split('='))) for x in linesplit]
                for entry in entrysplit:
                    if len(entry) == 2:
                        if entry[0] == 'Score':
                            allignInfo.update({k:allignInfo[k].replace('\n','')
                                              for k in ['Qseq','Sseq','Match']})
                            info = [allignInfo[k] for k in allignKeys]
                            fout.write('\t'.join(info) + '\n')
                            allignInfo.update({k:'-' for k in allignKeys})
                            allignInfo.update({k:'' for k in ['Qseq','Sseq','Match']})
                            allignInfo.update({k:currentInfo[k] for k in currentInfo})
                            allignInfo["Score"] = entry[1].split()[0]
                        elif entry[0] == 'Expect':
                            allignInfo["Expect"] = entry[1]
                        elif entry[0] in ['Identities', 'Positives', 'Gaps']:
                            allignInfo[entry[0]] = entry[1].split('/')[0].strip()
            elif stripline.startswith('Query'):
                linesplit = stripline.split()
                if allignInfo["Qstart"] == '-':
                    allignInfo["Qstart"] = linesplit[1]
                allignInfo["Qend"] = linesplit[3]
                allignInfo["Qseq"] = allignInfo["Qseq"] + linesplit[2]
                allignLen = len(linesplit[2])
                previous = 'query seq'
            elif stripline.startswith('Sbjct'):
                linesplit = stripline.split()
                if allignInfo["Sstart"] == '-':
                    allignInfo["Sstart"] = linesplit[1]
                allignInfo["Send"] = linesplit[3]
                allignInfo["Sseq"] = allignInfo["Sseq"] + linesplit[2]
                previous = 'subject seq'
            elif (len(line) > 0) and (previous == 'query seq'):
                allignInfo["Match"] = allignInfo["Match"] + line[ - allignLen : ]
        info = [allignInfo[k] for k in allignKeys]
        fout.write('\t'.join(info) + '\n')
    print('\t%d lines parsed from ' % (i + 1) + str(inPath) +  
          'and written to file ' + str(outPath))
