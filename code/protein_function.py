#----------------------------------------------------------------------------------------
# Modules for computations on protein function.
#----------------------------------------------------------------------------------------

import os
import io
import sys
import pickle
import pandas as pd
import numpy as np
import subprocess
from scipy.stats.stats import pearsonr
from simple_tools import valid_uniprot_id, is_numeric, hamming_dist
from geo_tools import read_gds_softfile

def partner_sim (p1, p2, partners):
    """Calculate the fraction of interaction partners shared by two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        partners (dict): dictionary containing set of partners for each protein.
    
    Returns:
        float: fraction of shared partners if at least one protein has a partner, otherwise NaN.

    """
    partners_1 = (partners[p1] if p1 in partners else set()) - {p2}
    partners_2 = (partners[p2] if p2 in partners else set()) - {p1}
    total = len(partners_1 | partners_2)
    if total > 0:
        return len(partners_1 & partners_2) / total
    else:
        return np.nan

def go_sim (p1, p2, goAssoc):
    """Calculate the fraction of gene ontology (GO) terms shared by two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        goAssoc (dict): dictionary containing list of GO terms for each protein.
    
    Returns:
        float: fraction of shared GO terms if at least one protein is associated with a 
                GO term, otherwise NaN.
    
    """
    go1 = (goAssoc[p1] if p1 in goAssoc else set())
    go2 = (goAssoc[p2] if p2 in goAssoc else set())
    if (len(go1) > 0) or (len(go2) > 0):
        return len(go1 & go2) / len(go1 | go2)
    else:
        return np.nan

def coexpr (p1, p2, expr, minPts = 3, method = 'pearson_corr', singleExp = True):
    """Calculate tissue co-expression for two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        expr (dict): dictionary containing tissue expression values for each protein.
        minPts (int): minimum required number of samples with expression levels 
                            defined for both proteins.
        method (str): method used to calculate coexpression. Set to 'pearson_corr' to return 
                        Pearson's correlation coefficient, or 'hamming_dist' to return 
                        1 - hamming_distance / length of valid columns.
        singleExp (bool): set to False if protein expression dictionary has multiple 
                            expression datasets per protein. In this case, the mean of 
                            coexpression across all datasets is returned.
    
    Returns:
        float
    
    """
    if singleExp:
        if (p1 in expr) and (p2 in expr):
            e1, e2 = expr[p1], expr[p2]
            not_nan = (np.isnan(e1) == False) & (np.isnan(e2) == False)
            numPts = sum(not_nan)
            if numPts >= minPts:
                e1, e2 = e1[not_nan], e2[not_nan]
                if method is 'pearson_corr':
                    if (len(set(e1)) > 1) and (len(set(e2)) > 1):
                        corr, p = pearsonr(e1, e2)
                        return corr
                elif method is 'hamming_dist':
                    return 1 - hamming_dist (e1, e2) / numPts
    else:
        all = []
        if (p1 in expr) and (p2 in expr):
            for dataset in expr[p1]:
                if dataset in expr[p2]:
                    subExpr = {p1:expr[p1][dataset], p2:expr[p2][dataset]}
                    all.append(coexpr (p1,
                                       p2,
                                       subExpr,
                                       minPts = minPts,
                                       method = method,
                                       singleExp = True))
            all = [c for c in all if not np.isnan(c)]
            if all:
                return np.mean(all)
    return np.nan

def expr_log_diff (p1,
                   p2,
                   expr,
                   minPts = 1,
                   logBase = 10,
                   method = 'mean',
                   singleExp = True,
                   average = False):

    if singleExp:
        if (p1 in expr) and (p2 in expr):
            e1, e2 = expr[p1], expr[p2]
            not_nan = (np.isnan(e1) == False) & (np.isnan(e2) == False)
            e1, e2 = e1[not_nan], e2[not_nan]
            if logBase:
                not_zero = (e1 > 0) & (e2 > 0)
                e1, e2 = e1[not_zero], e2[not_zero]
                e1 = np.log10(e1) / np.log10(logBase)
                e2 = np.log10(e2) / np.log10(logBase)
            if len(e1) >= minPts:
                if method is 'mean':
                    return np.mean(np.abs(e1 - e2))
        return np.nan
    else:
        all = []
        if (p1 in expr) and (p2 in expr):
            for dataset in expr[p1]:
                if dataset in expr[p2]:
                    subExpr = {p1:expr[p1][dataset], p2:expr[p2][dataset]}
                    all.append(expr_log_diff (p1,
                                              p2,
                                              subExpr,
                                              minPts = minPts,
                                              logBase = logBase,
                                              method = method,
                                              singleExp = True))
            return np.nanmean(all) if average else all
        return np.nan if average else all

def expr_ratio (p1, p2, expr, minPts = 1, method = 'mean'):
    """Calculate the log of expression ratio for two proteins across multiple experiments.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        expr (dict): dictionary containing expression values for each protein.
        minPts (int): 
        method (str): method used to merge expression ratio from all experiments.
    
    Returns:
        float
    
    """
    if (p1 in expr) and (p2 in expr):
        e1, e2 = expr[p1], expr[p2]
        not_nan = (np.isnan(e1) == False) & (np.isnan(e2) == False)
        e1, e2 = e1[not_nan], e2[not_nan]
        not_zero = (e1 > 0) & (e2 > 0)
        if sum(not_zero) >= minPts:
            e1, e2 = e1[not_zero], e2[not_zero]
            if method is 'mean':
                ratio = np.true_divide(e1, e2)
                ratio = [max(r, 1/r) for r in ratio] 
                logRatio = list(map(np.log10, ratio))
                return np.mean(logRatio)
    return np.nan

def expr_ratio_symbolic (p1,
                         p2,
                         expr,
                         minPts = 1,
                         values = None,
                         high = None,
                         method = 'majority'):
    
    if not values:
        values = ['Low', 'Medium', 'High']
    if not high:
        high = [{'Low', 'High'}]
    
    if (p1 in expr) and (p2 in expr):
        expr1, expr2 = np.array(expr[p1]), np.array(expr[p2])
        valid = [(e1 in values) and (e2 in values) for e1, e2 in zip(expr1, expr2)]
        if sum(valid) >= minPts:
            expr1, expr2 = expr1[valid], expr2[valid]
            if method is 'majority':
                ratio = []
                for e1, e2 in zip(expr1, expr2):
                    ratio.append('high' if {e1, e2} in high else 'low')
                numHigh, numLow = ratio.count('high'), ratio.count('low')
                if (numHigh + numLow) >= minPts:
                    return 'high' if numHigh >= numLow else 'low'
    return '-'

def produce_protein_go_dictionaries (inPath,
                                     GO_outPath,
                                     MF_outPath,
                                     BP_outPath,
                                     CC_outPath):
    """Make dictionaries of protein gene ontology (GO) associations, with each root GO 
        in a separate dictionary.

    Args:
        inPath (Path): path to file containing all protein GO associations.
        GO_outPath (Path): file path to save dict of all GO terms.
        MF_outPath (Path): file path to save dict of molecular function (F) terms.
        BP_outPath (Path): file path to save dict of biological process (P) terms.
        CC_outPath (Path): file path to save dict of cellular component (C) terms.

    """
    produce_protein_go_dict (inPath, GO_outPath)
    produce_protein_go_dict (inPath, MF_outPath, aspect = 'F')
    produce_protein_go_dict (inPath, BP_outPath, aspect = 'P')
    produce_protein_go_dict (inPath, CC_outPath, aspect = 'C')

def produce_protein_go_dict (inPath, outPath, aspect = None):
    """Make a dictionary of protein gene ontology (GO) associations.

    Args:
        inPath (Path): path to file containing protein GO associations.
        outPath (Path): file path to save output dict to.
        aspect (str): GO aspects to select: P, F, C, or all if not provided.

    """
    goa = pd.read_table(inPath, header=None, sep="\t")
    goa.columns = ["db",
                   "id",
                   "gene_symbol",
                   "qualifier",
                   "go_id",
                   "db_reference",
                   "evid_code",
                   "with_or_from",
                   "go_aspect",
                   "gene_name",
                   "synonym",
                   "product_type",
                   "taxon",
                   "date",
                   "source",
                   "extension",
                   "variant"]
    goa = goa.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    goa = goa[goa["qualifier"].apply(lambda x: ('NOT' not in x) if isinstance(x, str) else True) &
              (goa["taxon"] == 'taxon:9606') &
              goa["id"].apply(valid_uniprot_id) &
              (goa["product_type"] == 'protein')]
    if aspect is not None:
        goa = goa[goa["go_aspect"] == aspect]
    go = {}
    for _, row in goa.iterrows():
        if row.id in go:
            go[row.id].add(row.go_id)
        else:
            go[row.id] = {row.go_id}
    with open(outPath, 'wb') as fOut:
        pickle.dump(go, fOut)

def get_all_go_terms (inPath, aspect = None):
    """Return a list of all GO terms.

    Args:
        inPath (Path): path to file containing protein GO associations.
        aspect (str): GO aspects to select: P, F, C, or all if not provided.

    Returns:
        list

    """
    goa = pd.read_table(inPath, header=None, sep="\t")
    goa.columns = ["db",
                   "id",
                   "gene_symbol",
                   "qualifier",
                   "go_id",
                   "db_reference",
                   "evid_code",
                   "with_or_from",
                   "go_aspect",
                   "gene_name",
                   "synonym",
                   "product_type",
                   "taxon",
                   "date",
                   "source",
                   "extension",
                   "variant"]
    goa = goa.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    if aspect is None:
        return list(goa["go_id"])
    else:
        return list(goa.loc[goa["go_aspect"] == aspect, "go_id"])

def produce_fastsemsim_protein_gosim_dict (inPath,
                                           outPath,
                                           task = 'SS',
                                           ont_type = 'GeneOntology',
                                           sim_measure = 'SimGIC',
                                           mix_method = 'BMA',
                                           cutoff = -1,
                                           remove_nan = True,
                                           query_mode = 'pairs',
                                           query_ss_type = 'obj',
                                           ac_species = 'human',
                                           ont_root = 'biological_process',
                                           ont_ignore = None,
                                           ec_ignore = None,
                                           verbosity = '-vv',
                                           ontologyFile = 'go-basic.obo',
                                           annotationFile = None,
                                           paramOutFile = 'fastsemsim_parameters',
                                           fastsemsimOutFile = 'fastsemsim_output'):
    """Calculate gene ontology (GO) similarity using fastsemsim library.
        See https://pythonhosted.org/fastsemsim/
    Args:
        inPath (Path): path to file containing list of proteins (pairs or singles).
        outPath (Path): path to file where protein GO similarity dictionary is saved to.
        task (str): calculation to perform.
        ont_type (str): type of ontology.
        sim_measure (str): measure used to calculate GO similarity.
        mix_method (str): mixing strategy used by pairwise semantic similarity measures.
        cutoff (numeric): filter out GO similarity results below this cutoff.
        remove_nan (bool): remove NaN values from results.
        query_mode (str): input file format; either pairs or list of proteins.
        query_ss_type (str): query type; ex objects annotated with ontology terms.
        ac_species (str): species name.
        ont_root (str): root ontology; biological_process, molecular function or cellular component.
        ont_ignore (list): relationships to ignore.
        ec_ignore (list): evidence codes to ignore.
        verbosity (str): verbosity level.
        ontologyFile (Path): path to gene ontology file.
        annotationFile (Path): path to gene ontology associations.
        paramOutFile (Path): path to save Fastsemsim input parameters.
        fastsemsimOutFile (Path): path to save Fastsemsim output.
    
    """
    if annotationFile:
        cmd = ['fastsemsim', '--ac_file', annotationFile]
    else:
        cmd = ['fastsemsim', '--ac_species', ac_species]
    cmd += [verbosity, '--task', task, '-o', ont_type, '--ontology_file', ontologyFile, 
            '--query_input', 'file', '--query_mode', query_mode, '--query_ss_type', query_ss_type, 
            '--query_file', inPath, '--tss', sim_measure, '--tmix', mix_method, '--root', 
            ont_root, '--cut', str(cutoff), '--remove_nan', '--save_params', paramOutFile, 
            '--output_file', fastsemsimOutFile]
    if ec_ignore:
        for ec in ec_ignore:
            cmd += ['--ignore_EC', ec]
    if ont_ignore:
        for ont in ont_ignore:
            cmd += ['--ontology_ignore', ont]
    if remove_nan:
        cmd.append('--remove_nan')
    print(' '.join(map(str, cmd)))
    subprocess.run(cmd)
    table = pd.read_table(fastsemsimOutFile, sep='\t')
    gosim = {tuple(sorted((p1, p2))):sim for p1, p2, sim in table.values}
    with open(outPath, 'wb') as fOut:
        pickle.dump(gosim, fOut)

def produce_illumina_expr_dict (inPath,
                                uniprotIDmapFile,
                                outPath,
                                normalize = False):
    """Make a dictionary of protein tissue expression from Illumina Body Map dataset.

    Args:
        inPath (Path): path to file containing Illumina Body Map tissue expression data.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.
        normalize (bool): if True, normalize gene expression across tissues.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    expr = pd.read_table(inPath, comment='#', sep="\t")
    
    e = {'Tissues':list(expr.columns[2:])}
    expr['Protein'] = expr['Gene Name'].apply(lambda x: uniprotID[x] if x in uniprotID else '-')
    expr = expr[expr['Protein'] != '-'].reset_index(drop=True)
    
    for _, row in expr.iterrows():
       if row['Gene Name'] in uniprotID:
            values = np.array(row[e['Tissues']].values, dtype=float)
            if normalize:
                e[row.Protein] = (values - values.mean()) / values.std()
            else:
                e[row.Protein] = values
    with open(outPath, 'wb') as fOut:
        pickle.dump(e, fOut)

def produce_gtex_expr_dict (inDir,
                            uniprotIDmapFile,
                            outPath,
                            uniprotIDlistFile = None):
    """Make a dictionary of protein tissue expression from the GTEx dataset.

    Args:
        inDir (Path): file directory containing GTEx tissue expression data files.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.
        uniprotIDlistFile (Path): path to file containing list of UniProt IDs.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    if uniprotIDlistFile:
        with open(uniprotIDlistFile, 'r') as f:
            uniprotIDlist = list(set(f.read().split()))
    else:
        uniprotIDlist = list(uniprotID.values())
    tissueExpr = {k:[] for k in uniprotIDlist}
    filenames = [f for f in os.listdir(inDir) if f.endswith('.bed') and not f.startswith('.')]
    for i, filename in enumerate(filenames):
        print('processing file %d of %d: %s' % (i + 1, len(filenames), filename))
        expr = {}
        inPath = inDir / filename
        with io.open(inPath, "r", errors = 'replace') as f:
            next(f)
            for j, line in enumerate(f):
                linesplit = list(map(str.strip, line.strip().split('\t')))
                if len(linesplit) > 4:
                    chr, start, end, gene_id = linesplit[:4]
                    gene_id = gene_id.split('.')[0]
                    if gene_id in uniprotID:
                        id = uniprotID[gene_id]
                        expr[id] = np.nanmean([float(e) for e in linesplit[4:] if is_numeric(e)])
        for id in uniprotIDlist:
            tissueExpr[id].append(expr[id] if id in expr else np.nan)
        print('%d lines processed' % (j + 1))
    for k, v in tissueExpr.items():
        tissueExpr[k] = np.array(v)
    with open(outPath, 'wb') as fOut:
        pickle.dump(tissueExpr, fOut)

def produce_hpa_expr_dict (inPath, uniprotIDmapFile, outPath):
    """Make a dictionary of protein tissue expression from the HPA dataset.

    Args:
        inPath (Path): path to file containing HPA tissue expression.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    matrix = pd.read_table(inPath, sep='\t')
    matrix = matrix [matrix['Reliability'] != 'Uncertain']
    matrix["unique_tissue"] = matrix["Tissue"] + '_' + matrix["Cell type"]
    matrix.rename(columns={"Gene name":"Gene_name"}, inplace=True)
    
    tissue_expr = {}
    IDs = set()
    for _, row in matrix.iterrows():
        if row.Gene in uniprotID:
            id = uniprotID[row.Gene]
        elif row.Gene_name in uniprotID:
            id = uniprotID[row.Gene_name]
        else:
            id = row.Gene_name
        IDs.add(id)
        tissue_expr[(id, row.unique_tissue)] = row.Level
    
    unique_tissues = list(set(matrix["unique_tissue"]))
    expr_vectors = {}
    for id in IDs:
        expr = []
        for tissue in unique_tissues:
            k = id, tissue
            expr.append(tissue_expr[k] if k in tissue_expr else '-')
        expr_vectors[id] = expr
        
    with open(outPath, 'wb') as fOut:
        pickle.dump(expr_vectors, fOut)

def produce_fantom5_expr_dict (inPath,
                               uniprotIDmapFile,
                               outPath,
                               normalize = False,
                               sampleTypes = None,
                               sampleTypeFile = None,
                               uniprotIDlistFile = None):
    """Make a dictionary of protein tissue expression from the Fantom5 dataset.

    Args:
        inDir (Path): file directory containing Fantom5 tissue expression data files.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.
        normalize (bool): if True, normalize gene expression across tissues.
        sampleTypes (str): type of sample; ex tissues.
        sampleTypeFile (Path): path to Fantom5 sample type spreadsheet.
        uniprotIDlistFile (Path): path to file containing list of UniProt IDs for output.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    if uniprotIDlistFile:
        with open(uniprotIDlistFile, 'r') as f:
            uniprotIDlist = list(set(f.read().split()))
    else:
        uniprotIDlist = list(uniprotID.values())
    
    if sampleTypes:
        sampleCategory = pd.read_excel(sampleTypeFile)
        if isinstance(sampleTypes, str):
            sampleTypes = [sampleTypes]
        selcols = sampleCategory.loc[sampleCategory["Sample category"].apply(lambda x: x in sampleTypes), "FF ontology id"]
        selcols = list(set(selcols.values))
    else:
        selcols = None
    
    with io.open(inPath, "r") as f:
        while True:
            line = f.readline()
            if line.startswith('##ParemeterValue[genome_assembly]='):
                geneAssembly = line.strip().split('=')[1]
                break
    
    with io.open(inPath, "r") as f:
        i, hgncIndex, tpmIndex, tissues = -1, -1, [], []
        while True:
            line = f.readline()
            if (line == '') or (line[:2] != '##'):
                break
            i += line.startswith('##ColumnVariables')
            if line.startswith('##ColumnVariables[tpm'):
                if selcols:
                    for col in selcols:
                        if '.%s.%s' % (col, geneAssembly) in line:
                            tpmIndex.append(i)
                            tissues.append(col)
                            break
                else:
                    tpmIndex.append(i)
                    tissues.append(line.strip()[22:].split(']', maxsplit=1)[0])
            elif line.startswith('##ColumnVariables[hgnc_id]'):
                hgncIndex = i
    print('%d columns selected: %s' % (len(tissues), tissues))
    
    print('Processing expression values')
    with io.open(inPath, "r") as f:
        for n, line in enumerate(f):
            pass
    n += 1
    with io.open(inPath, "r", errors = 'replace') as f:
        expr = {}
        for i, line in enumerate(f):
            sys.stdout.write('  Line %d out of %d (%.2f%%) \r' % (i+1, n, 100*(i+1)/n))
            sys.stdout.flush()
            if not line.startswith('##'):
                linesplit = list(map(str.strip, line.strip().split('\t')))
                if len(linesplit) > 7:
                    hgncID = linesplit[hgncIndex]
                    if hgncID in uniprotID:
                        id = uniprotID[hgncID]
                        if id in uniprotIDlist:
                            tpms = [tpm for ind, tpm in enumerate(linesplit) if ind in tpmIndex]
                            tpms =  list(map(lambda x: float(x) if is_numeric(x) else np.nan, tpms))
                            if id in expr:
                                expr[id].append(tpms)
                            else:
                                expr[id] = [tpms]
    print()
    for k, v in expr.items():
        values = np.nanmean(v, axis=0)
        if normalize:
            expr[k] = (values - values.mean()) / values.std()
        else:
            expr[k] = values
    expr['Tissues'] = tissues
    
    with open(outPath, 'wb') as fOut:
        pickle.dump(expr, fOut)

def produce_geo_expr_dict (inPath,
                           uniprotIDmapFile,
                           gdsDir,
                           outPath,
                           numPoints = 1,
                           avg = 'none'):
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    gdsData = pd.read_table (inPath, sep='\t')
    gdsData = gdsData [(gdsData["Time_point_count"] >= numPoints) &
                       (gdsData["Time_measure"] != '-')].reset_index(drop=True)
    
    expr, n = {}, len(gdsData)
    for i, (id, tm) in enumerate(gdsData[["GDS_ID", "Time_measure"]].values):
        print('*****************************************************')
        print('Dataset %d out of %d' % (i+1, n))
        print('*****************************************************')
        
        gds = read_gds_softfile (id, inDir = gdsDir)
        if gds:
            sampleIDs = []
            if avg == 'none':
                for subset in gds.subsets.values():
                    if tm in subset.metadata['type']:
                        sampleIDs.extend(subset.metadata['sample_id'][0].split(','))
                sampleIDs = list(set(sampleIDs))
                exprTable = gds.table[["IDENTIFIER"] + sampleIDs]
            elif avg == 'all':
                k = 0
                exprTable = pd.DataFrame(data = {"IDENTIFIER":gds.table["IDENTIFIER"].values})
                for subset in gds.subsets.values():
                    if tm in subset.metadata['type']:
                        k += 1
                        meanSampleID = 'Mean_' + str(k)
                        sampleIDs.append(meanSampleID)
                        samples = subset.metadata['sample_id'][0].split(',')
                        exprTable[meanSampleID] = gds.table[samples].mean(axis=1).values
        
            genes = list(set(exprTable["IDENTIFIER"].values))
            for gene in genes:
                e = exprTable.loc[exprTable["IDENTIFIER"] == gene, sampleIDs]
                g = gene.upper()
                if g in uniprotID:
                    p = uniprotID[g]
                    if p not in expr:
                        expr[p] = {}
                    expr[p][id] = e.mean().values if e.shape[0] > 1 else e.iloc[0].values
    
    with open(outPath, 'wb') as fOut:
        pickle.dump(expr, fOut)

def is_transient (p1,
                  p2,
                  expr,
                  minPts = 3,
                  maxCoexpr = 0.05,
                  singleExp = True):
    
    if singleExp:
        c = coexpr (p1, p2, expr, minPts = minPts)
        if not np.isnan(c):
            return 'transient' if c < maxCoexpr else 'permanent'
        else:
            return '-'
    else:
        all = []
        if (p1 in expr) and (p2 in expr):
            for dataset in expr[p1]:
                if dataset in expr[p2]:
                    subExpr = {p1:expr[p1][dataset], p2:expr[p2][dataset]}
                    all.append(is_transient (p1,
                                             p2,
                                             subExpr,
                                             minPts = minPts,
                                             maxCoexpr = maxCoexpr,
                                             singleExp = True))
        numTran = all.count('transient')
        numPerm = all.count('permanent')
        if numTran > numPerm:
            return 'transient'
        elif (numPerm > 0) and (numTran <= numPerm):
            return 'permanent'
        else:
            return '-'

def is_weak (p1, p2, energy, minEnergy = -20):
    
    ppi = tuple(sorted([p1, p2]))
    if ppi in energy:
        return 'weak' if energy[ppi] >= minEnergy else 'strong'
    else:
        return '-'

def is_unbalanced (p1,
                   p2,
                   expr,
                   minPts = 1,
                   logBase = 10,
                   minDiff = 1,
                   singleExp = True):

    if singleExp:
        meanDiff = expr_log_diff (p1,
                                  p2,
                                  expr,
                                  minPts = minPts,
                                  logBase = logBase,
                                  method = 'mean')
    
        if np.isnan(meanDiff):
            return '-'
        else:
            return 'unbalanced' if meanDiff >= minDiff else 'balanced'
    else:
        meanDiff = expr_log_diff (p1,
                                  p2,
                                  expr,
                                  minPts = 1,
                                  logBase = 10,
                                  method = 'mean',
                                  singleExp = False,
                                  average = False)
        
        all = []
        for d in meanDiff:
            if not np.isnan(d):
                all.append('unbalanced' if d >= minDiff else 'balanced')
        
        unbalanced = all.count('unbalanced')
        balanced = all.count('balanced')
        if unbalanced > balanced:
            return 'unbalanced'
        elif balanced > 0:
            return 'balanced'
        else:
            return '-'
