#----------------------------------------------------------------------------------------
# Write PPI properties to file
#----------------------------------------------------------------------------------------

import os
import pickle
from pathlib import Path
from energy_tools import read_ppi_energy
from interactome_tools import (read_single_interface_annotated_interactome,
                               mutExcl_simult_partners,
                               max_num_mutExcl)
from protein_function import (produce_illumina_expr_dict,
                              produce_fantom5_expr_dict,
                              produce_geo_expr_dict,
                              is_weak,
                              is_transient,
                              is_unbalanced)

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # structure to use for PPI energy
    # options: template, model
    energyStructure = 'template'
    
    # use the median of PPI binding free energy as a cutoff instead of fixed cutoff
    medEnergy = False
    
    # median interaction free energy for all PPIs in structural interactome
    medianEnergy = {'HuRI':     {'template': -20, 'model': -14},
                    'IntAct':   {'template': -15, 'model': -9}}
    
    # minumum interaction free energy required for weak PPIs
    if medEnergy:
        minEnergy = medianEnergy[interactome_name][energyStructure]
    else:
        minEnergy = -25
    
    # minimum number of expression point values required for protein pair tissue
    # co-expression to be considered
    minCoexprPoints = 5
    
    # median of PPI co-expression distribution in structural interactome
    medianCoexpr = {'HuRI':    {'Illumina':0.39, 'GTEx':0.80, 'HPA':0.22, 'Fantom5':0.16, 'GEO':0.10},
                    'IntAct':  {'Illumina':0.45, 'GTEx':0.83, 'HPA':0.26, 'Fantom5':0.22, 'GEO':0.11}}
    
    # minimum number of expression point values required for protein pair 
    # expression ratios to be considered
    minRatioPoints = 1
    
    # log base used for expression ratio
    logBase = 10
    
    # median of log difference in expression for PPIs in structural interactome
    medianDiff = {'HuRI':    {'Illumina':0.63, 'Fantom5':0.68, 'GEO':0.38},
                  'IntAct':  {'Illumina':0.56, 'Fantom5':0.59, 'GEO':0.40}}
    
    # max fraction of interface overlap allowed for simultaneous PPIs
    maxSimultOverlap = 0.1
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513-query-results.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    geoDir = extDir / 'GEO' / 'datasets'
    gdsTypeFile = procDir / 'gds_subset_type.txt'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    interactomeFile = interactomeDir / 'structural_interactome.txt'
    energyFile = interactomeDir / ('ppi_%s_energy_foldx.txt' % energyStructure)
    
    # output data files
    illuminaDictFile = procDir / 'protein_expr_Illumina.pkl'
    fantomDictFile = procDir / 'protein_expr_Fantom5.pkl'
    geoDictFile = procDir / 'protein_expr_GEO.pkl'
    outFile = interactomeDir / 'ppi_properties.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    #------------------------------------------------------------------------------------
    # Produce expression dictionaries
    #------------------------------------------------------------------------------------
    
    if not illuminaDictFile.is_file():
        produce_illumina_expr_dict (illuminaExprFile,
                                    uniprotIDmapFile,
                                    illuminaDictFile)
    
    if not fantomDictFile.is_file():
        produce_fantom5_expr_dict (fantomExprFile,
                                   uniprotIDmapFile,
                                   fantomDictFile,
                                   sampleTypes = 'tissues',
                                   sampleTypeFile = fantomSampleTypeFile,
                                   uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    if not geoDictFile.is_file():
        produce_geo_expr_dict (gdsTypeFile,
                               uniprotIDmapFile,
                               geoDir,
                               geoDictFile,
                               numPoints = 5,
                               avg = 'all')
    
    #------------------------------------------------------------------------------------
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    interactome = read_single_interface_annotated_interactome (interactomeFile)
    
    print()
    print('Interactome: %s' % interactome_name)
    
    #------------------------------------------------------------------------------------
    # Label weak and strong PPIs
    #------------------------------------------------------------------------------------
    
    energy = read_ppi_energy (energyFile)
    interactome ["Strength"] = interactome.apply (lambda x: is_weak (x["Protein_1"],
                                                                     x["Protein_2"],
                                                                     energy,
                                                                     minEnergy = minEnergy), axis=1)
    
    print()
    print('Weak PPIs: %d' % sum(interactome ["Strength"] == 'weak'))
    print('Strong PPIs: %d' % sum(interactome ["Strength"] == 'strong'))
    
    #------------------------------------------------------------------------------------
    # Label transient in time PPIs using GEO data
    #------------------------------------------------------------------------------------
    
    with open(geoDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Transient_in_time"] = interactome.apply (lambda x: 
                                        is_transient (x["Protein_1"],
                                                      x["Protein_2"],
                                                      expr,
                                                      minPts = minCoexprPoints,
                                                      maxCoexpr = medianCoexpr[interactome_name]['GEO'],
                                                      singleExp = False), axis=1)
    
    print()
    print('Transient in time PPIs: %d' % sum(interactome ["Transient_in_time"] == 'transient'))
    print('Permanent in time PPIs: %d' % sum(interactome ["Transient_in_time"] == 'permanent'))
    
    #------------------------------------------------------------------------------------
    # Label transient in space PPIs using Illumina data
    #------------------------------------------------------------------------------------
    
    with open(illuminaDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Transient_in_space_(Illumina)"] = interactome.apply (lambda x: 
                                        is_transient (x["Protein_1"],
                                                      x["Protein_2"],
                                                      expr,
                                                      minPts = minCoexprPoints,
                                                      maxCoexpr = medianCoexpr[interactome_name]['Illumina'],
                                                      singleExp = True), axis=1)
    
    print()
    print('Transient in space PPIs (Illumina): %d' % sum(interactome ["Transient_in_space_(Illumina)"] == 'transient'))
    print('Permanent in spcae PPIs (Illumina): %d' % sum(interactome ["Transient_in_space_(Illumina)"] == 'permanent'))
    
    #------------------------------------------------------------------------------------
    # Label transient in space PPIs using Fantom5 data
    #------------------------------------------------------------------------------------
    
    with open(fantomDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Transient_in_space_(Fantom5)"] = interactome.apply (lambda x: 
                                        is_transient (x["Protein_1"],
                                                      x["Protein_2"],
                                                      expr,
                                                      minPts = minCoexprPoints,
                                                      maxCoexpr = medianCoexpr[interactome_name]['Fantom5'],
                                                      singleExp = True), axis=1)
    
    print()
    print('Transient in space PPIs (Fantom5): %d' % sum(interactome ["Transient_in_space_(Fantom5)"] == 'transient'))
    print('Permanent in spcae PPIs (Fantom5): %d' % sum(interactome ["Transient_in_space_(Fantom5)"] == 'permanent'))
    
    #------------------------------------------------------------------------------------
    # Label unbalanced over time PPIs using GEO data
    #------------------------------------------------------------------------------------
    
    with open(geoDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Balance_over_time"] = interactome.apply (lambda x: 
                                        is_unbalanced (x["Protein_1"],
                                                       x["Protein_2"],
                                                       expr,
                                                       minPts = minRatioPoints,
                                                       logBase = logBase,
                                                       minDiff = medianDiff[interactome_name]['GEO'],
                                                       singleExp = False), axis=1)
    
    print()
    print('Unbalanced over time PPIs: %d' % sum(interactome ["Balance_over_time"] == 'unbalanced'))
    print('Balanced over time PPIs: %d' % sum(interactome ["Balance_over_time"] == 'balanced'))
    
    #------------------------------------------------------------------------------------
    # Label unbalanced over space PPIs using Illumina data
    #------------------------------------------------------------------------------------
    
    with open(illuminaDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Balance_over_space_(Illumina)"] = interactome.apply (lambda x: 
                                        is_unbalanced (x["Protein_1"],
                                                       x["Protein_2"],
                                                       expr,
                                                       minPts = minRatioPoints,
                                                       logBase = logBase,
                                                       minDiff = medianDiff[interactome_name]['Illumina'],
                                                       singleExp = True), axis=1)
    
    print()
    print('Unbalanced over space PPIs (Illumina): %d' % sum(interactome ["Balance_over_space_(Illumina)"] == 'unbalanced'))
    print('Balanced over space PPIs (Illumina): %d' % sum(interactome ["Balance_over_space_(Illumina)"] == 'balanced'))
    
    #------------------------------------------------------------------------------------
    # Label unbalanced over space PPIs using Fantom5 data
    #------------------------------------------------------------------------------------
    
    with open(fantomDictFile, 'rb') as f:
        expr = pickle.load(f)
    
    interactome ["Balance_over_space_(Fantom5)"] = interactome.apply (lambda x: 
                                        is_unbalanced (x["Protein_1"],
                                                       x["Protein_2"],
                                                       expr,
                                                       minPts = minRatioPoints,
                                                       logBase = logBase,
                                                       minDiff = medianDiff[interactome_name]['Fantom5'],
                                                       singleExp = True), axis=1)
    
    print()
    print('Unbalanced over space PPIs (Fantom5): %d' % sum(interactome ["Balance_over_space_(Fantom5)"] == 'unbalanced'))
    print('Balanced over space PPIs (Fantom5): %d' % sum(interactome ["Balance_over_space_(Fantom5)"] == 'balanced'))
    print()
    
    #------------------------------------------------------------------------------------
    # Count the maximum number of mutually exclusive PPIs for each PPI
    #------------------------------------------------------------------------------------
    
    mutExclusive, simultaneous = mutExcl_simult_partners (interactome, maxSimultOverlap = maxSimultOverlap)
    
    interactome["Mutually_exclusive_PPIs"] = interactome.apply (lambda x: 
                                                                max_num_mutExcl (x["Protein_1"],
                                                                                 x["Protein_2"],
                                                                                 mutExclusive), axis=1)
    
    print()
    print('0 mutually exclusive PPIs: %d' % sum(interactome ["Mutually_exclusive_PPIs"] == 0))
    print('1-4 mutually exclusive PPIs: %d' % sum((interactome ["Mutually_exclusive_PPIs"] > 0) & 
                                                  (interactome ["Mutually_exclusive_PPIs"] < 5)))
    print('≥5 mutually exclusive PPIs: %d' % sum(interactome ["Mutually_exclusive_PPIs"] >= 5))
    
    #------------------------------------------------------------------------------------
    # Write output file
    #------------------------------------------------------------------------------------
    
    interactome.to_csv (outFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
