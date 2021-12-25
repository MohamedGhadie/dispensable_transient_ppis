#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

from pathlib import Path
from text_tools import reduce_fasta_headers
from id_mapping import (produce_geneName_dict,
                        produce_uniqueGene_swissProtIDs,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict)

def main():
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # input data files
    uniprotRefSeqFile = extDir / 'UP000000589_10090.fasta'
    idmapFile = extDir / 'MOUSE_10090_idmapping.dat'
    proteomeListFile = extDir / 'uniprot_reviewed_mouse_proteome.list'
    
    # output data files
    refSeqFile = procDir / 'mouse_reference_sequences.fasta'
    geneMapFile = procDir / 'to_mouse_geneName_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_mouse_proteome.list'
    proteinSeqFile = procDir / 'mouse_reference_sequences.pkl'
    uniprotIDmapFile = procDir / 'to_mouse_uniprotID_map.pkl'
    
    # create output directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    
    #------------------------------------------------------------------------------------
    # Load interactome perturbations
    #------------------------------------------------------------------------------------
    
    if not refSeqFile.is_file():
        print('Reducing headers in protein sequence fasta file')
        reduce_fasta_headers (uniprotRefSeqFile, '|', 2, 2, refSeqFile)
    
    if not proteinSeqFile.is_file():
        print('producing protein sequence dictionary')
        produce_proteinSeq_dict (refSeqFile, proteinSeqFile)
    
    if not geneMapFile.is_file():
        print('producing UniProtID-to-geneName dictionary')
        produce_geneName_dict (idmapFile, proteomeListFile, geneMapFile)
    
    if not uniqueGeneSwissProtIDFile.is_file():
        print('producing list of unique-gene UniProt IDs')
        produce_uniqueGene_swissProtIDs (proteomeListFile, geneMapFile, uniqueGeneSwissProtIDFile)
    
    if not uniprotIDmapFile.is_file():
        print('producing to-UniProt-ID dictionary')
        produce_uniprotID_dict (idmapFile, uniqueGeneSwissProtIDFile, uniprotIDmapFile)

if __name__ == "__main__":
    main()
