#----------------------------------------------------------------------------------------
# Map mutations that are edgetic by geometry (interfacial) onto PPI structural models 
# to be submitted for binding ∆∆G calculations.
#----------------------------------------------------------------------------------------

from pathlib import Path
from text_tools import read_list_table
from threeD_structure_tools import write_interfacial_mutation_structure_maps

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # allow PDB structure files to be downloaded
    download_pdbs = False
    
    # suppress PDB warnings
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of processed model-related data files specific to interactome
    modellingDir = Path('../../edgotype_fitness_effect/data/processed/') / interactome_name / 'model_based'
    
    # directory of PPI structural models
    modelDir = modellingDir / 'ppi_models'
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    chainSeqFile = modellingDir / 'ppi_chain_sequences.pkl'
    chainStrucResFile = modellingDir / 'ppi_chain_strucRes.pkl'
    chainInterfaceFile = modellingDir / 'model_interfaces.txt'
    chainMapFile = modellingDir / 'struc_interactome_chain_map.txt'
    natMutEdgotypeFile = interactomeDir / 'geometry' / 'nondisease_mutation_edgetics.txt'
    disMutEdgotypeFile = interactomeDir / 'geometry' / 'disease_mutation_edgetics.txt'
    
    # output data files
    natural_mutations_onchains_file = interactomeDir / ('nondis_mut_binding_ddg_%s.txt' % ddg_method)
    disease_mutations_onchains_file = interactomeDir / ('dis_mut_binding_ddg_%s.txt' % ddg_method)
    
    #------------------------------------------------------------------------------------
    # write edgetic mutations mapped onto PPI structural models to file
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (natMutEdgotypeFile,
                                        cols = ["partners", "perturbations"],
                                        dtyp = [str, int])
    diseaseMutations = read_list_table (disMutEdgotypeFile,
                                        cols = ["partners", "perturbations"],
                                        dtyp = [str, int])
    
    naturalMutations = naturalMutations[naturalMutations["edgotype"] == 'edgetic'].reset_index(drop=True)
    diseaseMutations = diseaseMutations[diseaseMutations["edgotype"] == 'edgetic'].reset_index(drop=True)

    # write edgetic non-disease mutations
    if not natural_mutations_onchains_file.is_file():
        print( '\n' + 'Writing edgetic non-disease mutations with structure mapping to file ' 
                + str( natural_mutations_onchains_file ) )
        write_interfacial_mutation_structure_maps (naturalMutations,
                                                   structuralInteractomeFile,
                                                   chainMapFile,
                                                   chainSeqFile,
                                                   chainStrucResFile,
                                                   modelDir,
                                                   natural_mutations_onchains_file,
                                                   chainInterfaceFile = chainInterfaceFile,
                                                   downloadPDB = download_pdbs,
                                                   suppressWarnings = suppress_pdb_warnings)
    
    # write edgetic disease mutations
    if not disease_mutations_onchains_file.is_file():
        print( '\n' + 'Writing edgetic disease mutations with structure mapping to file ' 
                + str( disease_mutations_onchains_file ) )
        write_interfacial_mutation_structure_maps (diseaseMutations,
                                                   structuralInteractomeFile,
                                                   chainMapFile,
                                                   chainSeqFile,
                                                   chainStrucResFile,
                                                   modelDir,
                                                   disease_mutations_onchains_file,
                                                   chainInterfaceFile = chainInterfaceFile,
                                                   downloadPDB = download_pdbs,
                                                   suppressWarnings = suppress_pdb_warnings)

if __name__ == "__main__":
    main()
