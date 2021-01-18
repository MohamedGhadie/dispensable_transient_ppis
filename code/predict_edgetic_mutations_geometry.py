#----------------------------------------------------------------------------------------
# Predict PPI perturbations based on geometry.
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from text_tools import write_list_table
from interactome_tools import read_interface_annotated_interactome
from mutation_processing_tools import remove_mutation_overlaps
from mutation_interface_edgotype import (mutation_PPI_interface_perturbations, 
                                         assign_edgotypes,
                                         create_perturbed_network)
from stat_tools import sderror_on_fraction, fisher_test
from plot_tools import network_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # predict mono-edgetic mutations instead of edgetic mutations
    mono_edgetic = False
    
    # plot perturbed interactome and produce files for use by Cytoscape
    plot_perturbations = False
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of calculation method
    methodDir = interactomeDir / 'geometry'
    
    # directory of network perturbation output data files for use by Cytoscape
    cytoscapeDir = methodDir / 'cytoscape'
        
    # figure directory
    figDir = Path('../figures') / interactome_name / 'geometry'
    
    # input data files
    naturalMutationsFile = procDir / 'dbsnp_mutations4.txt'
    diseaseMutationsFile = procDir / 'clinvar_mutations6.txt'
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    
    # output data files
    natMutEdgotypeFile = methodDir / 'nondisease_mutation_edgetics.txt'
    disMutEdgotypeFile = methodDir / 'disease_mutation_edgetics.txt'
    naturalMutEdgeFile = cytoscapeDir / 'nondiseaseMut_perturbed_edges'
    naturalMutNodeFile = cytoscapeDir / 'nondiseaseMut_node_colors'
    diseaseMutEdgeFile = cytoscapeDir / 'diseaseMut_perturbed_edges'
    diseaseMutNodeFile = cytoscapeDir / 'diseaseMut_node_colors'
    
    # create output directories if not existing
    if not methodDir.exists():
        os.makedirs(methodDir)
    if not cytoscapeDir.exists():
        os.makedirs(cytoscapeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # further process mutations
    #------------------------------------------------------------------------------------
    
    naturalMutations, diseaseMutations = remove_mutation_overlaps (naturalMutationsFile,
                                                                   diseaseMutationsFile)
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations
    #------------------------------------------------------------------------------------
    
    # Consider PPI perturbations only for PPIs with this maximum number of interfaces 
    # mapped from distinct models.
    maxInterfaces = np.inf
    
    # Predict PPI perturbation if mutation is this number of residues away in sequence 
    # from an interface residue. Set to 0 if mutation must be exactly at interface.
    numResFromInterface = 0
    
    # Consider PPI perturbations only for PPIs with this minimum number of partners
    minPartners = 1
    
    structuralInteractome = read_interface_annotated_interactome (structuralInteractomeFile)
    
    print( '\n' + 'Predicting PPI perturbations by non-disease mutations based on geometry' )
    naturalMutation_perturbs = mutation_PPI_interface_perturbations (naturalMutations,
                                                                     structuralInteractome,
                                                                     maxInterfaces = maxInterfaces,
                                                                     dist = numResFromInterface)
    print( '\n' + 'Predicting PPI perturbations by disease mutations based on geometry' )
    diseaseMutation_perturbs = mutation_PPI_interface_perturbations (diseaseMutations,
                                                                     structuralInteractome,
                                                                     maxInterfaces,
                                                                     dist = numResFromInterface)
    
    naturalMutations["partners"], naturalMutations["perturbations"] = zip(* naturalMutation_perturbs)
    naturalMutations = naturalMutations[naturalMutations["partners"].apply(len) >= minPartners]
    
    diseaseMutations["partners"], diseaseMutations["perturbations"] = zip(* diseaseMutation_perturbs)
    diseaseMutations = diseaseMutations[diseaseMutations["partners"].apply(len) >= minPartners]
    
    # drop duplicate mutations based on location, regardless of residue type 
    naturalMutations = naturalMutations.drop_duplicates(subset=["protein", "mut_position"]).reset_index(drop=True)
    diseaseMutations = diseaseMutations.drop_duplicates(subset=["protein", "mut_position"]).reset_index(drop=True)
    
    print( '\n' + 'Number of mutations with PPI perturbation predictions after removing duplicates by position' )
    print( 'non-disease: %d' % len(naturalMutations) )
    print( 'disease: %d' % len(diseaseMutations) )
    
    #------------------------------------------------------------------------------------
    # Assign mutation edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations["edgotype"] = assign_edgotypes (naturalMutations["perturbations"].tolist(),
                                                     mono_edgetic = False)
    diseaseMutations["edgotype"] = assign_edgotypes (diseaseMutations["perturbations"].tolist(),
                                                     mono_edgetic = False)
    
    nat_mono_edgotype = assign_edgotypes (naturalMutations["perturbations"].tolist(), mono_edgetic = True)
    dis_mono_edgotype = assign_edgotypes (diseaseMutations["perturbations"].tolist(), mono_edgetic = True)
    
    if mono_edgetic:
        print( '\n' + 'Labeling mono-edgetic mutations' )
        naturalMutations["mono-edgotype"] = nat_mono_edgotype
        diseaseMutations["mono-edgotype"] = dis_mono_edgotype
    
    naturalMutations = naturalMutations [naturalMutations["edgotype"] != '-'].reset_index(drop=True)
    diseaseMutations = diseaseMutations [diseaseMutations["edgotype"] != '-'].reset_index(drop=True)
    
    if mono_edgetic:
        numNaturalMut_edgetic = sum(naturalMutations["mono-edgotype"] == 'mono-edgetic')
        numNaturalMut_nonedgetic = sum(naturalMutations["mono-edgotype"].apply(lambda x: 
                                                        x in ('non-edgetic', 'edgetic')))
        numDiseaseMut_edgetic = sum(diseaseMutations["mono-edgotype"] == 'mono-edgetic')
        numDiseaseMut_nonedgetic = sum(diseaseMutations["mono-edgotype"].apply(lambda x: 
                                                        x in ('non-edgetic', 'edgetic')))
    else:
        numNaturalMut_edgetic = sum(naturalMutations["edgotype"] == 'edgetic')
        numNaturalMut_nonedgetic = sum(naturalMutations["edgotype"] == 'non-edgetic')
        numDiseaseMut_edgetic = sum(diseaseMutations["edgotype"] == 'edgetic')
        numDiseaseMut_nonedgetic = sum(diseaseMutations["edgotype"] == 'non-edgetic')
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    label = 'monoedgetic' if mono_edgetic else 'edgetic'
    print( '\n' + 'Fraction of predicted %s mutations:' % label )
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_edgetic / numNaturalMut_considered,
               sderror_on_fraction (numNaturalMut_edgetic, numNaturalMut_considered),
               numNaturalMut_edgetic,
               numNaturalMut_considered))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_edgetic / numDiseaseMut_considered,
               sderror_on_fraction (numDiseaseMut_edgetic, numDiseaseMut_considered),
               numDiseaseMut_edgetic,
               numDiseaseMut_considered))
    
    fisher_test ([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                 [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    # write predicted mutation edgotypes to file
    write_list_table (naturalMutations, ["partners", "perturbations"], natMutEdgotypeFile)
    write_list_table (diseaseMutations, ["partners", "perturbations"], disMutEdgotypeFile)
        
    #------------------------------------------------------------------------------------
    # plot network perturbations
    #------------------------------------------------------------------------------------
    
    if plot_perturbations:
        print( '\n' + 'Creating network perturbed by non-disease mutations' )
        nodes, edges, nodeColors, edgeColors = create_perturbed_network (structuralInteractome,
                                                                         naturalMutations,
                                                                         naturalMutEdgeFile,
                                                                         naturalMutNodeFile)
        network_plot (edges,
                      nodes = nodes,
                      nodeSizes = [20] * len(nodes),
                      edgeWidth = 1,
                      nodeColors = nodeColors,
                      edgeColors = edgeColors,
                      show = showFigs,
                      figdir = figDir,
                      figname = 'nondisease_mut_perturbed_interactome_geometry')
    
        print( '\n' + 'Creating network perturbed by disease mutations' )
        nodes, edges, nodeColors, edgeColors = create_perturbed_network (structuralInteractome,
                                                                         diseaseMutations,
                                                                         diseaseMutEdgeFile,
                                                                         diseaseMutNodeFile)
        network_plot (edges,
                      nodes = nodes,
                      nodeSizes = [20] * len(nodes),
                      edgeWidth = 1,
                      nodeColors = nodeColors,
                      edgeColors = edgeColors,
                      show = showFigs,
                      figdir = figDir,
                      figname = 'disease_mut_perturbed_interactome_geometry')

if __name__ == "__main__":
    main()
