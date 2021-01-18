#----------------------------------------------------------------------------------------
# Predict edgetic mutations based on physics.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from text_tools import read_list_table, write_list_table
from interactome_tools import read_single_interface_annotated_interactome
from energy_tools import read_protein_mutation_ddg
from mutation_interface_edgotype import (energy_based_perturbation,
                                         assign_edgotypes,
                                         create_perturbed_network)
from stat_tools import sderror_on_fraction, fisher_test
from plot_tools import network_plot

def main():
    
    # reference interactome name: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # minimum reduction in binding free energy ∆∆G required for PPI disruption
    ddgCutoff = 0.5
    
    # if True predict mono-edgetic mutations instead of edgetic mutations
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
    
    # directory of edgetic mutation calculation method
    edgeticDir = interactomeDir / 'physics' / (ddg_method + '_edgetics')
    
    # directory of network perturbation output data files for use by Cytoscape
    cytoscapeDir = edgeticDir / 'cytoscape'
    
    # figure directory
    figDir = Path('../figures') / interactome_name / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    structuralInteractomeFile = interactomeDir / 'structural_interactome.txt'
    natMutGeomEdgotypeFile = interactomeDir / 'geometry' / 'nondisease_mutation_edgetics.txt'
    disMutGeomEdgotypeFile = interactomeDir / 'geometry' / 'disease_mutation_edgetics.txt'
    natMutDDGFile = interactomeDir / ('nondis_mut_binding_ddg_%s.txt' % ddg_method)
    disMutDDGFile = interactomeDir / ('dis_mut_binding_ddg_%s.txt' % ddg_method)
    
    # output data files
    natMutEdgotypeFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    disMutEdgotypeFile = edgeticDir / 'disease_mutation_edgetics.txt'
    naturalMutEdgeFile = cytoscapeDir / 'nondiseaseMut_perturbed_edges'
    naturalMutNodeFile = cytoscapeDir / 'nondiseaseMut_node_colors'
    diseaseMutEdgeFile = cytoscapeDir / 'diseaseMut_perturbed_edges'
    diseaseMutNodeFile = cytoscapeDir / 'diseaseMut_node_colors'
    
    # create output directories if not existing
    if not edgeticDir.exists():
        os.makedirs(edgeticDir)
    if not cytoscapeDir.exists():
        os.makedirs(cytoscapeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # Read mutation geometry-based edgetic perturbations and mutation ∆∆G
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (natMutGeomEdgotypeFile,
                                        ["partners", "perturbations"],
                                        [str, int])
    diseaseMutations = read_list_table (disMutGeomEdgotypeFile,
                                        ["partners", "perturbations"],
                                        [str, int])
    
    # read change in binding free energy for interfacial mutations
    naturalMutationsDDG = read_protein_mutation_ddg (natMutDDGFile, 'binding')
    diseaseMutationsDDG = read_protein_mutation_ddg (disMutDDGFile, 'binding')
    
    numNatMutGeomEdgetic = sum(naturalMutations["edgotype"] == 'edgetic')
    numDisMutGeomEdgetic = sum(diseaseMutations["edgotype"] == 'edgetic')
    
    print()
    print('Total number of mutations')
    print('Non-disease mutations: %d' % len(naturalMutations))
    print('Disease mutations: %d' % len(diseaseMutations))
    
    print()
    print('Number of known edgetic mutations from geometry calculations:')
    print('Non-disease mutations: %d' % numNatMutGeomEdgetic)
    print('Disease mutations: %d' % numDisMutGeomEdgetic)
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations based on physics
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Performing physics-based edgotype prediction for non-disease mutations' )
    naturalMutations["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (naturalMutations,
                                                                                         naturalMutationsDDG,
                                                                                         ddgCutoff)
    print( '\n' + 'Performing physics-based edgotype prediction for disease mutations' )
    diseaseMutations["perturbations"], knownDDG, unknownDDG = energy_based_perturbation (diseaseMutations,
                                                                                         diseaseMutationsDDG,
                                                                                         ddgCutoff)
    
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
        print('\n' + 'Labeling mono-edgetic mutations')
        naturalMutations["mono-edgotype"] = nat_mono_edgotype
        diseaseMutations["mono-edgotype"] = dis_mono_edgotype
    else:
        if "mono-edgotype" in naturalMutations.columns.values:
            naturalMutations = naturalMutations.drop("mono-edgotype", axis=1)
        if "mono-edgotype" in diseaseMutations.columns.values:
            diseaseMutations = diseaseMutations.drop("mono-edgotype", axis=1)
    
    unknown_nat = naturalMutations["edgotype"] == '-'
    unknown_dis = diseaseMutations["edgotype"] == '-'
    
    print()
    print('Interfacial mutations with known edgetic predictions:')
    print('Non-disease mutations: %d out of %d' % (numNatMutGeomEdgetic - sum(unknown_nat), numNatMutGeomEdgetic))
    print('Disease mutations: %d out of %d' % (numDisMutGeomEdgetic - sum(unknown_dis), numDisMutGeomEdgetic))
    
    naturalMutations = naturalMutations [unknown_nat == False].reset_index(drop=True)
    diseaseMutations = diseaseMutations [unknown_dis == False].reset_index(drop=True)
    
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
    
    # write predicted mutation edgotypes to tab-delimited file
    write_list_table (naturalMutations, ["partners", "perturbations"], natMutEdgotypeFile)
    write_list_table (diseaseMutations, ["partners", "perturbations"], disMutEdgotypeFile)
    
    #------------------------------------------------------------------------------------
    # plot network perturbations
    #------------------------------------------------------------------------------------
    
    if plot_perturbations:
        structuralInteractome = read_single_interface_annotated_interactome (structuralInteractomeFile)
        
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
                      figname = 'nondisease_mut_perturbed_interactome_%s' % ddg_method)
    
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
                      figname = 'disease_mut_perturbed_interactome_%s' % ddg_method)

if __name__ == "__main__":
    main()
