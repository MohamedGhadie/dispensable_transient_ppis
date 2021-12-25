#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

from pathlib import Path
from text_tools import read_list_table
from stat_tools import fisher_test, sderror_on_fraction

def main():
    
    # reference interactome name
    # options: HuRI, IntAct
    interactome_name = 'HuRI'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'foldx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of edgetic mutation calculation method
    edgeticDir = interactomeDir / 'physics' / (ddg_method + '_edgetics')
    
    # input data files
    naturalMutationsFile = edgeticDir / 'nondisease_mutation_edgetics.txt'
    diseaseMutationsFile = edgeticDir / 'disease_mutation_edgetics.txt'
    
    #------------------------------------------------------------------------------------
    # Load mutation edgotypes
    #------------------------------------------------------------------------------------
    
    naturalMutations = read_list_table (naturalMutationsFile,
                                        ["partners", "perturbations"],
                                        [str, int])
    diseaseMutations = read_list_table (diseaseMutationsFile,
                                        ["partners", "perturbations"],
                                        [str, int])
    
    naturalMutations = naturalMutations [(naturalMutations["edgotype"] == 'edgetic') |
                                         (naturalMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
    diseaseMutations = diseaseMutations [(diseaseMutations["edgotype"] == 'edgetic') |
                                         (diseaseMutations["edgotype"] == 'non-edgetic')].reset_index(drop=True)
        
    #------------------------------------------------------------------------------------
    # Print mutation numbers
    #------------------------------------------------------------------------------------
    
    print('\n' + 'Total number of mutations:')
    print('non-disease mutations: %d' % len(naturalMutations))
    print('disease mutations: %d' % len(diseaseMutations))
    
    natMutationProteins = set(naturalMutations["protein"].tolist())
    disMutationProteins = set(diseaseMutations["protein"].tolist())
    mutationProteins = natMutationProteins | disMutationProteins
    
    print('\n' + 'Total number of proteins carrying mutations:')
    print('non-disease mutations: %d' % len(natMutationProteins))
    print('disease mutations: %d' % len(disMutationProteins))
    print('all mutations: %d' % len(mutationProteins))
    
    #------------------------------------------------------------------------------------
    # Dispensable content among all PPIs
    #------------------------------------------------------------------------------------
    
    numNaturalMut_edgetic = sum(naturalMutations["edgotype"] == 'edgetic')
    numDiseaseMut_edgetic = sum(diseaseMutations["edgotype"] == 'edgetic')
    
    numNaturalMut_monoedgetic = sum((naturalMutations["edgotype"] == 'edgetic') & 
                                    (naturalMutations["perturbations"].apply(
                                     lambda x: sum([i > 0 for i in x])) == 1))

    numDiseaseMut_monoedgetic = sum((diseaseMutations["edgotype"] == 'edgetic') & 
                                    (diseaseMutations["perturbations"].apply(
                                     lambda x: sum([i > 0 for i in x])) == 1))
    
    print('Fraction of predicted monoedgetic mutations among edgetic mutations:')
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numNaturalMut_monoedgetic / numNaturalMut_edgetic,
               sderror_on_fraction (numNaturalMut_monoedgetic, numNaturalMut_edgetic),
               numNaturalMut_monoedgetic,
               numNaturalMut_edgetic))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' 
            % (numDiseaseMut_monoedgetic / numDiseaseMut_edgetic,
               sderror_on_fraction (numDiseaseMut_monoedgetic, numDiseaseMut_edgetic),
               numDiseaseMut_monoedgetic,
               numDiseaseMut_edgetic))
    
    fisher_test ([numNaturalMut_monoedgetic, numNaturalMut_edgetic - numNaturalMut_monoedgetic],
                 [numDiseaseMut_monoedgetic, numDiseaseMut_edgetic - numDiseaseMut_monoedgetic])

if __name__ == "__main__":
    main()
