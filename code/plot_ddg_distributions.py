#----------------------------------------------------------------------------------------
# Plot binding ∆∆G distributions.
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from energy_tools import read_protein_mutation_ddg
from plot_tools import multi_histogram_plot

def main():
    
    # reference interactome names: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
    
    # method used for calculating PPI binding ∆∆G
    # options: foldx, mCSM
    ddg_method = 'foldx'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    natMutFile = interactomeDir / ('nondis_mut_binding_ddg_%s.txt' % ddg_method)
    disMutFile = interactomeDir / ('dis_mut_binding_ddg_%s.txt' % ddg_method)
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load mutation ∆∆G
    #------------------------------------------------------------------------------------
    
    natMutDDG = read_protein_mutation_ddg (natMutFile, type = 'binding')
    disMutDDG = read_protein_mutation_ddg (disMutFile, type = 'binding')
    
    natMutDDG = [ddg[-1] for ddg in natMutDDG.values() if not np.isnan(ddg[-1])]
    disMutDDG = [ddg[-1] for ddg in disMutDDG.values() if not np.isnan(ddg[-1])]
    
    print('Average ∆∆G for nondisease mutations = %f' % np.mean(natMutDDG))
    print('Average ∆∆G for disease mutations = %f' % np.mean(disMutDDG))
    
    multi_histogram_plot (natMutDDG,
                          'g',
                          xlabel = 'Change in PPI binding free energy\n∆∆G (kcal/mol)',
                          ylabel = 'Number of mutations',
                          fontsize = 22,
                          edgecolor = 'k',
                          bins = 40,
                          xlim = [-12, 15],
                          ylim = [0, 64],
                          ylabels = [0, 16, 32, 48, 64],
                          #ylim = [0, 176],
                          #ylabels = [0, 44, 88, 132, 176],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'nondisease_mut_binding_ddg_distribution_%s' % ddg_method)
    
    multi_histogram_plot (disMutDDG,
                          'm',
                          xlabel = 'Change in PPI binding free energy\n∆∆G (kcal/mol)',
                          ylabel = 'Number of mutations',
                          fontsize = 22,
                          edgecolor = 'k',
                          bins = 40,
                          xlim = [-12, 15],
                          ylim = [0, 64],
                          ylabels = [0, 16, 32, 48, 64],
                          #ylim = [0, 124],
                          #ylabels = [0, 31, 62, 93, 124],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'disease_mut_binding_ddg_distribution_%s' % ddg_method)
    
if __name__ == "__main__":
    main()
