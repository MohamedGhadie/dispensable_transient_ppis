
import numpy as np
import pandas as pd
from pathlib import Path
from plot_tools import bar_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, HuRI, IntAct
    interactome_name = 'HuRI'
        
    # parent directory of all data files
    dataDir = Path('../data')
        
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    ppiFile = interactomeDir / 'ppi_properties.txt'
    
    #------------------------------------------------------------------------------------
    # Load structural interactome
    #------------------------------------------------------------------------------------
    
    ppis = pd.read_table (ppiFile, sep='\t')
    
    print()
    print('Interactome: %s' % interactome_name)
    print('PPIs: %d' % len(ppis))
    
    #------------------------------------------------------------------------------------
    # Label weak and strong PPIs
    #------------------------------------------------------------------------------------
    
    transient_labels = {'weak', 'transient', 'unbalanced'}
    permanent_labels = {'strong', 'permanent', 'balanced', 0}
    
    numTrans = []
    for _, ppi in ppis.iterrows():
        labels = ppi[["Strength",
                      "Transient_in_time",
                      "Transient_in_space_(Illumina)",
                      "Transient_in_space_(Fantom5)",
                      "Balance_over_time",
                      "Balance_over_space_(Illumina)",
                      "Balance_over_space_(Fantom5)",
                      "Mutually_exclusive_PPIs"]].values
        s = sum([p in transient_labels for p in labels if p != '-']) + (int(labels[-1]) > 0)
        numTrans.append(s)
            
    pdf = [numTrans.count(i) for i in range(9)]
    total = sum(pdf)
    
    print('Average number of times classified as transient = %f' % np.mean(numTrans))
    print('PPIs classified as tranient in 5 or more datasets = %f%% (%d out %d)' 
                % (100 * sum(pdf[5:]) / total, sum(pdf[5:]), total))
    
    bar_plot (pdf,
              xlabels = range(9),
              ylabels = None,
              xlabel = 'Number of times classified as transient',
              ylabel = 'Number of PPIs',
              barwidth = 0.5,
              colors = 'magenta',
              ewidth = 1,
              edgecolor = 'black',
              fontsize = 16,
              opacity = None,
              xlim = None,
              ylim = [0, 450],
              xticks = None,
              yMinorTicks = False,
              bottom = 0,
              adjustBottom = False,
              shiftBottomAxis = None,
              xbounds = None,
              show = True,
              figdir = figDir,
              figname = 'ppi_classification_distribution')

if __name__ == "__main__":
    main()
