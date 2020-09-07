#----------------------------------------------------------------------------------------
# Plot graph legend.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from plot_tools import multi_bar_plot

def main():
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
        
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
        
    alldisp = [[1,2],[3,4],[5,6],[7,8]]
    alldisp = [[1,2],[3,4],[5,6]]
    multi_bar_plot (alldisp,
                    #colors = ['lightsteelblue', 'mediumseagreen', 'orange', 'red'],
                    colors = ['lightsteelblue', 'mediumseagreen', 'red'],
                    fmt = '.k',
                    msize = 18,
                    ewidth = 1.25,
                    edgecolor = 'black',
                    ecolors = 'k',
                    fontsize = 20,
                    ylim = [0, 10],
                    yMinorTicks = 4,
                    leg = ['a', 'b','c'],
                    show = showFigs,
                    figdir = figDir,
                    figname = 'legend')

if __name__ == "__main__":
    main()
