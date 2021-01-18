#----------------------------------------------------------------------------------------
# Plot example pie chart.
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from plot_tools import pie_plot

def main():
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of processed data files
    procDir = dataDir / 'processed'
        
    # figure directory
    figDir = Path('../figures')
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    pie_plot ([65, 23, 12],
              angle = 90,
              colors = ['lightsteelblue', 'red', 'mediumseagreen'],
              edgewidth = 2,
              show = showFigs,
              figdir = figDir,
              figname = 'example_pie_chart')

if __name__ == "__main__":
    main()
