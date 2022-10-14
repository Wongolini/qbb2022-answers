from bdg_loader import *
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

if __name__ == "__main__":
    dir = sys.argv[1]

    #X = location on chr
    #Y = peak
    bdgs = []
    for dirpath,_,filenames in os.walk(dir):
        num_plots = len(filenames)
        fig,axs = plt.subplots(nrows=num_plots,
                                figsize=[5,10],
                                sharex=True,sharey=True)

        for i,f in enumerate(filenames):
            bdg = os.path.abspath(os.path.join(dirpath, f))
            name = f.strip('.bdg')
            d = load_data(bdg)
            axs[i].set_title(name)
            axs[i].plot(d['X'],d['Y'],lw=.1)
            axs[i].fill_between(d['X'],d['Y'],facecolor='blue',alpha=.5)

    #fig.text(0.5, 0.04, 'Chr location', ha='center')
    fig.supylabel('Peaks')
    fig.supxlabel('Chr Position')
    fig.suptitle('Sox2 Overlapping Peaks')
    plt.tight_layout()
    plt.savefig('Sox2_overlapping.png')

    #num_plots = len(args)
    #for i,arg in enumerate(args):


    #print(d1)
