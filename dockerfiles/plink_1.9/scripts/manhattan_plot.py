#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import seaborn as sns


def manhattan_plot(chromnum, infile, outfile):
    # Data
    df = pd.read_csv(infile)
    df = df.dropna()
    
    # Figure
    fig, ax = plt.subplots(figsize=(15, 10))
    #fig.patch.set_facecolor('#f5f5f5')
    #ax.set(frame_on=False)

    # Colors
    c1 = '#0095ff' # Blue
    c2 = '#ffbf00'

    # Plot
    sns.stripplot(np.array([chromnum]*len(df.BP)), -np.log10(df.P), c=c1, ax=ax)
    #ax.hlines(-np.log10(5e-8), df.BP.min(), df.BP.max(), color=c2)

    # Labels
    ax.set_xlabel('Chromosome', fontsize=28)
    ax.set_ylabel('-log10(P)', fontsize=28)

    # Ticks
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    # Show
    fig.savefig(outfile)
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromnum', help="Chromosome number", required=True)
    parser.add_argument('-i', '--infile', help="Assoc file", required=True)
    parser.add_argument('-o', '--outfile', help="Output file of image", required=True)
    
    args = parser.parse_args()
    manhattan_plot(args.chromnum, args.infile, args.outfile)
