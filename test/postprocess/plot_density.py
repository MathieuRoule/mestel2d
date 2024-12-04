import argparse
import warnings
import numpy as np
import h5py
import matplotlib.pyplot as plt
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap

def plot_density(file, nsnap, vmin=0, vmax=15):
    tabxy = file['full/snap'+f"{nsnap:04}"+'/xy']
    x, y = tabxy[::2], tabxy[1::2]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ax.scatter_density(x, y, vmin=vmin, vmax=vmax)
    ax.set_aspect('equal')
    ax.tick_params(axis='both', which='both', bottom=False, left=False, 
        labelbottom=False, labelleft=False)
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    time = round(file['full/tabt'][nsnap],1)
    plt.title("$t =$"+f"{time}"+"$t_{\\rm dyn}$", loc="left", x=0.4)
    return fig

if __name__ == "__main__":
    # Parsing the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vmin', type=int, default=0)
    parser.add_argument('-vmax', type=int, default=1)
    parser.add_argument('filename')
    parser.add_argument('outpath')
    args = parser.parse_args()
    # Setting the arguments
    vmin = args.vmin
    vmax = args.vmax
    filename = args.filename
    outpath = args.outpath

    # Open the file
    file = h5py.File(filename, 'r')
    # Read the number of snapshots 
    snapshots = len(file['full'].keys()) - 1 # -1 because of the field 'tabt' 
    print("Plotting", snapshots, "snapshots.")
    for i in range(snapshots):
        print(f"\r {i+1} / {snapshots}",end="")
        pl = plot_density(file, i, vmin, vmax)
        pl.savefig(outpath+'frame_'+f"{i:03}"+'.png')
        plt.close(pl)
    print("") # newline
    # Close the file
    file.close()