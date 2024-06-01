import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_mask():
    markerFile = open('../marker.dat', 'r')
    lines = markerFile.readlines()
    Nr = 500
    Nz = 1000
    Nphi = 1
    marker = np.zeros([Nr, Nz, Nphi])
    count = 0;
    for i in range(Nr):
        for j in range(Nz):
            for k in range(Nphi):
                line = lines[count]
                count = count + 1
                marker[i][j][k] = float(line.split()[0])

    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)

    im2 = ax.imshow(marker[:,:,0].T, origin='upper', aspect='equal', interpolation= None)

    plt.savefig('marker.png', bbox_inches='tight')
    #plt.show()
    plt.close()