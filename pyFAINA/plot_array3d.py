import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_array3d():
    radiationFile = open('../area.dat','r')
    lines = radiationFile.readlines()
    N = len(lines)

    Nrho = 50
    Nz = 50
    Nphi = 5


    radiation = np.zeros([Nrho,Nz,Nphi])
    count = 0
    for i in range(Nrho):
        for j in range(Nz):
            for k in range(Nphi):
                s = float(lines[count].split()[0])
                count=count+1
                radiation[i][j][k] = s



    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)

    rad = np.linspace(0, 1, Nrho)
    azm = np.linspace(0, 2*np.pi, Nphi)
    r, th = np.meshgrid(rad, azm)

    ax = plt.subplot(projection="polar")
    ax.axis("off")

    zpoint = int(Nz/2)

    minradiation = amin(radiation[:,zpoint,:])
    maxradiation = amax(radiation[:,zpoint,:])
    if(minradiation == 0):
        radiation2 = np.copy(radiation[:,zpoint,:])
        for i in range(Nrho):
            for j in range(Nphi):
                if(radiation2[i][j] == 0):
                    radiation2[i][j] = 2*maxradiation
        minradiation = amin(radiation2)

        for i in range(Nrho):
            for j in range(Nphi):
                if(radiation[i][zpoint][j] == 0):
                    radiation[i][zpoint][j] = 0.1*minradiation

    im2 = plt.pcolormesh(th, r, radiation[:,zpoint,:].T, norm = colors.LogNorm(vmin = minradiation, vmax = amax(radiation)))

    #cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])
    #plt.colorbar(im2, cax=cax2, orientation='horizontal')

    plt.savefig('array.png', bbox_inches='tight')
    #plt.show()
    plt.close()