import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_image(filename, name):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    Nr = len(lines)

    Nphi = len(lines[0].split())

    radiation = np.zeros([Nr,Nphi])
    for i in range(Nr):
        s = lines[i].split()
        for j in range(Nphi):
            radiation[i][j] = float(s[j])


    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)

    rad = np.linspace(0, 1, Nr)
    azm = np.linspace(0, 2*np.pi, Nphi)
    r, th = np.meshgrid(rad, azm)

    ax = plt.subplot(projection="polar")
    ax.axis("off")

    minradiation = amin(radiation)
    maxradiation = amax(radiation)
    if(minradiation == 0):
        radiation2 = np.copy(radiation)
        for i in range(Nr):
            for j in range(Nphi):
                if(radiation2[i][j] == 0):
                    radiation2[i][j] = 2*maxradiation
        minradiation = amin(radiation2)

        for i in range(Nr):
            for j in range(Nphi):
                if(radiation[i][j] == 0):
                    radiation[i][j] = 0.1*minradiation

    im2 = plt.pcolormesh(th, r, radiation.T, norm = colors.LogNorm(vmin = minradiation, vmax = amax(radiation)), shading='gouraud')

    cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])
    plt.colorbar(im2, cax=cax2, orientation='horizontal')

    plt.savefig(name + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()