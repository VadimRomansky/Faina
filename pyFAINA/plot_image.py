import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_image(filename, name):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    Nr = len(lines)

    Nphi = len(lines[0].split())

    h = 6.626E-27

    radiation = np.zeros([Nr,Nphi])
    for i in range(Nr):
        s = lines[i].split()
        for j in range(Nphi):
            radiation[i][j] = float(s[j])*h


    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    #plt.style.use('dark_background')
    #plt.figure(facecolor="black")

    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)

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

    Nghost = 2
    radiation3 = np.zeros([Nr+Nghost, Nphi])

    for i in range(Nr):
        for j in range(Nphi):
            radiation3[i][j]=radiation[i][j]

    for i in range(Nghost):
        for j in range(Nphi):
            radiation3[Nr+i][j]=0.1*minradiation

    rad = np.linspace(0, 1, Nr+Nghost)
    azm = np.linspace(0, 2 * np.pi, Nphi)
    r, th = np.meshgrid(rad, azm)

    im2 = plt.pcolormesh(th, r, radiation3.T, norm = colors.LogNorm(vmin = minradiation, vmax = amax(radiation)), shading='auto')

    #cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])
    #plt.colorbar(im2, cax=cax2, orientation='horizontal', label = '$cm^{-2}s^{-1}$')
    plt.colorbar(im2, orientation='horizontal', label = '$erg~cm^{-2}s^{-1}Hz^{-1}sr^{-1}$')
    #ax.set_facecolor("black")

    plt.savefig(name + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()