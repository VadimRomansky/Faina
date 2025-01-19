import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_rectangular_array3d(filename, name):
    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams['axes.linewidth'] = 0.1
    plt.rcParams['text.usetex'] = True
    f1 = plt.figure(figsize=[10, 8])
    ax = f1.add_subplot(111)

    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    minX = -2.4688e+20
    maxX = 2.4688e+20
    minY = 0
    maxY = 1.2344e+20

    Nrho = 200
    Nz = 50
    Nphi = 25


    radiation = np.zeros([Nrho,Nz,Nphi])
    count = 0
    for i in range(Nrho):
        for j in range(Nz):
            for k in range(Nphi):
                s = float(lines[count].split()[0])
                count=count+1
                radiation[i][j][k] = s


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

    minradiation = 1E-4 * maxradiation

    radiation3 = np.copy(radiation[:,zpoint,:])

    im2 = plt.imshow(radiation3.T, origin = 'lower', norm = colors.Normalize(vmin = minradiation, vmax = amax(radiation)), aspect = 'equal', extent = [minX, maxX, minY, maxY])

    cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])

    plt.colorbar(im2, cax=cax2, orientation='horizontal')  # vertical colorbar for fluid data.
    ax.minorticks_on()

    plt.savefig(name + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()