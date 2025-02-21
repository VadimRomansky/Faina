import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_rectangular_image_with_axis(filename, xfilename, yfilename, name, aspect = 'equal'):
    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams['axes.linewidth'] = 0.1
    plt.rcParams['text.usetex'] = True
    f1 = plt.figure(figsize=[10, 8])
    ax = f1.add_subplot(111)
    #plt.style.use('dark_background')
    #plt.figure(facecolor="black")
    radiation = np.loadtxt(filename)
    xdata = np.loadtxt(xfilename)
    ydata = np.loadtxt(yfilename)
    Nx = xdata.shape[0]
    Ny = ydata.shape[0]


    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()

    h = 6.626E-27

    minradiation = amin(radiation)
    maxradiation = amax(radiation)
    if(minradiation == 0):
        radiation2 = np.copy(radiation)
        for i in range(Nx):
            for j in range(Ny):
                if(radiation2[i][j] == 0):
                    radiation2[i][j] = 2*maxradiation
        minradiation = 1E-3*maxradiation

    minradiation = 1E-4 * maxradiation
        #for i in range(Nx):
            #for j in range(Ny):
                #if(radiation[i][j] == 0):
                    #radiation[i][j] = 0.1*minradiation


    im2 = plt.pcolormesh(xdata, ydata, radiation.T, norm = colors.LogNorm(vmin = minradiation, vmax = amax(radiation)), shading='auto')
    #im2 = plt.imshow(radiation.T, origin = 'lower', norm = colors.LogNorm(vmin = minradiation, vmax = amax(radiation)), aspect = aspect, extent = [minX, maxX, minY, maxY])
    #cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])
    #plt.colorbar(im2, cax=cax2, orientation='horizontal', label = '$cm^{-2}s^{-1}$')
    #plt.colorbar(im2, orientation='horizontal', label = '$erg~cm^{-2}s^{-1}sr^{-1}$')
    #ax.set_facecolor("black")
    cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])
    #ax.set_yscale("log")

    plt.colorbar(im2, cax=cax2, orientation='horizontal')  # vertical colorbar for fluid data.
    ax.minorticks_on()

    plt.savefig(name + '.png', bbox_inches='tight')
    #plt.show()
    plt.close()