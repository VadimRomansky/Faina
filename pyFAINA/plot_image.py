import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_image():
    radiationFile = open('../image.dat','r')
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

    im2 = plt.pcolormesh(th, r, radiation.T, norm = colors.LogNorm(vmin = amin(radiation), vmax = amax(radiation)))

    plt.savefig('image.png', bbox_inches='tight')
    #plt.show()
    plt.close()