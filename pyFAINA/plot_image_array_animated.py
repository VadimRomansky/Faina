import matplotlib.pyplot as plt
from matplotlib import animation, colors
from matplotlib.animation import FuncAnimation
from pylab import *
import numpy as np

def plot_image_array_animated(filename, name):
    f1 = plt.figure(figsize=[8, 8])

    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    Nrho = 50
    Nphi = len(lines[0].split())

    Nnu = 20


    radiation = np.zeros([Nnu,Nrho,Nphi])
    count = 0
    for i in range(Nnu):
        for j in range(Nrho):
            s =  lines[count].split()
            count = count + 1
            for k in range(Nphi):
                radiation[i][j][k] = float(s[k])



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


    minradiation = amin(radiation[:,:,:])
    maxradiation = amax(radiation[:,:,:])

    if (minradiation == 0):
        radiation2 = np.copy(radiation[:, :, :])
        for k in range(Nnu):
            for i in range(Nrho):
                for j in range(Nphi):
                    if(radiation2[k][i][j] == 0):
                        radiation2[k][i][j] = 2*maxradiation

        minradiation = amin(radiation2)

        for i in range(Nrho):
            for j in range(Nphi):
                if(radiation[k][i][j] == 0):
                    radiation[k][i][j] = 0.1*minradiation

    def update(frame_number):
        #f1 = plt.figure(figsize=[6, 6])
        f1.clear()
        f1.set_figheight(8)
        f1.set_figwidth(8)
        ax = plt.subplot(projection="polar")
        im2 = ax.pcolormesh(th, r, radiation[frame_number,:,:].T, norm = colors.Normalize(vmin = minradiation, vmax = amax(radiation)))
        plt.colorbar(im2, orientation='horizontal')
        return im2
    #cax2 = f1.add_axes([0.125, 0.92, 0.775, 0.03])
    #plt.colorbar(im2, cax=cax2, orientation='horizontal')

    anim = FuncAnimation(f1, update, interval=10, frames=Nnu)

    f = name +".gif"
    writergif = animation.PillowWriter(fps=4)
    anim.save(f, writer=writergif)
    plt.close()