import matplotlib.pyplot as plt
from matplotlib import animation, colors
from pylab import *
import numpy as np

def plot_error_profile(Np1, Np2):

    parameterFile1 = open('../parameter' + str(Np1) + '.dat','r')
    parameterFile2 = open('../parameter' + str(Np2) + '.dat','r')
    errorFile = open('../error_optimized_'+str(Np1)+'_'+str(Np2)+'.dat','r')

    lines1 = parameterFile1.readlines()
    lines2 = parameterFile2.readlines()
    linesE = errorFile.readlines()

    N = len(lines1)

    parameter1 = np.zeros([N])
    parameter2 = np.zeros([N])
    error = np.zeros([N,N])

    for i in range(N):
        parameter1[i] = float(lines1[i])
        parameter2[i] = float(lines2[i])
        s = linesE[i].split()
        for j in range(N):
            error[i][j] = float(s[j])

    minErr = np.amin(error)
    maxErr = np.amax(error)

    plt.rcParams.update({'font.size': 20})
    plt.rcParams['text.usetex'] = True
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    cax2 = f1.add_axes([0.125, 0.97, 0.775, 0.03])
    #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_xlabel(r'$n~cm^{-3}$', fontsize=30)
    ax.set_ylabel(r'$\sigma$', fontsize=30)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.minorticks_on()
    #x,y = np.meshgrid(parameter1, parameter2)
    #x = np.log10(parameter1)
    #y = np.log10(parameter2)
    # plt.axis([0.0,1.0,0.0,1.0])
    #im2 = ax.imshow(error,origin='lower', norm=colors.LogNorm(vmin=5, vmax=1E5), aspect='auto',
    #          extent=[parameter1[0], parameter1[N-1],parameter2[0],parameter2[N-1]], cmap='jet', interpolation='gaussian')

    im2 = ax.pcolormesh(parameter2, parameter1, error, cmap='jet', norm=colors.LogNorm(vmin=5, vmax=1E5))
    plt.colorbar(im2, cax=cax2, orientation='horizontal')
    plt.show()
    #plt.savefig('error.png')