import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_profile_Brinkmann(filename, name):
    distance = 5500*3E18
    minutetoradian = (360.0/(24.0*60.0))*np.pi/(180)
    shift = -14.2
    xmmdata = np.loadtxt("../examples_data/W50/profileBrinkmann.dat")


    Nxmm = xmmdata.shape[0]
    xmmprofile = np.zeros([2,Nxmm])
    maxBrinkmann = np.amax(xmmdata[:,1])

    for i in range(Nxmm):
        xmmprofile[0,i] = (-xmmdata[i,0] - shift)*minutetoradian*distance
        xmmprofile[1,i] = xmmdata[i,1]/maxBrinkmann

    modelxmm = np.loadtxt("../output/xmmprofile_Brinkmann.dat")
    modelxmm = modelxmm.T
    maxflux = np.amax(modelxmm[1,:])
    N1 = modelxmm.shape[1]
    for i in range(N1):
        modelxmm[1,i] = modelxmm[1,i]/maxflux

    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": 'Times New Roman'
    })

    f1 = plt.figure(figsize=[5, 4])
    ax = f1.add_subplot(111)

    ax.set_ylabel(r'$F/F_{max}$', fontsize=20,fontweight='bold')

    ax.set_xlabel(r'z [pc]', fontsize=20, fontweight='bold')

    #ax1.set_xscale("log")
    #ax.set_yscale("log")
    #ax1.set_xlim([1E12, 1E16])
    #ax1.set_xlim([1E11, 1E15])
    #ax1.set_xlim([1E3, 5E4])
    #ax1.set_ylim([6E-15, 5E-10])
    ax.tick_params(axis='x', size=5, width=1)
    ax.tick_params(axis='y', size=5, width=1)
    ax.minorticks_on()

    #plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    #plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4])
    #plt.yticks([6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])

   # ax1.plot(xmmprofile[0], xmmprofile[1], 'r', linewidth=4)
    #plt.plot(lhaaso[0], lhaaso[1],'b', linewidth=4)
    ax.plot(modelxmm[0]/3E18, modelxmm[1], 'b', linewidth = 2, label = 'model 0.3-10 keV')
    ax.plot(xmmprofile[0]/3E18, xmmprofile[1], 'r', linewidth=2, label = 'XMM')

    ax.legend(fontsize = "10")

    #ax.set_xlim([-5E19, 1E19])

    #ax3.set_yscale("log")
    #ax3.set_ylim([9E-7, 6E-5])


    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax1.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')