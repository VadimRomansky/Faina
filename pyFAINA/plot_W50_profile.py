import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_profile(filename, name):
    distance = 5500*3E18
    secondtoradian = np.pi/(180*3600)
    shift = 550
    xmmdata = np.loadtxt("../examples_data/W50/profileXMM.dat")
    nustardata = np.loadtxt("../examples_data/W50/profileNuSTAR.dat")

    Bturb = np.loadtxt("../output/Bturb.dat")
    Bturb = Bturb.T

    Nxmm = xmmdata.shape[0]
    Nnustar = nustardata.shape[0]

    xmmprofile = np.zeros([4,Nxmm])
    nustarprofile = np.zeros([4,Nnustar])

    for i in range(Nxmm):
        xmmprofile[0,i] = (xmmdata[i,0] - shift)*secondtoradian*distance
        xmmprofile[1,i] = xmmdata[i,1]
        xmmprofile[2,i] = xmmdata[i,2] - xmmprofile[1,i]
        xmmprofile[3,i] = xmmprofile[1,i] - xmmdata[i,3]

    for i in range(Nnustar):
        nustarprofile[0,i] = (nustardata[i,0] - shift)*secondtoradian*distance
        nustarprofile[1,i] = nustardata[i,1]
        nustarprofile[2,i] = nustardata[i,2] - nustarprofile[1,i]
        nustarprofile[3,i] = nustarprofile[1,i] - nustardata[i,3]

    modelxmm = np.loadtxt("../output/xmmprofile.dat")
    modelxmm = modelxmm.T
    maxflux = np.amax(modelxmm[1,:])
    N1 = modelxmm.shape[1]
    for i in range(N1):
        modelxmm[1,i] = modelxmm[1,i]/maxflux

    modelnustar = np.loadtxt("../output/nustarprofile.dat")
    modelnustar = modelnustar.T
    maxflux = np.amax(modelnustar[1,:])
    N2 = modelnustar.shape[1]
    for i in range(N2):
        modelnustar[1,i] = 1.1*modelnustar[1,i]/maxflux

    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": 'Times New Roman'
    })
    f1, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col')
    f1.set_figheight(7)
    f1.set_figwidth(5)
    plt.subplots_adjust(hspace=.0)
    #ax1.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax1.set_ylabel(r'$F/F_{max}$', fontsize=14,fontweight='bold')
    ax2.set_ylabel(r'$F/F_{max}$', fontsize=14, fontweight='bold')
    #ax1.set_xlabel(r'$\nu~Ггц$', fontsize=40,fontweight='bold')
    #ax1.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax1.set_xlabel(r'z [pc]', fontsize=14, fontweight='bold')
    ax2.set_xlabel(r'z [pc]', fontsize=14,fontweight='bold')
    ax3.set_xlabel(r'z [pc]', fontsize=14, fontweight='bold')

    ax3.set_ylabel(r'B [$\mu$G]', fontsize=14, fontweight='bold')
    #ax1.set_ylabel(r'$F_{\nu}~мЯн$', fontsize=40,fontweight='bold')
    #ax1.set_yscale("log")
    #ax1.set_xscale("log")
    #ax1.set_xlim([1E12, 1E16])
    #ax1.set_xlim([1E11, 1E15])
    #ax1.set_xlim([1E3, 5E4])
    #ax1.set_ylim([6E-15, 5E-10])
    ax1.tick_params(axis='x', size=5, width=1)
    ax1.tick_params(axis='y', size=5, width=1)
    ax1.minorticks_on()
    ax2.tick_params(axis='x', size=5, width=1)
    ax2.tick_params(axis='y', size=5, width=1)
    ax2.minorticks_on()
    #plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    #plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4])
    #plt.yticks([6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])

   # ax1.plot(xmmprofile[0], xmmprofile[1], 'r', linewidth=4)
    #plt.plot(lhaaso[0], lhaaso[1],'b', linewidth=4)
    ax1.plot(modelxmm[0]/3E18, modelxmm[1], 'b', linewidth = 2, label = 'model 0.3-10 keV')
    ax1.errorbar(xmmprofile[0, :]/3E18, xmmprofile[1, :], yerr = [xmmprofile[3, :], xmmprofile[2, :]], ecolor='r', elinewidth=1.5, linewidth=0, capsize=2.5, capthick=1.5, label = 'XMM')
    ax2.plot(modelnustar[0]/3E18, modelnustar[1],'b',linewidth = 2, label = 'model 10-20 keV')
    ax2.errorbar(nustarprofile[0,:]/3E18, nustarprofile[1,:], yerr = [nustarprofile[3,:], nustarprofile[2,:]], ecolor = 'r', elinewidth=1.5, linewidth=0, capsize=2.5, capthick=1.5, label = 'NuSTAR')

    ax3.plot(Bturb[0,:]/3E18, Bturb[1,:]*1E6,'b', linewidth = 2, label = r'B$_{\|}$')
    ax3.plot(Bturb[0,:]/3E18, Bturb[2,:]*1E6,'r', linewidth = 2, label = r'B$_{\perp}$')

    ax1.legend(fontsize = "10")
    ax2.legend(fontsize="10")
    ax3.legend(fontsize="10")

    #ax1.set_xlim([-5E19, 1E19])
    #ax2.set_xlim([-5E19, 1E19])
    #ax3.set_xlim([-5E19, 1E19])
    ax1.set_xlim([-15, 3])
    ax2.set_xlim([-15, 3])
    ax3.set_xlim([-15, 3])
    ax3.set_yscale("linear")
    #ax3.set_ylim([9E-7, 6E-5])


    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax1.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')