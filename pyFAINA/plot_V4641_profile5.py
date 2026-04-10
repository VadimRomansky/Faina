import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_V4641_profile5(filename1, filename2, filename3, filename4, filename5, name):
    distance = (20200/3.26)*3E18
    secondtoradian = np.pi/(180*3600)



    xmmwidth = 10*secondtoradian*distance
    nustarwidth = 15*secondtoradian*distance

    modelxmm = np.loadtxt(filename1)
    modelxmm2 = np.loadtxt(filename2)
    modelxmm3 = np.loadtxt(filename3)
    modelxmm4 = np.loadtxt(filename4)
    modelxmm5 = np.loadtxt(filename5)
    modelxmm = modelxmm.T
    modelxmm2 = modelxmm2.T
    modelxmm3 = modelxmm3.T
    modelxmm4 = modelxmm4.T
    modelxmm5 = modelxmm5.T
    N1 = modelxmm.shape[1]

    averagedFluxXmm = 0
    lXmm = 0
    for i in range(N1):
        if((modelxmm[0,i] > -xmmwidth*9.5/10.0) and (modelxmm[0,i] < xmmwidth*0.5/10.0)):
            averagedFluxXmm = averagedFluxXmm + modelxmm[1,i]*(modelxmm[0,i] - modelxmm[0,i-1])
            lXmm = lXmm + (modelxmm[0,i] - modelxmm[0,i-1])
    averagedFluxXmm = averagedFluxXmm/lXmm
    for i in range(N1):
        if ((modelxmm[0, i] > -xmmwidth*9.5/10.0) and (modelxmm[0, i] < xmmwidth*0.5/10.0)):
            modelxmm[1,i] = averagedFluxXmm

    averagedFluxXmm = 0
    lXmm = 0
    for i in range(N1):
        if((modelxmm2[0,i] > -xmmwidth*9.5/10.0) and (modelxmm2[0,i] < xmmwidth*0.5/10.0)):
            averagedFluxXmm = averagedFluxXmm + modelxmm2[1,i]*(modelxmm2[0,i] - modelxmm2[0,i-1])
            lXmm = lXmm + (modelxmm2[0,i] - modelxmm2[0,i-1])
    averagedFluxXmm = averagedFluxXmm/lXmm
    for i in range(N1):
        if ((modelxmm2[0, i] > -xmmwidth*9.5/10.0) and (modelxmm2[0, i] < xmmwidth*0.5/10.0)):
            modelxmm2[1,i] = averagedFluxXmm

    averagedFluxXmm = 0
    lXmm = 0
    for i in range(N1):
        if((modelxmm3[0,i] > -xmmwidth*9.5/10.0) and (modelxmm3[0,i] < xmmwidth*0.5/10.0)):
            averagedFluxXmm = averagedFluxXmm + modelxmm3[1,i]*(modelxmm3[0,i] - modelxmm3[0,i-1])
            lXmm = lXmm + (modelxmm3[0,i] - modelxmm3[0,i-1])
    averagedFluxXmm = averagedFluxXmm/lXmm
    for i in range(N1):
        if ((modelxmm3[0, i] > -xmmwidth*9.5/10.0) and (modelxmm3[0, i] < xmmwidth*0.5/10.0)):
            modelxmm3[1,i] = averagedFluxXmm

    averagedFluxXmm = 0
    lXmm = 0
    for i in range(N1):
        if((modelxmm4[0,i] > -xmmwidth*9.5/10.0) and (modelxmm4[0,i] < xmmwidth*0.5/10.0)):
            averagedFluxXmm = averagedFluxXmm + modelxmm4[1,i]*(modelxmm4[0,i] - modelxmm4[0,i-1])
            lXmm = lXmm + (modelxmm4[0,i] - modelxmm4[0,i-1])
    averagedFluxXmm = averagedFluxXmm/lXmm
    for i in range(N1):
        if ((modelxmm4[0, i] > -xmmwidth*9.5/10.0) and (modelxmm4[0, i] < xmmwidth*0.5/10.0)):
            modelxmm4[1,i] = averagedFluxXmm

    averagedFluxXmm = 0
    lXmm = 0
    for i in range(N1):
        if((modelxmm5[0,i] > -xmmwidth*9.5/10.0) and (modelxmm5[0,i] < xmmwidth*0.5/10.0)):
            averagedFluxXmm = averagedFluxXmm + modelxmm5[1,i]*(modelxmm5[0,i] - modelxmm5[0,i-1])
            lXmm = lXmm + (modelxmm5[0,i] - modelxmm5[0,i-1])
    averagedFluxXmm = averagedFluxXmm/lXmm
    for i in range(N1):
        if ((modelxmm5[0, i] > -xmmwidth*9.5/10.0) and (modelxmm5[0, i] < xmmwidth*0.5/10.0)):
            modelxmm5[1,i] = averagedFluxXmm

    maxflux = np.amax(modelxmm[1, :])
    for i in range(N1):
        modelxmm[1,i] = modelxmm[1,i]/maxflux
        modelxmm2[1,i] = modelxmm2[1,i]/maxflux
        modelxmm3[1,i] = modelxmm3[1,i]/maxflux
        modelxmm4[1,i] = modelxmm4[1,i]/maxflux
        modelxmm5[1,i] = modelxmm5[1,i]/maxflux

    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": 'Times New Roman'
    })
    f1 = plt.figure(figsize=[5, 4])
    ax1 = f1.add_subplot(111)
    plt.subplots_adjust(hspace=.0)
    #ax1.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax1.set_ylabel(r'$F/F_{max}$', fontsize=14,fontweight='bold')
    #ax1.set_xlabel(r'$\nu~Ггц$', fontsize=40,fontweight='bold')
    #ax1.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax1.set_xlabel(r'z [pc]', fontsize=14, fontweight='bold')

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
    #plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    #plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4])
    #plt.yticks([6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])

   # ax1.plot(xmmprofile[0], xmmprofile[1], 'r', linewidth=4)
    #plt.plot(lhaaso[0], lhaaso[1],'b', linewidth=4)
    ax1.plot(modelxmm[0]/3E18, modelxmm[1], 'r', linewidth = 2, label = 'L0 = 1E17')
    ax1.plot(modelxmm2[0]/3E18, modelxmm2[1], 'g', linewidth = 2, label = 'L0 = 3E18')
    ax1.plot(modelxmm3[0]/3E18, modelxmm3[1], 'b', linewidth = 2, label = 'L0 = 5E18')
    ax1.plot(modelxmm4[0]/3E18, modelxmm4[1], 'orange', linewidth = 2, label = 'L0 = 1E18')
    ax1.plot(modelxmm5[0]/3E18, modelxmm5[1], 'salmon', linewidth = 2, label = 'L0 = 2E18')


    ax1.legend(fontsize = "10")


    #ax1.set_xlim([-5E19, 1E19])
    #ax2.set_xlim([-5E19, 1E19])
    #ax3.set_xlim([-5E19, 1E19])
    ax1.set_xlim([-15, 0])
    #ax3.set_ylim([9E-7, 6E-5])


    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax1.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')