import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_W50_EFE_synchrotron(filename, name, factor = 1.0):
    radiationFile = open(filename,'r')
    lines = radiationFile.readlines()
    N = len(lines)

    radiation = np.zeros([4,N])
    for i in range(N):
        s = lines[i].split()
        radiation[0,i] = float(s[0])
        radiation[1,i] = factor*float(s[1])
        radiation[2,i] = factor*float(s[2])
        radiation[3,i] = factor*float(s[3])

    outFile = open('Wsynchandcompt.dat', 'w')
    for i in range(N):
        print(radiation[0,i], radiation[1,i], radiation[2,i], radiation[3,i], sep = ' ', file = outFile)

    outFile.close()


    plt.rcParams.update({'font.size': 15})
    plt.rcParams["figure.dpi"] = 500
    plt.rcParams['axes.linewidth'] = 1.0
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": 'Times New Roman'
    })

    f1 = plt.figure(figsize=[5, 4])
    ax = f1.add_subplot(111)
    #ax.set_xlabel(r'$\nu~GHz$', fontsize=40,fontweight='bold')
    ax.set_xlabel(r'E [eV]', fontsize=14,fontweight='bold')
    #ax.set_xlabel(r'$\nu~Ггц$', fontsize=40,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~mJy$', fontsize=40,fontweight='bold')
    ax.set_ylabel(r'EF(E) [erg cm$^{-2}$ s$^{-1}$]', fontsize=14,fontweight='bold')
    #ax.set_ylabel(r'$F_{\nu}~мЯн$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #ax.set_xlim([1E12, 1E16])
    #ax.set_xlim([1E11, 1E15])
    ax.set_xlim([1E-7, 5E5])
    ax.set_ylim([1E-18, 0.3E-11])
    ax.tick_params(axis='x', size=5, width=1)
    ax.tick_params(axis='y', size=5, width=1)
    ax.minorticks_on()
    #plt.xticks([1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14])
    #plt.xticks([1E2, 2E2, 3E2, 4E2, 5E2, 6E2, 7E2, 8E2, 9E2, 1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4, 6E4, 7E4, 8E4, 9E4, 1E5, 2E5, 3E5, 4E5, 6E5, 7E5, 8E5, 9E5, 1E6])
    #plt.yticks([6E-11, 7E-11, 8E-11, 9E-11, 1E-10, 2E-10])

    plt.plot(radiation[0], radiation[2], 'r', linewidth=2, label = 'Head model')
    plt.plot(radiation[0], radiation[3], 'b', linewidth=2, label = 'Cone model')

    Ne = 10
    E1 = np.linspace(500, 10000, Ne)
    F1h = np.zeros([Ne])
    F1hl = np.zeros([Ne])
    F1hr = np.zeros([Ne])
    F1c = np.zeros([Ne])
    F1cl = np.zeros([Ne])
    F1cr = np.zeros([Ne])
    F1d = np.zeros([Ne])
    Ah = 1.985E-9
    Ahl = 6.079978E-10
    Ahr = 6.463794443E-9
    Ac = 5.2E-10
    Acl = 1.931678841E-10
    Acr = 1.406750397E-9
    Ad = 5.89E-11
    for i in range(Ne):
        currentE = E1[i]*1.6E-12
        F1h[i] = Ah*currentE**(2-1.58)
        F1hl[i] = Ahl * currentE ** (2 - (1.58+0.06))
        F1hr[i] = Ahr * currentE ** (2 - (1.58-0.06))
        F1c[i] = Ac*currentE**(2-1.65)
        F1cl[i] = Acl * currentE ** (2 - (1.65+0.05))
        F1cr[i] = Acr * currentE ** (2 - (1.65-0.05))
        F1d[i] = Ad*currentE**(2-1.75)

    E2 = np.linspace(3000,30000, Ne)
    F2h = np.zeros([Ne])
    F2hl = np.zeros([Ne])
    F2hr = np.zeros([Ne])
    F2c = np.zeros([Ne])
    F2cl = np.zeros([Ne])
    F2cr = np.zeros([Ne])
    F2d = np.zeros([Ne])
    Bh = 1.125E-9
    Bhl = 1.792936566E-10
    Bhr = 7.009018705E-9
    Bc = 6.86E-13
    Bcl = 1.530579984E-13
    Bcr = 2.570021148E-12
    Bd = 4.34E-13
    for i in range(Ne):
        currentE = E2[i]*1.6E-12
        F2h[i] = Bh*currentE**(2-1.6)
        F2hl[i] = Bhl * currentE ** (2 - (1.6+0.1))
        F2hr[i] = Bhr * currentE ** (2 - (1.6-0.1))
        F2c[i] = Bc*currentE**(2.0 - 2.0)
        F2cl[i] = Bcl * currentE ** (2.0 - (2.0+0.08))
        F2cr[i] = Bcr * currentE ** (2.0 - (2.0-0.07))
        F2d[i] = Bd*currentE**(2.0 - 2.0)

    EminXmm = 500*1.6E-12
    EmaxXmm = 10*1000*1.6E-12
    headIndexXmm = 1.58
    headMaxIndexXmm = 1.58+0.06
    headMinIndexXmm = 1.58-0.06
    headFluxXmm = 1.8E-12
    headMaxFluxXmm = (1.8 + 0.06)*1E-12
    headMinFluxXmm = (1.8 - 0.06)*1E-12
    coneIndexXmm = 1.65
    coneMaxIndexXmm = 1.65 + 0.05
    coneMinIndexXmm = 1.65 - 0.05
    coneFluxXmm = 1.81E-12
    coneMaxFluxXmm = (1.81 + 0.06) * 1E-12
    coneMinFluxXmm = (1.81 - 0.06) * 1E-12

    lowHeadFluxXmm = np.zeros([Ne])
    upperHeadFluxXmm = np.zeros([Ne])
    indexes = np.linspace(headMinIndexXmm, headMaxIndexXmm, 100)
    tempFlux = np.zeros([100])

    for i in range(Ne):
        for j in range(100):
            tempFlux[j] = headMaxFluxXmm*(2-indexes[j])*(E1[i]*1.6E-12)**(2-indexes[j])/(EmaxXmm**(2-indexes[j]) - EminXmm**(2-indexes[j]))
        upperHeadFluxXmm[i] = np.amax(tempFlux)
        for j in range(100):
            tempFlux[j] = headMinFluxXmm*(2-indexes[j])*(E1[i]*1.6E-12)**(2-indexes[j])/(EmaxXmm**(2-indexes[j]) - EminXmm**(2-indexes[j]))
        lowHeadFluxXmm[i] = np.amin(tempFlux)

    lowConeFluxXmm = np.zeros([Ne])
    upperConeFluxXmm = np.zeros([Ne])
    indexes = np.linspace(coneMinIndexXmm, coneMaxIndexXmm, 100)
    tempFlux = np.zeros([100])

    for i in range(Ne):
        for j in range(100):
            tempFlux[j] = coneMaxFluxXmm*(2-indexes[j])*(E1[i]*1.6E-12)**(2-indexes[j])/(EmaxXmm**(2-indexes[j]) - EminXmm**(2-indexes[j]))
        upperConeFluxXmm[i] = np.amax(tempFlux)
        for j in range(100):
            tempFlux[j] = coneMinFluxXmm*(2-indexes[j])*(E1[i]*1.6E-12)**(2-indexes[j])/(EmaxXmm**(2-indexes[j]) - EminXmm**(2-indexes[j]))
        lowConeFluxXmm[i] = np.amin(tempFlux)



    EminNustar = 3*1000*1.6E-12
    EmaxNustar = 30*1000*1.6E-12
    headIndexNustar = 1.6
    headMaxIndexNustar = 1.6+0.1
    headMinIndexNustar = 1.6-0.1
    headFluxNustar = 2.0E-12
    headMaxFluxNustar = (2.0 + 0.1)*1E-12
    headMinFluxNustar = (2.0 - 0.1)*1E-12
    coneIndexNustar = 2.0
    coneMaxIndexNustar = 2.0 + 0.08
    coneMinIndexNustar = 2.0 - 0.07
    coneFluxNustar = 1.58E-12
    coneMaxFluxNustar = (1.58 + 0.1) * 1E-12
    coneMinFluxNustar = (1.58 - 0.09) * 1E-12

    lowHeadFluxNustar = np.zeros([Ne])
    upperHeadFluxNustar = np.zeros([Ne])
    indexes = np.linspace(headMinIndexNustar, headMaxIndexNustar, 100)
    tempFlux = np.zeros([100])

    for i in range(Ne):
        for j in range(100):
            tempFlux[j] = headMaxFluxNustar*(2-indexes[j])*(E2[i]*1.6E-12)**(2-indexes[j])/(EmaxNustar**(2-indexes[j]) - EminNustar**(2-indexes[j]))
        upperHeadFluxNustar[i] = np.amax(tempFlux)
        for j in range(100):
            tempFlux[j] = headMinFluxNustar*(2-indexes[j])*(E2[i]*1.6E-12)**(2-indexes[j])/(EmaxNustar**(2-indexes[j]) - EminNustar**(2-indexes[j]))
        lowHeadFluxNustar[i] = np.amin(tempFlux)

    lowConeFluxNustar = np.zeros([Ne])
    upperConeFluxNustar = np.zeros([Ne])
    indexes = np.linspace(coneMinIndexNustar, coneMaxIndexNustar, 100)
    tempFlux = np.zeros([100])

    for i in range(Ne):
        for j in range(100):
            tempFlux[j] = coneMaxFluxNustar*(2-indexes[j])*(E2[i]*1.6E-12)**(2-indexes[j])/(EmaxNustar**(2-indexes[j]) - EminNustar**(2-indexes[j]))
        upperConeFluxNustar[i] = np.amax(tempFlux)
        for j in range(100):
            tempFlux[j] = coneMinFluxNustar*(2-indexes[j])*(E2[i]*1.6E-12)**(2-indexes[j])/(EmaxNustar**(2-indexes[j]) - EminNustar**(2-indexes[j]))
        lowConeFluxNustar[i] = np.amin(tempFlux)

    Ne1 = Ne-3

    plt.plot(E1, F1h, linewidth = 1, label = 'Head XMM', color = 'm', linestyle = 'dashed')
    plt.fill_between(E1, lowHeadFluxXmm, upperHeadFluxXmm, color = 'm', alpha = 0.2)
    #plt.plot(E1, F1hl, linewidth=1, color = 'm', linestyle = 'dashed')
    #plt.plot(E1, F1hr, linewidth=1, color = 'm', linestyle = 'dashed')
    plt.plot(E2[0:Ne1-1], F2h[0:Ne1-1], linewidth = 1, label = 'Head NuSTAR', color = 'm', linestyle = 'dotted')
    plt.fill_between(E2[0:Ne1-1], lowHeadFluxNustar[0:Ne1-1], upperHeadFluxNustar[0:Ne1-1], color = 'm', alpha = 0.2)
    #plt.plot(E2, F2hl, linewidth=1, color = 'm', linestyle = 'dotted')
    #plt.plot(E2, F2hr, linewidth=1, color = 'm', linestyle = 'dotted')

    plt.plot(E1, F1c, linewidth = 1, label = 'Cone XMM', color = 'c', linestyle = 'dashed')
    plt.fill_between(E1, lowConeFluxXmm, upperConeFluxXmm, color = 'c', alpha = 0.4)
    #plt.plot(E1, F1cl, linewidth=1, color = 'c', linestyle = 'dashed')
    #plt.plot(E1, F1cr, linewidth=1, color = 'c', linestyle = 'dashed')
    plt.plot(E2[0:Ne1-1], F2c[0:Ne1-1], linewidth = 1, label = 'Cone NuSTAR', color = 'c', linestyle = 'dotted')
    plt.fill_between(E2[0:Ne1-1], lowConeFluxNustar[0:Ne1-1], upperConeFluxNustar[0:Ne1-1], color = 'c', alpha = 0.4)
    #plt.plot(E2, F2cl, linewidth=1, color = 'c', linestyle = 'dotted')
    #plt.plot(E2, F2cr, linewidth=1, color = 'c', linestyle = 'dotted')

    ax.legend(fontsize = "10")
    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')