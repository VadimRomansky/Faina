import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def plot_V4641_EFE_synchrotron(filename1, name, factor = 1.0):

    radiationFile1 = open(filename1,'r')
    lines1 = radiationFile1.readlines()
    N1 = len(lines1)

    radiation1 = np.zeros([3,N1])
    #radiation5 = np.zeros([2, N1])
    for i in range(N1):
        s = lines1[i].split()
        radiation1[0,i] = float(s[0])
        radiation1[1,i] = factor*float(s[1])
        radiation1[2,i] = factor*float(s[2])


    erositaEnergy = 0
    for i in range(N1):
        if((radiation1[0,i] > 500) and (radiation1[0,i] < 6000)):
            erositaEnergy = erositaEnergy + radiation1[2,i]*(radiation1[0,i] - radiation1[0,i-1])/radiation1[0,i]

    erositaFlux = erositaEnergy
    print('flux 0.5-6 kev = ', erositaFlux)

    integralEnergy = 0
    for i in range(N1):
        if((radiation1[0,i] > 17000) and (radiation1[0,i] < 60000)):
            integralEnergy = integralEnergy + radiation1[1,i]*(radiation1[0,i] - radiation1[0,i-1])/radiation1[0,i]

    integralFlux = integralEnergy
    print('flux 17-60 kev = ', integralFlux)

    lhaasoFile = open("../examples_data/V4641/LHAASO.dat",'r')
    lhaasoLines = lhaasoFile.readlines()
    lhaasoN = len(lhaasoLines)

    lhaaso = np.zeros([4, lhaasoN])
    lhaasoLimits = np.zeros([lhaasoN])
    for i in range(lhaasoN):
        s = lhaasoLines[i].split()
        lhaaso[0, i] = float(s[0])
        lhaaso[1, i] = float(s[1])
        lhaaso[2, i] = float(s[2]) - lhaaso[1,i]
        lhaaso[3, i] = lhaaso[1, i] - float(s[3])
        lhaasoLimits[i] = False

    hessFile = open("../examples_data/V4641/HESS.dat",'r')
    hessLines = hessFile.readlines()
    hessN = len(hessLines)

    hess = np.zeros([4, hessN])
    hessLimits = np.zeros([hessN])
    for i in range(hessN):
        s = hessLines[i].split()
        hess[0, i] = float(s[0])
        hess[1, i] = float(s[1])
        hess[2, i] = float(s[2]) - hess[1,i]
        hess[3, i] = hess[1, i] - float(s[3])
        hessLimits[i] = False

    hawcFile = open("../examples_data/V4641/HAWC.dat", 'r')
    hawcLines = hawcFile.readlines()
    hawcN = len(hawcLines)

    hawc = np.zeros([4, hawcN])
    hawcLimits = np.zeros([hawcN])
    for i in range(hawcN):
        s = hawcLines[i].split()
        hawc[0, i] = float(s[0])*1E12
        hawc[1, i] = float(s[1])*1E12*(1.6*1E-12)
        hawc[2, i] = float(s[2])*1E12*(1.6*1E-12)
        hawc[3, i] = -float(s[3])*1E12*(1.6*1E-12)
        hawcLimits[i] = False

    #fermiFile = open("../examples_data/W50/Fermi.dat", 'r')
    #fermiLines = fermiFile.readlines()
    #fermiN = len(fermiLines)

    #fermi = np.zeros([5, fermiN])
    #fermiLimits = np.zeros([fermiN])
    #for i in range(fermiN):
    #    s = fermiLines[i].split()
    #    fermi[0, i] = float(s[0])
    #    fermi[1, i] = float(s[1])
    #    fermi[2, i] = float(s[2]) - fermi[0, i]
    #    fermi[3, i] = fermi[0, i] - float(s[3])
    #    fermi[4, i] = 0.5*fermi[1, i]
    #    fermiLimits[i] = True

    meerkatFile = open("../examples_data/V4641/Meerkat.dat", 'r')
    meerkatLines = meerkatFile.readlines()
    meerkatN = len(meerkatLines)

    meerkat = np.zeros([4, meerkatN])
    meerkatLimits = np.zeros([meerkatN])
    for i in range(meerkatN):
        s = meerkatLines[i].split()
        meerkat[0, i] = float(s[0])
        meerkat[1, i] = float(s[1])
        meerkat[2, i] = float(s[2]) - meerkat[1, i]
        meerkat[3, i] = meerkat[1, i] - float(s[3])
        meerkatLimits[i] = False

    xrismFile = open("../examples_data/V4641/Xrism.dat", 'r')
    xrismLines = xrismFile.readlines()
    xrismN = len(xrismLines)

    xrism = np.zeros([4, xrismN])
    xrismLimits = np.zeros([xrismN])
    for i in range(xrismN):
        s = xrismLines[i].split()
        xrism[0, i] = float(s[0])
        xrism[1, i] = float(s[1])
        xrism[2, i] = float(s[2]) - xrism[1, i]
        xrism[3, i] = xrism[1, i] - float(s[3])
        xrismLimits[i] = False

    erosita = np.zeros([2,1])
    erosita[0,0] = 4000
    erosita[1,0] = 2E-12
    erositalimits = np.zeros([1])
    erositalimits[0] = True

    Ne = 10
    E1 = np.linspace(500, 6000, Ne)
    F1 = np.zeros([Ne])
    B1 = 2.5E-12
    A1= 0.5*B1/(np.sqrt(6000* 1.6E-12) - np.sqrt(500* 1.6E-12))

    for i in range(Ne):
        currentE = E1[i] * 1.6E-12
        F1[i] = A1 * np.sqrt(currentE)

    E2 = np.linspace(17000, 60000, Ne)
    F2 = np.zeros([Ne])
    B2 = 6.4E-12
    A2= 0.5*B2/(np.sqrt(60000* 1.6E-12) - np.sqrt(17000* 1.6E-12))

    for i in range(Ne):
        currentE = E2[i] * 1.6E-12
        F2[i] = A2 * np.sqrt(currentE)

    flux = 0
    for i in range(Ne-1):
        flux = flux + F2[i+1]*(E2[i+1] - E2[i])/E2[i+1]
    print('flux 17-60 kev 2 = ', flux)


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
    #ax.set_xlim([100, 1E15])
    #ax.set_xlim([5E10, 1E15])
    ax.set_xlim([1E2, 1E7])
    ax.set_ylim([1E-13, 1E-10])
    ax.tick_params(axis='x', size=5, width=1)
    ax.tick_params(axis='y', size=5, width=1)
    ax.minorticks_on()
    #plt.xticks([1E9, 2E9, 3E9, 4E9, 5E9, 6E9, 7E9, 8E9, 9E9, 1E10, 2E10, 3E10, 4E10, 5E10, 6E10, 7E10, 8E10, 9E10, 1E11,2E11, 3E11, 4E11,5E11, 6E11, 7E11, 8E11, 9E11,1E12,2E12, 3E12, 4E12,5E12, 6E12, 7E12, 8E12, 9E12,1E13,2E13, 3E13, 4E13,5E13,6E13,7E13,8E13,9E13,1E14,2E14,3E14,4E14,5E14,6E14,7E14,8E14,9E14,1E15])
    #plt.xticks(
    #    [5E10, 6E10, 7E10, 8E10, 9E10, 1E11, 2E11,
   #      3E11, 4E11, 5E11, 6E11, 7E11, 8E11, 9E11, 1E12, 2E12, 3E12, 4E12, 5E12, 6E12, 7E12, 8E12, 9E12, 1E13, 2E13,
   #      3E13, 4E13, 5E13, 6E13, 7E13, 8E13, 9E13, 1E14, 2E14, 3E14, 4E14, 5E14, 6E14, 7E14, 8E14, 9E14, 1E15])
    #plt.xticks([1E3, 2E3, 3E3, 4E3, 5E3, 6E3, 7E3, 8E3, 9E3, 1E4, 2E4, 3E4, 4E4, 5E4])
    #plt.yticks([2E-16, 3E-16, 4E-16, 5E-16, 6E-16, 7E-16, 8E-16, 9E-16, 1E-15, 2E-15, 3E-15, 4E-15, 5E-15, 6E-15, 7E-15, 8E-15, 9E-15, 1E-14, 2E-14, 3E-14, 4E-14, 5E-14, 6E-14, 7E-14, 8E-14, 9E-14, 1E-13, 2E-13, 3E-13, 4E-13, 5E-13, 6E-13, 7E-13, 8E-13, 9E-13, 1E-12, 2E-12, 3E-12, 4E-12, 5E-12])

    plt.plot(radiation1[0], radiation1[1], 'r', linewidth=2, label = 'full')
    plt.plot(radiation1[0], radiation1[2], 'pink', linewidth=2, label = 'in 10\'')
    plt.plot(E1, F1, 'orange', linewidth=2, label = 'eROSITA')
    plt.plot(E2, F2, 'salmon', linewidth=2, label='INTEGRAL')
    #plt.plot(radiation3[0], hessmodel[1], 'salmon', linewidth=2, label='HESS model')
    #plt.plot(radiation3[0], lhaasomodel[1], 'c', linewidth=2, label = 'LHAASO model')
    #plt.plot(radiationPion[0], radiationPion[1], 'orange', linewidth = 2, label = 'pion')

    #plt.errorbar(meerkat[0, :], meerkat[1, :], yerr=[meerkat[3, :], meerkat[2, :]], marker='s',
    #             markerfacecolor='green', markeredgecolor='green', markersize=3.5, uplims=meerkatLimits,
    #             ecolor='green', elinewidth=1.5, linewidth=0, capsize=2.5, capthick=1.5, label="MeerKAT")
    #plt.errorbar(xrism[0, :], xrism[1, :], yerr=[xrism[3, :], xrism[2, :]], marker='s',
    #             markerfacecolor='c', markeredgecolor='c', markersize=3.5, uplims=xrismLimits,
    #             ecolor='c', elinewidth=1.5, linewidth=0, capsize=2.5, capthick=1.5, label="XRISM")
    #plt.errorbar(erosita[0, :], erosita[1, :] , marker='s',
    #             markerfacecolor='c', markeredgecolor='c', markersize=3.5, uplims=erositalimits,
    #             ecolor='c', elinewidth=1.5, linewidth=0, capsize=2.5, capthick=1.5, label="e-ROSITA")
    #plt.errorbar(hess[0, :], hess[1, :], yerr=[hess[3, :], hess[2, :]], marker='s', markerfacecolor='b',
    #             markeredgecolor='b', markersize=3.5, uplims=hessLimits, ecolor='b', elinewidth=1.5, linewidth=0,
    #             capsize=2.5, capthick=1.5, label='HESS')
    #plt.errorbar(hawc[0, :], hawc[1, :], yerr=[hawc[3, :], hawc[2, :]], marker='s', markerfacecolor='purple',
    #             markeredgecolor='purple', markersize=3.5, uplims=hawcLimits, ecolor='purple', elinewidth=1.5, linewidth=0,
    #             capsize=2.5, capthick=1.5, label='HAWC')
    #plt.errorbar(lhaaso[0, :], lhaaso[1, :], yerr = [lhaaso[3, :], lhaaso[2, :]],marker='s',markerfacecolor='magenta',markeredgecolor='magenta', markersize = 3.5, uplims = lhaasoLimits, ecolor='magenta', elinewidth=1.5, linewidth=0, capsize=2.5, capthick=1.5, label = 'LHAASO')
    ax.legend(fontsize = "10", loc = 'upper left')
    #plt.plot(radiation[0], radiation[2], 'b', linewidth=4)
    #plt.errorbar(cssx1, cssy1, cssError1, ecolor = 'b', elinewidth = 4, linewidth=0, capsize = 5, capthick = 4)
    #ax.legend([r'BremsstrahlungThermalEvaluator', r'BremsstrahlungEbaluator'], fontsize="20")
    #plt.show()
    plt.savefig(name + '.png', bbox_inches='tight')