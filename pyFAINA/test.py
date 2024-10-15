import matplotlib.pyplot as plt
from matplotlib import animation
from pylab import *
import numpy as np

def test():
    fig, ax = plt.subplots()

    ax.set_xlabel('E, GeV', fontsize=11)
    ax.set_ylabel(r'Flux, erg cm$^{-2}$ s$^{-1}$', fontsize=11, loc='bottom')

    aerr1 = [data1[:, 1] - data1lerr[:, 1], data1herr[:, 1] - data1[:, 1]]
    aerr2 = [data2[:, 1] - data2lerr[:, 1], data2herr[:, 1] - data2[:, 1]]
    ax.loglog(d1[100:999, 0], d1[100:999, 1] / 0.624, label=r"$p-p$ $\gamma$-ray radiation")
    ax.loglog(d2[1:800, 0], d4[1:800, 1] / 0.624, 'k-', linewidth=0.9, label=r"SYNCH - primary $e$")
    ax.loglog(d2[1:800, 0], d2[1:800, 1] / 0.624, 'k--', linewidth=0.9, label=r"SYNCH - secondary $e$")
    ax.loglog(d4[1:800, 0], d4[1:800, 1] / 0.624 + d2[1:800, 1] / 0.624, label=r"SYNCH total")
    ax.loglog(d5[1:100, 0] * 1e-9, d5[1:100, 1] / 0.624, '-.', label=r"IC from primary $e$", color='tab:red')
    ax.loglog(En[:], En[:] * En[:] * Rentg[:] / 624, '--', label=r"power-law: $\Gamma=2.4$", linewidth=2.9)
    ax.errorbar(1000 * data1[:, 0], data1[:, 1], yerr=aerr1, fmt=".g", markersize=6, label=r"Fermi")
    ax.errorbar(1000 * data2[:, 0], data2[:, 1], yerr=aerr2, fmt=".b", markersize=6, label=r"H.E.S.S.")
    ax.loglog(wd2radio[0] * 1.0e-9, wd2radio[1], 'r+', markersize=9, label=r"radio")  # , Eg[:], 1.0e4*Eg[:]**(-3))

    # ax.loglog(En2[:], 1.0e-5*En2[:]*En2[:]*IC[:]/624, '--', label=r"IC")
    ax.set_ylim([1.0e-17, 3.0e-11])
    ax.set_xlim([1.0e-14, 3.0e4])
    # ax.errorbar(1.0e-7, 0.65e-14, yerr=0.5e-14, uplims=True, color='black')

    # ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.tick_params(axis='y', which='both', direction='in', right=True, labelsize=11)
    ax.tick_params(axis='x', which='both', direction='in', top=True, labelsize=11)
    # ax.legend(loc=4, fontsize=11)
    ax.legend(loc='upper right', bbox_to_anchor=(0.65, 1.25),
              ncol=2, fancybox=True, shadow=True)
    ax.set_box_aspect(0.5)
    plt.show()
    fig.savefig('Synch_ru.png', dpi=400, pad_inches=3)