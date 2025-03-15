from matplotlib import animation
from pylab import *
import numpy as np

def plot_distributions_MC2():
    x_data1 = np.loadtxt('../examples_data/W50/lowfield0.4/MC/x_grid.dat')
    p_data1 = np.loadtxt('../examples_data/W50/lowfield0.4/MC/p_grid.dat')
    dist1 = np.loadtxt('../examples_data/W50/lowfield0.4/MC/pdf.dat')

    x_data2 = np.loadtxt('../examples_data/W50/lowfield0.4/advection/x_grid.dat')
    p_data2 = np.loadtxt('../examples_data/W50/lowfield0.4/advection/p_grid.dat')
    dist2 = np.loadtxt('../examples_data/W50/lowfield0.4/advection/pdf.dat')

    Nx1 = x_data1.shape[0]
    Nx2 = x_data2.shape[0]
    N1 = p_data1.shape[0]
    N2 = p_data2.shape[0]

    diff = Nx1 - Nx2

    distribution1 = np.zeros([Nx1, N1])
    distribution2 = np.zeros([Nx2, N2])

    count = 0
    for i in range(Nx1):
        for j in range(N1):
            distribution1[i][j] = dist1[count][1]
            count = count + 1

    count = 0
    for i in range(Nx2):
        for j in range(N2):
            distribution2[i][j] = dist2[count][1]
            count = count + 1


    plt.rcParams.update({'font.size': 40})
    plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.linewidth'] = 4
    f1 = plt.figure(figsize=[10, 10])
    ax = f1.add_subplot(111)
    #plt.subplots_adjust(bottom=0.12)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    #ax.tick_params(axis='both', which='minor', labelsize=8)
    #ax.set_xlabel('$\\rm E_{kin}/m_p c^2$', fontsize=40,fontweight='bold')
    ax.set_xlabel('$\\rm p/m_p c$', fontsize=40,fontweight='bold')
    #ax.set_ylabel('$\\rm f(E) E^2$', fontsize=40, fontweight='bold')
    ax.set_ylabel('$\\rm f(p/m_p c) (p/m_p c)^4$', fontsize=40,fontweight='bold')
    ax.set_yscale("log")
    ax.set_xscale("log")
    #extraticks=[1,100]
    #plt.xticks(list(plt.xticks()[0]))
    ax.tick_params(axis='x', size=10, width=4)
    ax.tick_params(axis='y', size=10, width=4)
    #ax.minorticks_on()
    # plt.axis([0.0,1.0,0.0,1.0])

    n = 1
    plt.plot(p_data1, distribution1[diff + n,0:N1],'r', linewidth = 4)
    plt.plot(p_data2, 0.2*distribution2[n, 0:N2],'r', linewidth = 4, linestyle = 'dashed')

    n = 10
    plt.plot(p_data1, distribution1[diff + n,0:N1],'g', linewidth = 4)
    plt.plot(p_data2, 0.2*distribution2[n, 0:N2],'g', linewidth = 4, linestyle = 'dashed')

    n = 20
    plt.plot(p_data1, distribution1[diff + n,0:N1],'b', linewidth = 4)
    plt.plot(p_data2, 0.2*distribution2[n, 0:N2],'b', linewidth = 4, linestyle = 'dashed')

    n = 50
    plt.plot(p_data1, distribution1[diff + n, 0:N1], 'magenta', linewidth=4)
    plt.plot(p_data2, 0.2*distribution2[n, 0:N2], 'magenta', linewidth=4, linestyle='dashed')

    ax.set_xlim(1E1, 1E7)
    ax.set_ylim(5E-7,2E-5)

    plt.savefig('distribution_diffusion and advection.png', bbox_inches='tight')
