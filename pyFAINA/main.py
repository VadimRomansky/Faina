from matplotlib import animation
from pylab import *
import numpy as np

from plot_bremsstrahlung import plot_bremsstrahlung
from plot_compton_radiation import plot_compton_radiation
from plot_distributions import plot_distributions
from plot_distributions2 import plot_distributions2
from plot_distributions_protons import plot_distributions_protons
from plot_error_profile import plot_error_profile
from plot_image import plot_image
from plot_long_radiation import plot_long_radiation
from plot_mask import plot_mask
from plot_pion import plot_pion
from plot_radiation import plot_radiation
from plot_radiation2 import plot_radiation2
from plot_radiation4 import plot_radiation4
from plot_synchrotron import plot_synchrotron

if __name__ == '__main__':
    plt.rcParams['image.cmap'] = 'jet'
    #plot_dummy()
    #plot_distributions()
    #plot_distributions_protons()
    #plot_distributions2()
    #plot_radiation()
    #plot_synchrotron()
    #plot_bremsstrahlung()
    #plot_pion()
    #plot_radiation4()
    #plot_compton_radiation()
    #plot_long_radiation()
    #plot_radiation2()
    #plot_error_profile(1,2)
    plot_mask()
    plot_image()