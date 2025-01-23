from matplotlib import animation
from pylab import *
import numpy as np

from plot_array3d import plot_array3d
from plot_array3d_animated import plot_array3d_animated
from plot_bremsstrahlung import plot_bremsstrahlung
from plot_compton_radiation import plot_compton_radiation
from plot_distributions import plot_distributions
from plot_distributions2 import plot_distributions2
from plot_distributions9 import plot_distributions9
from plot_distributions_4 import plot_distributions_4
from plot_distributions_mc import plot_distributions_mc
from plot_distributions_protons import plot_distributions_protons
from plot_error_profile import plot_error_profile
from plot_image import plot_image
from plot_image_array_animated import plot_image_array_animated
from plot_long_radiation import plot_long_radiation
from plot_mask import plot_mask
from plot_pion import plot_pion
from plot_radiation import plot_radiation
from plot_radiation10 import plot_radiation10
from plot_radiation2 import plot_radiation2
from plot_radiation3 import plot_radiation3
from plot_radiation4 import plot_radiation4
from plot_radiation_EFE import plot_radiation_EFE
from plot_rectangular_array3d import plot_rectangular_array3d
from plot_rectangular_image import plot_rectangular_image
from plot_synchrotron import plot_synchrotron
from plot_synchrotron3 import plot_synchrotron3

if __name__ == '__main__':
    plt.rcParams['image.cmap'] = 'jet'
    #plt.rcParams['image.cmap'] = 'hot'
    #plot_dummy()
    #plot_distributions()
    #plot_distributions9()
    #plot_distributions_protons()
    #plot_distributions2()
    #plot_distributions_4()
    #plot_distributions_mc()
    #plot_radiation("../W50bremsstrahlung.dat")
    #plot_radiation_EFE("../W50bremsstrahlung.dat", "W50bremsstrahlung")
    #plot_radiation_EFE("../W50pion.dat", "W50pion")
    plot_radiation_EFE("../W50compton.dat", "W50compton")
    plot_radiation_EFE("../W50synchrotron.dat", "W50synchrotron")
    #plot_synchrotron()
    #plot_synchrotron3()
    #plot_bremsstrahlung()
    #plot_pion()
    #plot_radiation4()
    #plot_radiation3()
    #plot_radiation10()
    #plot_compton_radiation()
    #plot_long_radiation()
    #plot_radiation2()
    #plot_error_profile(1,2)
    #plot_mask()
    #plot_image("../image.dat", "image")
    #plot_image("../image1.dat", "image1")
    #plot_rectangular_image("../W50bremsstrahlungImageeV.dat","W50bremsstrahlungeV")
    #plot_rectangular_image("../W50bremsstrahlungImageKeV.dat", "W50bremsstrahlungKeV")
    #plot_rectangular_image("../W50bremsstrahlungImageMeV.dat", "W50bremsstrahlungMeV")
    #plot_rectangular_image("../W50pionImageGeV.dat", "W50pionGeV")
    #plot_rectangular_image("../W50pionImageTeV.dat", "W50pionTeV")
    #plot_rectangular_image("../W50pionImagePeV.dat", "W50pionPeV")
    #plot_image_array_animated("../image_array.dat","image_array")
    #plot_array3d()
    #plot_rectangular_array3d("../concentration.dat", "concentration")
    #plot_rectangular_array3d("../distribution.dat", "distribution")
    #plot_array3d_animated("../area.dat","area")
    #plot_array3d_animated("../length.dat","length")
    #plot_array3d_animated("../volume.dat","volume")