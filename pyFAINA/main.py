from matplotlib import animation
from pylab import *
import numpy as np

from plot_W50_EFE import plot_W50_EFE
from plot_array3d import plot_array3d
from plot_array3d_animated import plot_array3d_animated
from plot_bremsstrahlung import plot_bremsstrahlung
from plot_compton_radiation import plot_compton_radiation
from plot_data import plot_data
from plot_data3 import plot_data3
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
from plot_rectangular_image_with_axis import plot_rectangular_image_with_axis
from plot_synchrotron import plot_synchrotron
from plot_synchrotron3 import plot_synchrotron3

if __name__ == '__main__':
    plt.rcParams['image.cmap'] = 'jet'
    #plt.rcParams['image.cmap'] = 'hot'
    #plot_dummy()
    #plot_data("../nishina_losses.dat","nishina_losses", 3, 'log', 'log')
    #plot_data("../Bturb.dat", "B", 2)
    plot_data3("../distributionRight.dat", "../distributionMiddle.dat", "../distributionLeft.dat", "distributions", 1)
    #plot_data3("../W50_Bconst_2E-5_CMB+140.dat","../W50_Bprofile_CMB+140.dat","../W50_Bprofile_CMB.dat", "W50",1, "Bconst CMB+140", "Bprofile CMB+140", "Bprofile+CMB", r'$E~eV$', r'$EF(E)~erg~cm^{-2} s^{-1}$')
    #plot_distributions()
    #plot_distributions9()
    #plot_distributions_protons()
    #plot_distributions2()
    #plot_distributions_4()
    #plot_distributions_mc()
    #plot_radiation("../W50bremsstrahlung.dat")
    #plot_radiation_EFE("../W50bremsstrahlung.dat", "W50bremsstrahlung")
    #plot_radiation_EFE("../W50pion.dat", "W50pion")
    #plot_radiation_EFE("../W50compton.dat", "W50compton")
    #plot_radiation_EFE("../W50synchrotron.dat", "W50synchrotron")
    #plot_radiation_EFE("../W50synchandcompt.dat", "W50synchandcompt")
    plot_W50_EFE("../W50synchandcompt.dat", "W50synchandcompt")
    #plot_radiation_EFE("../W50synchandcompt2.dat", "W50synchandcompt2")
    #plot_radiation_EFE("../W50highenergy.dat", "W50highenergy")
    #plot_radiation_EFE("../W50highenergy2.dat", "W50highenergy2")
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
    #plot_rectangular_image("../W50bremsstrahlungImageeV.dat","W50bremsstrahlungeV", -2.4688e+20, 2.4688e+20, 0, 1.2344e+20)
    #plot_rectangular_image("../W50bremsstrahlungImageKeV.dat", "W50bremsstrahlungKeV", -2.4688e+20, 2.4688e+20, 0, 1.2344e+20)
    #plot_rectangular_image("../W50bremsstrahlungImageMeV.dat", "W50bremsstrahlungMeV", -2.4688e+20, 2.4688e+20, 0, 1.2344e+20)
    #plot_rectangular_image("../W50pionImageGeV.dat", "W50pionGeV", -2.4688e+20, 2.4688e+20, 0, 1.2344e+20)
    #plot_rectangular_image("../W50pionImageTeV.dat", "W50pionTeV", -2.4688e+20, 2.4688e+20, 0, 1.2344e+20)
    #plot_rectangular_image("../W50pionImagePeV.dat", "W50pionPeV", -2.4688e+20, 2.4688e+20, 0, 1.2344e+20)
    #plot_rectangular_image("../W50scImageeV.dat", "W50sceV", 0, 5E20, 1E18, 2E18, aspect = 'auto')
    #plot_rectangular_image("../W50scImageKeV.dat", "W50scKeV", 0, 5E20, 1E18, 2E18, aspect = 'auto')
    #plot_rectangular_image("../W50scImageMeV.dat", "W50scMeV", 0, 5E20, 1E18, 2E18, aspect = 'auto')
    #plot_rectangular_image("../W50scImageGeV.dat", "W50scGeV", 0, 5E20, 1E18, 2E18, aspect = 'auto')
    #plot_rectangular_image("../W50scImageTeV.dat", "W50scTeV", 0, 5E20, 1E18, 2E18, aspect = 'auto')
    #plot_rectangular_image("../W50scImagePeV.dat", "W50scPeV", 0, 5E20, 1E18, 2E18, aspect = 'auto')

    #plot_rectangular_image("../xE.dat", "xE", 0, 5E29, 1, 16, aspect = 'auto')
    plot_rectangular_image_with_axis("../xE.dat", "../x_grid.dat", "../E_grid.dat", "xE", 0, 5E29, 1, 16, aspect = 'auto')
    #plot_image_array_animated("../image_array.dat","image_array")
    #plot_array3d()
    #plot_rectangular_array3d("../concentration.dat", "concentration")
    #plot_rectangular_array3d("../distribution.dat", "distribution")
    #plot_array3d_animated("../area.dat","area")
    #plot_array3d_animated("../length.dat","length")
    #plot_array3d_animated("../volume.dat","volume")