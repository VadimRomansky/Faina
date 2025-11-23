from matplotlib import animation
from pylab import *
import numpy as np

from plot_W50_EFE import plot_W50_EFE
from plot_W50_EFE_higenergy import plot_W50_EFE_highenergy
from plot_W50_EFE_higenergy2 import plot_W50_EFE_highenergy2
from plot_W50_EFE_higenergy3 import plot_W50_EFE_highenergy3
from plot_W50_EFE_synchrotron import plot_W50_EFE_synchrotron
from plot_W50_compton2 import plot_W50_compton2
from plot_W50_compton3 import plot_W50_compton3
from plot_W50_profile import plot_W50_profile
from plot_W50_profile_Brinkmann import plot_W50_profile_Brinkmann
from plot_array3d import plot_array3d
from plot_array3d_animated import plot_array3d_animated
from plot_bremsstrahlung import plot_bremsstrahlung
from plot_compton_radiation import plot_compton_radiation
from plot_data import plot_data
from plot_data2 import plot_data2
from plot_data3 import plot_data3
from plot_data4 import plot_data4
from plot_data5 import plot_data5
from plot_data6 import plot_data6
from plot_diffusionConvection import plot_diffusionConvection
from plot_distributions import plot_distributions
from plot_distributions2 import plot_distributions2
from plot_distributions9 import plot_distributions9
from plot_distributions_4 import plot_distributions_4
from plot_distributions_MC2 import plot_distributions_MC2
from plot_distributions_MC3 import plot_distributions_MC3
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
    #plot_diffusionConvection("../output/diffusionConvection.dat", "D", 6, yscale = "log", xscale = "symlog")

    plot_data("../output/Bturb.dat", "B", 2, yscale = "linear", xscale = "linear")

    plot_data2("../output/electrons.dat","../output/protons.dat", "distributions", 1, label1 = 'electrons', label2 = 'protons', ylim1 = None, ylim2 = None)
    plot_data2("../examples_data/W50/newPeV/electrons.dat", "../examples_data/W50/newPeV/protons.dat", "distributions_p", 1, label1='electrons',
               label2='protons', ylim1=None, ylim2=None)
    plot_data2("../examples_data/W50/newPeV/electrons.dat", "../examples_data/W50/newdistribution/electrons.dat","electrons", 1, label1 = 'maxE = 5 PeV', label2 = 'maxE = 0.5 PeV', ylim1 = 3E-2, ylim2 = 5E-1)
    plot_data2("../examples_data/W50/newPeV/protons.dat", "../examples_data/W50/newdistribution/protons.dat",
               "protons", 1, label1='maxE = 5 PeV', label2='maxE = 0.5 PeV', ylim1 = 1E-10, ylim2 = 1E-6)
    #plot_data2("../output/thinDistribution.dat","../output/thickDistribution.dat","distributions", 1, "thin", "thick", "$p/mc$", "$F(p)p^4$")
    plot_W50_profile("","profile")
    plot_W50_profile_Brinkmann("", "profile_Brinkmann")
    #plot_W50_compton2("../output/W50comptonBigSource2.dat","../output/W50thickCompton.dat","W50compton",1,"numerical","thick","E eV", "$EF(E)~erg~cm^{-2} s^{-1}$")
    #plot_W50_compton3("../output/W50comptonBigSource.dat","../output/W50comptonBigSource2.dat","../output/W50thickCompton.dat","W50compton3",1,"numerical", "numerical thicker","thick","E eV", "$EF(E)~erg~cm^{-2} s^{-1}$")
    #plot_data3("../output/distributionRight.dat", "../output/distributionMiddle.dat", "../output/distributionLeft.dat", "distributions", 1)
    #plot_data3("../W501.dat","../W502.dat","../W503.dat", "W50",1, "MC function B profile", "MC function B 20mkG", "3*advection function B 20 mkG", r'$E~eV$', r'$EF(E)~erg~cm^{-2} s^{-1}$')
    #plot_data4("../W501.dat","../W502.dat","../W503.dat", "../W504.dat", "W50_multi",1, "MC function B profile", "MC function B 20mkG", "3*advection function B 20 mkG", "MC function above 10 TeV B 20 mkG", r'$E~eV$', r'$EF(E)~erg~cm^{-2} s^{-1}$')
    #plot_data5("../W501.dat","../W502.dat","../W503.dat", "../W504.dat", "../W505.dat", "W50_multi",1, "MC function B profile", "MC function B 20mkG", "3*advection function B 20 mkG", "MC function above 10 TeV B 20 mkG", "advection above 1 TeV B 20 mkG", '$E~eV$', '$EF(E)~erg~cm^{-2} s^{-1}$')
    #plot_data6("../W501.dat","../W502.dat","../W503.dat", "../W504.dat", "../W505.dat", "../W506.dat", "W50_multi",1, "MC function B profile + upstream", "MC function B profile", "MC function B 20 mkG", "advection B profile", "advection B 20 mkG", "advection B 20 mkG above 10 TeV", '$E~eV$', '$EF(E)~erg~cm^{-2} s^{-1}$')
    #plot_data6("../W501.dat","../W502.dat","../W503.dat", "../W504.dat", "../W505.dat", "../W506.dat", "W50_multi",1, "1 GeV", "10 GeV", "100 GeV", "1 TeV", "10 TeV", "100 TeV", '$E~eV$', '$EF(E)~erg~cm^{-2} s^{-1}$')
    #plot_distributions()
    #plot_distributions9()
    #plot_distributions_proton
    #plot_distributions_MC2()
    #plot_distributions_MC3()
    #plot_distributions2()
    #plot_distributions_4()
    #plot_distributions_mc()
    #plot_radiation("../W50bremsstrahlung.dat")
    #plot_radiation_EFE("../W50bremsstrahlung.dat", "W50bremsstrahlung")
    #plot_radiation_EFE("../W50pion.dat", "W50pion")
    #plot_radiation_EFE("../W50compton.dat", "W50compton")
    #plot_radiation_EFE("../W50synchrotron.dat", "W50synchrotron")
    #plot_radiation_EFE("../W50synchandcompt.dat", "W50synchandcompt")
    factor = 5E-3
    plot_W50_EFE("../output/W50synchandcompt.dat", "W50synchandcompt", factor)
    plot_W50_EFE_highenergy("../output/W50synchandcompt.dat", "W50сompton", factor)
    plot_W50_EFE_highenergy2("../output/W50synchandcompt.dat","../output/W50thickcompton.dat","W50compton2", factor)
    plot_W50_EFE_synchrotron("../output/W50synchandcompt.dat", "W50synch", factor)
    #plot_W50_EFE_highenergy3("../output/W50synchandcompt.dat","../output/W50thickcompton.dat","../output/W50thickcompton2.dat", "W50compton3", factor)
    #plot_radiation_EFE("../W50synchandcompt2.dat", "W50synchandcompt2")
    #plot_radiation_EFE("../output/W50highenergy.dat", "W50highenergy", 1E12, 1E16)
    #plot_radiation_EFE("../output/W50synchandcompt.dat", "W50kev", 1000, 200000, 2E-12, 2E-10)
    #plot_synchrotron()
    #plot_synchrotron3()
    #plot_bremsstrahlung()
    plot_pion()
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
    #plot_rectangular_image_with_axis("../output/xE.dat", "../output/x_grid.dat", "../output/E_grid.dat", "xE", aspect = 'auto')
    #plot_image_array_animated("../image_array.dat","image_array")
    #plot_array3d()
    #plot_rectangular_array3d("../concentration.dat", "concentration")
    #plot_rectangular_array3d("../distribution.dat", "distribution")
    #plot_array3d_animated("../area.dat","area")
    #plot_array3d_animated("../length.dat","length")
    #plot_array3d_animated("../volume.dat","volume")