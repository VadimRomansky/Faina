#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "electronDistribution.h"
#include "photonDistribution.h"

#include "inverseCompton.h"

double evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, ElectronDistribution* electronDistribution, const double& volume, const double& distance, const double& Emin, const double& Emax, const int Ne, const int Nmu, const int Nphi)
{
	double* cosTheta = new double[Nmu];
	double* cosThetaLeft = new double[Nmu];
	double* dcosTheta = new double[Nmu];

	double thetamin = 0.001 / (Emax/me_c2);
	double dlogtheta = log(pi / thetamin) / (Nmu - 1);

	cosThetaLeft[0] = 1.0;
	cosTheta[0] = cos(thetamin);
	for (int i = 1; i < Nmu; ++i) {
		cosTheta[i] = cos(thetamin * exp(dlogtheta * (i)));
		cosThetaLeft[i] = (cosTheta[i] + cosTheta[i - 1]) / 2.0;
	}
	for (int i = 0; i < Nmu - 1; ++i) {
		dcosTheta[i] = -(cosThetaLeft[i + 1] - cosThetaLeft[i]);
	}
	dcosTheta[Nmu - 1] = 1.0 + cosThetaLeft[Nmu - 1];

	double factor = pow(Emax / Emin, 1.0 / (Ne - 1));

	double* Ee = new double[Ne];
	Ee[0] = Emin;
	for (int i = 1; i < Ne; ++i) {
		Ee[i] = Ee[i - 1] * factor;
	}

	double I = 0;

	double photonFinalCosTheta = cos(photonFinalTheta);
	for (int k = 0; k < Ne; ++k) {
		double electronInitialEnergy = Ee[k];
		double delectronEnergy;
		double electronInitialGamma = electronInitialEnergy / (massElectron * speed_of_light2);
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
		if (k == 0) {
			delectronEnergy = Ee[1] - Ee[0];
		}
		else {
			delectronEnergy = Ee[k] - Ee[k - 1];
		}		

		for (int imue = 0; imue < Nmu; ++imue) {
			double mu_e = cosTheta[imue];
			for (int iphie = 0; iphie < Nphi; ++iphie) {			
				double phi_e = 2 * pi * (iphie + 0.5) / Nphi;
				double dphi_e = 1.0 / Nphi;

				double electronDist = electronDistribution->distribution(electronInitialEnergy, mu_e, phi_e);
				//rotation
				double photonFinalCosThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(mu_e, phi_e, photonFinalCosTheta, photonFinalPhi, photonFinalCosThetaRotated, photonFinalPhiRotated);
				//double photonFinalEnergyPrimed;
				//double photonFinalCosThetaPrimed;
				//LorentzTransformationPhotonZ(electronInitialBeta, photonFinalEnergy, photonFinalCosThetaRotated, photonFinalEnergyPrimed, photonFinalCosThetaPrimed);
				double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalCosThetaPrimed = (photonFinalCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonFinalCosThetaRotated);
				double photonFinalSinThetaPrimed = sqrt(1.0 - photonFinalCosThetaPrimed * photonFinalCosThetaPrimed);
				for (int imuph = 0; imuph < Nmu; ++imuph) {
					double mu_ph = -cosTheta[imuph];
					for (int iphiph = 0; iphiph < Nphi; ++iphiph) {
						double phi_ph = 2 * pi * (iphiph + 0.5) / Nphi;
						double dphi_ph = 1.0 / Nphi;
						double photonInitialCosThetaPrimed = mu_ph;
						double photonInitialPhiRotated = phi_ph;
						//inverseRotationSphericalCoordinates(mu_e, phi_e, mu_ph, phi_ph, photonInitialCosThetaRotated, photonInitialPhiRotated);
						//double photonInitialCosThetaPrimed = (photonInitialCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonInitialCosThetaRotated);

						double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed * photonInitialCosThetaPrimed);

						double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (me_c2)) * (1.0 - cosXiPrimed));
						double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
						double photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1 + photonInitialCosThetaPrimed * electronInitialBeta);

						double photonInitialCosTheta;
						double photonInitialPhi;
						inverseRotationSphericalCoordinates(mu_e, phi_e, photonInitialCosThetaRotated, photonInitialPhiRotated, photonInitialCosTheta, photonInitialPhi);

						double dI = volume * 0.5 * re2 * speed_of_light * electronDist *
							(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / (1.0 - photonFinalCosThetaRotated * electronInitialBeta)) *
							(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / me_c2) * sqr(1 - cosXiPrimed) / (1 - (photonFinalEnergyPrimed / me_c2) * (1 - cosXiPrimed))) *
							photonDistribution->distribution(photonInitialEnergy, photonInitialCosTheta, photonInitialPhi) *
							dphi_e * dphi_ph * dcosTheta[imue] * dcosTheta[imuph] * delectronEnergy;

						if (dI < 0) {
							printf("dI[i] <  0\n");
							printLog("dI[i] < 0\n");
							exit(0);
						}

						if (dI != dI) {
							printf("I[i] = NaN\n");
							printLog("I[i] = NaN\n");
							exit(0);
						}

						I += dI / (distance * distance);
					}
				}
			}
		}
	}

	delete[] cosTheta;
	delete[] cosThetaLeft;
	delete[] dcosTheta;

	delete[] Ee;

	return I;
}