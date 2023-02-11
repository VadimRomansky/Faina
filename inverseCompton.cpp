#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "inverseCompton.h"

InverseComptonEvaluator::InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, PhotonIsotropicDistribution* photonDistribution) : RadiationEvaluator(Ne, Emin, Emax){
	my_Nmu = Nmu;
	my_Nphi = Nphi;

    my_photonDistribution = photonDistribution;

	my_cosTheta = new double[my_Nmu];
	my_cosThetaLeft = new double[my_Nmu];
	my_dcosTheta = new double[my_Nmu];

	//todo what if not electrons
	double thetamin = 0.1 / (my_Emax / me_c2);
	double dlogtheta = log(pi / thetamin) / (Nmu - 1);

	my_cosThetaLeft[0] = 1.0;
	my_cosTheta[0] = cos(thetamin);
	for (int i = 1; i < my_Nmu; ++i) {
		my_cosTheta[i] = cos(thetamin * exp(dlogtheta * (i)));
		my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
	}
	for (int i = 0; i < my_Nmu - 1; ++i) {
		my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
	}
	my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];
}

InverseComptonEvaluator::~InverseComptonEvaluator() {
	delete[] my_cosTheta;
	delete[] my_cosThetaLeft;
	delete[] my_dcosTheta;
}

double InverseComptonEvaluator::evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance) {
	double I = 0;
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double photonFinalCosTheta = cos(photonFinalTheta);
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double delectronEnergy;
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}

		for (int imue = 0; imue < my_Nmu; ++imue) {
			double mu_e = my_cosTheta[imue];
			for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				double dphi_e = 1.0 / my_Nphi;

				double electronDist = electronDistribution->distribution(electronInitialEnergy, mu_e, phi_e);
				//rotation
				double photonFinalCosThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(mu_e, phi_e, photonFinalCosTheta, photonFinalPhi, photonFinalCosThetaRotated, photonFinalPhiRotated);

				double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalCosThetaPrimed = (photonFinalCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonFinalCosThetaRotated);
				double photonFinalSinThetaPrimed = sqrt(1.0 - photonFinalCosThetaPrimed * photonFinalCosThetaPrimed);
				for (int imuph = 0; imuph < my_Nmu; ++imuph) {
					double mu_ph = -my_cosTheta[imuph];
					double photonInitialCosThetaPrimed = mu_ph;
					double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed * photonInitialCosThetaPrimed);
					for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
						double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
						double dphi_ph = 1.0 / my_Nphi;
						double photonInitialPhiRotated = phi_ph;

						double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * (1.0 - cosXiPrimed));
						double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
						double photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1 + photonInitialCosThetaPrimed * electronInitialBeta);

						double photonInitialCosTheta;
						double photonInitialPhi;
						inverseRotationSphericalCoordinates(mu_e, phi_e, photonInitialCosThetaRotated, photonInitialPhiRotated, photonInitialCosTheta, photonInitialPhi);

						double dI = volume * 0.5 * r2 * speed_of_light * electronDist *
							(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / (1.0 - photonFinalCosThetaRotated * electronInitialBeta)) *
							(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(1 - cosXiPrimed) / (1 - (photonFinalEnergyPrimed / m_c2) * (1 - cosXiPrimed))) *
							photonDistribution->distribution(photonInitialEnergy, photonInitialCosTheta, photonInitialPhi) *
							dphi_e * dphi_ph * my_dcosTheta[imue] * my_dcosTheta[imuph] * delectronEnergy;

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

	return I;
}

void InverseComptonEvaluator::resetParameters(const double *parameters, const double *normalizationUnits){
//todo change photon distribution
}

double InverseComptonEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);
	
	double I = 0;

	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double delectronEnergy;
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		double electronDist = electronDistribution->distribution(electronInitialEnergy);
		for (int imue = 0; imue < my_Nmu; ++imue) {
			double mu_e = my_cosTheta[imue];
			//for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				//double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				//double dphi_e = 1.0 / my_Nphi;
				//rotation
				double photonFinalCosThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(mu_e, 0.0, photonFinalCosTheta, photonFinalPhi, photonFinalCosThetaRotated, photonFinalPhiRotated);

				double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalCosThetaPrimed = (photonFinalCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonFinalCosThetaRotated);
				double photonFinalSinThetaPrimed = sqrt(1.0 - photonFinalCosThetaPrimed * photonFinalCosThetaPrimed);
				for (int imuph = 0; imuph < my_Nmu; ++imuph) {
					double mu_ph = -my_cosTheta[imuph];
					double photonInitialCosThetaPrimed = mu_ph;
					double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed * photonInitialCosThetaPrimed);
					for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
						double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
						double dphi_ph = 1.0 / my_Nphi;
						double photonInitialPhiRotated = phi_ph;

						double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * (1.0 - cosXiPrimed));
						double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
						double photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1 + photonInitialCosThetaPrimed * electronInitialBeta);

						//double photonInitialCosTheta;
						//double photonInitialPhi;
						//inverseRotationSphericalCoordinates(mu_e, phi_e, photonInitialCosThetaRotated, photonInitialPhiRotated, photonInitialCosTheta, photonInitialPhi);

                        double dI = photonFinalEnergy*volume * 0.5 * r2 * speed_of_light * electronDist *
							(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / (1.0 - photonFinalCosThetaRotated * electronInitialBeta)) *
							(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(1 - cosXiPrimed) / (1 - (photonFinalEnergyPrimed / m_c2) * (1 - cosXiPrimed))) *
                            my_photonDistribution->distribution(photonInitialEnergy) *
							2*pi * dphi_ph * my_dcosTheta[imue] * my_dcosTheta[imuph] * delectronEnergy;

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

	return I;
}

double InverseComptonEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource * source) {
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
                result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	return result;
}
