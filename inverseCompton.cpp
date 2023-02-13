#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "inverseCompton.h"

InverseComptonEvaluator::InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, PhotonIsotropicDistribution* photonDistribution, COMPTON_SOLVER_TYPE solverType) : RadiationEvaluator(Ne, Emin, Emax){
	my_solverType = solverType;
	
	my_Nmu = Nmu;
	my_Nphi = Nphi;

    my_photonDistribution = photonDistribution;

	my_cosTheta = new double[my_Nmu];
	my_cosThetaLeft = new double[my_Nmu];
	my_dcosTheta = new double[my_Nmu];

	//todo what if not electrons
	double thetamin = 0.01 / (my_Emax / me_c2);
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

double InverseComptonEvaluator::evaluateComptonLuminocityKleinNisinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance) {
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

double InverseComptonEvaluator::evaluateComptonLuminocityThompsonIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	MassiveParticlePowerLawDistribution* powerlawDistribution = dynamic_cast<MassiveParticlePowerLawDistribution*>(electronDistribution);
	if (powerlawDistribution == NULL) {
		printf("Thompson solver works only with powerlaw distribution\n");
		printLog("Thompson solver works only with powerlaw distribution\n");
		exit(0);
	}
	//double photonConcentration = 3.536638846E7; //for CMB
	//double photonConcentration = 400; //for CMB
	double photonConcentration = photonDistribution->getConcentration(); //for CMB
	double electronConcentration = electronDistribution->getConcentration();
	double sigmat = 8 * pi * re2 / 3.0;
	double index = powerlawDistribution->getIndex();
	double E0 = my_Emin;
	double K = (index - 1) * pow(E0, index - 1);
	double I = photonConcentration * electronConcentration * volume * 0.5 * sigmat * pow(massElectron * speed_of_light2, 1.0 - index) * pow(photonFinalEnergy, -(index - 1.0) / 2.0) * pow((4.0 / 3.0) * 2.7012 * kBoltzman * 2.7, (index - 1.0) / 2.0) * K * speed_of_light / (4 * pi);
	return I/(4*pi*sqr(distance));
}
double InverseComptonEvaluator::evaluateComptonLuminocityKangJonesIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	int Nph = 1000;
	double* Eph = new double[Nph];
	double Ephmin = 0.001*kBoltzman*2.7;
	double Ephmax = 10000 * Ephmin;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	Eph[0] = Ephmin;
	for (int i = 1; i < Nph; ++i) {
		Eph[i] = Eph[i - 1] * factor;
	}

	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double delectronEnergy;
		double electronInitialGamma = electronInitialEnergy / (massElectron * speed_of_light2);
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		for (int l = 0; l < Nph; ++l) {
			double photonInitialEnergy = Eph[l];
			double dphotonInitialEnergy;
			if (l == 0) {
				dphotonInitialEnergy = Eph[1] - Eph[0];
			}
			else {
				dphotonInitialEnergy = Eph[l] - Eph[l - 1];
			}

			if (electronInitialGamma > photonFinalEnergy / (massElectron * speed_of_light2)) {
				double G = 4 * electronInitialGamma * photonInitialEnergy / (massElectron * speed_of_light2);
				double q = (photonFinalEnergy / (massElectron * speed_of_light2)) / ((electronInitialGamma - photonFinalEnergy / (massElectron * speed_of_light2)) * G);
				if (q <= 1.0) {
					double sigma = 2 * pi * re2 * (2 * q * log(q) + 1 + q - 2 * q * q + 0.5 * q * q * (1 - q) * G * G / (1 + q * G)) / (electronInitialGamma * electronInitialGamma * photonInitialEnergy / (massElectron * speed_of_light2));
					//divide by energy to get number of photons
					//I += electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy)  * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / (massElectron * speed_of_light2);
					I += photonFinalEnergy*electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy) * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / (massElectron * speed_of_light2);
					if (I < 0) {
						printf("I < 0\n");
						printLog("I < 0\n");
						exit(0);
					}

					if (I != I) {
						printf("I = NaN\n");
						printLog("I = NaN\n");
						exit(0);
					}
				}
			}
		}

	}

	delete[] Eph;
	return I / (distance * distance);
}

double InverseComptonEvaluator::evaluateComptonLuminocityKleinNisinaIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
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
		double I1 = 0;
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
                            photonDistribution->distribution(photonInitialEnergy) *
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

						I += dI;
						I1 += dI;
					}
				}
		}
		//printf("%g\n", I1);
	}

	return I / (distance * distance);
}

double InverseComptonEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	if (my_solverType == COMPTON_SOLVER_TYPE::ISOTROPIC_THOMPSON) {
		return evaluateComptonLuminocityThompsonIsotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else if (my_solverType == COMPTON_SOLVER_TYPE::ISOTROPIC_KANG_JONES) {
		return evaluateComptonLuminocityKangJonesIsotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else if (my_solverType == COMPTON_SOLVER_TYPE::KLEIN_NISINA) {
		return evaluateComptonLuminocityKleinNisinaIsotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else {
		printf("unknown compton solver type\n");
		printLog("unknown compton solver type\n");
		exit(0);
	}
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
