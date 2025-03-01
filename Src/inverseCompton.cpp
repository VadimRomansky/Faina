#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "inverseCompton.h"

InverseComptonEvaluator::InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, int Nph, double Ephmin, double Ephmax, PhotonDistribution* photonDistribution, const double& photonConcentration, ComptonSolverType solverType, bool absorption, bool doppler) : RadiationEvaluator(Ne, Emin, Emax, absorption, doppler){
	my_solverType = solverType;
	
	my_Ephmin = Ephmin;
	my_Ephmax = Ephmax;

	my_Nmu = Nmu;
	my_Nphi = Nphi;

	my_Nph = Nph;

    my_photonDistribution = photonDistribution;
	my_photonConcentration = photonConcentration;

	my_theta = new double[my_Nmu];
	my_cosTheta = new double[my_Nmu];
	my_cosThetaLeft = new double[my_Nmu];
	my_dcosTheta = new double[my_Nmu];

	//todo what if not electrons
	double fraction = 0.5 / my_Nmu;
	double thetamin = min(0.1/ (my_Emax / me_c2), pi/(2*my_Nmu));
	double dlogtheta = log((pi + (1-fraction)*thetamin) / (thetamin)) / (Nmu - 1);

	my_theta[0] = 0;
	my_cosThetaLeft[0] = 1.0;
	my_cosTheta[0] = cos(fraction*thetamin);
	for (int i = 1; i < my_Nmu; ++i) {
		my_theta[i] = thetamin * exp(dlogtheta * (i)) - (1 - fraction) * thetamin;
		my_cosTheta[i] = cos(my_theta[i]);
		my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
	}
	for (int i = 0; i < my_Nmu - 1; ++i) {
		my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
	}
	my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];

	/*double dtheta = pi / Nmu;

	my_cosThetaLeft[0] = 1.0;
	my_cosTheta[0] = cos(dtheta/2);
	for (int i = 1; i < my_Nmu; ++i) {
		my_cosTheta[i] = cos(dtheta*(i+0.5));
		my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
	}
	for (int i = 0; i < my_Nmu - 1; ++i) {
		my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
	}
	my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];*/
}

InverseComptonEvaluator::~InverseComptonEvaluator() {
	delete[] my_theta;
	delete[] my_cosTheta;
	delete[] my_cosThetaLeft;
	delete[] my_dcosTheta;
}

PhotonDistribution* InverseComptonEvaluator::getPhotonDistribution(const double& rho, const double& z, const double& phi) {
	return my_photonDistribution;
}

double InverseComptonEvaluator::getPhotonConcentration(const double& rho, const double& z, const double& phi)
{
	return my_photonConcentration;
}

void InverseComptonEvaluator::resetParameters(const double *parameters, const double *normalizationUnits){
//todo change photon distribution
}

double InverseComptonEvaluator::evaluateComptonEmissivityThomsonIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	MassiveParticlePowerLawDistribution* powerlawDistribution = dynamic_cast<MassiveParticlePowerLawDistribution*>(electronDistribution);
	if (powerlawDistribution == NULL) {
		printf("Thompson solver works only with powerlaw distribution\n");
		printLog("Thompson solver works only with powerlaw distribution\n");
		exit(0);
	}
	PhotonPlankDistribution* plankDistribution = dynamic_cast<PhotonPlankDistribution*>(photonDistribution);
	if (plankDistribution == NULL) {
		printf("Thompson solver works only with plank photon distribution\n");
		printLog("Thompson solver works only with plank photon distribution\n");
		exit(0);
	}
	//double photonConcentration = 3.536638846E7; //for CMB
	//double photonConcentration = 400; //for CMB
	double temperature = plankDistribution->getTemperature();
	double sigmat = 8 * pi * re2 / 3.0;
	double index = powerlawDistribution->getIndex();
	double meanPlank = 2.701178034;
	double E0 = my_Emin;
	double K = (index - 1) * pow(E0, index - 1);
	double I = photonConcentration * electronConcentration * 0.5 * sigmat * pow(massElectron * speed_of_light2, 1.0 - index) * pow(photonFinalEnergy, -(index - 1.0) / 2.0) * pow((4.0 / 3.0) * meanPlank * kBoltzman * temperature, (index - 1.0) / 2.0) * K * speed_of_light / (4 * pi);
	return I;
}

//returns fluxes in cm^-2 s^-1 = d energy flux/ d energy
double InverseComptonEvaluator::evaluateComptonEmissivityJonesIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	MassiveParticleIsotropicDistribution* isotropicElectronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(electronDistribution);
	if (isotropicElectronDistribution == NULL) {
		printf("Jones solver works only with isotropic electrons distribution\n");
		printLog("Jones solver works only with isotropic electrons distribution\n");
		exit(0);
	}

	MassiveParticleTabulatedIsotropicDistribution* tabulatedIsotropicDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);

	PhotonIsotropicDistribution* isotropicPhotonDistribution = dynamic_cast<PhotonIsotropicDistribution*>(photonDistribution);
	if (isotropicPhotonDistribution == NULL) {
		printf("Jones solver works only with isotropic photons distribution\n");
		printLog("Jones solver works only with isotropic photons distribution\n");
		exit(0);
	}

	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double* Eph = new double[my_Nph];
	double Ephmin = my_Ephmin;
	double Ephmax = my_Ephmax;
	double factor = pow(Ephmax / Ephmin, 1.0 / (my_Nph - 1));
	Eph[0] = Ephmin;
	for (int i = 1; i < my_Nph; ++i) {
		Eph[i] = Eph[i - 1] * factor;
	}


	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	double* Ee;
	int localNe = my_Ne;
	if (tabulatedIsotropicDistribution == NULL) {
		Ee = new double[my_Ne];
	}
	else {
		localNe = tabulatedIsotropicDistribution->getN();
		Ee = tabulatedIsotropicDistribution->getEnergyArray();
	}

	for (int l = 0; l < my_Nph; ++l) {
		double photonInitialEnergy = Eph[l];
		double dphotonInitialEnergy;
		if (l == 0) {
			dphotonInitialEnergy = Eph[1] - Eph[0];
		}
		else {
			dphotonInitialEnergy = Eph[l] - Eph[l - 1];
		}

		if (tabulatedIsotropicDistribution == NULL) {
			double Eemax = my_Emax;
			double tempEemax = electronDistribution->maxEnergy();
			if (tempEemax > 0) {
				if (tempEemax < Eemax) {
					Eemax = tempEemax;
				}
			}

			double Eemin = max(0.5 * photonFinalEnergy * (1.0 + sqrt(1 + m_c2 * m_c2 / (photonFinalEnergy * photonInitialEnergy))), m_c2);
			if (Eemin < my_Emin) {
				Eemin = my_Emin;
			}
			if (Eemin < electronDistribution->minEnergy()) {
				Eemin = electronDistribution->minEnergy();
			}
			if (Eemin > Eemax) {
				continue;
			}
			double Eetran = 2 * Eemin;
			if (Eemax < Eemin) {
				Eemax = 4 * Eemin;
			}
			if (Eemax <= Eetran) {
				Eetran = (Eemax + Eemin) / 2;
			}
			int linearN = my_Ne / 2;
			int logN = my_Ne - linearN;
			double delta = (Eetran - Eemin) / linearN;
			Ee[0] = Eemin;
			for (int i = 0; i < linearN; ++i) {
				Ee[i] = Eemin + delta * i;
			}
			double factor = pow(Eemax / Eetran, 1.0 / (logN - 1));
			Ee[linearN] = Eetran;
			for (int i = linearN + 1; i < my_Ne; ++i) {
				Ee[i] = Ee[i - 1] * factor;
			}
		}

		


		/*if (Eemin < 20 * m_c2) {
			printf("ok\n");
		}*/
		
		for (int k = 0; k < localNe; ++k) {
			double electronInitialEnergy = Ee[k];
			if (k > 0) {
				electronInitialEnergy = 0.5 * (Ee[k] + Ee[k - 1]);
			}
			double electronInitialGamma = electronInitialEnergy / m_c2;
			if (electronInitialGamma < photonFinalEnergy / m_c2) {
				continue;
			}
			double delectronEnergy;
			double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
			if (k == 0) {
				delectronEnergy = Ee[1] - Ee[0];
			}
			else {
				delectronEnergy = Ee[k] - Ee[k - 1];
			}
			if (delectronEnergy < 0) {
				omp_set_lock(&my_lock);
				printf("evaluateComptonFluxJonesIsotropic\n");
				printLog("evaluateComptonFluxJonesIsotropic\n");
				printf("delectron energy = %g < 0\n", delectronEnergy);
				printLog("delectron energy = %g < 0\n", delectronEnergy);
				printf("l photon energy = %d\n", l);
				printLog("l photon energy = %d\n", l);
				printf("k electron energy = %d\n", k);
				printLog("k electron energy = %d\n", k);
				printf("photon energy = %g\n", photonFinalEnergy);
				printLog("photon energy = %g\n", photonFinalEnergy);
				//printf("electron energy min = %g\n", Eemin);
				//printLog("electron energy min = %g\n", Eemin);
				//printf("electron energy tran = %g\n", Eetran);
				//printLog("electron energy tran = %g\n", Eetran);
				//printf("electron energy max = %g\n", Eemax);
				//printLog("electron energy max = %g\n", Eemax);
				for (int i = 0; i < my_Ne; ++i) {
					printf("electron energy[%d] = %g\n", i, my_Ee[i]);
					printLog("electron energy[%d] = %g\n", i, my_Ee[i]);
				}
				omp_unset_lock(&my_lock);
				exit(0);
			}

			if (electronInitialGamma > photonFinalEnergy / (m_c2)) {
				//if (electronInitialGamma < 0.1 * m_c2 / (4 * photonFinalEnergy)) {

				//} else {
				double G = 4 * electronInitialGamma * photonInitialEnergy / (m_c2);
				double q = (photonFinalEnergy / (m_c2)) / ((electronInitialGamma - photonFinalEnergy / (m_c2)) * G);
				if (q <= 1.0) {
					double sigma = 2 * pi * r2 * (2 * q * log(q) + 1 + q - 2 * q * q + 0.5 * q * q * (1 - q) * G * G / (1 + q * G)) / (electronInitialGamma * electronInitialGamma * photonInitialEnergy / m_c2);
					//divide by energy to get number of photons
					//I += electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy)  * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / (massElectron * speed_of_light2);
					double electronDist = isotropicElectronDistribution->distributionNormalized(electronInitialEnergy);
					I += photonFinalEnergy * 4 * pi * electronConcentration*electronDist * sigma * photonConcentration*isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy) * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / m_c2;
					if (I < 0) {
						omp_set_lock(&my_lock);
						printf("evaluateComptonFluxJonesIsotropic\n");
						printLog("evaluateComptonFluxJonesIsotropic\n");
						printf("I < 0\n");
						printLog("I < 0\n");
						printf("l photon energy = %d\n", l);
						printLog("l photon energy = %d\n", l);
						printf("k electron energy = %d\n", k);
						printLog("k electron energy = %d\n", k);
						printf("sigma = %g\n", sigma);
						printLog("sigma = %g\n", sigma);
						printf("photonFinalEnergy = %g\n", photonFinalEnergy);
						printLog("photonFinalEnergy = %g\n", photonFinalEnergy);
						printf("electronInitialEnergy = %g\n", electronInitialEnergy);
						printLog("electronInitialEnergy = %g\n", electronInitialEnergy);
						printf("electron concentration = %g\n", electronConcentration);
						printLog("electron concentration = %g\n", electronConcentration);
						printf("electron distribution = %g\n", isotropicElectronDistribution->distributionNormalized(electronInitialEnergy));
						printLog("electron distribution = %g\n", isotropicElectronDistribution->distributionNormalized(electronInitialEnergy));
						printf("photon initial energy = %g\n", photonInitialEnergy);
						printLog("photon initial energy = %g\n", photonInitialEnergy);
						printf("photon concentration = %g\n", photonConcentration);
						printLog("photon concentration = %g\n", photonConcentration);
						printf("photon distrinution = %g\n", isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy));
						printLog("photon distrinution = %g\n", isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy));
						printf("electron initial beta = %g\n", electronInitialBeta);
						printLog("electron initial beta = %g\n", electronInitialBeta);
						printf("delectronEnergy = %g\n", delectronEnergy);
						printLog("delectronEnergy = %g\n", delectronEnergy);
						printf("dphotonInitialEnergy = %g\n", dphotonInitialEnergy);
						printLog("dphotonInitialEnergy = %g\n", dphotonInitialEnergy);
						omp_unset_lock(&my_lock);
						exit(0);
					}

					if (I != I) {
						omp_set_lock(&my_lock);
						printf("I = NaN\n");
						printLog("I = NaN\n");
						printf("sigma = %g\n", sigma);
						printLog("sigma = %g\n", sigma);
						printf("photonFinalEnergy = %g\n", photonFinalEnergy);
						printLog("photonFinalEnergy = %g\n", photonFinalEnergy);
						printf("electronInitialEnergy = %g\n", electronInitialEnergy);
						printLog("electronInitialEnergy = %g\n", electronInitialEnergy);
						printf("electron concentration = %g\n", electronConcentration);
						printLog("electron concentration = %g\n", electronConcentration);
						printf("electron distribution = %g\n", isotropicElectronDistribution->distributionNormalized(electronInitialEnergy));
						printLog("electron distribution = %g\n", isotropicElectronDistribution->distributionNormalized(electronInitialEnergy));
						printf("photon initial energy = %g\n", photonInitialEnergy);
						printLog("photon initial energy = %g\n", photonInitialEnergy);
						printf("photon concentration = %g\n", photonConcentration);
						printLog("photon concentration = %g\n", photonConcentration);
						printf("photon distrinution = %g\n", isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy));
						printLog("photon distrinution = %g\n", isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy));
						printf("electron initial beta = %g\n", electronInitialBeta);
						printLog("electron initial beta = %g\n", electronInitialBeta);
						printf("delectronEnergy = %g\n", delectronEnergy);
						printLog("delectronEnergy = %g\n", delectronEnergy);
						printf("dphotonInitialEnergy = %g\n", dphotonInitialEnergy);
						printLog("dphotonInitialEnergy = %g\n", dphotonInitialEnergy);
						omp_unset_lock(&my_lock);
						exit(0);
					}
				}
			}
		}

	}

	delete[] Ee;
	delete[] Eph;
	return I;
}

//uses primed integration and formula 4.30
double InverseComptonEvaluator::evaluateComptonFluxKleinNishinaIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	MassiveParticleIsotropicDistribution* isotropicElectronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(electronDistribution);
	if (isotropicElectronDistribution == NULL) {
		printf("Isotropic Klein Nishina solver works only with isotropic distribution\n");
		printLog("Isotropic Klein Nishina solver works only with isotropic distribution\n");
		exit(0);
	}
	PhotonIsotropicDistribution* isotropicPhotonDistribution = dynamic_cast<PhotonIsotropicDistribution*>(photonDistribution);
	if (isotropicPhotonDistribution == NULL) {
		printf("Isotropic Klein Nishina solver works only with isotropic photons distribution\n");
		printLog("Isotropic Klein Nishina solver works only with isotropic photons distribution\n");
		exit(0);
	}
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double* cosThetaLeft = new double[my_Nmu];
	double* theta = new double[my_Nmu];
	double* cosTheta = new double[my_Nmu];
	double* dcosTheta = new double[my_Nmu];

	double photonFinalTheta = 0;
	double photonFinalPhi = 0;

	double Emin = my_Emin;
	if (Emin < electronDistribution->minEnergy()) {
		Emin = electronDistribution->minEnergy();
	}
	double Emax = my_Emax;
	double tempEemax = electronDistribution->maxEnergy();
	if (tempEemax > 0) {
		if (tempEemax < Emax) {
			Emax = tempEemax;
		}
	}

	if (Emin > Emax) {
		printf("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		printLog("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}

	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		//double thetamin = min(0.1 * m_c2/photonFinalEnergy, pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 2);

		cosThetaLeft[0] = 1.0;
		theta[0] = 0;
		cosTheta[0] = cos(theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			theta[i] = thetamin * exp(dlogtheta * (i - 1)) - (1 - fraction) * thetamin;
			cosTheta[i] = cos(theta[i]);
			cosThetaLeft[i] = (cosTheta[i] + cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			dcosTheta[i] = -(cosThetaLeft[i + 1] - cosThetaLeft[i]);
			dcosTheta[i] = versin(theta[i + 1]) - versin(theta[i]);
			/*if (theta[i] < 1E-7) {
				if (i == 0) {
					dcosTheta[i] = sin(0.5 * (theta[i] + theta[i + 1])) * (theta[i + 1] - theta[i]);
				}
				else {
					dcosTheta[i] = sin(theta[i]) * 0.5 * (theta[i + 1] - theta[i - 1]);
				}
			}*/
		}
		theta[my_Nmu - 1] = (pi+ theta[my_Nmu-2])/2;
		dcosTheta[my_Nmu - 1] = 1.0 + cosThetaLeft[my_Nmu - 1];

		double delectronEnergy;
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
		//double electronInitialDelta = 0.5 / (electronInitialGamma * electronInitialGamma) + 0.125/(electronInitialGamma*electronInitialGamma*electronInitialGamma*electronInitialGamma); //diff from 1 for large gamma
		double electronInitialDelta = relativisticDelta(electronInitialGamma); //diff from 1 for large gamma
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		double I1 = 0;
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double dthetae;
			if (imue == 0) {
				dthetae = theta[1] - theta[0];
			}
			else {
				dthetae = theta[imue] - theta[imue - 1];
			}
			/*if (imue == 0) {
				dthetae = 0.5 * (theta[1] - theta[0]);
			}
			else if (imue == my_Nmu - 1) {
				dthetae = 0.5 * (theta[my_Nmu - 1] - theta[my_Nmu - 2]);
			}
			else {
				dthetae = 0.5 * (theta[imue + 1] - theta[imue - 1]);
			}*/
			double theta_e = theta[imue];
			double sintheta_e = sin(theta_e);
			//for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				//double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				//double dphi_e = 2 * pi / my_Nphi;
			double phi_e = 0;
				double electronDist = isotropicElectronDistribution->distributionNormalized(electronInitialEnergy);
				//double phi_e = 0;
				//rotation
				double photonFinalThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
				double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
				double photonFinalVersinThetaRotated = versin(photonFinalThetaRotated);
				//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalEnergyPrimed;
				double photonFinalThetaPrimed;
				//if ((1.0 - electronInitialBeta * photonFinalCosThetaRotated) < 1E16) {
				//	photonFinalCosThetaPrimed = 1.0;
				//}
				//else {
					//photonFinalCosThetaPrimed = (photonFinalCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonFinalCosThetaRotated);
				//}
				LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
				double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
				double photonFinalVersinThetaPrimed = versin(photonFinalThetaPrimed);
				double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
				for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
					double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
					double dphi_ph = 2 * pi / my_Nphi;
					double photonInitialPhiRotated = phi_ph;
					if (photonFinalEnergyPrimed / m_c2 > 0.5) {
						//break;
						/*double a = photonFinalCosThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						double b = 1.0 -  m_c2 / photonFinalEnergyPrimed;
						double D = a * a * a * a - b * b * a * a + a * a * photonFinalCosThetaPrimed * photonFinalCosThetaPrimed;
						if (D < 0) {
							printf("D < 0\n");
						}
						else {
							double x1 = (b * photonFinalCosThetaPrimed - sqrt(D)) / (a * a + photonFinalCosThetaPrimed * photonFinalCosThetaPrimed);
							double x2 = (b * photonFinalCosThetaPrimed + sqrt(D)) / (a * a + photonFinalCosThetaPrimed * photonFinalCosThetaPrimed);
							double f1 = photonFinalCosThetaPrimed * x1 - a * sqrt(1.0 - x1 * x1);
							double f2 = photonFinalCosThetaPrimed * x2 - a * sqrt(1.0 - x2 * x2);
							double x3 = (x1 + x2) / 2.0;
							double f3 = photonFinalCosThetaPrimed * x3 - a * sqrt(1.0 - x3 * x3);
							//printf("x1 = %g x2 =  %g f1 = %g f2 = %g f3 = %g b = %g\n", x1, x2, f1, f2, f3, b);
						}*/
					}

					for (int imuph = 0; imuph < my_Nmu; ++imuph) {
						double photonInitialCosThetaPrimed = -cosTheta[imuph];
						double photonInitialThetaPrimed = pi - theta[imuph];
						//alpha = pi - theta
						double photonInitialAlphaPrimed = theta[imuph];
						/*double photonInitialThetaPrimed = pi - (imuph + 0.5) * pi / my_Nmu;
						double photonInitialAlphaPrimed = (imuph + 0.5) * pi / my_Nmu;
						double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);*/
						double dthetaph = pi / my_Nmu;
						if (imuph == 0) {
							dthetaph = theta[1] - theta[0];
						}
						else {
							dthetaph = theta[imuph] - theta[imuph - 1];
						}
						/*if (imuph == 0) {
							dthetaph = 0.5 * (theta[1] - theta[0]);
						}
						else if (imuph == my_Nmu - 1) {
							dthetaph = 0.5 * (theta[my_Nmu - 1] - theta[my_Nmu - 2]);
						}
						else {
							dthetaph = 0.5 * (theta[imuph + 1] - theta[imuph - 1]);
						}*/
						//double photonInitialVersinAlphaPrimed = versin(photonInitialAlphaPrimed);
						//double photonInitialThetaPrimed = pi - imuph*(pi/my_Nmu);
						//double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);

						double photonInitialSinThetaPrimed = sin(photonInitialAlphaPrimed);
						//double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);
						double photonInitialVersinThetaPrimed = versin(photonInitialThetaPrimed);

						//double photonInitialPhiRotated = 0;

						//double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						//double Xidelta = 1.0 - cosXiPrimed;
						double Xidelta = photonInitialVersinThetaPrimed + photonFinalVersinThetaPrimed - photonInitialVersinThetaPrimed * photonFinalVersinThetaPrimed - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						checkAndFixVersin(Xidelta);
						//double Xiepsilon = photonInitialVersinAlphaPrimed - photonFinalVersinThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						//double Xidelta = 2.0 - Xiepsilon;
						/*if (cosXiPrimed >= (1.0 - 1E-7) && photonInitialSinThetaPrimed > 0 && photonFinalSinThetaPrimed > 0) {
							Xidelta = 0.5*(photonInitialSinThetaPrimed* photonInitialSinThetaPrimed + photonFinalSinThetaPrimed* photonFinalSinThetaPrimed) - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						}*/

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * Xidelta);
						if (photonInitialEnergyPrimed <= 0) {
							continue;
						}
						//double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
						double photonInitialEnergy;
						double photonInitialThetaRotated;
						double photonInitialAlphaRotated;
						//if (1.0 + photonInitialCosThetaPrimed * electronInitialBeta < 1E16) {
						//	photonInitialCosThetaRotated = -1.0;
						//}
						//else {
							//photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1.0 + photonInitialCosThetaPrimed * electronInitialBeta);
						//}
						//LorentzTransformationPhotonReverseZ(electronInitialGamma, photonFinalEnergyPrimed, photonInitialThetaPrimed, photonInitialEnergy, photonInitialThetaRotated);
						LorentzTransformationPhotonZ(electronInitialGamma, photonInitialEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
						photonInitialThetaRotated = pi - photonInitialAlphaRotated;
						//LorentzTransformationPhotonReverseZalpha(electronInitialGamma, photonFinalEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
						//double photonInitialCosThetaRotated = cos(photonInitialThetaRotated);
						//double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
						double photonInitialVersinThetaRotated = 2.0 - versin(photonInitialAlphaRotated);
						double photonInitialTheta;
						double photonInitialPhi;
						inverseRotationSphericalCoordinates(theta_e, phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);
						//rotationSphericalCoordinates(pi - theta_e, -phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);


						double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta * photonFinalVersinThetaRotated;
						double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;
						/*if (1.0 - photonFinalCosThetaRotated * electronInitialBeta < 1E-7) {
							denom = electronInitialDelta;
						}
						else {
							denom = (1.0 - photonFinalCosThetaRotated * electronInitialBeta);
						}*/
						if (photonInitialSinThetaPrimed < 0 && photonInitialSinThetaPrimed > -1E-14) {
							printf("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
							printLog("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
						}
						double photonsN = photonConcentration * isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy);
						double dI = 0.5 * r2 * speed_of_light * electronConcentration* electronDist*
							//(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / denom) *
							(sqr(numenator) / denom) *
							//(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
							(2 - 2 * Xidelta + Xidelta * Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
							photonsN *
							//photonFinalEnergy * 2*pi * dphi_ph * dcosTheta[imue] * dcosTheta[imuph] * delectronEnergy;
						photonFinalEnergy * 2*pi * dphi_ph * photonInitialSinThetaPrimed * dthetaph * sintheta_e * dthetae * delectronEnergy;
						if (dI < 0) {
							omp_set_lock(&my_lock);
							printf("dI[i] <  0\n");
							printLog("dI[i] < 0\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						if (dI != dI) {
							omp_set_lock(&my_lock);
							printf("I[i] = NaN\n");
							printLog("I[i] = NaN\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						I += dI;
						I1 += dI;
					}
				//}
			}
		}
		//printf("%g\n", I1);
	}

	delete[] cosTheta;
	delete[] theta;
	delete[] cosThetaLeft;
	delete[] dcosTheta;

	return I;
}

//uses lab rotated integration and formula 4.31
double InverseComptonEvaluator::evaluateComptonEmissivityKleinNishinaIsotropic2(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	MassiveParticleIsotropicDistribution* isotropicElectronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(electronDistribution);
	if (isotropicElectronDistribution == NULL) {
		printf("Isotropic Klein Nishina solver works only with isotropic distribution\n");
		printLog("Isotropic Klein Nishina solver works only with isotropic distribution\n");
		exit(0);
	}
	PhotonIsotropicDistribution* isotropicPhotonDistribution = dynamic_cast<PhotonIsotropicDistribution*>(photonDistribution);
	if (isotropicPhotonDistribution == NULL) {
		printf("Isotropic Klein Nishina solver works only with isotropic photons distribution\n");
		printLog("Isotropic Klein Nishina solver works only with isotropic photons distribution\n");
		exit(0);
	}

	double photonFinalTheta = 0;
	double photonFinalPhi = 0;

	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double* cosThetaLeft = new double[my_Nmu];
	double* theta = new double[my_Nmu];
	double* cosTheta = new double[my_Nmu];
	double* dcosTheta = new double[my_Nmu];

	double Emin = my_Emin;
	if (Emin < electronDistribution->minEnergy()) {
		Emin = electronDistribution->minEnergy();
	}
	double Emax = my_Emax;
	double tempEemax = electronDistribution->maxEnergy();
	if (tempEemax > 0) {
		if (tempEemax < Emax) {
			Emax = tempEemax;
		}
	}

	if (Emin > Emax) {
		printf("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		printLog("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}

	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		//double thetamin = min(0.1 * m_c2 / photonFinalEnergy, pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 2);

		cosThetaLeft[0] = 1.0;
		theta[0] = 0;
		cosTheta[0] = cos(theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			theta[i] = thetamin * exp(dlogtheta * (i - 1)) - (1 - fraction) * thetamin;
			cosTheta[i] = cos(theta[i]);
			cosThetaLeft[i] = (cosTheta[i] + cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			dcosTheta[i] = -(cosThetaLeft[i + 1] - cosThetaLeft[i]);
			dcosTheta[i] = versin(theta[i + 1]) - versin(theta[i]);
			/*if (theta[i] < 1E-7) {
				if (i == 0) {
					dcosTheta[i] = sin(0.5 * (theta[i] + theta[i + 1])) * (theta[i + 1] - theta[i]);
				}
				else {
					dcosTheta[i] = sin(theta[i]) * 0.5 * (theta[i + 1] - theta[i - 1]);
				}
			}*/
		}
		//theta[my_Nmu - 1] = pi;
		theta[my_Nmu - 1] = (pi + theta[my_Nmu - 2]) / 2;
		dcosTheta[my_Nmu - 1] = 1.0 + cosThetaLeft[my_Nmu - 1];

		double delectronEnergy;
		//double electronInitialDelta = 0.5 / (electronInitialGamma * electronInitialGamma) + 0.125/(electronInitialGamma*electronInitialGamma*electronInitialGamma*electronInitialGamma); //diff from 1 for large gamma
		double electronInitialDelta = relativisticDelta(electronInitialGamma); //diff from 1 for large gamma
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		double I1 = 0;
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double dthetae;
			if (imue == 0) {
				dthetae = theta[1] - theta[0];
			}
			else {
				dthetae = theta[imue] - theta[imue - 1];
			}
			/*if (imue == 0) {
				dthetae = 0.5 * (theta[1] - theta[0]);
			}
			else if (imue == my_Nmu - 1) {
				dthetae = 0.5 * (theta[my_Nmu - 1] - theta[my_Nmu - 2]);
			}
			else {
				dthetae = 0.5 * (theta[imue + 1] - theta[imue - 1]);
			}*/
			double theta_e = theta[imue];
			double sintheta_e = sin(theta_e);
			//for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				double phi_e = 0;
				double dphi_e = 2 * pi;
				double electronDist = isotropicElectronDistribution->distributionNormalized(electronInitialEnergy);
				//double phi_e = 0;
				//rotation
				double photonFinalThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
				double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
				double photonFinalVersinThetaRotated = versin(photonFinalThetaRotated);
				//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalEnergyPrimed;
				double photonFinalThetaPrimed;

				LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
				double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
				double photonFinalVersinThetaPrimed = versin(photonFinalThetaPrimed);
				double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
				for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
					double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
					double dphi_ph = 2 * pi / my_Nphi;
					double photonInitialPhiRotated = phi_ph;

					for (int imuph = 0; imuph < my_Nmu; ++imuph) {
						double photonInitialThetaRotated = pi - theta[imuph];
						//alpha = pi - theta
						double photonInitialAlphaRotated = theta[imuph];
						photonInitialThetaRotated = pi - (imuph + 0.5) * pi / my_Nmu;
						photonInitialAlphaRotated = pi - photonInitialThetaRotated;
						double dthetaph = pi / my_Nmu;
						/*if (imuph == 0) {
							dthetaph = theta[1] - theta[0];
						}
						else {
							dthetaph = theta[imuph] - theta[imuph - 1];
						}*/

						double photonInitialSinThetaRotated = sin(photonInitialThetaRotated);
						double photonInitialTheta;
						double photonInitialPhi;
						inverseRotationSphericalCoordinates(theta_e, phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);

						double dummy;
						double photonInitialThetaPrimed;
						//rotationSphericalCoordinates(theta_e, phi_e, photonInitialTheta, photonInitialPhi, photonInitialThetaRotated, photonInitialPhiRotated);
						LorentzTransformationPhotonZ(electronInitialGamma, 1.0, photonInitialThetaRotated, dummy, photonInitialThetaPrimed);

						double photonInitialVersinThetaPrimed = versin(photonInitialThetaPrimed);
						double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);

						double Xidelta = photonInitialVersinThetaPrimed + photonFinalVersinThetaPrimed - photonInitialVersinThetaPrimed * photonFinalVersinThetaPrimed - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						checkAndFixVersin(Xidelta);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * Xidelta);
						if (photonInitialEnergyPrimed <= 0) {
							continue;
						}

						double photonInitialEnergy;

						double photonInitialAlphaPrimed = pi - photonInitialThetaPrimed;
						LorentzTransformationPhotonZ(electronInitialGamma, photonInitialEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);


						double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
						double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta * photonFinalVersinThetaRotated;
						double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;

						if (photonInitialSinThetaPrimed < 0 && photonInitialSinThetaPrimed > -1E-14) {
							printf("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
							printLog("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
						}
						double photonsN = photonConcentration*isotropicPhotonDistribution->distributionNormalized(photonInitialEnergy);
						double sigmaKN = (2 - 2 * Xidelta + Xidelta * Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta)));
						int index = floor(photonInitialTheta * my_Nmu / pi);
						if (index == my_Nmu) {
							index = my_Nmu - 1;
						}

						double dI = 0.5 * r2 * speed_of_light * electronConcentration*electronDist *
							(1.0 / (electronInitialGamma * electronInitialGamma * denom)) *
							sigmaKN *
							photonsN *
							//photonFinalEnergy* dphi_e * dphi_ph * dcosTheta[imue] * dcosTheta[imuph] * delectronEnergy;
							photonFinalEnergy * dphi_e * dphi_ph * photonInitialSinThetaRotated * dthetaph * sintheta_e * dthetae * delectronEnergy;
						if (dI < 0) {
							omp_set_lock(&my_lock);
							printf("dI[i] <  0\n");
							printLog("dI[i] < 0\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						if (dI != dI) {
							omp_set_lock(&my_lock);
							printf("I[i] = NaN\n");
							printLog("I[i] = NaN\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						I += dI;
						I1 += dI;
					}
				}
		}
		//printf("%g\n", I);
	}

	delete[] cosTheta;
	delete[] theta;
	delete[] cosThetaLeft;
	delete[] dcosTheta;

	return I;
}


//uses primed integration and formula 4.30 
double InverseComptonEvaluator::evaluateComptonEmissivityKleinNishinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double* cosThetaLeft = new double[my_Nmu];
	double* theta = new double[my_Nmu];
	double* cosTheta = new double[my_Nmu];
	double* dcosTheta = new double[my_Nmu];

	double* debugFlux = new double[my_Nmu];
	for (int i = 0; i < my_Nmu; ++i) {
		debugFlux[i] = 0;
	}

	double Emin = my_Emin;
	if (Emin < electronDistribution->minEnergy()) {
		Emin = electronDistribution->minEnergy();
	}
	double Emax = my_Emax;
	double tempEemax = electronDistribution->maxEnergy();
	if (tempEemax > 0) {
		if (tempEemax < Emax) {
			Emax = tempEemax;
		}
	}

	if (Emin > Emax) {
		printf("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		printLog("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}

	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		//double thetamin = min(0.1 * m_c2/photonFinalEnergy, pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 2);

		cosThetaLeft[0] = 1.0;
		theta[0] = 0;
		cosTheta[0] = cos(theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			theta[i] = thetamin * exp(dlogtheta * (i-1)) - (1 - fraction) * thetamin;
			cosTheta[i] = cos(theta[i]);
			cosThetaLeft[i] = (cosTheta[i] + cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			dcosTheta[i] = -(cosThetaLeft[i + 1] - cosThetaLeft[i]);
			dcosTheta[i] = versin(theta[i + 1]) - versin(theta[i]);
		}
		//theta[my_Nmu - 1] = pi;
		theta[my_Nmu - 1] = (pi + theta[my_Nmu - 2]) / 2;
		dcosTheta[my_Nmu - 1] = 1.0 + cosThetaLeft[my_Nmu - 1];

		double delectronEnergy;
		//double electronInitialDelta = 0.5 / (electronInitialGamma * electronInitialGamma) + 0.125/(electronInitialGamma*electronInitialGamma*electronInitialGamma*electronInitialGamma); //diff from 1 for large gamma
		double electronInitialDelta = relativisticDelta(electronInitialGamma); //diff from 1 for large gamma
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		double I1 = 0;
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double dthetae;
			if (imue == 0) {
				dthetae = theta[1] - theta[0];
			}
			else {
				dthetae = theta[imue] - theta[imue - 1];
			}

			double theta_e = theta[imue];
			double sintheta_e = sin(theta_e);
			for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				double dphi_e = 2*pi / my_Nphi;
				double electronDist = electronDistribution->distributionNormalized(electronInitialEnergy, cos(theta_e), phi_e);
			//double phi_e = 0;
			//rotation
			double photonFinalThetaRotated;
			double photonFinalPhiRotated;
			rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
			double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
			double photonFinalVersinThetaRotated = versin(photonFinalThetaRotated);
			//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
			double photonFinalEnergyPrimed;
			double photonFinalThetaPrimed;

			LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
			double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
			double photonFinalVersinThetaPrimed = versin(photonFinalThetaPrimed);
			double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
			for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
				double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
				double dphi_ph = 2 * pi / my_Nphi;
				double photonInitialPhiRotated = phi_ph;


				for (int imuph = 0; imuph < my_Nmu; ++imuph) {
					double photonInitialCosThetaPrimed = -cosTheta[imuph];
					double photonInitialThetaPrimed = pi - theta[imuph];
					//alpha = pi - theta
					double photonInitialAlphaPrimed = theta[imuph];
					/*double photonInitialThetaPrimed = pi - (imuph + 0.5) * pi / my_Nmu;
					double photonInitialAlphaPrimed = (imuph + 0.5) * pi / my_Nmu;
					double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);*/
					double dthetaph = pi / my_Nmu;
					if (imuph == 0) {
						dthetaph = theta[1] - theta[0];
					}
					else {
						dthetaph = theta[imuph] - theta[imuph - 1];
					}
					/*if (imuph == 0) {
						dthetaph = 0.5 * (theta[1] - theta[0]);
					}
					else if (imuph == my_Nmu - 1) {
						dthetaph = 0.5 * (theta[my_Nmu - 1] - theta[my_Nmu - 2]);
					}
					else {
						dthetaph = 0.5 * (theta[imuph + 1] - theta[imuph - 1]);
					}*/
					//double photonInitialVersinAlphaPrimed = versin(photonInitialAlphaPrimed);
					//double photonInitialThetaPrimed = pi - imuph*(pi/my_Nmu);
					//double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);

					double photonInitialSinThetaPrimed = sin(photonInitialAlphaPrimed);
					//double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);
					double photonInitialVersinThetaPrimed = versin(photonInitialThetaPrimed);

					//double photonInitialPhiRotated = 0;

					//double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
					//double Xidelta = 1.0 - cosXiPrimed;
					double Xidelta = photonInitialVersinThetaPrimed + photonFinalVersinThetaPrimed - photonInitialVersinThetaPrimed * photonFinalVersinThetaPrimed - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
					checkAndFixVersin(Xidelta);
					//double Xiepsilon = photonInitialVersinAlphaPrimed - photonFinalVersinThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
					//double Xidelta = 2.0 - Xiepsilon;
					/*if (cosXiPrimed >= (1.0 - 1E-7) && photonInitialSinThetaPrimed > 0 && photonFinalSinThetaPrimed > 0) {
						Xidelta = 0.5*(photonInitialSinThetaPrimed* photonInitialSinThetaPrimed + photonFinalSinThetaPrimed* photonFinalSinThetaPrimed) - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
					}*/

					double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * Xidelta);
					if (photonInitialEnergyPrimed <= 0) {
						continue;
					}
					//double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
					double photonInitialEnergy;
					double photonInitialThetaRotated;
					double photonInitialAlphaRotated;
					//if (1.0 + photonInitialCosThetaPrimed * electronInitialBeta < 1E16) {
					//	photonInitialCosThetaRotated = -1.0;
					//}
					//else {
						//photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1.0 + photonInitialCosThetaPrimed * electronInitialBeta);
					//}
					//LorentzTransformationPhotonReverseZ(electronInitialGamma, photonFinalEnergyPrimed, photonInitialThetaPrimed, photonInitialEnergy, photonInitialThetaRotated);
					LorentzTransformationPhotonZ(electronInitialGamma, photonInitialEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
					photonInitialThetaRotated = pi - photonInitialAlphaRotated;
					//LorentzTransformationPhotonReverseZalpha(electronInitialGamma, photonFinalEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
					//double photonInitialCosThetaRotated = cos(photonInitialThetaRotated);
					//double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
					double photonInitialVersinThetaRotated = 2.0 - versin(photonInitialAlphaRotated);
					double photonInitialTheta;
					double photonInitialPhi;
					inverseRotationSphericalCoordinates(theta_e, phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);
					//rotationSphericalCoordinates(pi - theta_e, -phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);


					double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta * photonFinalVersinThetaRotated;
					double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;
					/*if (1.0 - photonFinalCosThetaRotated * electronInitialBeta < 1E-7) {
						denom = electronInitialDelta;
					}
					else {
						denom = (1.0 - photonFinalCosThetaRotated * electronInitialBeta);
					}*/
					if (photonInitialSinThetaPrimed < 0 && photonInitialSinThetaPrimed > -1E-14) {
						printf("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
						printLog("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
					}
					//double photonsN = photonConcentration*photonDistribution->distributionNormalized(photonInitialEnergy, cos(photonInitialTheta), photonInitialPhi);
					double sigmaKN = (2 - 2 * Xidelta + Xidelta * Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta)));
					int index = floor(photonInitialTheta * my_Nmu / pi);
					if (index == my_Nmu) {
						index = my_Nmu - 1;
					}
					if (photonInitialEnergy > 4.8E-16 && photonInitialEnergy < 4.9E-16) {
						debugFlux[index] += (sqr(numenator) / denom) * sigmaKN;
					}
					double dI = 0.5 * r2 * speed_of_light * electronConcentration* electronDist *
						//(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / denom) *
						(sqr(numenator) / denom) *
						//(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
						sigmaKN *
						photonConcentration*photonDistribution->distributionNormalized(photonInitialEnergy, cos(photonInitialTheta), photonInitialPhi) *
						//photonFinalEnergy* dphi_e * dphi_ph * dcosTheta[imue] * dcosTheta[imuph] * delectronEnergy;
						photonFinalEnergy * dphi_e * dphi_ph * photonInitialSinThetaPrimed * dthetaph * sintheta_e * dthetae * delectronEnergy;
					if (dI < 0) {
						omp_set_lock(&my_lock);
						printf("dI[i] <  0\n");
						printLog("dI[i] < 0\n");
						printf("electron concentration = %g\n", electronConcentration);
						printLog("electron concentration = %g\n", electronConcentration);
						printf("electron distribution = %g\n", electronDist);
						printLog("electron distribution = %g\n", electronDist);
						printf("numenator = %g\n", numenator);
						printLog("numenator = %g\n", numenator);
						printf("denom = %g\n", denom);
						printLog("denom = %g\n", denom);
						printf("Xidelta = %g\n", Xidelta);
						printLog("Xidelta = %g\n", Xidelta);
						printf("photon concentration = %g\n", photonConcentration);
						printLog("photon concentration = %g\n", photonConcentration);
						printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
						printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
						printf("photon initial energy = %g\n", photonInitialEnergy);
						printLog("photon initial energy = %g\n", photonInitialEnergy);
						printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
						printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
						printf("sintheta_e = %g\n", sintheta_e);
						printLog("sintheta_e = %g\n", sintheta_e);
						omp_unset_lock(&my_lock);
						exit(0);
					}

					if (dI != dI) {
						omp_set_lock(&my_lock);
						printf("I[i] = NaN\n");
						printLog("I[i] = NaN\n");
						printf("electron concentration = %g\n", electronConcentration);
						printLog("electron concentration = %g\n", electronConcentration);
						printf("electron distribution = %g\n", electronDist);
						printLog("electron distribution = %g\n", electronDist);
						printf("numenator = %g\n", numenator);
						printLog("numenator = %g\n", numenator);
						printf("denom = %g\n", denom);
						printLog("denom = %g\n", denom);
						printf("Xidelta = %g\n", Xidelta);
						printLog("Xidelta = %g\n", Xidelta);
						printf("photon concentration = %g\n", photonConcentration);
						printLog("photon concentration = %g\n", photonConcentration);
						printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
						printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
						printf("photon initial energy = %g\n", photonInitialEnergy);
						printLog("photon initial energy = %g\n", photonInitialEnergy);
						printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
						printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
						printf("sintheta_e = %g\n", sintheta_e);
						printLog("sintheta_e = %g\n", sintheta_e);
						omp_unset_lock(&my_lock);
						exit(0);
					}

					I += dI;
					I1 += dI;
				}
			}
			}
		}
		//printf("%g\n", I1);
	}

	FILE* debugFile = fopen("debugFLux.dat", "w");
	for (int i = 0; i < my_Nmu; ++i) {
		fprintf(debugFile, "%g %g\n", i* pi / my_Nmu, debugFlux[i]);
	}
	fclose(debugFile);

	delete[] debugFlux;

	delete[] cosTheta;
	delete[] theta;
	delete[] cosThetaLeft;
	delete[] dcosTheta;

	return I;
}

//uses lab rotated integration and formula 4.31
double InverseComptonEvaluator::evaluateComptonEmissivityKleinNishinaAnisotropic2(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double* cosThetaLeft = new double[my_Nmu];
	double* theta = new double[my_Nmu];
	double* cosTheta = new double[my_Nmu];
	double* dcosTheta = new double[my_Nmu];

	double Emin = my_Emin;
	if (Emin < electronDistribution->minEnergy()) {
		Emin = electronDistribution->minEnergy();
	}
	double Emax = my_Emax;
	double tempEemax = electronDistribution->maxEnergy();
	if (tempEemax > 0) {
		if (tempEemax < Emax) {
			Emax = tempEemax;
		}
	}

	if (Emin > Emax) {
		printf("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		printLog("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}

	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		//double thetamin = min(0.1 * m_c2 / photonFinalEnergy, pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 2);

		cosThetaLeft[0] = 1.0;
		theta[0] = 0;
		cosTheta[0] = cos(theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			theta[i] = thetamin * exp(dlogtheta * (i - 1)) - (1 - fraction) * thetamin;
			cosTheta[i] = cos(theta[i]);
			cosThetaLeft[i] = (cosTheta[i] + cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			dcosTheta[i] = -(cosThetaLeft[i + 1] - cosThetaLeft[i]);
			dcosTheta[i] = versin(theta[i + 1]) - versin(theta[i]);
			/*if (theta[i] < 1E-7) {
				if (i == 0) {
					dcosTheta[i] = sin(0.5 * (theta[i] + theta[i + 1])) * (theta[i + 1] - theta[i]);
				}
				else {
					dcosTheta[i] = sin(theta[i]) * 0.5 * (theta[i + 1] - theta[i - 1]);
				}
			}*/
		}
		//theta[my_Nmu - 1] = pi;
		theta[my_Nmu - 1] = (pi + theta[my_Nmu - 2]) / 2;
		dcosTheta[my_Nmu - 1] = 1.0 + cosThetaLeft[my_Nmu - 1];

		double delectronEnergy;
		//double electronInitialDelta = 0.5 / (electronInitialGamma * electronInitialGamma) + 0.125/(electronInitialGamma*electronInitialGamma*electronInitialGamma*electronInitialGamma); //diff from 1 for large gamma
		double electronInitialDelta = relativisticDelta(electronInitialGamma); //diff from 1 for large gamma
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		double I1 = 0;
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double dthetae;
			if (imue == 0) {
				dthetae = theta[1] - theta[0];
			}
			else {
				dthetae = theta[imue] - theta[imue - 1];
			}
			/*if (imue == 0) {
				dthetae = 0.5 * (theta[1] - theta[0]);
			}
			else if (imue == my_Nmu - 1) {
				dthetae = 0.5 * (theta[my_Nmu - 1] - theta[my_Nmu - 2]);
			}
			else {
				dthetae = 0.5 * (theta[imue + 1] - theta[imue - 1]);
			}*/
			double theta_e = theta[imue];
			double costheta_e = cos(theta_e);
			double sintheta_e = sin(theta_e);
			for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				double dphi_e = 2 * pi / my_Nphi;
				double electronDist = electronDistribution->distributionNormalized(electronInitialEnergy, costheta_e, phi_e);
				//double phi_e = 0;
				//rotation
				double photonFinalThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
				double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
				double photonFinalVersinThetaRotated = versin(photonFinalThetaRotated);
				//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalEnergyPrimed;
				double photonFinalThetaPrimed;

				LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
				double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
				double photonFinalVersinThetaPrimed = versin(photonFinalThetaPrimed);
				double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
				for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
					double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
					double dphi_ph = 2 * pi / my_Nphi;
					double photonInitialPhiRotated = phi_ph;

					for (int imuph = 0; imuph < my_Nmu; ++imuph) {
						double photonInitialThetaRotated = pi - theta[imuph];
						//alpha = pi - theta
						double photonInitialAlphaRotated = theta[imuph];
						photonInitialThetaRotated = pi - (imuph + 0.5) * pi / my_Nmu;
						photonInitialAlphaRotated = pi - photonInitialThetaRotated;
						double dthetaph = pi / my_Nmu;
						/*if (imuph == 0) {
							dthetaph = theta[1] - theta[0];
						}
						else {
							dthetaph = theta[imuph] - theta[imuph - 1];
						}*/

						double photonInitialSinThetaRotated = sin(photonInitialThetaRotated);
						double photonInitialTheta;
						double photonInitialPhi;
						inverseRotationSphericalCoordinates(theta_e, phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);

						double dummy;
						double photonInitialThetaPrimed;
						//rotationSphericalCoordinates(theta_e, phi_e, photonInitialTheta, photonInitialPhi, photonInitialThetaRotated, photonInitialPhiRotated);
						LorentzTransformationPhotonZ(electronInitialGamma, 1.0, photonInitialThetaRotated, dummy, photonInitialThetaPrimed);

						double photonInitialVersinThetaPrimed = versin(photonInitialThetaPrimed);
						double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);

						double Xidelta = photonInitialVersinThetaPrimed + photonFinalVersinThetaPrimed - photonInitialVersinThetaPrimed * photonFinalVersinThetaPrimed - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						checkAndFixVersin(Xidelta);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * Xidelta);
						if (photonInitialEnergyPrimed <= 0) {
							continue;
						}

						double photonInitialEnergy;

						double photonInitialAlphaPrimed = pi - photonInitialThetaPrimed;
						LorentzTransformationPhotonZ(electronInitialGamma, photonInitialEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);


						double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
						double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta * photonFinalVersinThetaRotated;
						double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;

						if (photonInitialSinThetaPrimed < 0 && photonInitialSinThetaPrimed > -1E-14) {
							printf("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
							printLog("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
						}
						double photonsN = photonConcentration*photonDistribution->distributionNormalized(photonInitialEnergy, cos(photonInitialTheta), photonInitialPhi);
						double sigmaKN = (2 - 2 * Xidelta + Xidelta * Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta)));
						/*int index = floor(photonInitialTheta * my_Nmu / pi);
						if (index == my_Nmu) {
							index = my_Nmu - 1;
						}*/

						double dI = 0.5 * r2 * speed_of_light * electronConcentration*electronDist *
							(1.0 / (electronInitialGamma*electronInitialGamma*denom)) *
							sigmaKN *
							photonsN *
							//photonFinalEnergy* dphi_e * dphi_ph * dcosTheta[imue] * dcosTheta[imuph] * delectronEnergy;
							photonFinalEnergy * dphi_e * dphi_ph * photonInitialSinThetaRotated * dthetaph * sintheta_e * dthetae * delectronEnergy;
						if (dI < 0) {
							omp_set_lock(&my_lock);
							printf("dI[i] <  0\n");
							printLog("dI[i] < 0\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						if (dI != dI) {
							omp_set_lock(&my_lock);
							printf("I[i] = NaN\n");
							printLog("I[i] = NaN\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						I += dI;
						I1 += dI;
					}
				}
			}
		}
		//printf("%g\n", I);
	}

	delete[] cosTheta;
	delete[] theta;
	delete[] cosThetaLeft;
	delete[] dcosTheta;

	return I;
}

//uses lab not rotated integration and formula 4.31
double InverseComptonEvaluator::evaluateComptonEmissivityKleinNishinaAnisotropic3(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double* cosThetaLeft = new double[my_Nmu];
	double* theta = new double[my_Nmu];
	double* cosTheta = new double[my_Nmu];
	double* dcosTheta = new double[my_Nmu];

	double Emin = my_Emin;
	if (Emin < electronDistribution->minEnergy()) {
		Emin = electronDistribution->minEnergy();
	}
	double Emax = my_Emax;
	double tempEemax = electronDistribution->maxEnergy();
	if (tempEemax > 0) {
		if (tempEemax < Emax) {
			Emax = tempEemax;
		}
	}

	if (Emin > Emax) {
		printf("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		printLog("Emin > Emax in evaluate flux Klein Nishina Isotropic\n");
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}

	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		double electronInitialGamma = electronInitialEnergy / m_c2;
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		//double thetamin = min(0.1 * m_c2 / photonFinalEnergy, pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 2);

		cosThetaLeft[0] = 1.0;
		theta[0] = 0;
		cosTheta[0] = cos(theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			theta[i] = thetamin * exp(dlogtheta * (i - 1)) - (1 - fraction) * thetamin;
			cosTheta[i] = cos(theta[i]);
			cosThetaLeft[i] = (cosTheta[i] + cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			dcosTheta[i] = -(cosThetaLeft[i + 1] - cosThetaLeft[i]);
			dcosTheta[i] = versin(theta[i + 1]) - versin(theta[i]);
			/*if (theta[i] < 1E-7) {
				if (i == 0) {
					dcosTheta[i] = sin(0.5 * (theta[i] + theta[i + 1])) * (theta[i + 1] - theta[i]);
				}
				else {
					dcosTheta[i] = sin(theta[i]) * 0.5 * (theta[i + 1] - theta[i - 1]);
				}
			}*/
		}
		//theta[my_Nmu - 1] = pi;
		theta[my_Nmu - 1] = (pi + theta[my_Nmu - 2]) / 2;
		dcosTheta[my_Nmu - 1] = 1.0 + cosThetaLeft[my_Nmu - 1];

		double delectronEnergy;
		//double electronInitialDelta = 0.5 / (electronInitialGamma * electronInitialGamma) + 0.125/(electronInitialGamma*electronInitialGamma*electronInitialGamma*electronInitialGamma); //diff from 1 for large gamma
		double electronInitialDelta = relativisticDelta(electronInitialGamma); //diff from 1 for large gamma
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		double I1 = 0;
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double dthetae;
			if (imue == 0) {
				dthetae = theta[1] - theta[0];
			}
			else {
				dthetae = theta[imue] - theta[imue - 1];
			}
			/*if (imue == 0) {
				dthetae = 0.5 * (theta[1] - theta[0]);
			}
			else if (imue == my_Nmu - 1) {
				dthetae = 0.5 * (theta[my_Nmu - 1] - theta[my_Nmu - 2]);
			}
			else {
				dthetae = 0.5 * (theta[imue + 1] - theta[imue - 1]);
			}*/
			double theta_e = theta[imue];
			double sintheta_e = sin(theta_e);
			for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				double dphi_e = 2 * pi / my_Nphi;
				double electronDist = electronDistribution->distributionNormalized(electronInitialEnergy, cos(theta_e), phi_e);
				//double phi_e = 0;
				//rotation
				double photonFinalThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
				double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
				double photonFinalVersinThetaRotated = versin(photonFinalThetaRotated);
				//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				double photonFinalEnergyPrimed;
				double photonFinalThetaPrimed;

				LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
				double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
				double photonFinalVersinThetaPrimed = versin(photonFinalThetaPrimed);
				double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
				for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
					double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
					double dphi_ph = 2 * pi / my_Nphi;
					double photonInitialPhi = phi_ph;

					for (int imuph = 0; imuph < my_Nmu; ++imuph) {
						double photonInitialTheta = pi - theta[imuph];
						//alpha = pi - theta
						double photonInitialAlpha = theta[imuph];
						photonInitialTheta = pi - (imuph + 0.5) * pi / my_Nmu;
						photonInitialAlpha = pi - photonInitialTheta;
						double dthetaph = pi / my_Nmu;
						/*if (imuph == 0) {
							dthetaph = theta[1] - theta[0];
						}
						else {
							dthetaph = theta[imuph] - theta[imuph - 1];
						}*/

						double photonInitialSinTheta = sin(photonInitialTheta);
						double photonInitialThetaRotated;
						double photonInitialPhiRotated;
						rotationSphericalCoordinates(theta_e, phi_e, photonInitialTheta, photonInitialPhi, photonInitialThetaRotated, photonInitialPhiRotated);

						double dummy;
						double photonInitialThetaPrimed;
						//rotationSphericalCoordinates(theta_e, phi_e, photonInitialTheta, photonInitialPhi, photonInitialThetaRotated, photonInitialPhiRotated);
						LorentzTransformationPhotonZ(electronInitialGamma, 1.0, photonInitialThetaRotated, dummy, photonInitialThetaPrimed);

						double photonInitialVersinThetaPrimed = versin(photonInitialThetaPrimed);
						double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);

						double Xidelta = photonInitialVersinThetaPrimed + photonFinalVersinThetaPrimed - photonInitialVersinThetaPrimed * photonFinalVersinThetaPrimed - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
						checkAndFixVersin(Xidelta);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * Xidelta);
						if (photonInitialEnergyPrimed <= 0) {
							continue;
						}

						double photonInitialEnergy;

						double photonInitialAlphaPrimed = pi - photonInitialThetaPrimed;
						double photonInitialAlphaRotated;
						LorentzTransformationPhotonZ(electronInitialGamma, photonInitialEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);


						double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
						double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta * photonFinalVersinThetaRotated;
						double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;

						if (photonInitialSinThetaPrimed < 0 && photonInitialSinThetaPrimed > -1E-14) {
							printf("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
							printLog("photon Initial Sin Theta Primed = %g reduced to 0\n", photonInitialSinThetaPrimed);
						}
						double photonsN = photonConcentration*photonDistribution->distributionNormalized(photonInitialEnergy, cos(photonInitialTheta), photonInitialPhi);
						double sigmaKN = (2 - 2 * Xidelta + Xidelta * Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta)));
						int index = floor(photonInitialTheta * my_Nmu / pi);
						if (index == my_Nmu) {
							index = my_Nmu - 1;
						}

						double dI = 0.5 * r2 * speed_of_light * electronConcentration*electronDist *
							(1.0 / (electronInitialGamma * electronInitialGamma * denom)) *
							sigmaKN *
							photonsN *
							//photonFinalEnergy* dphi_e * dphi_ph * dcosTheta[imue] * dcosTheta[imuph] * delectronEnergy;
							photonFinalEnergy * dphi_e * dphi_ph * photonInitialSinTheta * dthetaph * sintheta_e * dthetae * delectronEnergy;
						if (dI < 0) {
							omp_set_lock(&my_lock);
							printf("dI[i] <  0\n");
							printLog("dI[i] < 0\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						if (dI != dI) {
							omp_set_lock(&my_lock);
							printf("I[i] = NaN\n");
							printLog("I[i] = NaN\n");
							printf("electron concentration = %g\n", electronConcentration);
							printLog("electron concentration = %g\n", electronConcentration);
							printf("electron distribution = %g\n", electronDist);
							printLog("electron distribution = %g\n", electronDist);
							printf("numenator = %g\n", numenator);
							printLog("numenator = %g\n", numenator);
							printf("denom = %g\n", denom);
							printLog("denom = %g\n", denom);
							printf("Xidelta = %g\n", Xidelta);
							printLog("Xidelta = %g\n", Xidelta);
							printf("photon concentration = %g\n", photonConcentration);
							printLog("photon concentration = %g\n", photonConcentration);
							printf("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printLog("photon final energy primed = %g\n", photonFinalEnergyPrimed);
							printf("photon initial energy = %g\n", photonInitialEnergy);
							printLog("photon initial energy = %g\n", photonInitialEnergy);
							printf("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printLog("photon initial sin theta primed = %g\n", photonInitialSinThetaPrimed);
							printf("sintheta_e = %g\n", sintheta_e);
							printLog("sintheta_e = %g\n", sintheta_e);
							omp_unset_lock(&my_lock);
							exit(0);
						}

						I += dI;
						I1 += dI;
					}
				}
			}
		}
		//printf("%g\n", I);
	}

	delete[] cosTheta;
	delete[] theta;
	delete[] cosThetaLeft;
	delete[] dcosTheta;

	return I;
}

double InverseComptonEvaluator::evaluateComptonEmissivityKleinNishinaAnisotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration) {
	//return evaluateComptonFluxKleinNishinaAnisotropic(photonFinalEnergy, 0, 0, photonDistribution, electronDistribution, volume, distance);
	return evaluateComptonEmissivityKleinNishinaAnisotropic2(photonFinalEnergy, 0, 0, photonDistribution, electronDistribution, photonConcentration, electronConcentration);
}

/*double InverseComptonEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi) {
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double rho = source->getRho(irho);
	double phi = source->getPhi(iphi);

	double result = 0;
	for (int iz = 0; iz < Nz; ++iz) {
		double z = source->getZ(iz);
		result += evaluateFluxFromFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance(), getPhotonDistribution(rho, z, phi));
	}
	return result;
}*/

double InverseComptonEvaluator::evaluateFluxFromSourceAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, const double& photonConcentration, RadiationSourceInCylindrical* source, ComptonSolverType solverType) {
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;
	int irho = 0;

	omp_init_lock(&my_lock);

#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nrho, Nz, Nphi) reduction(+:result)

	for (irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double electronConcentration = source->getConcentration(irho, iz, iphi);
				if (solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA) {
					result += evaluateComptonEmissivityKleinNishinaAnisotropic2(photonFinalEnergy, photonFinalTheta, photonFinalPhi, photonDistribution, source->getParticleDistribution(irho, iz, iphi), photonConcentration, electronConcentration)*source->getVolume(irho, iz, iphi);
				}
				//else if (solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA2) {
				//	result += evaluateComptonEmissivityKleinNishinaAnisotropic2(photonFinalEnergy, photonFinalTheta, photonFinalPhi, photonDistribution, source->getParticleDistribution(irho, iz, iphi), photonConcentration, electronConcentration) * source->getVolume(irho, iz, iphi);
				//}
				//else if (solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA3) {
				//	result += evaluateComptonEmissivityKleinNishinaAnisotropic3(photonFinalEnergy, photonFinalTheta, photonFinalPhi, photonDistribution, source->getParticleDistribution(irho, iz, iphi), photonConcentration, electronConcentration) * source->getVolume(irho, iz, iphi);
				//}
				else {
					printf("wrong Compton solver type for anisotropic flux\n");
					printLog("wrong Compton solver type for anisotropic flux\n");
					exit(0);
				}
			}
		}
	}

	omp_destroy_lock(&my_lock);

	return result/sqr(source->getDistance());
}

double InverseComptonEvaluator::evaluateTotalFluxInEnergyRangeAnisotropic(const double& Ephmin, const double& Ephmax, const double& photonFinalTheta, const double& photonFinalPhi, int Nph, PhotonDistribution* photonDistribution, const double& photonConcentration, RadiationSourceInCylindrical* source, ComptonSolverType solverType ) {
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	double flux = 0;
	for (int i = 0; i < Nph; ++i) {
		printf("%d\n", i);
		double dE = currentE * (factor - 1.0);
		flux += evaluateFluxFromSourceAnisotropic(currentE, photonFinalTheta, photonFinalPhi, photonDistribution, photonConcentration, source, solverType) * dE;
		currentE = currentE * factor;
	}
	return flux;
}

double InverseComptonEvaluator::evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	double x1 = source->getX1(ix1);
	double z = source->getZ(iz);
	double x2 = source->getX2(ix2);
	double photonConcentration = getPhotonConcentration(x1, z, x2);
	double electronConcentration = source->getConcentration(ix1, iz, ix2);
	if (my_solverType == ComptonSolverType::ISOTROPIC_THOMSON) {
		return evaluateComptonEmissivityThomsonIsotropic(photonFinalEnergy, getPhotonDistribution(x1, z, x2), source->getParticleDistribution(ix1, iz, ix2), photonConcentration, electronConcentration);
	}
	else if (my_solverType == ComptonSolverType::ISOTROPIC_JONES) {
		return evaluateComptonEmissivityJonesIsotropic(photonFinalEnergy, getPhotonDistribution(x1, z, x2), source->getParticleDistribution(ix1, iz, ix2), photonConcentration, electronConcentration);
	}
	else if (my_solverType == ComptonSolverType::ISOTROPIC_KLEIN_NISHINA) {
		return evaluateComptonEmissivityKleinNishinaIsotropic2(photonFinalEnergy, getPhotonDistribution(x1, z, x2), source->getParticleDistribution(ix1, iz, ix2), photonConcentration, electronConcentration);
	}
	//else if (my_solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA) {
	//	return evaluateComptonEmissivityKleinNishinaAnisotropic(photonFinalEnergy, 0, 0, getPhotonDistribution(x1, z, x2), source->getParticleDistribution(ix1, iz, ix2), photonConcentration, electronConcentration);
	//}
	else if (my_solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA) {
		return evaluateComptonEmissivityKleinNishinaAnisotropic2(photonFinalEnergy, 0, 0, getPhotonDistribution(x1, z, x2), source->getParticleDistribution(ix1, iz, ix2), photonConcentration, electronConcentration);
	}
	//else if (my_solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA3) {
	//	return evaluateComptonEmissivityKleinNishinaAnisotropic3(photonFinalEnergy, 0, 0, getPhotonDistribution(x1, z, x2), source->getParticleDistribution(ix1, iz, ix2), photonConcentration, electronConcentration);
	//}
	else {
		printf("unknown compton solver type\n");
		printLog("unknown compton solver type\n");
		exit(0);
	}
}

double InverseComptonEvaluator::evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	return 0.0;
}

InverseComptonEvaluatorWithSource::InverseComptonEvaluatorWithSource(int Ne, int Nmu, int Nphi, double Emin, double Emax, int Nph, double Ephmin, double Ephmax, PhotonDistribution* photonDistribution, const double& photonConcentration, ComptonSolverType solverType, const double& sourceR, const double& sourceZ, const double& sourcePhi, bool absorption, bool doppler) : InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution, photonConcentration, solverType, absorption, doppler) {
	my_sourceR = sourceR;
	my_sourceZ = sourceZ;
	my_sourcePhi = sourcePhi;
	my_defaultPhotonConcentration = photonConcentration;
}
InverseComptonEvaluatorWithSource::~InverseComptonEvaluatorWithSource() {

}
PhotonDistribution* InverseComptonEvaluatorWithSource::getPhotonDistribution(const double& rho, const double& z, const double& phi) {
	return my_photonDistribution;
}

double InverseComptonEvaluatorWithSource::getPhotonConcentration(const double& rho, const double& z, const double& phi)
{
	double x = rho * cos(phi);
	double y = rho * sin(phi);
	double sourceX = my_sourceR * cos(my_sourcePhi);
	double sourceY = my_sourceR * sin(my_sourcePhi);
	double r = sqrt(sqr(x - sourceX) + sqr(y - sourceY) + sqr(z - my_sourceZ));
	double r0 = sqrt(my_sourceR * my_sourceR + my_sourceZ * my_sourceZ);

	double newConcentration = my_defaultPhotonConcentration * sqr(r0 / r);
	return newConcentration;
}
