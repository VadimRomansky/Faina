#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "inverseCompton.h"

InverseComptonEvaluator::InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, double Ephmin, double Ephmax, PhotonIsotropicDistribution* photonDistribution, ComptonSolverType solverType) : RadiationEvaluator(Ne, Emin, Emax){
	my_solverType = solverType;
	
	my_Ephmin = Ephmin;
	my_Ephmax = Ephmax;

	my_Nmu = Nmu;
	my_Nphi = Nphi;

    my_photonDistribution = photonDistribution;

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

void InverseComptonEvaluator::outputDifferentialFlux(const char* fileName) {
	FILE* outFile = fopen(fileName, "w");

	double photonFinalEnergy = 1E15*1.6E-12;
	double photonFinalTheta = 0;
	double photonFinalPhi = 0;

	double electronInitialEnergy = 1E10 * massElectron * speed_of_light2;
	double electronInitialTheta = 1E-10;
	double electronInitialPhi = 0;

	double photonInitialTheta = pi - 1E-10;
	double photonInitialPhi = 0;

	for (int i = 0; i < my_Nmu; ++i) {
	//for (int i = 0; i < my_Nphi; ++i) {
	//for (int i = 0; i < my_Ne; ++i) {
		photonInitialTheta = pi-my_theta[i];
		//photonFinalTheta = my_theta[i];
		//electronInitialTheta = my_theta[i];
		//photonFinalPhi = (i + 0.5) * 2 * pi / my_Nphi;
		//photonInitialPhi = (i + 0.5) * 2 * pi / my_Nphi;
		//electronInitialPhi = (i + 0.5) * 2 * pi / my_Nphi;
		//photonFinalEnergy = my_Ee[i];
		//electronInitialEnergy = my_Ee[i];
		fprintf(outFile, "%15.10g %15.10g  %15.10g\n", photonInitialTheta, my_dcosTheta[i]*evaluateDifferentialFlux(photonFinalEnergy, photonFinalTheta, photonFinalPhi, electronInitialEnergy, electronInitialTheta, electronInitialPhi, photonInitialTheta, photonInitialPhi), my_dcosTheta[i]);
	}

	fclose(outFile);
}

void InverseComptonEvaluator::outputDifferentialFluxJones(const char* fileName, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution) {
	FILE* outFile = fopen(fileName, "w");

	double photonFinalEnergy = 1E18 * 1.6E-12;
	double photonFinalTheta = 0;
	double photonFinalPhi = 0;

	double electronInitialEnergy = 1E10 * massElectron * speed_of_light2;
	double electronInitialGamma = electronInitialEnergy / me_c2;
	double photonInitialEnergy = kBoltzman * 10000;


	int Nph = 100;
	double* Eph = new double[Nph];
	double Ephmin = my_Ephmin;
	double Ephmax = my_Ephmax;
	Ephmin = photonFinalEnergy / (4 * electronInitialGamma * (electronInitialGamma - photonFinalEnergy / me_c2));
	if (Ephmin <= 0) {
		Ephmin = photonFinalEnergy / (4 * electronInitialGamma * electronInitialGamma);
	}
	Ephmax = 2 * photonFinalEnergy;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	Eph[0] = Ephmin;
	for (int j = 1; j < Nph; ++j) {
		Eph[j] = Eph[j - 1] * factor;
	}

	for (int i = 0; i < Nph; ++i) {
	//for (int i = 0; i < my_Ne; ++i) {
		photonFinalEnergy = Eph[i];
		//photonInitialEnergy = Eph[i];
		//electronInitialEnergy = my_Ee[i];


		fprintf(outFile, "%15.10g  %15.10g\n", photonFinalEnergy, evaluateDifferentialFluxJones(photonFinalEnergy, electronInitialEnergy, photonInitialEnergy, photonDistribution, electronDistribution));
		//fprintf(outFile, "%15.10g  %15.10g\n", photonInitialEnergy, evaluateDifferentialFluxJones(photonFinalEnergy, electronInitialEnergy, photonInitialEnergy, photonDistribution, electronDistribution));
		//fprintf(outFile, "%15.10g  %15.10g\n", electronInitialEnergy, evaluateDifferentialFluxJones(photonFinalEnergy, electronInitialEnergy, photonInitialEnergy, photonDistribution, electronDistribution));
	}


	delete[] Eph;
	fclose(outFile);

	FILE* file1 = fopen("output3.dat", "w");
	for (int i = 0; i < my_Ne; ++i) {
		double electronInitialEnergy = my_Ee[i];
		if (i > 0) {
			electronInitialEnergy = 0.5 * (my_Ee[i] + my_Ee[i - 1]);
		}
		double delectronEnergy;
		if (i == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[i] - my_Ee[i - 1];
		}
		fprintf(file1, "%g %g %g %g\n", electronInitialEnergy, delectronEnergy, electronDistribution->distributionNormalized(electronInitialEnergy), electronInitialEnergy * delectronEnergy * electronDistribution->distributionNormalized(electronInitialEnergy));
	}
	fclose(file1);
}

double InverseComptonEvaluator::evaluateDifferentialFlux(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, const double& electronInitialEnergy, const double& theta_e, const double& phi_e, const double& photonInitialThetaPrimed, const double& photonInitialPhiRotated) {
	double m = massElectron;
	double m_c2 = m * speed_of_light2;
	double fraction = 0.5 / my_Nmu;
	double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
	double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 1);

	my_cosThetaLeft[0] = 1.0;
	my_theta[0] = 0;
	my_cosTheta[0] = cos(fraction * thetamin);
	for (int i = 1; i < my_Nmu; ++i) {
		my_theta[i] = thetamin * exp(dlogtheta * (i)) - (1 - fraction) * thetamin;
		my_cosTheta[i] = cos(my_theta[i]);
		my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
	}
	for (int i = 0; i < my_Nmu - 1; ++i) {
		my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
		if (my_theta[i] < 1E-7) {
			if (i == 0) {
				my_dcosTheta[i] = sin(my_theta[i]) * (my_theta[i + 1] - my_theta[i]);
			}
			else {
				my_dcosTheta[i] = sin(my_theta[i]) * 0.5 * (my_theta[i + 1] - my_theta[i - 1]);
			}
		}
	}
	my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];

	double electronInitialGamma = electronInitialEnergy / m_c2;
	double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
	double electronInitialDelta = 0.5 / (electronInitialGamma * electronInitialGamma); //diff from 1 for large gamma

	//rotation
	double photonFinalThetaRotated;
	double photonFinalPhiRotated;
	rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
	double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
	//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
	double photonFinalEnergyPrimed;
	double photonFinalThetaPrimed;

	LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
	double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
	double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
	double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);

	double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);


	double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
	double Xidelta = 1.0 - cosXiPrimed;
	if (cosXiPrimed >= (1.0 - 1E-10) && photonInitialSinThetaPrimed > 0 && photonFinalSinThetaPrimed > 0) {
		Xidelta = 0.5 * (photonInitialSinThetaPrimed * photonInitialSinThetaPrimed + photonFinalSinThetaPrimed * photonFinalSinThetaPrimed) - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
	}

	double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * Xidelta);
	if (photonInitialEnergyPrimed <= 0) {
		return 0;
	}
				
	double photonInitialEnergy;
	double photonInitialThetaRotated;

	LorentzTransformationPhotonReverseZ(electronInitialGamma, photonFinalEnergyPrimed, photonInitialThetaPrimed, photonInitialEnergy, photonInitialThetaRotated);
	double photonInitialCosThetaRotated = cos(photonInitialThetaRotated);
				

	double denom;
	if (1.0 - photonFinalCosThetaRotated * electronInitialBeta == 0) {
		denom = electronInitialDelta;
	}
	else {
		denom = (1.0 - photonFinalCosThetaRotated * electronInitialBeta);
	}
	double result = (sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / denom) *
		(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
		photonFinalEnergy;
	if (result < 0) {
		printf("dI[i] <  0\n");
		printLog("dI[i] < 0\n");
		exit(0);
	}

	if (result != result) {
		printf("I[i] = NaN\n");
		printLog("I[i] = NaN\n");
		exit(0);
	}




	return result;
}

double InverseComptonEvaluator::evaluateDifferentialFluxJones(const double& photonFinalEnergy, const double& electronInitialEnergy, const double& photonInitialEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution)
{
	double electronInitialGamma = electronInitialEnergy / me_c2;
	if (electronInitialGamma < photonFinalEnergy / me_c2) {
		return 0;
	}
	double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
	double G = 4 * electronInitialGamma * photonInitialEnergy / (me_c2);
	double q = (photonFinalEnergy / (me_c2)) / ((electronInitialGamma - photonFinalEnergy / (me_c2)) * G);
	if (q <= 1.0) {
		double sigma = 2 * pi * re2 * (2 * q * log(q) + 1 + q - 2 * q * q + 0.5 * q * q * (1 - q) * G * G / (1 + q * G)) / (electronInitialGamma * electronInitialGamma * photonInitialEnergy / me_c2);
		//divide by energy to get number of photons
		//I += electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy)  * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / (massElectron * speed_of_light2);
		double I = electronDistribution->distributionNormalized(electronInitialEnergy)*photonDistribution->distributionNormalized(photonInitialEnergy)*photonFinalEnergy * 4 * pi * sigma * speed_of_light * electronInitialBeta / me_c2;
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
		return I;
	}
	return 0;
}

void InverseComptonEvaluator::resetParameters(const double *parameters, const double *normalizationUnits){
//todo change photon distribution
}

double InverseComptonEvaluator::evaluateComptonFluxThomsonIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
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
	double photonConcentration = photonDistribution->getConcentration(); //for CMB
	double electronConcentration = electronDistribution->getConcentration();
	double sigmat = 8 * pi * re2 / 3.0;
	double index = powerlawDistribution->getIndex();
	double meanPlank = 2.701178034;
	double E0 = my_Emin;
	double K = (index - 1) * pow(E0, index - 1);
	double I = photonConcentration * electronConcentration * volume * 0.5 * sigmat * pow(massElectron * speed_of_light2, 1.0 - index) * pow(photonFinalEnergy, -(index - 1.0) / 2.0) * pow((4.0 / 3.0) * meanPlank * kBoltzman * temperature, (index - 1.0) / 2.0) * K * speed_of_light / (4 * pi);
	return I/(sqr(distance));
}
/*double InverseComptonEvaluator::evaluateComptonFluxJonesIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	int Nph = 100;
	double* Eph = new double[Nph];
	double Ephmin = my_Ephmin;
	double Ephmax = my_Ephmax;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	Eph[0] = Ephmin;
	for (int i = 1; i < Nph; ++i) {
		Eph[i] = Eph[i - 1] * factor;
	}


	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		if (k > 0) {
			electronInitialEnergy = 0.5 * (my_Ee[k] + my_Ee[k - 1]);
		}
		double electronInitialGamma = electronInitialEnergy / m_c2;
		if (electronInitialGamma < photonFinalEnergy / m_c2) {
			continue;
		}

		//Ephmin = photonFinalEnergy - (electronInitialEnergy - m_c2);
		Ephmin = photonFinalEnergy / (4 * electronInitialGamma * (electronInitialGamma - photonFinalEnergy/m_c2));
		if (Ephmin <= 0) {
			Ephmin = photonFinalEnergy/ (4*electronInitialGamma*electronInitialGamma);
		}
		Ephmax = 2 * photonFinalEnergy;
		factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		Eph[0] = Ephmin;
		for (int i = 1; i < Nph; ++i) {
			Eph[i] = Eph[i - 1] * factor;
		}

		double delectronEnergy;
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

			if (electronInitialGamma > photonFinalEnergy / (m_c2)) {
				double G = 4 * electronInitialGamma * photonInitialEnergy / (m_c2);
				double q = (photonFinalEnergy / (m_c2)) / ((electronInitialGamma - photonFinalEnergy / (m_c2)) * G);
				if (q <= 1.0) {
					double sigma = 2 * pi * r2 * (2 * q * log(q) + 1 + q - 2 * q * q + 0.5 * q * q * (1 - q) * G * G / (1 + q * G)) / (electronInitialGamma * electronInitialGamma * photonInitialEnergy / m_c2);
					//divide by energy to get number of photons
					//I += electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy)  * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / (massElectron * speed_of_light2);
					I += photonFinalEnergy*4*pi*electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy) * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / m_c2;
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
}*/

double InverseComptonEvaluator::evaluateComptonFluxJonesIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	int Nph = 100;
	double* Eph = new double[Nph];
	double Ephmin = my_Ephmin;
	double Ephmax = my_Ephmax;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	Eph[0] = Ephmin;
	for (int i = 1; i < Nph; ++i) {
		Eph[i] = Eph[i - 1] * factor;
	}


	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	for (int l = 0; l < Nph; ++l) {
		double photonInitialEnergy = Eph[l];
		double dphotonInitialEnergy;
		if (l == 0) {
			dphotonInitialEnergy = Eph[1] - Eph[0];
		}
		else {
			dphotonInitialEnergy = Eph[l] - Eph[l - 1];
		}

		double Eemin = max(0.5 * photonFinalEnergy * (1.0 + sqrt(1 + m_c2 * m_c2 / (photonFinalEnergy * photonInitialEnergy))), m_c2);
		if (Eemin < my_Emin) {
			Eemin = my_Emin;
		}
		double Eemax = my_Emax;
		double Eetran = 2 * Eemin;
		if (Eemax < Eemin) {
			Eemax = 2 * Eemin;
		}
		if (Eemax < Eetran) {
			Eetran = (Eemax + Eemin) / 2;
		}
		int linearN = my_Ne / 2;
		int logN = my_Ne - linearN;
		double delta = (Eetran - Eemin) / linearN;
		my_Ee[0] = Eemin;
		for (int i = 0; i < linearN; ++i) {
			my_Ee[i] = Eemin + delta * i;
		}
		double factor = pow(Eemax / Eetran, 1.0 / (logN - 1));
		my_Ee[linearN] = Eetran;
		for (int i = linearN+1; i < my_Ne; ++i) {
			my_Ee[i] = my_Ee[i - 1] * factor;
		}



		
		for (int k = 0; k < my_Ne; ++k) {
			double electronInitialEnergy = my_Ee[k];
			if (k > 0) {
				electronInitialEnergy = 0.5 * (my_Ee[k] + my_Ee[k - 1]);
			}
			double electronInitialGamma = electronInitialEnergy / m_c2;
			if (electronInitialGamma < photonFinalEnergy / m_c2) {
				continue;
			}
			double delectronEnergy;
			double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
			if (k == 0) {
				delectronEnergy = my_Ee[1] - my_Ee[0];
			}
			else {
				delectronEnergy = my_Ee[k] - my_Ee[k - 1];
			}

			if (electronInitialGamma > photonFinalEnergy / (m_c2)) {
				double G = 4 * electronInitialGamma * photonInitialEnergy / (m_c2);
				double q = (photonFinalEnergy / (m_c2)) / ((electronInitialGamma - photonFinalEnergy / (m_c2)) * G);
				if (q <= 1.0) {
					double sigma = 2 * pi * r2 * (2 * q * log(q) + 1 + q - 2 * q * q + 0.5 * q * q * (1 - q) * G * G / (1 + q * G)) / (electronInitialGamma * electronInitialGamma * photonInitialEnergy / m_c2);
					//divide by energy to get number of photons
					//I += electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy)  * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / (massElectron * speed_of_light2);
					I += photonFinalEnergy * 4 * pi * electronDistribution->distribution(electronInitialEnergy) * volume * sigma * photonDistribution->distribution(photonInitialEnergy) * speed_of_light * electronInitialBeta * delectronEnergy * dphotonInitialEnergy / m_c2;
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

double InverseComptonEvaluator::evaluateComptonFluxKleinNishinaIsotropic1(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);
	
	double I = 0;

	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	double photonFinalTheta = 0;
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 1);

		my_cosThetaLeft[0] = 1.0;
		my_theta[0] = 0;
		my_cosTheta[0] = cos(my_theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			my_theta[i] = thetamin * exp(dlogtheta * (i)) - (1 - fraction) * thetamin;
			my_cosTheta[i] = cos(my_theta[i]);
			my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
			if (my_theta[i] < 1E-7) {
				if (i == 0) {
					my_dcosTheta[i] = sin(0.5*(my_theta[i] + my_theta[i+1])) * (my_theta[i + 1] - my_theta[i]);
				}
				else {
					my_dcosTheta[i] = sin(my_theta[i]) * 0.5* (my_theta[i + 1] - my_theta[i-1]);
				}
			}
		}
		my_theta[my_Nmu - 1] = pi;
		my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];
		
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
		double electronDist = electronDistribution->distribution(electronInitialEnergy);
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double dthetae;
			if (imue == 0) {
				dthetae = my_theta[1] - my_theta[0];
			}
			else {
				dthetae = my_theta[imue] - my_theta[imue - 1];
			}
			double theta_e = my_theta[imue];
			double sintheta_e = sin(theta_e);
			//for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				//double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				//double dphi_e = 2*pi / my_Nphi;
				double phi_e = 0;
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
				for (int imuph = 0; imuph < my_Nmu; ++imuph) {
					double photonInitialCosThetaPrimed = -my_cosTheta[imuph];
					double photonInitialThetaPrimed = pi - my_theta[imuph];
					//alpha = pi - theta
					double photonInitialAlphaPrimed = my_theta[imuph];
					double dthetaph;
					if (imuph == 0) {
						dthetaph = my_theta[1] - my_theta[0];
					}
					else {
						dthetaph = my_theta[imuph] - my_theta[imuph-1];
					}
					//double photonInitialVersinAlphaPrimed = versin(photonInitialAlphaPrimed);
					//double photonInitialThetaPrimed = pi - imuph*(pi/my_Nmu);
					//double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);
					
					double photonInitialSinThetaPrimed = sin(photonInitialAlphaPrimed);
					//double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);
					double photonInitialVersinThetaPrimed = versin(photonInitialThetaPrimed);
					for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
						double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
						double dphi_ph = 2*pi / my_Nphi;
						double photonInitialPhiRotated = phi_ph;
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
						LorentzTransformationPhotonReverseZ(electronInitialGamma, photonFinalEnergyPrimed, photonInitialThetaPrimed, photonInitialEnergy, photonInitialThetaRotated);
						LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
						//LorentzTransformationPhotonReverseZalpha(electronInitialGamma, photonFinalEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
						//double photonInitialCosThetaRotated = cos(photonInitialThetaRotated);
						//double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
						double photonInitialVersinThetaRotated = 2.0 - versin(photonInitialAlphaRotated);
						//double photonInitialCosTheta;
						//double photonInitialPhi;
						//inverseRotationSphericalCoordinates(mu_e, phi_e, photonInitialCosThetaRotated, photonInitialPhiRotated, photonInitialCosTheta, photonInitialPhi);

						
						double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta*photonFinalVersinThetaRotated;
						double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;
						/*if (1.0 - photonFinalCosThetaRotated * electronInitialBeta < 1E-7) {
							denom = electronInitialDelta;
						}
						else {
							denom = (1.0 - photonFinalCosThetaRotated * electronInitialBeta);
						}*/
                        double dI = volume * 0.5 * r2 * speed_of_light * electronDist *
							//(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / denom) *
							(sqr(numenator) / denom) *
							//(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
							(2 - 2*Xidelta + Xidelta*Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
                            photonDistribution->distribution(photonInitialEnergy) *
							//photonFinalEnergy*2*pi * dphi_ph * my_dcosTheta[imue] * my_dcosTheta[imuph] * delectronEnergy;
							photonFinalEnergy*2*pi * dphi_ph * photonInitialSinThetaPrimed * dthetaph * sintheta_e * dthetae * delectronEnergy;
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

double InverseComptonEvaluator::evaluateComptonFluxKleinNishinaIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double I = 0;

	double photonFinalCosTheta = 1.0;
	double photonFinalPhi = 0.0;
	double photonFinalTheta = 0;
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];

		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 1);

		my_cosThetaLeft[0] = 1.0;
		my_theta[0] = 0;
		my_cosTheta[0] = cos(my_theta[0]);
		for (int i = 1; i < my_Nmu; ++i) {
			my_theta[i] = thetamin * exp(dlogtheta * (i)) - (1 - fraction) * thetamin;
			my_cosTheta[i] = cos(my_theta[i]);
			my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
			if (my_theta[i] < 1E-7) {
				if (i == 0) {
					my_dcosTheta[i] = sin(0.5 * (my_theta[i] + my_theta[i + 1])) * (my_theta[i + 1] - my_theta[i]);
				}
				else {
					my_dcosTheta[i] = sin(my_theta[i]) * 0.5 * (my_theta[i + 1] - my_theta[i - 1]);
				}
			}
		}
		my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];

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
		double electronDist = electronDistribution->distribution(electronInitialEnergy);
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}
		for (int imue = 0; imue < my_Nmu; ++imue) {
			//double mu_e = my_cosTheta[imue];
			double theta_e = my_theta[imue];
			//for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				//double phi_e = 2 * pi * (iphie + 0.5) / my_Nphi;
				//double dphi_e = 2*pi / my_Nphi;
			double phi_e = 0;
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
			for (int imuph = 0; imuph < my_Nmu; ++imuph) {
				double photonInitialCosThetaRotated = -my_cosTheta[imuph];
				double photonInitialThetaRotated = pi - my_theta[imuph];
				//double photonInitialThetaPrimed = pi - my_theta[imuph];
				//alpha = pi - theta
				//double photonInitialAlphaRotated = my_theta[imuph];
				//double photonInitialVersinAlphaPrimed = versin(photonInitialAlphaRotated);
				double photonInitialEnergy = 1.0;
				double photonInitialEnergyPrimed = 1.0;
				double photonInitialThetaPrimed;
				//double photonInitialThetaPrimed = pi - imuph*(pi/my_Nmu);
				//double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);
				LorentzTransformationPhotonZ(electronInitialGamma, photonInitialEnergy, photonInitialThetaRotated, photonInitialEnergyPrimed, photonInitialThetaPrimed);
				double photonInitialCosThetaPrimed = cos(photonInitialThetaPrimed);
				double photonInitialSinThetaPrimed = sin(photonInitialThetaPrimed);
				//double photonInitialSinTheta = sin(photonInitialAlpha);
				for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
					double phi_ph = 2 * pi * (iphiph + 0.5) / my_Nphi;
					double dphi_ph = 2 * pi / my_Nphi;
					double photonInitialPhiRotated = phi_ph;
					//double photonInitialPhiRotated = 0;

					double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
					double Xidelta = 1.0 - cosXiPrimed;
					//double Xidelta = photonInitialVersinThetaPrimed + photonFinalVersinThetaPrimed - photonInitialVersinThetaPrimed * photonFinalVersinThetaPrimed - photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);
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
					//double photonInitialThetaRotated;
					double photonInitialAlphaRotated;
					//if (1.0 + photonInitialCosThetaPrimed * electronInitialBeta < 1E16) {
					//	photonInitialCosThetaRotated = -1.0;
					//}
					//else {
						//photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1.0 + photonInitialCosThetaPrimed * electronInitialBeta);
					//}
					//LorentzTransformationPhotonReverseZ(electronInitialGamma, photonFinalEnergyPrimed, photonInitialThetaPrimed, photonInitialEnergy, photonInitialThetaRotated);
					double epsilon = 1 + photonInitialCosThetaPrimed;
					photonInitialEnergy = photonInitialEnergyPrimed * electronInitialGamma * (epsilon + electronInitialDelta - epsilon * electronInitialDelta);
					//LorentzTransformationPhotonReverseZalpha(electronInitialGamma, photonFinalEnergyPrimed, photonInitialAlphaPrimed, photonInitialEnergy, photonInitialAlphaRotated);
					//double photonInitialCosThetaRotated = cos(photonInitialThetaRotated);
					double photonInitialVersinThetaRotated = versin(photonInitialThetaRotated);
					//double photonInitialVersinThetaRotated = 2.0 - versin(photonInitialAlphaRotated);
					//double photonInitialCosTheta;
					//double photonInitialPhi;
					//inverseRotationSphericalCoordinates(mu_e, phi_e, photonInitialCosThetaRotated, photonInitialPhiRotated, photonInitialCosTheta, photonInitialPhi);


					double denom = electronInitialDelta + photonFinalVersinThetaRotated - electronInitialDelta * photonFinalVersinThetaRotated;
					double numenator = electronInitialDelta + photonInitialVersinThetaRotated - electronInitialDelta * photonInitialVersinThetaRotated;
					/*if (1.0 - photonFinalCosThetaRotated * electronInitialBeta < 1E-7) {
						denom = electronInitialDelta;
					}
					else {
						denom = (1.0 - photonFinalCosThetaRotated * electronInitialBeta);
					}*/
					double dI = volume * 0.5 * r2 * speed_of_light * electronDist *
						//(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / denom) *
						(numenator / (electronInitialGamma* electronInitialGamma * electronInitialGamma * electronInitialGamma *denom)) *
						//(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
						(2 - 2 * Xidelta + Xidelta * Xidelta + sqr(photonFinalEnergyPrimed / m_c2) * sqr(Xidelta) / (1 - (photonFinalEnergyPrimed / m_c2) * (Xidelta))) *
						photonDistribution->distribution(photonInitialEnergy) *
						photonFinalEnergy * 2 * pi * dphi_ph * my_dcosTheta[imue] * my_dcosTheta[imuph] * delectronEnergy;
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

double InverseComptonEvaluator::evaluateComptonFluxKleinNishinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance) {
	double I = 0;
	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double r2 = sqr(electron_charge * electron_charge / m_c2);

	double photonFinalCosTheta = cos(photonFinalTheta);
	for (int k = 0; k < my_Ne; ++k) {
		double electronInitialEnergy = my_Ee[k];
		//experimental
		double fraction = 0.5 / my_Nmu;
		double thetamin = min(0.1 / (electronInitialEnergy / m_c2), pi / (2 * my_Nmu));
		double dlogtheta = log((pi + (1 - fraction) * thetamin) / (thetamin)) / (my_Nmu - 1);

		my_cosThetaLeft[0] = 1.0;
		my_cosTheta[0] = cos(fraction * thetamin);
		for (int i = 1; i < my_Nmu; ++i) {
			my_cosTheta[i] = cos(thetamin * exp(dlogtheta * (i)) - (1 - fraction) * thetamin);
			my_cosThetaLeft[i] = (my_cosTheta[i] + my_cosTheta[i - 1]) / 2.0;
		}
		for (int i = 0; i < my_Nmu - 1; ++i) {
			my_dcosTheta[i] = -(my_cosThetaLeft[i + 1] - my_cosThetaLeft[i]);
		}
		my_dcosTheta[my_Nmu - 1] = 1.0 + my_cosThetaLeft[my_Nmu - 1];
		double delectronEnergy;
		double electronInitialGamma = electronInitialEnergy / (m_c2);
		double electronInitialBeta = sqrt(1.0 - 1.0 / (electronInitialGamma * electronInitialGamma));
		if (k == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[k] - my_Ee[k - 1];
		}
		//if (photonFinalEnergy > electronInitialEnergy) {
		//	printf("aaa\n");
		//}

		for (int imue = 0; imue < my_Nmu; ++imue) {
			double mu_e = my_cosTheta[imue];
			double theta_e = my_theta[imue];
			for (int iphie = 0; iphie < my_Nphi; ++iphie) {
				double phi_e = 2 * pi * (iphie + 0.0) / my_Nphi;
				double electronDist = electronDistribution->distribution(electronInitialEnergy, mu_e, phi_e);
				double dphi_e = 2*pi / my_Nphi;
				//rotation
				double photonFinalThetaRotated;
				double photonFinalPhiRotated;
				rotationSphericalCoordinates(theta_e, phi_e, photonFinalTheta, photonFinalPhi, photonFinalThetaRotated, photonFinalPhiRotated);
				double photonFinalCosThetaRotated = cos(photonFinalThetaRotated);
				double photonFinalEnergyPrimed;
				double photonFinalThetaPrimed;
				LorentzTransformationPhotonZ(electronInitialGamma, photonFinalEnergy, photonFinalThetaRotated, photonFinalEnergyPrimed, photonFinalThetaPrimed);
				//double photonFinalEnergyPrimed = electronInitialGamma * (1 - photonFinalCosThetaRotated * electronInitialBeta) * photonFinalEnergy;
				//double photonFinalCosThetaPrimed = (photonFinalCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonFinalCosThetaRotated);
				double photonFinalCosThetaPrimed = cos(photonFinalThetaPrimed);
				double photonFinalSinThetaPrimed = sin(photonFinalThetaPrimed);
				for (int imuph = 0; imuph < my_Nmu; ++imuph) {
					double mu_ph = -my_cosTheta[imuph];
					for (int iphiph = 0; iphiph < my_Nphi; ++iphiph) {
						double phi_ph = 2 * pi * (iphiph + 0.0) / my_Nphi;
						double dphi_ph = 2*pi / my_Nphi;
						double photonInitialCosThetaPrimed = mu_ph;
						double photonInitialPhiRotated = phi_ph;
						//inverseRotationSphericalCoordinates(mu_e, phi_e, mu_ph, phi_ph, photonInitialCosThetaRotated, photonInitialPhiRotated);
						//double photonInitialCosThetaPrimed = (photonInitialCosThetaRotated - electronInitialBeta) / (1.0 - electronInitialBeta * photonInitialCosThetaRotated);

						double photonInitialSinThetaPrimed = sqrt(1.0 - photonInitialCosThetaPrimed * photonInitialCosThetaPrimed);

						double cosXiPrimed = photonInitialCosThetaPrimed * photonFinalCosThetaPrimed + photonInitialSinThetaPrimed * photonFinalSinThetaPrimed * cos(photonFinalPhiRotated - photonInitialPhiRotated);

						double photonInitialEnergyPrimed = photonFinalEnergyPrimed / (1.0 - (photonFinalEnergyPrimed / (m_c2)) * (1.0 - cosXiPrimed));
						if (photonInitialEnergyPrimed <= 0) {
							continue;
						}

						double photonInitialEnergy;
						double photonInitialCosThetaRotated;
						LorentzTransformationPhotonReverseZ(electronInitialGamma, photonInitialEnergyPrimed, photonInitialCosThetaPrimed, photonInitialEnergy, photonInitialCosThetaRotated);
						//double photonInitialEnergy = electronInitialGamma * photonInitialEnergyPrimed + electronInitialBeta * electronInitialGamma * photonInitialEnergyPrimed * photonInitialCosThetaPrimed;
						//double photonInitialCosThetaRotated = (photonInitialCosThetaPrimed + electronInitialBeta) / (1 + photonInitialCosThetaPrimed * electronInitialBeta);
						double photonInitialThetaRotated = acos(photonInitialCosThetaRotated);

						double photonInitialTheta;
						double photonInitialPhi;
						inverseRotationSphericalCoordinates(theta_e, phi_e, photonInitialThetaRotated, photonInitialPhiRotated, photonInitialTheta, photonInitialPhi);
						double photonInitialCosTheta = cos(photonInitialTheta);

						double dI = electronDist * volume * 0.5 * r2 * speed_of_light *
							(sqr(1 - photonInitialCosThetaRotated * electronInitialBeta) / (1.0 - photonFinalCosThetaRotated * electronInitialBeta)) *
							(1 + cosXiPrimed * cosXiPrimed + sqr(photonFinalEnergyPrimed / m_c2) * sqr(1 - cosXiPrimed) / (1 - (photonFinalEnergyPrimed / m_c2) * (1 - cosXiPrimed))) *
							photonFinalEnergy * photonDistribution->distribution(photonInitialEnergy, photonInitialCosTheta, photonInitialPhi) *
							dphi_e * dphi_ph * my_dcosTheta[imue] * my_dcosTheta[imuph] * delectronEnergy;

						if (dI < 0) {
							printf("dI <  0\n");
							printLog("dI < 0\n");
							exit(0);
						}

						I += dI;
						if (I != I) {
							printf("I = NaN\n");
							printLog("I = NaN\n");
							exit(0);
						}
					}
				}
			}
		}
	}
	return I / sqr(distance);
}

double InverseComptonEvaluator::evaluateComptonFluxKleinNishinaAnisotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	return evaluateComptonFluxKleinNishinaAnisotropic(photonFinalEnergy, 0, 0, photonDistribution, electronDistribution, volume, distance);
}

double InverseComptonEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) {
	if (my_solverType == ComptonSolverType::ISOTROPIC_THOMSON) {
		return evaluateComptonFluxThomsonIsotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else if (my_solverType == ComptonSolverType::ISOTROPIC_JONES) {
		return evaluateComptonFluxJonesIsotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else if (my_solverType == ComptonSolverType::ISOTROPIC_KLEIN_NISHINA) {
		return evaluateComptonFluxKleinNishinaIsotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else if (my_solverType == ComptonSolverType::ISOTROPIC_KLEIN_NISHINA1) {
		return evaluateComptonFluxKleinNishinaIsotropic1(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
	}
	else if (my_solverType == ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA) {
		return evaluateComptonFluxKleinNishinaAnisotropic(photonFinalEnergy, my_photonDistribution, electronDistribution, volume, distance);
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

double InverseComptonEvaluator::evaluateFluxFromSourceAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, RadiationSource* source) {
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				result += evaluateComptonFluxKleinNishinaAnisotropic(photonFinalEnergy, photonFinalTheta, photonFinalPhi, photonDistribution, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	return result;
}
