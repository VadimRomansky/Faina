#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "massiveParticleDistribution.h"
#include "coordinateTransform.h"

#include "synchrotron.h"

double evaluateMcDonaldIntegral(const double& nu) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		//printf("x < UvarovX[0]\n");
		return 0;
	}
	if (nu > UvarovX[Napprox - 1]) {
		//printf("x > UvarovX[Napprox - 1]\n");
		return sqrt(nu * pi / 2) * exp(-nu);
	}
	int leftIndex = 0;
	int rightIndex = Napprox - 1;
	while (rightIndex - leftIndex > 1) {
		int currentIndex = (rightIndex + leftIndex) / 2;
		if (UvarovX[currentIndex] > nu) {
			rightIndex = currentIndex;
		}
		else {
			leftIndex = currentIndex;
		}
	}

	//double result = (UvarovValue[rightIndex]*(nu - UvarovX[leftIndex]) + UvarovValue[leftIndex]*(UvarovX[rightIndex] - nu))/(UvarovX[rightIndex] - UvarovX[leftIndex]);
	double result = UvarovValue[leftIndex] * exp(log(UvarovValue[rightIndex] / UvarovValue[leftIndex]) * ((nu - UvarovX[leftIndex]) / (UvarovX[rightIndex] - UvarovX[leftIndex])));
	if (result < 0) {
		printf("evaluateMcDonaldIntegral\n");
		printLog("evaluateMcDonaldIntegral\n");
		printf("result < 0\n");
		printLog("result < 0\n");
	}
	return result;
}

double evaluateMcDonaldFunction5_3(const double& nu) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		return McDonaldValue5_3[0];
	}
	if (nu > UvarovX[Napprox - 1]) {
		return 0;
	}
	int leftIndex = 0;
	int rightIndex = Napprox - 1;
	while (rightIndex - leftIndex > 1) {
		int currentIndex = (rightIndex + leftIndex) / 2;
		if (UvarovX[currentIndex] > nu) {
			rightIndex = currentIndex;
		}
		else {
			leftIndex = currentIndex;
		}
	}

	double result = (McDonaldValue5_3[rightIndex] * (nu - UvarovX[leftIndex]) + McDonaldValue5_3[leftIndex] * (UvarovX[rightIndex] - nu)) / (UvarovX[rightIndex] - UvarovX[leftIndex]);
	//double result = McDonaldValue[curIndex - 1] * exp(
	//	log(McDonaldValue[curIndex] / McDonaldValue[curIndex - 1]) * ((nu - UvarovX[curIndex - 1]) / (UvarovX[curIndex] - UvarovX[curIndex - 1])));
	if (result < 0) {
		printf("evaluateMcDonaldFunction5_3\n");
		printLog("evaluateMcDonaldFunction5_3\n");
		printf("result < 0\n");
		printLog("result < 0\n");
	}
	return result;
}

void evaluateMcDonaldFunctions(const double& nu, double& K1_3, double& K2_3, double& K4_3, double& K5_3) {
	int curIndex = 0;
	if (nu < UvarovX[0]) {
		K1_3 = McDonaldValue1_3[0];
		K2_3 = McDonaldValue2_3[0];
		K4_3 = McDonaldValue4_3[0];
		K5_3 = McDonaldValue5_3[0];
		return;
	}
	if (nu > UvarovX[Napprox - 1]) {
		K1_3 = 0;
		K2_3 = 0;
		K4_3 = 0;
		K5_3 = 0;
		return;
	}
	int leftIndex = 0;
	int rightIndex = Napprox - 1;
	while (rightIndex - leftIndex > 1) {
		int currentIndex = (rightIndex + leftIndex) / 2;
		if (UvarovX[currentIndex] > nu) {
			rightIndex = currentIndex;
		}
		else {
			leftIndex = currentIndex;
		}
	}

	K1_3 = (McDonaldValue1_3[rightIndex] * (nu - UvarovX[leftIndex]) + McDonaldValue1_3[leftIndex] * (UvarovX[rightIndex] - nu)) / (UvarovX[rightIndex] - UvarovX[leftIndex]);
	K2_3 = (McDonaldValue2_3[rightIndex] * (nu - UvarovX[leftIndex]) + McDonaldValue2_3[leftIndex] * (UvarovX[rightIndex] - nu)) / (UvarovX[rightIndex] - UvarovX[leftIndex]);
	K4_3 = (McDonaldValue4_3[rightIndex] * (nu - UvarovX[leftIndex]) + McDonaldValue4_3[leftIndex] * (UvarovX[rightIndex] - nu)) / (UvarovX[rightIndex] - UvarovX[leftIndex]);
	K5_3 = (McDonaldValue5_3[rightIndex] * (nu - UvarovX[leftIndex]) + McDonaldValue5_3[leftIndex] * (UvarovX[rightIndex] - nu)) / (UvarovX[rightIndex] - UvarovX[leftIndex]);
}

double criticalNu(const double& E, const double& sinhi, const double& H, const double& coef) {
	return coef * H * sinhi * E * E;
}

SynchrotronEvaluator::SynchrotronEvaluator(int Ne, double Emin, double Emax, bool absorption, bool doppler):RadiationEvaluator(Ne, Emin, Emax, absorption, doppler)
{
}

SynchrotronEvaluator::~SynchrotronEvaluator()
{
}

void SynchrotronEvaluator::resetParameters(const double *parameters, const double *normalizationUnits){
    //do nothing?
}

void SynchrotronEvaluator::evaluateEmissivityAndAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source, double& I, double& A)
{
	double photonFinalFrequency = photonFinalEnergy / hplank;
	double B = source->getB(ix1, iz, ix2);
	if (B != B) {
		printf("B = NaN\n");
		printLog("B = NaN\n");
		exit(0);
	}
	if (B == 0) {
		I = 0;
		A = 0;
		return;
	}
	double sinTheta = source->getSinTheta(ix1, iz, ix2);
	double concentration = source->getConcentration(ix1, iz, ix2);

	double v, vtheta, vphi;
	source->getVelocity(ix1, iz, ix2, v, vtheta, vphi);
	double beta = v / speed_of_light;
	double gamma = 1.0 / sqrt(1.0 - beta * beta);
	double mu = cos(vtheta);
	double D = gamma * (1 - mu * beta);

	if (my_doppler) {
		double Btheta = source->getBTheta(ix1, iz, ix2);
		double Bphi = source->getBPhi(ix1, iz, ix2);
		double Bpar = B * (cos(Btheta) * cos(vtheta) + sin(Btheta) * sin(vtheta) * (cos(Bphi) * cos(vphi) + sin(Bphi) * sin(vphi)));
		//from invarian B^2 - E^2
		double Bprimed = sqrt(B*B / (gamma*gamma) + Bpar*Bpar*beta*beta);
		if (Bprimed == 0) {
			I = 0;
			A = 0;
			return;
		}
		double Bnorm = sqrt(Bprimed * Bprimed - Bpar * Bpar);
		double concentrationPrimed = concentration / gamma;

		double thetar;
		double phir;
		rotationSphericalCoordinates(vtheta, vphi, 0, 0, thetar, phir);
		double thetaBr;
		double phiBr;
		rotationSphericalCoordinates(vtheta, vphi, Btheta, Bphi, thetaBr, phiBr);
		double photonFinalEnergyPrimed;
		double thetaPrimed;
		LorentzTransformationPhotonZ(gamma, photonFinalEnergy, thetar, photonFinalEnergyPrimed, thetaPrimed);
		double cosbtheta = Bpar / Bprimed;
		double BthetaPrimed;
		if (cosbtheta > 1.0){
			if (cosbtheta < 1.000000001) {
				BthetaPrimed = 0.0;
			}
			else {
				printf("cos B theta primed %g > 1\n", cosbtheta);
				printLog("cos B theta primed %g > 1\n", cosbtheta);
				exit(0);
			}
		}
		else if (cosbtheta < -1.0) {
			if (cosbtheta > -1.000000001) {
				BthetaPrimed = pi;
			}
			else {
				printf("cos B theta primed %g < -1\n", cosbtheta);
				printLog("cos B theta primed %g < -1\n", cosbtheta);
				exit(0);
			}
		}
		else {
			BthetaPrimed = acos(cosbtheta);
		}
		if (thetaBr > pi / 2) {
			BthetaPrimed = pi - acos(Bpar/Bprimed);
		}
		if (BthetaPrimed != BthetaPrimed) {
			printf("BthetaPrimed = NaN\n");
			printLog("BthetaPrimed = NaN\n");
			printLog("B = %g, Btheta = %g, Bphi = %g, Bpar = %g, Bprimed = %g, gamma = %g, vtheta = %g, vphi = %g, thetaBr = %g, phiBr = %g\n", B, Btheta, Bphi, Bpar, Bprimed, gamma, vtheta, vphi, thetaBr, phiBr);
			exit(0);
		}

		double Bpar1 = Bprimed * (cos(BthetaPrimed) * cos(thetaPrimed) + sin(BthetaPrimed) * sin(thetaPrimed) * (cos(phiBr) * cos(phir) + sin(phiBr) * sin(phir)));
		double Bnorm1 = sqrt(Bprimed * Bprimed - Bpar1 * Bpar1);

		photonFinalFrequency = photonFinalEnergyPrimed / hplank;
		B = Bprimed;
		if (Bprimed > 0) {
			sinTheta = Bnorm1 / Bprimed;
		}
		else {
			sinTheta = 0;
		}
		if (sinTheta != sinTheta) {
			printf("sinTheta = NaN\n");
			printLog("sinTheta = NaN\n");
			exit(0);
		}
		concentration = concentrationPrimed;
	}
	MassiveParticleDistribution* distribution = source->getParticleDistribution(ix1, iz, ix2);
	MassiveParticleIsotropicDistribution* electronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(distribution);
	if (electronDistribution == NULL) {
		printf("Synchrotron evaluator works only with isotropic distributions\n");
		printLog("Synchrotron evaluator works only with isotropic distributions\n");
		exit(0);
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
		printf("Emin > Emax in synchrotron evaluator\n");
		printLog("Emin > Emax in synchrotron evaluator\n");
		I = 0;
		A = 0;
		return;
	}

	if ((B == 0) || (sinTheta == 0)) {
		I = 0;
		A = 0;
		return;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	double* Ee = new double[my_Ne];
	Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		Ee[i] = Ee[i - 1] * factor;
	}
	//Anu from ghiselini simple
	I = 0;
	A = 0;

	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / m_c2;
	double criticalNuCoef = 3 * electron_charge / (4 * pi * m * m * m * speed_of_light * speed_of_light4);


	if (sinTheta == 0.0) {
		return;
	}

	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinTheta * sinTheta);

	for (int j = 0; j < my_Ne; ++j) {
		double delectronEnergy = 0;
		if (j == 0) {
			delectronEnergy = Ee[j + 1] - Ee[j];
			if (delectronEnergy < 0) {
				printf("delectronEnergy %g < 0\n", delectronEnergy);
				printLog("delectronEnergy %g < 0\n", delectronEnergy);
				printf("j = %d, Ee[j+1] = %g, Ee[j] = %g\n", j, Ee[j+1], Ee[j]);
				printLog("j = %d, Ee[j+1] = %g, Ee[j] = %g\n", j, Ee[j+1], Ee[j]);
				exit(0);
			}
		}
		else {
			delectronEnergy = Ee[j] - Ee[j - 1];
			if (delectronEnergy < 0) {
				printf("delectronEnergy %g < 0\n", delectronEnergy);
				printLog("delectronEnergy %g < 0\n", delectronEnergy);
				printf("j = %d, Ee[j] = %g, Ee[j-1] = %g\n", j, Ee[j], Ee[j - 1]);
				printLog("j = %d, Ee[j] = %g, Ee[j-1] = %g\n", j, Ee[j], Ee[j - 1]);
				exit(0);
			}
		}
		double electronDist = concentration*electronDistribution->distributionNormalized(Ee[j]);
		if (electronDist < 0) {
			printf("Fe = %g < 0\n", electronDist);
			printLog("Fe = %g < 0\n", electronDist);
			exit(0);
		}
		if (electronDist > 0) {
			double dFe = electronDist * delectronEnergy;
			double nuc = criticalNu(Ee[j], sinTheta, B, criticalNuCoef);
			double gamma = Ee[j] / m_c2;

			double mcDonaldIntegral = evaluateMcDonaldIntegral(photonFinalFrequency / nuc);

            //todo h?
            I = I + emissivityCoef * dFe * B * sinTheta * mcDonaldIntegral/hplank;
			if (I < 0) {
				printf("evaluateSynchrotronIandA\n");
				printLog("evaluateSynchrotronIandA\n");
				printf("I < 0\n");
				printLog("I < 0\n");
				printf("dFe[j] = %g\n", dFe);
				printLog("dFe[j] = %g\n", dFe);
				printf("B = %g\n", B);
				printLog("B = %g\n", B);
				printf("sinhi = %g\n", sinTheta);
				printLog("sinhi = %g\n", sinTheta);
				printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				exit(0);
			}
			if (my_absorption) {
				////todo! F(g +dg)
				double tempP = gamma * gamma * emissivityCoef * B * sinTheta * mcDonaldIntegral;
				double dg = 0.1 * delectronEnergy / (m_c2);

				double tempGamma = gamma + dg;
				double tempnuc = criticalNu(m_c2 * tempGamma, sinTheta, B, criticalNuCoef);
				double tempP2 = tempGamma * tempGamma * emissivityCoef * B * sinTheta * evaluateMcDonaldIntegral(photonFinalFrequency / tempnuc);
				double Pder = (tempP2 - tempP) / dg;



				A = A + (1.0 / (2 * m * photonFinalFrequency * photonFinalFrequency)) * dFe * Pder / (gamma * gamma);
				if (A != A) {
					printf("evaluateSynchrotronIandA\n");
					printLog("evaluateSynchrotronIandA\n");
					printf("A = NaN\n");
					printLog("A = NaN\n");
					printf("dFe[j] = %g\n", dFe);
					printLog("dFe[j] = %g\n", dFe);
					printf("B = %g\n", B);
					printLog("B = %g\n", B);
					printf("sinhi = %g\n", sinTheta);
					printLog("sinhi = %g\n", sinTheta);
					printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
					printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
					printf("Pder = %g\n", Pder);
					printLog("Pder = %g\n", Pder);
					exit(0);
				}
			}
			if (I != I) {
				printf("evaluateSynchrotronIandA\n");
				printLog("evaluateSynchrotronIandA\n");
				printf("I = NaN\n");
				printLog("I = NaN\n");
				printf("dFe[j] = %g\n", dFe);
				printLog("dFe[j] = %g\n", dFe);
				printf("B = %g\n", B);
				printLog("B = %g\n", B);
				printf("sinhi = %g\n", sinTheta);
				printLog("sinhi = %g\n", sinTheta);
				printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				exit(0);
			}
		}
    }

	delete[] Ee;

	if (my_doppler) {
		I = I / (D * D);
		A = A * D;
	}
}

double SynchrotronEvaluator::evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	double photonFinalFrequency = photonFinalEnergy / hplank;
	double B = source->getB(ix1, iz, ix2);
	double sinTheta = source->getSinTheta(ix1, iz, ix2);
	double concentration = source->getConcentration(ix1, iz, ix2);
	MassiveParticleIsotropicDistribution* electronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(ix1, iz, ix2));
	if (electronDistribution == NULL) {
		printf("Synchrotron evaluator works only with isotropic distributions\n");
		printLog("Synchrotron evaluator works only with isotropic distributions\n");
		exit(0);
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
		printf("Emin > Emax in synchrotron evaluator\n");
		printLog("Emin > Emax in synchrotron evaluator\n");
		return 0;
	}

	if ((B == 0) || (sinTheta == 0)) {
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}
	//Anu from ghiselini simple

	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / m_c2;
	double criticalNuCoef = 3 * electron_charge / (4 * pi * m * m * m * speed_of_light * speed_of_light4);

	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double cosTheta = sqrt(1.0 - sinTheta * sinTheta);

	double I = 0;

	for (int j = 0; j < my_Ne; ++j) {
		double delectronEnergy = 0;
		if (j == 0) {
			delectronEnergy = my_Ee[j + 1] - my_Ee[j];
		}
		else {
			delectronEnergy = my_Ee[j] - my_Ee[j - 1];
		}
		double electronDist = concentration * electronDistribution->distributionNormalized(my_Ee[j]);
		if (electronDist > 0) {
			double dFe = electronDist * delectronEnergy;
			double nuc = criticalNu(my_Ee[j], sinTheta, B, criticalNuCoef);
			double gamma = my_Ee[j] / m_c2;

			double mcDonaldIntegral = evaluateMcDonaldIntegral(photonFinalFrequency / nuc);

			//todo h?
			I = I + emissivityCoef * dFe * B * sinTheta * mcDonaldIntegral / hplank;
			if (I < 0) {
				printf("evaluateSynchrotronIandA\n");
				printLog("evaluateSynchrotronIandA\n");
				printf("I < 0\n");
				printLog("I < 0\n");
				printf("dFe[j] = %g\n", dFe);
				printLog("dFe[j] = %g\n", dFe);
				printf("B = %g\n", B);
				printLog("B = %g\n", B);
				printf("sinhi = %g\n", sinTheta);
				printLog("sinhi = %g\n", sinTheta);
				printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				exit(0);
			}

			if (I != I) {
				printf("evaluateSynchrotronIandA\n");
				printLog("evaluateSynchrotronIandA\n");
				printf("I = NaN\n");
				printLog("I = NaN\n");
				printf("dFe[j] = %g\n", dFe);
				printLog("dFe[j] = %g\n", dFe);
				printf("B = %g\n", B);
				printLog("B = %g\n", B);
				printf("sinhi = %g\n", sinTheta);
				printLog("sinhi = %g\n", sinTheta);
				printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				exit(0);
			}
		}
	}

	return I;
}

double SynchrotronEvaluator::evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	double photonFinalFrequency = photonFinalEnergy / hplank;

	double B = source->getB(ix1, iz, ix2);
	double sinTheta = source->getSinTheta(ix1, iz, ix2);
	double concentration = source->getConcentration(ix1, iz, ix2);
	MassiveParticleIsotropicDistribution* electronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(ix1, iz, ix2));
	if (electronDistribution == NULL) {
		printf("Synchrotron evaluator works only with isotropic distributions\n");
		printLog("Synchrotron evaluator works only with isotropic distributions\n");
		exit(0);
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
		printf("Emin > Emax in synchrotron evaluator\n");
		printLog("Emin > Emax in synchrotron evaluator\n");
		return 0;
	}

	if ((B == 0) || (sinTheta == 0)) {
		return 0;
	}

	double factor = pow(Emax / Emin, 1.0 / (my_Ne - 1));

	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}
	//Anu from ghiselini simple


	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / m_c2;
	double criticalNuCoef = 3 * electron_charge / (4 * pi * m * m * m * speed_of_light * speed_of_light4);


	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double cosTheta = sqrt(1.0 - sinTheta * sinTheta);

	double A = 0;

	for (int j = 0; j < my_Ne; ++j) {
		double delectronEnergy = 0;
		if (j == 0) {
			delectronEnergy = my_Ee[j + 1] - my_Ee[j];
		}
		else {
			delectronEnergy = my_Ee[j] - my_Ee[j - 1];
		}
		double electronDist = concentration * electronDistribution->distributionNormalized(my_Ee[j]);
		if (electronDist > 0) {
			double dFe = electronDist * delectronEnergy;
			double nuc = criticalNu(my_Ee[j], sinTheta, B, criticalNuCoef);
			double gamma = my_Ee[j] / m_c2;

			double mcDonaldIntegral = evaluateMcDonaldIntegral(photonFinalFrequency / nuc);

			if (my_absorption) {
				////todo! F(g +dg)
				double tempP = gamma * gamma * emissivityCoef * B * sinTheta * mcDonaldIntegral;
				double dg = 0.1 * delectronEnergy / (m_c2);

				double tempGamma = gamma + dg;
				double tempnuc = criticalNu(m_c2 * tempGamma, sinTheta, B, criticalNuCoef);
				double tempP2 = tempGamma * tempGamma * emissivityCoef * B * sinTheta * evaluateMcDonaldIntegral(photonFinalFrequency / tempnuc);
				double Pder = (tempP2 - tempP) / dg;



				A = A + (1.0 / (2 * m * photonFinalFrequency * photonFinalFrequency)) * dFe * Pder / (gamma * gamma);
				if (A != A) {
					printf("evaluateSynchrotronIandA\n");
					printLog("evaluateSynchrotronIandA\n");
					printf("A = NaN\n");
					printLog("A = NaN\n");
					printf("dFe[j] = %g\n", dFe);
					printLog("dFe[j] = %g\n", dFe);
					printf("B = %g\n", B);
					printLog("B = %g\n", B);
					printf("sinhi = %g\n", sinTheta);
					printLog("sinhi = %g\n", sinTheta);
					printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
					printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
					printf("Pder = %g\n", Pder);
					printLog("Pder = %g\n", Pder);
					exit(0);
				}
			}
		}
	}

	return A;
}