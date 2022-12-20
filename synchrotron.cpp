#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "electronDistribution.h"

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
		printf("result < 0\n");
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
		printf("result < 0\n");
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

double criticalNu(const double& E, const double& sinhi, const double& H) {
	return criticalNuCoef * H * sinhi * E * E;
}

double evaluateSynchrotron(const double& photonFinalFrequency, const double& photonFinalTheta, const double& photonFinalPhi, const double& B, const double& sinhi, ElectronDistribution* electronDistribution, const double& volume, const double& distance, const double& Emin, const double& Emax, const int Ne, const int Nmu, const int Nphi) {
	double factor = pow(Emax / Emin, 1.0 / (Ne - 1));

	double* Ee = new double[Ne];
	Ee[0] = Emin;
	for (int i = 1; i < Ne; ++i) {
		Ee[i] = Ee[i - 1] * factor;
	}
	
	//Anu from ghiselini simple
	double Inu = 0;
	double Anu = 0;

	if (sinhi == 0.0) {
		return;
	}

	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi * sinhi);


	double oldA = 0;
	for (int j = 1; j < Ne; ++j) {
		double delectronEnergy = Ee[j] - Ee[j - 1];
		double electronDist = electronDistribution->distribution(Ee[j], cos(photonFinalTheta), photonFinalPhi);
		if (electronDist > 0) {
			double dFe = electronDist * delectronEnergy;
			double nuc = criticalNu(Ee[j], sinhi, B);
			double gamma = Ee[j] / (massElectron * speed_of_light2);
			double gamma4 = gamma * gamma * gamma * gamma;
			double sigmaCoef = dFe / (B * gamma4);

			Inu = Inu + emissivityCoef * dFe * B * sinhi * evaluateMcDonaldIntegral(photonFinalFrequency / nuc);
			if (Inu < 0) {
				printf("Inu < 0\n");
				printf("dFe[j] = %g\n", dFe);
				exit(0);
			}
			////todo! F(g +dg)
			double tempP = gamma * gamma * emissivityCoef * B * dFe * sinhi * evaluateMcDonaldIntegral(photonFinalFrequency / nuc);
			double dg = 0.1 * (Ee[j] - Ee[j - 1]) / (massElectron * speed_of_light2);

			double tempGamma = gamma + dg;
			double tempnuc = criticalNu(massElectron * speed_of_light2 * tempGamma, sinhi, B);
			double nextdFe = evaluateNextdFe(Ee, dFe, dg, j, Np);
			if (nextdFe < 0) {
				printf("nextdFe < 0\n");
			}
			double tempP2 = tempGamma * tempGamma * emissivityCoef * B * nextdFe * sinhi * evaluateMcDonaldIntegral(nu / tempnuc);
			double Pder = (tempP2 - tempP) / dg;



			Anu = Anu + (1.0 / (2 * massElectron * nu * nu)) * Pder / (gamma * gamma);

			if (Inu != Inu) {
				printf("Inu NaN\n");
				exit(0);
			}
			if (Anu != Anu) {
				printf("Anu Nan\n");
				exit(0);
			}
		}
	}
}