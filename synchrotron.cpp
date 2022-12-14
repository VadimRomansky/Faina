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

void evaluateSynchrotronIandA(const double& photonFinalFrequency, const double& photonFinalTheta, const double& photonFinalPhi, const double& B, const double& sinhi, const double& concentration, ElectronIsotropicDistribution* electronDistribution, const double& Emin, const double& Emax, const int Ne, double& I, double& A) {
	double factor = pow(Emax / Emin, 1.0 / (Ne - 1));

	double* Ee = new double[Ne];
	Ee[0] = Emin;
	for (int i = 1; i < Ne; ++i) {
		Ee[i] = Ee[i - 1] * factor;
	}
	
	//Anu from ghiselini simple
	I = 0;
	A = 0;

	if (sinhi == 0.0) {
		return;
	}

	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi * sinhi);


	double oldA = 0;
	for (int j = 1; j < Ne; ++j) {
		double delectronEnergy = Ee[j] - Ee[j - 1];
		double electronDist = concentration*electronDistribution->distributionNormalized(Ee[j]);
		if (electronDist > 0) {
			double dFe = electronDist * delectronEnergy;
			double nuc = criticalNu(Ee[j], sinhi, B);
			double gamma = Ee[j] / (massElectron * speed_of_light2);

			double mcDonaldIntegral = evaluateMcDonaldIntegral(photonFinalFrequency / nuc);

			I = I + emissivityCoef * dFe * B * sinhi * mcDonaldIntegral;
			if (I < 0) {
				printf("I < 0\n");
				printf("dFe[j] = %g\n", dFe);
				exit(0);
			}
			////todo! F(g +dg)
			double tempP = gamma * gamma * emissivityCoef * B * sinhi * mcDonaldIntegral;
			double dg = 0.1 * (Ee[j] - Ee[j - 1]) / (me_c2);

			double tempGamma = gamma + dg;
			double tempnuc = criticalNu(massElectron * speed_of_light2 * tempGamma, sinhi, B);
			double tempP2 = tempGamma * tempGamma * emissivityCoef * B * sinhi * evaluateMcDonaldIntegral(photonFinalFrequency / tempnuc);
			double Pder = (tempP2 - tempP) / dg;



			A = A + (1.0 / (2 * massElectron * photonFinalFrequency * photonFinalFrequency)) * dFe * Pder / (gamma * gamma);

			if (I != I) {
				printf("I NaN\n");
				exit(0);
			}
			if (A != A) {
				printf("A Nan\n");
				exit(0);
			}
		}
	}

	delete[] Ee;
}

double evaluateSynchrotronFluxFromSource(RadiationSource* source, const double& photonFinalFrequency, const double& Emin, const double& Emax, const int Ne) {
	int Nrho = source->getNrho();
	double maxRho = source->getMaxRho();
	int Nz = source->getNz();
	double minZ = source->getMinZ();
	double maxZ = source->getMaxZ();
	int Nphi = source->getNphi();

	double drho = maxRho / Nrho;
	double dz = (maxZ - minZ) / Nz;
	double dphi = 2 * pi / Nphi;

	double result = 0;

	for (int irho = 0; irho < Nrho; ++irho) {
		double area = source->getArea(irho);
		for (int iphi = 0; iphi < Nphi; ++iphi) {
			double localI = 0;
			for (int iz = 0; iz < Nz; ++iz) {
				double A = 0;
				double I = 0;
				evaluateSynchrotronIandA(photonFinalFrequency, 0, 0, source->getB(irho, iz, iphi), source->getSinTheta(irho, iz, iphi), source->getConcentration(irho, iz, iphi), source->getElectronDistribution(irho, iz, iphi), Emin, Emax, Ne, I, A);
				double length = source->getLength(irho, iz, iphi);
				if (length > 0) {
					double I0 = localI;
					double Q = I * area;
					double tau = A * length;
					double S = 0;
					if (Q > 0) {
						S = Q / A;
					}
					if (fabs(tau) < 1E-15) {
						localI = I0 * (1.0 - tau) + S * tau;
					}
					else {
						localI = S + (I0 - S) * exp(-tau);
					}
				}
			}
			result += localI;
		}
	}

	return result/sqr(source->getDistance());
}

SynchrotronEvaluator::SynchrotronEvaluator(int Ne, double Emin, double Emax)
{
	double factor = pow(Emax / Emin, 1.0 / (Ne - 1));
	my_Ne = Ne;
	my_Emin = Emin;
	my_Emax = Emax;

	my_Ee = new double[my_Ne];
	my_Ee[0] = Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}
}

SynchrotronEvaluator::~SynchrotronEvaluator()
{
	delete[] my_Ee;
}

void SynchrotronEvaluator::evaluateSynchrotronIandA(const double& photonFinalFrequency, const double& photonFinalTheta, const double& photonFinalPhi, const double& B, const double& sinhi, const double& concentration, ElectronIsotropicDistribution* electronDistribution, double& I, double& A)
{
	//Anu from ghiselini simple
	I = 0;
	A = 0;

	if (sinhi == 0.0) {
		return;
	}

	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi * sinhi);


	double oldA = 0;
	for (int j = 0; j < my_Ne; ++j) {
		double delectronEnergy = 0;
		if (j == 0) {
			delectronEnergy = my_Ee[j + 1] - my_Ee[j];
		}
		else {
			delectronEnergy = my_Ee[j] - my_Ee[j - 1];
		}
		double electronDist = concentration*electronDistribution->distributionNormalized(my_Ee[j]);
		if (electronDist > 0) {
			double dFe = electronDist * delectronEnergy;
			double nuc = criticalNu(my_Ee[j], sinhi, B);
			double gamma = my_Ee[j] / (massElectron * speed_of_light2);

			double mcDonaldIntegral = evaluateMcDonaldIntegral(photonFinalFrequency / nuc);

			I = I + emissivityCoef * dFe * B * sinhi * mcDonaldIntegral;
			if (I < 0) {
				printf("I < 0\n");
				printf("dFe[j] = %g\n", dFe);
				exit(0);
			}
			////todo! F(g +dg)
			double tempP = gamma * gamma * emissivityCoef * B * sinhi * mcDonaldIntegral;
			double dg = 0.1 * delectronEnergy / (me_c2);

			double tempGamma = gamma + dg;
			double tempnuc = criticalNu(massElectron * speed_of_light2 * tempGamma, sinhi, B);
			double tempP2 = tempGamma * tempGamma * emissivityCoef * B * sinhi * evaluateMcDonaldIntegral(photonFinalFrequency / tempnuc);
			double Pder = (tempP2 - tempP) / dg;



			A = A + (1.0 / (2 * massElectron * photonFinalFrequency * photonFinalFrequency)) * dFe * Pder / (gamma * gamma);

			if (I != I) {
				printf("I NaN\n");
				exit(0);
			}
			if (A != A) {
				printf("A Nan\n");
				exit(0);
			}
		}
	}
}

double SynchrotronEvaluator::evaluateSynchrotronFluxFromSource(RadiationSource* source, const double& photonFinalFrequency)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();


	double result = 0;

	for (int irho = 0; irho < Nrho; ++irho) {
		double area = source->getArea(irho);
		for (int iphi = 0; iphi < Nphi; ++iphi) {
			double localI = 0;
			for (int iz = 0; iz < Nz; ++iz) {
				double A = 0;
				double I = 0;
				evaluateSynchrotronIandA(photonFinalFrequency, 0, 0, source->getB(irho, iz, iphi), source->getSinTheta(irho, iz, iphi), source->getConcentration(irho, iz, iphi), source->getElectronDistribution(irho, iz, iphi), I, A);
				double length = source->getLength(irho, iz, iphi);
				if (length > 0) {
					double I0 = localI;
					double Q = I * area;
					double tau = A * length;
					double S = 0;
					if (Q > 0) {
						S = Q / A;
					}
					if (fabs(tau) < 1E-15) {
						localI = I0 * (1.0 - tau) + S * tau;
					}
					else {
						localI = S + (I0 - S) * exp(-tau);
					}
				}
			}
			result += localI;
		}
	}

	return result / sqr(source->getDistance());
}
