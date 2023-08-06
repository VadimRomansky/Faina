#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "massiveParticleDistribution.h"

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

SynchrotronEvaluator::SynchrotronEvaluator(int Ne, double Emin, double Emax, bool selfAbsorption, bool doppler, const double& defaultB, const double& defaultSinTheta, const double& defaultLength):RadiationEvaluator(Ne, Emin, Emax)
{
	my_selfAbsorption = selfAbsorption;
	my_defaultB = defaultB;
	my_defaultSinTheta = defaultSinTheta;
	my_defaultLength = defaultLength;
	my_doppler = doppler;
}

SynchrotronEvaluator::~SynchrotronEvaluator()
{
}

void SynchrotronEvaluator::resetParameters(const double *parameters, const double *normalizationUnits){
    //do nothing?
}

void SynchrotronEvaluator::evaluateSynchrotronIandA(const double& photonFinalFrequency, const double& photonFinalTheta, const double& photonFinalPhi, const double& B, const double& sinhi, const double& concentration, MassiveParticleIsotropicDistribution* electronDistribution, double& I, double& A)
{
	//Anu from ghiselini simple
	I = 0;
	A = 0;

	double m = electronDistribution->getMass();
	double m_c2 = m * speed_of_light2;
	double emissivityCoef = sqrt(3.0) * electron_charge * electron_charge * electron_charge / m_c2;
	double criticalNuCoef = 3 * electron_charge / (4 * pi * m * m * m * speed_of_light * speed_of_light4);


	if (sinhi == 0.0) {
		return;
	}

	//double coefAbsorb = concentration * absorpCoef/B;
	//todo what if < 0?
	double coshi = sqrt(1.0 - sinhi * sinhi);

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
			double nuc = criticalNu(my_Ee[j], sinhi, B, criticalNuCoef);
			double gamma = my_Ee[j] / m_c2;

			double mcDonaldIntegral = evaluateMcDonaldIntegral(photonFinalFrequency / nuc);

            //todo h?
            I = I + emissivityCoef * dFe * B * sinhi * mcDonaldIntegral/hplank;
			if (I < 0) {
				printf("evaluateSynchrotronIandA\n");
				printLog("evaluateSynchrotronIandA\n");
				printf("I < 0\n");
				printLog("I < 0\n");
				printf("dFe[j] = %g\n", dFe);
				printLog("dFe[j] = %g\n", dFe);
				printf("B = %g\n", B);
				printLog("B = %g\n", B);
				printf("sinhi = %g\n", sinhi);
				printLog("sinhi = %g\n", sinhi);
				printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				exit(0);
			}
			if (my_selfAbsorption) {
				////todo! F(g +dg)
				double tempP = gamma * gamma * emissivityCoef * B * sinhi * mcDonaldIntegral;
				double dg = 0.1 * delectronEnergy / (m_c2);

				double tempGamma = gamma + dg;
				double tempnuc = criticalNu(m_c2 * tempGamma, sinhi, B, criticalNuCoef);
				double tempP2 = tempGamma * tempGamma * emissivityCoef * B * sinhi * evaluateMcDonaldIntegral(photonFinalFrequency / tempnuc);
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
					printf("sinhi = %g\n", sinhi);
					printLog("sinhi = %g\n", sinhi);
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
				printf("sinhi = %g\n", sinhi);
				printLog("sinhi = %g\n", sinhi);
				printf("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				printLog("mcDonaldIntegral = %g\n", mcDonaldIntegral);
				exit(0);
			}
		}
    }
}

double SynchrotronEvaluator::evaluateFluxFromIsotropicFunction(const double &photonFinalEnergy, MassiveParticleIsotropicDistribution *electronDistribution, const double &volume, const double &distance)
{
	printf("don't use direct flux evaluation with distribution for synchrotron radiation\n");
	printLog("don't use direct flux evaluation with distribution for synchrotron radiation\n");
    double A = 0;
    double I = 0;
    double photonFinalFrequency = photonFinalEnergy/hplank;
    evaluateSynchrotronIandA(photonFinalFrequency, 0, 0, my_defaultB, my_defaultSinTheta, electronDistribution->getConcentration(), electronDistribution, I, A);

	double localI = 0;
	if (my_defaultLength > 0) {
		double area = volume / my_defaultLength;
		if (my_selfAbsorption) {
			double I0 = localI;
			double Q = I * area;
			double tau = A * my_defaultLength;
			double S = 0;
			if (A > 0) {
				S = Q / A;
			}
			if (fabs(tau) < 1E-15) {
				localI = I0 * (1.0 - tau) + S * tau;
			}
			else {
				localI = S + (I0 - S) * exp(-tau);
			}
		}
		else {
			localI = I * volume;
		}
	}

    return localI/sqr(distance);
}

/*double SynchrotronEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();
    double photonFinalFrequency = photonFinalEnergy/hplank;


	double result = 0;
	int irho = 0;

	omp_init_lock(&my_lock);
#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nrho, Nz, Nphi, photonFinalFrequency) reduction(+:result)

	for (irho = 0; irho < Nrho; ++irho) {
		for (int iphi = 0; iphi < Nphi; ++iphi) {
			double localI = 0;
			for (int iz = 0; iz < Nz; ++iz) {
				double area = source->getArea(irho, iz, iphi);
				double A = 0;
				double I = 0;
				double v;
				double theta;
				double phi;
				source->getVelocity(irho, iz, iphi, v, theta, phi);
				if (!my_doppler) {
					v = 0;
				}
				double beta = v / speed_of_light;
				double gamma = 1.0 / sqrt(1 - beta * beta);
				double mu = cos(theta);

				double D = gamma * (1.0 - beta * mu);
				double photonFinalFrequencyPrimed = photonFinalFrequency * D;
				//evaluateSynchrotronIandA(photonFinalFrequency, 0, 0, source->getB(irho, iz, iphi), source->getSinTheta(irho, iz, iphi), source->getConcentration(irho, iz, iphi), source->getParticleDistribution(irho, iz, iphi), I, A);
				evaluateSynchrotronIandA(photonFinalFrequencyPrimed, 0, 0, source->getB(irho, iz, iphi), source->getSinTheta(irho, iz, iphi), source->getConcentration(irho, iz, iphi), source->getParticleDistribution(irho, iz, iphi), I, A);
				double length = source->getLength(irho, iz, iphi);
				if (length > 0) {
					if (my_selfAbsorption) {
						double I0 = localI*D*D;
						double Q = I * area;
						//todo lorentz length
						double lnorm = fabs(length * sin(theta));
						double lpar = fabs(length * cos(theta));
						double lengthPrimed = sqrt(lnorm * lnorm + lpar * lpar * gamma * gamma);
						double tau = A * lengthPrimed;
						double S = 0;
						if (A > 0) {
							S = Q / A;
						}
						if (fabs(tau) < 1E-10) {
							localI = (I0 * (1.0 - tau) + S * tau)/(D*D);
						}
						else {
							localI = (S + (I0 - S) * exp(-tau))/(D*D);
						}
					}
					else {
						localI = localI + I * area * length/(D*D);
					}
				}
			}
			result += localI;
		}
	}

	omp_destroy_lock(&my_lock);

	return result / sqr(source->getDistance());
}*/

double SynchrotronEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi) {
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();
	double photonFinalFrequency = photonFinalEnergy / hplank;
	double localI = 0;
	for (int iz = 0; iz < Nz; ++iz) {
		double area = source->getArea(irho, iz, iphi);
		double A = 0;
		double I = 0;
		double v;
		double theta;
		double phi;
		source->getVelocity(irho, iz, iphi, v, theta, phi);
		if (!my_doppler) {
			v = 0;
		}
		double beta = v / speed_of_light;
		double gamma = 1.0 / sqrt(1 - beta * beta);
		double mu = cos(theta);

		double D = gamma * (1.0 - beta * mu);
		double photonFinalFrequencyPrimed = photonFinalFrequency * D;
		//evaluateSynchrotronIandA(photonFinalFrequency, 0, 0, source->getB(irho, iz, iphi), source->getSinTheta(irho, iz, iphi), source->getConcentration(irho, iz, iphi), source->getParticleDistribution(irho, iz, iphi), I, A);
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticlePowerLawDistribution*>(source->getParticleDistribution(irho, iz, iphi));
		if (distribution == NULL) {
			printf("Synchrotron evaluator works only with isotropic distributions\n");
			printLog("Synchrotron evaluator works only with isotropic distributions\n");
			exit(0);
		}
		evaluateSynchrotronIandA(photonFinalFrequencyPrimed, 0, 0, source->getB(irho, iz, iphi), source->getSinTheta(irho, iz, iphi), source->getConcentration(irho, iz, iphi), distribution, I, A);
		double length = source->getLength(irho, iz, iphi);
		if (length > 0) {
			if (my_selfAbsorption) {
				double I0 = localI * D * D;
				double Q = I * area;
				//todo lorentz length
				double lnorm = fabs(length * sin(theta));
				double lpar = fabs(length * cos(theta));
				double lengthPrimed = sqrt(lnorm * lnorm + lpar * lpar * gamma * gamma);
				double tau = A * lengthPrimed;
				double S = 0;
				if (A > 0) {
					S = Q / A;
				}
				if (fabs(tau) < 1E-10) {
					localI = (I0 * (1.0 - tau) + S * tau) / (D * D);
				}
				else {
					localI = (S + (I0 - S) * exp(-tau)) / (D * D);
				}
			}
			else {
				localI = localI + I * area * length / (D * D);
			}
		}
	}
	double distance = source->getDistance();
	return localI/ (distance * distance);
}