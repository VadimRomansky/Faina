#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "bremsstrahlung.h"


BremsstrahlungThermalEvaluator::BremsstrahlungThermalEvaluator() : RadiationEvaluator(2, massElectron*speed_of_light2, 2*massElectron*speed_of_light2)
{
}

BremsstrahlungThermalEvaluator::~BremsstrahlungThermalEvaluator()
{
}

void BremsstrahlungThermalEvaluator::resetParameters(const double* parameters, const double* normalizationUnits)
{
}

double BremsstrahlungThermalEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance)
{
	MassiveParticleMaxwellJuttnerDistribution* maxwellDistribution = dynamic_cast<MassiveParticleMaxwellJuttnerDistribution*>(electronDistribution);
	if (maxwellDistribution == NULL) {
		printf("primitive bremsstrahlung evaluator works only with maxwellian distribution\n");
		printLog("primitive bremsstrahlung evaluator works only with maxwellian distribution\n");
		exit(0);
	}
	double temperature = maxwellDistribution->getTemperature();
	double concentration = electronDistribution->getConcentration();
	double m = electronDistribution->getMass();
	double theta = photonFinalEnergy / (kBoltzman * temperature);
	double ritberg = m * electron_charge * electron_charge * electron_charge * electron_charge * 2 * pi *pi / (hplank*hplank);
	double eta = kBoltzman * temperature / ritberg;
	

	double gauntFactor = 1.0;
	if (eta >= 1.0) {
		if (theta >= 1.0) {
			gauntFactor = sqrt(3 /(pi * theta));
		}
		else {
			gauntFactor = (sqrt(3) / pi) * log(4 / (euler_mascheroni * theta));
		}
	}
	else {
		if (theta >= 1.0) {
			if (log(theta) >= -log(eta)) {
				double dzeta = eta * theta;
				gauntFactor = sqrt(12 /dzeta);
			} else {
				gauntFactor = 1.0;
			}
		}
		else {
			if (-log(theta) >= -0.5 * log(eta)) {
				//where is 4?
				gauntFactor = (sqrt(3) / pi) * log(4.0 / (pow(euler_mascheroni, 2.5) * theta) * sqrt(eta));
			}
			else {
				gauntFactor = 1.0;
			}
		}
	}

	double result = (1.0 / hplank) * (32 * pi * pow(electron_charge, 6) / (3 * m * speed_of_light * speed_of_light2)) * sqrt(2 * pi / (3 * kBoltzman * temperature * m)) * concentration * concentration * exp(-theta)*gauntFactor;
	result *= volume / (4*pi*sqr(distance));
	return result;
}

/*double BremsstrahlungThermalEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;
	int irho = 0;

	omp_init_lock(&my_lock);

#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nrho, Nz, Nphi) reduction(+:result)

	for (irho = 0; irho < Nrho; ++irho) {
		for (int iphi = 0; iphi < Nphi; ++iphi) {
			for (int iz = 0; iz < Nz; ++iz) {
				result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	omp_destroy_lock(&my_lock);

	return result;
}*/

double BremsstrahlungThermalEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi) {
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;

	for (int iz = 0; iz < Nz; ++iz) {
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(irho, iz, iphi));
		if (distribution == NULL) {
			printf("BremsstrahlungThermalEvaluator works only with MassiveParticleIsotropicDistribution\n");
			printLog("BremsstrahlungThermalEvaluator works only with MassiveParticleIsotropicDistribution\n");
			exit(0);
		}
		distribution->resetConcentration(source->getConcentration(irho, iz, iphi));
		result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, distribution, source->getVolume(irho, iz, iphi), source->getDistance());
	}

	return result;
}

double BremsstrahlungEvaluator::evaluateSigma1(const double& gammaE, const double& epsilonG)
{
	if (epsilonG >= gammaE - 1.0) {
		return 0;
	}
	return (4*re2*alpha/epsilonG)*(1.0 + (1.0/3.0 - epsilonG/(gammaE))*(1.0 - epsilonG / (gammaE)))*(log(2*(gammaE)*(gammaE-epsilonG)/epsilonG)-0.5);
}

double BremsstrahlungEvaluator::evaluateSigma2(const double& gammaE, const double& epsilonG)
{
	double coef = re2 * alpha / (3 * epsilonG);
	if (epsilonG < 0.5) {
		return coef * (16*(1.0-epsilonG+epsilonG*epsilonG)*log(gammaE/epsilonG) - 1.0/sqr(epsilonG) + 3.0/epsilonG - 4.0 + 4.0*epsilonG - 8.0*epsilonG*epsilonG - 2*(1 - 2*epsilonG)*log(1.0 - 2.0*epsilonG)*(1.0/(4.0*cube(epsilonG))-1/(2.0*sqr(epsilonG)+3.0/epsilonG-2.0+4.0*epsilonG)));
	}
	else {
		return coef * (2.0/epsilonG)*((4.0 - 1/epsilonG + 1/(4.0*sqr(epsilonG)))*log(2.0*gammaE)-2.0+2/epsilonG-5/(8.0*sqr(epsilonG)));
	}
}

double BremsstrahlungEvaluator::evaluateA(const double& gammaE, const double& epsilonG)
{
	double A = 1.0 - (10.0/3.0)*pow(gammaE - 1.0, 1.0/5.0)/(gammaE + 1.0)*pow(epsilonG/gammaE, 1.0/3.0);
	if (A < 0) {
		printf("BremsstrahlungEeEvaluator::evaluateA A < 0 gammaE = %g epsilonG = %g\n", gammaE, epsilonG);
		printLog("BremsstrahlungEeEvaluator::evaluateA A < 0 gammaE = %g epsilonG = %g\n", gammaE, epsilonG);
		exit(0);
	}
	return A;
}

double BremsstrahlungEvaluator::evaluateSigmaNR(const double& gammaE, const double& epsilonG) {

}

double BremsstrahlungEvaluator::evaluateSigmaee(const double& gammaE, const double& epsilonG) {

}
double BremsstrahlungEvaluator::evaluateSigmape(const double& gammaE, const double& epsilonG) {

}

double BremsstrahlungEvaluator::evaluateSigma(const double& gammaE, const double& epsilonG) {
	double sigma = 0;
	if (my_ionNumber > 0) {
		sigma += evaluateSigmape(gammaE, epsilonG) * my_effectiveProtonConcentration;
	}

	sigma += evaluateSigmaee(gammaE, epsilonG);

	return sigma;
}

BremsstrahlungEvaluator::BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax) : RadiationEvaluator(Ne, Emin, Emax) {
	my_ionNumber = 0;
	my_ionConcentrations = NULL;
	my_ionConcentrations = NULL;
	my_effectiveProtonConcentration = 0;
}

BremsstrahlungEvaluator::BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, double protonsRelativeConcentration) : RadiationEvaluator(Ne, Emin, Emax) {
	my_ionNumber = 1;
	my_ionConcentrations = new double[my_ionNumber];
	my_ionConcentrations[0] = protonsRelativeConcentration;
	my_ionCharges = new int[my_ionNumber];
	my_ionCharges[0] = 1;
	my_effectiveProtonConcentration = protonsRelativeConcentration;
}


BremsstrahlungEvaluator::BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, int ionNumber, double* ionConcentrations, int* ionCharges) : RadiationEvaluator(Ne, Emin, Emax) {
	my_ionNumber = ionNumber;
	my_ionConcentrations = new double[my_ionNumber];
	my_ionCharges = new int[my_ionNumber];
	my_effectiveProtonConcentration = 0;
	for (int i = 0; i < my_ionNumber; ++i) {
		my_ionConcentrations[i] = ionConcentrations[i];
		my_ionCharges[i] = ionCharges[i];
		my_effectiveProtonConcentration += my_ionConcentrations[i] * my_ionCharges[i] * my_ionCharges[i];
	}
}

BremsstrahlungEvaluator::~BremsstrahlungEvaluator()
{
	if (my_ionNumber > 0) {
		delete[] my_ionConcentrations;
		delete[] my_ionCharges;
	}
}

void BremsstrahlungEvaluator::resetParameters(const double* parameters, const double* normalizationUnits)
{
}

double BremsstrahlungEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance)
{
	double result = 0;

	for (int i = 0; i < my_Ne; ++i) {
		double electronEnergy = my_Ee[i];
		double delectronEnergy;
		if (i == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[i] - my_Ee[i - 1];
		}
		double electronGamma = electronEnergy / (massElectron * speed_of_light2);
		double electronBeta = sqrt(1.0 - 1.0 / (electronGamma * electronGamma));
		double electronKineticEnergy = massElectron * speed_of_light2 * (electronGamma - 1.0);
		double epsilonG = photonFinalEnergy / (massElectron * speed_of_light2);

		double sigma = evaluateSigma(electronGamma, epsilonG);

		result += photonFinalEnergy * (speed_of_light * electronBeta) * sigma * (4*Pi*electronDistribution->distribution(electronEnergy)) * volume * delectronEnergy / (4*pi*sqr(distance));

		if (result != result) {
			printf("result = NaN in bremsstrahlung\n");
			printLog("result = NaN in bremsstrahlung\n");
			exit(0);
		}
	}

	return result;
}

double BremsstrahlungEvaluator::evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;

	for (int iz = 0; iz < Nz; ++iz) {
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(irho, iz, iphi));
		if (distribution == NULL) {
			printf("Bremsstrahlung e-e evaluator works only with isotropic electrons distribution\n");
			printLog("Bremsstrahlung e-e evaluator works only with isotropic electrons distribution\n");
			exit(0);
		}
		result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, distribution, source->getVolume(irho, iz, iphi), source->getDistance());
	}

	return result;
}