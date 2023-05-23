#include "stdio.h"
#include "math.h"

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
	MassiveParticleMaxwellDistribution* maxwellDistribution = dynamic_cast<MassiveParticleMaxwellDistribution*>(electronDistribution);
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

double BremsstrahlungThermalEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;
#pragma omp parallel for shared(result, photonFinalEnergy, source, Nrho, Nz, Nphi, photonFinalFrequency)	

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	return result;
}
