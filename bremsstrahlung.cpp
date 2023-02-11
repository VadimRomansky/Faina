#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "bremsstrahlung.h"


BremsstrahlungPrimitiveEvaluator::BremsstrahlungPrimitiveEvaluator(const double& temperature) : RadiationEvaluator(2, massElectron*speed_of_light2, 2*massElectron*speed_of_light2)
{
    my_temperature = temperature;
}

BremsstrahlungPrimitiveEvaluator::~BremsstrahlungPrimitiveEvaluator()
{
}

void BremsstrahlungPrimitiveEvaluator::resetParameters(const double* parameters, const double* normalizationUnits)
{
}

double BremsstrahlungPrimitiveEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance)
{
	double concentration = electronDistribution->getConcentration();
	double m = electronDistribution->getMass();
	double theta = photonFinalEnergy / (kBoltzman * my_temperature);
	double ritberg = m * electron_charge * electron_charge * electron_charge * electron_charge * 2 * pi *pi / (hplank*hplank);
	double eta = kBoltzman * my_temperature / ritberg;
	

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
				gauntFactor = sqrt(12 * ritberg / photonFinalEnergy);
			} else {
				gauntFactor = 1.0;
			}
		}
		else {
			if (-log(theta) >= -0.5 * log(eta)) {
				gauntFactor = (sqrt(3) / pi) * log(1.0 / (4 * pow(euler_mascheroni, 2.5) * theta) * sqrt(eta));
			}
			else {
				gauntFactor = 1.0;
			}
		}
	}

	double result = (1.0 / hplank) * (32 * pi * pow(electron_charge, 6) / (3 * m * speed_of_light * speed_of_light2)) * sqrt(2 * pi / (3 * kBoltzman * my_temperature * m)) * concentration * concentration * exp(-theta)*gauntFactor;
	result *= volume / (4*pi*sqr(distance));
	return result;
}

double BremsstrahlungPrimitiveEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
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
