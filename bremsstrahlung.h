#ifndef BREMSSTRAHLUNG_H
#define BREMSSTRAHLUNG_H

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class BremsstrahlungPrimitiveEvaluator : public RadiationEvaluator {
protected:
	double my_temperature;
public:
	BremsstrahlungPrimitiveEvaluator(const double& temperature);
	~BremsstrahlungPrimitiveEvaluator();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
};

#endif
