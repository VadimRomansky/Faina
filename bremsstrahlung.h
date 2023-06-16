#ifndef BREMSSTRAHLUNG_H
#define BREMSSTRAHLUNG_H

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class BremsstrahlungThermalEvaluator : public RadiationEvaluator {
protected:
public:
	BremsstrahlungThermalEvaluator();
	~BremsstrahlungThermalEvaluator();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

#endif
