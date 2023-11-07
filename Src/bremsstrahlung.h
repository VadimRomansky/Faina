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

class BremsstrahlungEeEvaluator : public RadiationEvaluator {
protected:
	double my_ambientConcentration;

	double evaluateSigma1(const double& gammaE, const double& epsilonG);
	double evaluateSigma2(const double& gammaE, const double& epsilonG);
	double evaluateA(const double& gammaE, const double& epsilonG);
public:
	BremsstrahlungEeEvaluator(int Ne, const double& Emin, const double& Emax, const double& ambientConcentration);
	~BremsstrahlungEeEvaluator();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

class BremsstrahlungPeEvaluator : public RadiationEvaluator {
protected:
	double my_ambientConcentration;

	double evaluateSigma(const double& p1, const double& p2, const double& epsilonG);
public:
	BremsstrahlungPeEvaluator(int Ne, const double& Emin, const double& Emax, const double& ambientConcentration);
	~BremsstrahlungPeEvaluator();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

#endif
