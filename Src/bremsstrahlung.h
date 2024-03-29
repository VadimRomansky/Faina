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

class BremsstrahlungEvaluator : public RadiationEvaluator {
protected:
	int my_ionNumber;
	double* my_ionConcentrations;
	int* my_ionCharges;
	double my_effectiveProtonConcentration;

	double evaluateSigma1(const double& gammaE, const double& epsilonG);
	double evaluateSigma2(const double& gammaE, const double& epsilonG);
	double evaluateA(const double& gammaE, const double& epsilonG);

	double evaluateSigmaNR(const double& gammaE, const double& epsilonG);
	double evaluateB(const double& gammaE);
	double evaluateC(const double& gammaE, const double& x);

	double evaluateSigmaee(const double& gammaE, const double& epsilonG);
	double evaluateSigmape(const double& gammaE, const double& epsilonG);

public:
	BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax);
	BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, double protonsRelativeConcentration);
	BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, int ionNumber, double* ionConcentrations, int* ionCharges);
	~BremsstrahlungEvaluator();

	double evaluateSigma(const double& gammaE, const double& epsilonG);

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

#endif