#ifndef BREMSSTRAHLUNG_H
#define BREMSSTRAHLUNG_H

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

class BremsstrahlungThermalEvaluator : public RadiationEvaluator {
protected:
public:
	BremsstrahlungThermalEvaluator(bool absorption = false, bool doppler = false);
	~BremsstrahlungThermalEvaluator();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	double evaluateGauntFactor(double eta, double theta);
	virtual double evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);
	virtual double evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);
	//virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
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
	BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, bool absorption = false, bool doppler = false);
	BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, double protonsRelativeConcentration, bool absorption = false, bool doppler = false);
	BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, int ionNumber, double* ionConcentrations, int* ionCharges, bool absorption = false, bool doppler = false);
	~BremsstrahlungEvaluator();

	double evaluateSigma(const double& gammaE, const double& epsilonG);

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual void evaluateEmissivityAndAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source, double& I, double& A);
	virtual double evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);
	virtual double evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);
	//virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

#endif