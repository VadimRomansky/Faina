#ifndef KPI_EVALUATOR_H
#define KPI_EVALUATOR_H

#include "radiation.h"
#include "radiationSource.h"

class KPIevaluator {
public:
	virtual double evaluate(const double* vector, const double* maxParameters, RadiationSource* source, RadiationEvaluator* evaluator) = 0;
};

class SpectrumKPIevaluator : public KPIevaluator {
protected:
	int my_Ne;
	double* my_energy;
	double* my_observedFlux;
	double* my_observedError;
public:
	SpectrumKPIevaluator(double* energy, double* observedFlux, double* observedError, int Ne);
	~SpectrumKPIevaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationSource* source, RadiationEvaluator* evaluator);
};

class RadialProfileKPIevaluator : public KPIevaluator {
protected:
	int my_Nrho;
	double my_energy;
	double* my_observedFlux;
	double* my_observedError;
	double* my_rhoPoints;
public:
	RadialProfileKPIevaluator(double energy, double* observedFlux, double* observedError, double* rhoPoints, int Nrho);
	~RadialProfileKPIevaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationSource* source, RadiationEvaluator* evaluator);
};

#endif
