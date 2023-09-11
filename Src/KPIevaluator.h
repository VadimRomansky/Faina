#ifndef KPI_EVALUATOR_H
#define KPI_EVALUATOR_H

#include "radiation.h"
#include "radiationSource.h"

class KPIevaluator {
public:
	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator) = 0;
};

class SpectrumKPIevaluator : public KPIevaluator {
protected:
	RadiationSource* my_radiationSource;
	int my_Ne;
	double* my_energy;
	double* my_observedFlux;
	double* my_observedError;
public:
	SpectrumKPIevaluator(double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* radiatiornSource);
	~SpectrumKPIevaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator);
};

class TimeDependentSpectrumKPIevaluator : public KPIevaluator {
protected:
	RadiationTimeDependentSource* my_radiationSource;
	int my_Ntimes;
	int* my_Ne;
	double* my_times;
	double** my_energy;
	double** my_observedFlux;
	double** my_observedError;
public:
	TimeDependentSpectrumKPIevaluator(double** energy, double** observedFlux, double** observedError, int* Ne, double* times, int Ntimes, RadiationTimeDependentSource* radiationSource);
	~TimeDependentSpectrumKPIevaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator);
};

class RadialProfileKPIevaluator : public KPIevaluator {
protected:
	RadiationSource* my_radiationSource;
	int my_Nrho;
	double my_energy;
	double* my_observedFlux;
	double* my_observedError;
	double* my_rhoPoints;
public:
	RadialProfileKPIevaluator(double energy, double* observedFlux, double* observedError, double* rhoPoints, int Nrho, RadiationSource* radiaionSource);
	~RadialProfileKPIevaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator);
};

#endif
