#ifndef KPI_EVALUATOR_H
#define KPI_EVALUATOR_H

#include "radiation.h"
#include "radiationSource.h"

class LossEvaluator {
public:
	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator) = 0;
};

class SpectrumLossEvaluator : public LossEvaluator {
protected:
	RadiationSourceInCylindrical* my_radiationSource;
	int my_Ne;
	double* my_energy;
	double* my_observedFlux;
	double* my_observedError;
public:
	SpectrumLossEvaluator(double* energy, double* observedFlux, double* observedError, int Ne, RadiationSourceInCylindrical* radiatiornSource);
	~SpectrumLossEvaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator);
};

class TimeDependentSpectrumLossEvaluator : public LossEvaluator {
protected:
	RadiationTimeDependentSource* my_radiationSource;
	int my_Ntimes;
	int* my_Ne;
	double* my_times;
	double** my_energy;
	double** my_observedFlux;
	double** my_observedError;
public:
	TimeDependentSpectrumLossEvaluator(double** energy, double** observedFlux, double** observedError, int* Ne, double* times, int Ntimes, RadiationTimeDependentSource* radiationSource);
	~TimeDependentSpectrumLossEvaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator);
};

class RadialProfileLossEvaluator : public LossEvaluator {
protected:
	RadiationSourceInCylindrical* my_radiationSource;
	int my_Nrho;
	double my_energy;
	double* my_observedFlux;
	double* my_observedError;
	double* my_rhoPoints;
public:
	RadialProfileLossEvaluator(double energy, double* observedFlux, double* observedError, double* rhoPoints, int Nrho, RadiationSourceInCylindrical* radiaionSource);
	~RadialProfileLossEvaluator();

	virtual double evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator);
};

#endif
