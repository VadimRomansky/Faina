#include "stdio.h"
#include "math.h"

#include "util.h"
#include "constants.h"
#include "radiationSource.h"
#include "radiation.h"

#include "KPIevaluator.h"

SpectrumKPIevaluator::SpectrumKPIevaluator(double* energy, double* observedFlux, double* observedError, int Ne) {
	my_Ne = Ne;
	my_energy = new double[my_Ne];
	my_observedFlux = new double[my_Ne];
	my_observedError = new double[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = energy[i];
		my_observedFlux[i] = observedFlux[i];
		my_observedError[i] = observedError[i];
	}
}

SpectrumKPIevaluator::~SpectrumKPIevaluator()
{
	delete[] my_energy;
	delete[] my_observedFlux;
	delete[] my_observedError;
}

double SpectrumKPIevaluator::evaluate(const double* vector, const double* maxParameters, RadiationSource* source, RadiationEvaluator* evaluator)
{
	double* totalInu = new double[my_Ne];

	source->resetParameters(vector, maxParameters);
	evaluator->resetParameters(vector, maxParameters);

	for (int i = 0; i < my_Ne; ++i) {
		totalInu[i] = evaluator->evaluateFluxFromSource(my_energy[i], source);
	}
	double err = 0;
	//for debug only
	//err = sqr(vector[0] * my_maxParameters[0] - 1.4E17)/1E30 + sqr(vector[1] * my_maxParameters[1] - 0.6) + sqr(vector[2] * my_maxParameters[2] - 10) + sqr(vector[3] * my_maxParameters[3] - 0.5);
	//err = sqr(vector[0] * my_maxParameters[0] - 1.4E17) / 1E30 + sqr(vector[1] * my_maxParameters[1] - 0.6 + vector[2] * my_maxParameters[2] - 10) + sqr(vector[3] * my_maxParameters[3] - 0.5);
	for (int j = 0; j < my_Ne; ++j) {
		double err1 = 0;
		err1 = sqr(totalInu[j] - my_observedFlux[j]) / sqr(my_observedError[j]);

		err = err + err1;
	}

	delete[] totalInu;

	return err;
}

RadialProfileKPIevaluator::RadialProfileKPIevaluator(double energy, double* observedFlux, double* observedError, double* rhoPoints, int Nrho)
{
	my_energy = energy;
	my_Nrho = Nrho;
	my_observedFlux = new double[my_Nrho];
	my_observedError = new double[my_Nrho];
	my_rhoPoints = new double[my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_observedFlux[i] = observedFlux[i];
		my_observedError[i] = observedError[i];
		my_rhoPoints[i] = rhoPoints[i];
	}
}

RadialProfileKPIevaluator::~RadialProfileKPIevaluator()
{
	delete[] my_observedFlux;
	delete[] my_observedError;
	delete[] my_rhoPoints;
}

double RadialProfileKPIevaluator::evaluate(const double* vector, const double* maxParameters, RadiationSource* source, RadiationEvaluator* evaluator)
{
	double* totalInu = new double[my_Nrho];

	source->resetParameters(vector, maxParameters);
	evaluator->resetParameters(vector, maxParameters);

	for (int i = 0; i < my_Nrho; ++i) {
		if (my_rhoPoints[i] > source->getMaxRho() || my_rhoPoints[i] < source->getMinRho()) {
			totalInu[i] = 0;
		}
		else {
			int irho = source->getRhoIndex(my_rhoPoints[i]);
			int iphi = 0;
			totalInu[i] = evaluator->evaluateFluxFromSourceAtPoint(my_energy, source, irho, iphi) / source->getCrossSectionArea(irho, iphi);
		}
	}
	double err = 0;
	//for debug only
	//err = sqr(vector[0] * my_maxParameters[0] - 1.4E17)/1E30 + sqr(vector[1] * my_maxParameters[1] - 0.6) + sqr(vector[2] * my_maxParameters[2] - 10) + sqr(vector[3] * my_maxParameters[3] - 0.5);
	//err = sqr(vector[0] * my_maxParameters[0] - 1.4E17) / 1E30 + sqr(vector[1] * my_maxParameters[1] - 0.6 + vector[2] * my_maxParameters[2] - 10) + sqr(vector[3] * my_maxParameters[3] - 0.5);
	for (int j = 0; j < my_Nrho; ++j) {
		double err1 = 0;
		err1 = sqr(totalInu[j] - my_observedFlux[j]) / sqr(my_observedError[j]);

		err = err + err1;
	}

	delete[] totalInu;

	return err;
}
