#include "stdio.h"
#include "math.h"

#include "util.h"
#include "constants.h"
#include "radiationSource.h"
#include "radiation.h"

#include "KPIevaluator.h"

SpectrumLossEvaluator::SpectrumLossEvaluator(double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* radiationSource) {
	my_Ne = Ne;
	my_energy = new double[my_Ne];
	my_observedFlux = new double[my_Ne];
	my_observedError = new double[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = energy[i];
		my_observedFlux[i] = observedFlux[i];
		my_observedError[i] = observedError[i];
	}
	my_radiationSource = radiationSource;
}

SpectrumLossEvaluator::~SpectrumLossEvaluator()
{
	delete[] my_energy;
	delete[] my_observedFlux;
	delete[] my_observedError;
}

double SpectrumLossEvaluator::evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator)
{
	double* totalInu = new double[my_Ne];

	my_radiationSource->resetParameters(vector, maxParameters);
	evaluator->resetParameters(vector, maxParameters);

	for (int i = 0; i < my_Ne; ++i) {
		totalInu[i] = evaluator->evaluateFluxFromSource(my_energy[i], my_radiationSource);
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

RadialProfileLossEvaluator::RadialProfileLossEvaluator(double energy, double* observedFlux, double* observedError, double* rhoPoints, int Nrho, RadiationSource* radiationSource)
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
	my_radiationSource = radiationSource;
}

RadialProfileLossEvaluator::~RadialProfileLossEvaluator()
{
	delete[] my_observedFlux;
	delete[] my_observedError;
	delete[] my_rhoPoints;
}

double RadialProfileLossEvaluator::evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator)
{
	double* totalInu = new double[my_Nrho];

	my_radiationSource->resetParameters(vector, maxParameters);
	evaluator->resetParameters(vector, maxParameters);

	for (int i = 0; i < my_Nrho; ++i) {
		if (my_rhoPoints[i] > my_radiationSource->getMaxRho() || my_rhoPoints[i] < my_radiationSource->getMinRho()) {
			totalInu[i] = 0;
		}
		else {
			int irho = my_radiationSource->getRhoIndex(my_rhoPoints[i]);
			int iphi = 0;
			totalInu[i] = evaluator->evaluateFluxFromSourceAtPoint(my_energy, my_radiationSource, irho, iphi) / my_radiationSource->getCrossSectionArea(irho, iphi);
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

TimeDependentSpectrumLossEvaluator::TimeDependentSpectrumLossEvaluator(double** energy, double** observedFlux, double** observedError, int* Ne, double* times, int Ntimes, RadiationTimeDependentSource* radiationSource)
{
	my_Ntimes = Ntimes;

	my_Ne = new int[my_Ntimes];
	my_times = new double[my_Ntimes];
	my_energy = new double*[my_Ntimes];
	my_observedFlux = new double*[my_Ntimes];
	my_observedError = new double*[my_Ntimes];
	for (int i = 0; i < my_Ntimes; ++i) {
		my_times[i] = times[i];
		my_Ne[i] = Ne[i];
		my_energy[i] = new double[my_Ne[i]];
		my_observedFlux[i] = new double[my_Ne[i]];
		my_observedError[i] = new double[my_Ne[i]];
		for (int j = 0; j < my_Ne[i]; ++j) {
			my_energy[i][j] = energy[i][j];
			my_observedFlux[i][j] = observedFlux[i][j];
			my_observedError[i][j] = observedError[i][j];
		}
	}
	my_radiationSource = radiationSource;
}

TimeDependentSpectrumLossEvaluator::~TimeDependentSpectrumLossEvaluator()
{
	for (int i = 0; i < my_Ntimes; ++i) {
		delete[] my_energy[i];
		delete[] my_observedFlux[i];
		delete[] my_observedError[i];
	}
	delete[] my_energy;
	delete[] my_observedFlux;
	delete[] my_observedError;

	delete[] my_times;
	delete[] my_Ne;
}

double TimeDependentSpectrumLossEvaluator::evaluate(const double* vector, const double* maxParameters, RadiationEvaluator* evaluator)
{
	my_radiationSource->resetParameters(vector, maxParameters);
	evaluator->resetParameters(vector, maxParameters);
	double err = 0;
	for (int k = 0; k < my_Ntimes; ++k) {
		//double* totalInu = new double[Ne[k]];

		double totalInu = 0;
		RadiationSource* source1 = my_radiationSource->getRadiationSource(my_times[k], maxParameters);
		for (int j = 0; j < my_Ne[k]; ++j) {
			//1E26 from Jansky
			totalInu = evaluator->evaluateFluxFromSource(my_energy[k][j], source1);

			double err1 = 0;
			err1 = sqr(totalInu - my_observedFlux[k][j]) / sqr(my_observedError[k][j]);
			//todo hack!!!
			if (k == my_Ntimes - 1) {
				err1 *= 3;
			}

			err = err + err1;
		}

		//delete[] totalInu;
	}

	if (isinf(err)) {
		printf("error function is infinity\n");
		printLog("error function is infinity\n");
		exit(0);
	}

	if (err != err) {
		printf("error function is NaN\n");
		printLog("error function is NaN\n");
		exit(0);
	}

	return err;
}
