#include "stdio.h"
#include "math.h"

#include "util.h"
#include "constants.h"
#include "radiationSource.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "synchrotron.h"
#include "inverseCompton.h"

#include "optimization.h"


RadiationOptimizer::RadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams) {
	my_evaluator = evaluator;
	my_Nparams = Nparams;
	my_minParameters = new double[my_Nparams];
	my_maxParameters = new double[my_Nparams];
	my_minVector = new double[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_minParameters[i] = minParameters[i];
		my_maxParameters[i] = maxParameters[i];
		my_minVector[i] = my_minParameters[i] / my_maxParameters[i];
	}
}
RadiationOptimizer::~RadiationOptimizer() {
	delete[] my_minParameters;
	delete[] my_maxParameters;
	delete[] my_minVector;
}

double RadiationOptimizer::evaluateOptimizationFunction(const double* vector, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source) {
	double* totalInu = new double[Ne];
	
	source->resetParameters(vector, my_maxParameters);
    my_evaluator->resetParameters(vector, my_maxParameters);

	for (int i = 0; i < Ne; ++i) {
        totalInu[i] = my_evaluator->evaluateFluxFromSource(energy[i], source);
	}
	double err = 0;
	//for debug only
	//err = sqr(vector[0] * my_maxParameters[0] - 1.4E17)/1E30 + sqr(vector[1] * my_maxParameters[1] - 0.6) + sqr(vector[2] * my_maxParameters[2] - 10) + sqr(vector[3] * my_maxParameters[3] - 0.5);
	//err = sqr(vector[0] * my_maxParameters[0] - 1.4E17) / 1E30 + sqr(vector[1] * my_maxParameters[1] - 0.6 + vector[2] * my_maxParameters[2] - 10) + sqr(vector[3] * my_maxParameters[3] - 0.5);
	for (int j = 0; j < Ne; ++j) {
		double err1 = 0;
		err1 = sqr(totalInu[j] - observedFlux[j]) / sqr(observedError[j]);

		err = err + err1;
	}

	delete[] totalInu;

	return err;
}

void RadiationOptimizer::optimize(double* vector, bool* optPar, double* energy, double* observedFlux, int Ne, RadiationSource* source) {
	double* observedError = new double[Ne];
	for (int i = 0; i < Ne; ++i) {
		observedError[i] = 1.0;
	}

	optimize(vector, optPar, energy, observedFlux, observedError, Ne, source);

	delete[] observedError;
}

GradientDescentRadiationOptimizer::GradientDescentRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations) : RadiationOptimizer(evaluator, minParameters, maxParameters, Nparams) {
	my_Niterations = Niterations;

	tempVector = new double[my_Nparams];
	tempVector1 = new double[my_Nparams];
	tempVector2 = new double[my_Nparams];
	prevVector = new double[my_Nparams];
	currentVector = new double[my_Nparams];
	valley1 = new double[my_Nparams];
	valley2 = new double[my_Nparams];
	grad = new double[my_Nparams];
}

GradientDescentRadiationOptimizer::~GradientDescentRadiationOptimizer() {
	delete[] tempVector;
	delete[] tempVector1;
	delete[] tempVector2;
	delete[] prevVector;
	delete[] currentVector;
	delete[] valley1;
	delete[] valley2;
	delete[] grad;
}

void GradientDescentRadiationOptimizer::findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source, const double& currentF) {
	int Npar = 0;
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}

	double maxLambda = sqrt(1.0 * Npar);
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - my_minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}

	double step = maxLambda / 2;

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - grad[i] * step;
			if (tempVector1[i] != tempVector1[i]) {
				printf("tempvector[i] != tempVector[i]\n");
				exit(0);
			}
		}
	}

	double f = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, source);
	if (f > currentF) {
		int count = 0;
		while (f > currentF && count < 20) {
			count++;
			step = step / 2;

			for (int i = 0; i < my_Nparams; ++i) {
				tempVector1[i] = vector[i];
			}
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					tempVector1[i] = vector[i] - grad[i] * step;
					if (tempVector1[i] != tempVector1[i]) {
						printf("tempvector[i] != tempVector[i]\n");
						exit(0);
					}
				}
			}

			f = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, source);
		}

		if (f < currentF) {
			//todo!
			for (int i = 0; i < my_Nparams; ++i) {
				vector[i] = tempVector1[i];
			}
			//return;
		}
		else {
			return;
		}
	}

	f = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, source);

	maxLambda = sqrt(1.0 * Npar);
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - my_minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - maxLambda * grad[i];
			if (tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if (tempVector1[i] < my_minVector[i]) {
				tempVector1[i] = my_minVector[i];
			}
		}
	}
	double f3 = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, source);

	double leftLambda = 0;
	double rightLambda = maxLambda;


	for (int j = 0; j < 10; ++j) {
		double lambda1 = leftLambda + (rightLambda - leftLambda) / 3.0;
		double lambda2 = leftLambda + (rightLambda - leftLambda) * 2.0 / 3.0;

		for (int i = 0; i < my_Nparams; ++i) {
			tempVector1[i] = vector[i];
			tempVector2[i] = vector[i];
		}
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				tempVector1[i] = vector[i] - lambda1 * grad[i];
				tempVector2[i] = vector[i] - lambda2 * grad[i];
			}
		}

		double f1 = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, source);
		double f2 = evaluateOptimizationFunction(tempVector2, energy, observedFlux, observedError, Ne, source);

		if (f1 > f2) {
			leftLambda = lambda1;
		}
		else {
			rightLambda = lambda2;
		}
	}

	double lambda = (leftLambda + rightLambda) / 2.0;
	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - lambda * grad[i];
			if (tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if (tempVector1[i] < my_minVector[i]) {
				tempVector1[i] = my_minVector[i];
			}
		}
	}

	double tempF = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, source);

	if (f < tempF) {
		if (f < f3) {
			return;
		}
		else {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - maxLambda * grad[i];
				}
			}
		}
	}
	else {
		if (tempF < f3) {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - lambda * grad[i];
				}
			}
		}
		else {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - maxLambda * grad[i];
				}
			}
		}
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (vector[i] > 1.0) {
				vector[i] = 1.0;
			}
			if (vector[i] < my_minVector[i]) {
				vector[i] = my_minVector[i];
			}
		}
	}
}

void GradientDescentRadiationOptimizer::optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source) {
	int Npar = 0;
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}

	for (int i = 0; i < my_Nparams; ++i) {
		prevVector[i] = vector[i];
		currentVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, source);
	printf("optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}
	
	for (int k = 0; k < my_Niterations; ++k) {
		///randomization;
		for(int j = 0; j < 5; ++j){
			for(int i = 0; i < my_Nparams; ++i) {
				if(optPar[i]){
					//tempVector[i] = minVector[i] + (1.0 - minVector[i])*uniformDistribution();
					tempVector[i] = vector[i] + 0.2*my_minVector[i]*(0.5 - uniformDistribution());
					if(tempVector[i] > 1.0) {
						tempVector[i] = 1.0;
					}
					if(tempVector[i] < my_minVector[i]) {
						tempVector[i] = my_minVector[i];
					}
				}
			}

			double tempF = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, source);
			if(tempF < currentF){
				currentF = tempF;
				for(int i = 0; i < my_Nparams; ++i) {
					vector[i] = tempVector[i];
				}
				printf("optimization function = %g\n", currentF);
				printf("random search succeed\n");
				printLog("random search succeed\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		for (int i = 0; i < my_Nparams; ++i) {
			valley1[i] = vector[i];
		}

		printf("optimization k = %d\n", k);

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < my_Nparams; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i]) / 100;
				if (k > 0) {
					dx = min(dx,max(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2));
				}
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, source);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, source);

				grad[i] = (f - f1) / (2 * dx);
				if (grad[i] != grad[i]) {
					printf("grad[i] = NaN\n");
					exit(0);
				}
				if (grad[i] > 0) {
					if (vector[i] - my_minVector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
				if (grad[i] < 0) {
					if (1.0 - vector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
			}
		}

		double gradNorm = 0;
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm <= 0) {
			printf("gradNorm = %g\n", gradNorm);
			continue;
		}
		if (gradNorm != gradNorm) {
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = grad[i] / gradNorm;
			}
		}


		findMinParametersAtDirection(vector, optPar, grad, energy, observedFlux, observedError, Ne, source, currentF);

		currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, source);
		//valley second step

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < my_Nparams; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i]) / 100;
				if (k > 0) {
					dx = min(dx,max(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2));
				}
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, source);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, source);

				grad[i] = (f - f1) / (2 * dx);
				if (grad[i] != grad[i]) {
					printf("grad[i] = NaN\n");
					exit(0);
				}
				if (grad[i] > 0) {
					if (vector[i] - my_minVector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
				if (grad[i] < 0) {
					if (1.0 - vector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
			}
		}

		gradNorm = 0;
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm <= 0) {
			printf("gradNorm <= 0\n");
			continue;
		}
		if (gradNorm != gradNorm) {
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = grad[i] / gradNorm;
			}
		}

		findMinParametersAtDirection(vector, optPar, grad, energy, observedFlux, observedError, Ne, source, currentF);

		currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, source);

		for (int i = 0; i < my_Nparams; ++i) {
			valley2[i] = vector[i];
		}

		//// valley third step

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = valley1[i] - valley2[i];
			}
		}


		gradNorm = 0;
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm <= 0) {
			printf("gradNorm <= 0\n");
		}
		if (gradNorm != gradNorm) {
			printf("gradNorm = NaN\n");
			exit(0);
		}

		if (gradNorm > 0) {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					grad[i] = grad[i] / gradNorm;
				}
			}

			findMinParametersAtDirection(vector, optPar, grad, energy, observedFlux, observedError, Ne, source, currentF);

		}
		currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, source);

		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/
		for (int i = 0; i < my_Nparams; ++i) {
			prevVector[i] = currentVector[i];
			currentVector[i] = vector[i];
		}
		printf("optimization function = %g\n", currentF);
		for (int i = 0; i < my_Nparams; ++i) {
			printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
		}
	}
	printf("GradientDescentOptimizer: final optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}
	printf("finish optimization\n");
}

GridEnumRadiationOptimizer::GridEnumRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, const int* Npoints) : RadiationOptimizer(evaluator, minParameters, maxParameters, Nparams)
{
	my_Npoints = new int[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_Npoints[i] = Npoints[i];
	}

	my_points = new double* [my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_points[i] = new double[my_Npoints[i]];
		my_points[i][0] = my_minVector[i];
		my_points[i][my_Npoints[i] - 1] = 1.0;
		if (my_Npoints[i] > 1) {
			double delta = (my_points[i][my_Npoints[i] - 1] - my_points[i][0]) / (my_Npoints[i] - 1);
			double factor = pow(my_points[i][my_Npoints[i] - 1]/my_points[i][0], 1.0/(my_Npoints[i] - 1));
			for (int j = 1; j < my_Npoints[i] - 1; ++j) {
				my_points[i][j] = my_points[i][j - 1]*factor;
			}
		}
	}
}

GridEnumRadiationOptimizer::~GridEnumRadiationOptimizer()
{
	for (int i = 0; i < my_Nparams; ++i) {
		delete[] my_points[i];
	}
	delete[] my_points;
	delete[] my_Npoints;
}

void GridEnumRadiationOptimizer::optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source)
{
	double* tempVector = new double[my_Nparams];

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, source);
	int N = 1;
	for (int i = 0; i < my_Nparams; ++i) {
		N *= my_Npoints[i];
	}

	int* index = new int[my_Nparams];
	for (int i = 0; i < N; ++i) {
		int tempN = i;
		for (int j = 0; j < my_Nparams; ++j) {
			index[j] = tempN % my_Npoints[j];
			tempN = tempN / my_Npoints[j];
			if (optPar[j]) {
				tempVector[j] = my_points[j][index[j]];
			}
		}
		double tempF = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, source);
		if (tempF < currentF) {
			currentF = tempF;
			for (int j = 0; j < my_Nparams; ++j) {
				if (optPar[j]) {
					vector[j] = tempVector[j];
				}
			}
			printf("optimization function = %g\n", tempF);
		}
	}

	printf("GridENumOptimizer: final optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}

	delete[] index;
	delete[] tempVector;
}

RadiationTimeOptimizer::RadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams)
{
	my_Nparams = Nparams;
	my_evaluator = evaluator;
	my_minParameters = new double[my_Nparams];
	my_maxParameters = new double[my_Nparams];
	my_minVector = new double[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_minParameters[i] = minParameters[i];
		my_maxParameters[i] = maxParameters[i];
		my_minVector[i] = my_minParameters[i] / my_maxParameters[i];
	}
}

RadiationTimeOptimizer::~RadiationTimeOptimizer()
{
	delete[] my_minParameters;
	delete[] my_maxParameters;
	delete[] my_minVector;
}

double RadiationTimeOptimizer::evaluateOptimizationFunction(const double* vector, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source)
{
	source->resetParameters(vector, my_maxParameters);
    my_evaluator->resetParameters(vector, my_maxParameters);
	double err = 0;
	for (int k = 0; k < Ntimes; ++k) {
		double* totalInu = new double[Ne[k]];


		RadiationSource* source1 = source->getRadiationSource(times[k], my_maxParameters);
		for (int i = 0; i < Ne[k]; ++i) {
			//1E26 from Jansky
            totalInu[i] = my_evaluator->evaluateFluxFromSource(energy[k][i], source1);
		}
		for (int j = 0; j < Ne[k]; ++j) {
			double err1 = 0;
			err1 = sqr(totalInu[j] - observedFlux[k][j]) / sqr(observedError[k][j]);

			err = err + err1;
		}

		delete[] totalInu;
	}

	return err;
}

void RadiationTimeOptimizer::optimize(double* vector, bool* optPar, double** energy, double** observedFlux, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source)
{
	double** observedError = new double*[Ntimes];
	for (int i = 0; i < Ntimes; ++i) {
		observedError[i] = new double[Ne[i]];
		for (int j = 0; j < Ne[i]; ++j) {
			observedError[i][j] = 1.0;
		}
	}

	optimize(vector, optPar, energy, observedFlux, observedError, Ne, Ntimes, times, source);

	for (int i = 0; i < Ntimes; ++i) {
		delete[] observedError[i];
	}
	delete[] observedError;
}

GradientDescentRadiationTimeOptimizer::GradientDescentRadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations) : RadiationTimeOptimizer(evaluator, minParameters, maxParameters, Nparams)
{
	my_Niterations = Niterations;

	tempVector = new double[my_Nparams];
	tempVector1 = new double[my_Nparams];
	tempVector2 = new double[my_Nparams];
	prevVector = new double[my_Nparams];
	currentVector = new double[my_Nparams];
	valley1 = new double[my_Nparams];
	valley2 = new double[my_Nparams];
	grad = new double[my_Nparams];
}

GradientDescentRadiationTimeOptimizer::~GradientDescentRadiationTimeOptimizer()
{
	delete[] tempVector;
	delete[] tempVector1;
	delete[] tempVector2;
	delete[] prevVector;
	delete[] currentVector;
	delete[] valley1;
	delete[] valley2;
	delete[] grad;
}

void GradientDescentRadiationTimeOptimizer::findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source, const double& currentF)
{
	int Npar = 0;
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}

	double maxLambda = sqrt(1.0 * Npar);
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - my_minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}

	double step = maxLambda / 2;

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - grad[i] * step;
			if (tempVector1[i] != tempVector1[i]) {
				printf("tempvector[i] != tempVector[i]\n");
				exit(0);
			}
		}
	}

	double f = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, Ntimes, times, source);
	if (f > currentF) {
		int count = 0;
		while (f > currentF && count < 20) {
			count++;
			step = step / 2;

			for (int i = 0; i < my_Nparams; ++i) {
				tempVector1[i] = vector[i];
			}
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					tempVector1[i] = vector[i] - grad[i] * step;
					if (tempVector1[i] != tempVector1[i]) {
						printf("tempvector[i] != tempVector[i]\n");
						exit(0);
					}
				}
			}

			f = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, Ntimes, times, source);
		}

		if (f < currentF) {
			//todo!
			for (int i = 0; i < my_Nparams; ++i) {
				vector[i] = tempVector1[i];
			}
			//return;
		}
		else {
			return;
		}
	}

	f = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, Ntimes, times, source);

	maxLambda = sqrt(1.0 * Npar);
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - my_minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - maxLambda * grad[i];
			if (tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if (tempVector1[i] < my_minVector[i]) {
				tempVector1[i] = my_minVector[i];
			}
		}
	}
	double f3 = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, Ntimes, times, source);

	double leftLambda = 0;
	double rightLambda = maxLambda;


	for (int j = 0; j < 10; ++j) {
		double lambda1 = leftLambda + (rightLambda - leftLambda) / 3.0;
		double lambda2 = leftLambda + (rightLambda - leftLambda) * 2.0 / 3.0;

		for (int i = 0; i < my_Nparams; ++i) {
			tempVector1[i] = vector[i];
			tempVector2[i] = vector[i];
		}
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				tempVector1[i] = vector[i] - lambda1 * grad[i];
				tempVector2[i] = vector[i] - lambda2 * grad[i];
			}
		}

		double f1 = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, Ntimes, times, source);
		double f2 = evaluateOptimizationFunction(tempVector2, energy, observedFlux, observedError, Ne, Ntimes, times, source);

		if (f1 > f2) {
			leftLambda = lambda1;
		}
		else {
			rightLambda = lambda2;
		}
	}

	double lambda = (leftLambda + rightLambda) / 2.0;
	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			tempVector1[i] = vector[i] - lambda * grad[i];
			if (tempVector1[i] > 1.0) {
				tempVector1[i] = 1.0;
			}
			if (tempVector1[i] < my_minVector[i]) {
				tempVector1[i] = my_minVector[i];
			}
		}
	}

	double tempF = evaluateOptimizationFunction(tempVector1, energy, observedFlux, observedError, Ne, Ntimes, times, source);

	if (f < tempF) {
		if (f < f3) {
			return;
		}
		else {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - maxLambda * grad[i];
				}
			}
		}
	}
	else {
		if (tempF < f3) {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - lambda * grad[i];
				}
			}
		}
		else {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					vector[i] = vector[i] - maxLambda * grad[i];
				}
			}
		}
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (vector[i] > 1.0) {
				vector[i] = 1.0;
			}
			if (vector[i] < my_minVector[i]) {
				vector[i] = my_minVector[i];
			}
		}
	}
}

void GradientDescentRadiationTimeOptimizer::optimize(double* vector, bool* optPar, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source)
{
	int Npar = 0;
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			Npar++;
		}
	}

	for (int i = 0; i < my_Nparams; ++i) {
		prevVector[i] = vector[i];
		currentVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
	printf("optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}

	for (int k = 0; k < my_Niterations; ++k) {
		///randomization;
		for (int j = 0; j < 5; ++j) {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					//tempVector[i] = minVector[i] + (1.0 - minVector[i])*uniformDistribution();
					tempVector[i] = vector[i] + 0.2 * my_minVector[i] * (0.5 - uniformDistribution());
					if (tempVector[i] > 1.0) {
						tempVector[i] = 1.0;
					}
					if (tempVector[i] < my_minVector[i]) {
						tempVector[i] = my_minVector[i];
					}
				}
			}

			double tempF = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
			if (tempF < currentF) {
				currentF = tempF;
				for (int i = 0; i < my_Nparams; ++i) {
					vector[i] = tempVector[i];
				}
				printf("optimization function = %g\n", currentF);
				printf("random search succeed\n");
				printLog("random search succeed\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		for (int i = 0; i < my_Nparams; ++i) {
			valley1[i] = vector[i];
		}

		printf("optimization k = %d\n", k);

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < my_Nparams; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(my_minVector[i]) / 100;
				/*if (k > 0) {
					dx = min(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2);
				}*/
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, Ntimes, times, source);

				grad[i] = (f - f1) / (2 * dx);
				if (grad[i] != grad[i]) {
					printf("grad[i] = NaN\n");
					exit(0);
				}
				if (grad[i] > 0) {
					if (vector[i] - my_minVector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
				if (grad[i] < 0) {
					if (1.0 - vector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
			}
		}

		double gradNorm = 0;
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm <= 0) {
			printf("gradNorm = %g\n", gradNorm);
			continue;
		}
		if (gradNorm != gradNorm) {
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = grad[i] / gradNorm;
			}
		}


		findMinParametersAtDirection(vector, optPar, grad, energy, observedFlux, observedError, Ne, Ntimes, times, source, currentF);

		currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
		//valley second step

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < my_Nparams; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(my_minVector[i]) / 100;
				/*if (k > 0) {
					dx = min(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2);
				}*/
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, Ntimes, times, source);

				grad[i] = (f - f1) / (2 * dx);
				if (grad[i] != grad[i]) {
					printf("grad[i] = NaN\n");
					exit(0);
				}
				if (grad[i] > 0) {
					if (vector[i] - my_minVector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
				if (grad[i] < 0) {
					if (1.0 - vector[i] < 0.000001) {
						grad[i] = 0;
					}
				}
			}
		}

		gradNorm = 0;
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm <= 0) {
			printf("gradNorm <= 0\n");
			continue;
		}
		if (gradNorm != gradNorm) {
			printf("gradNorm = NaN\n");
			exit(0);
		}
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = grad[i] / gradNorm;
			}
		}

		findMinParametersAtDirection(vector, optPar, grad, energy, observedFlux, observedError, Ne, Ntimes, times, source, currentF);

		currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, Ntimes, times, source);

		for (int i = 0; i < my_Nparams; ++i) {
			valley2[i] = vector[i];
		}

		//// valley third step

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = valley1[i] - valley2[i];
			}
		}


		gradNorm = 0;
		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				gradNorm = gradNorm + grad[i] * grad[i];
			}
		}
		gradNorm = sqrt(gradNorm);
		if (gradNorm <= 0) {
			printf("gradNorm <= 0\n");
		}
		if (gradNorm != gradNorm) {
			printf("gradNorm = NaN\n");
			exit(0);
		}

		if (gradNorm > 0) {
			for (int i = 0; i < my_Nparams; ++i) {
				if (optPar[i]) {
					grad[i] = grad[i] / gradNorm;
				}
			}

			findMinParametersAtDirection(vector, optPar, grad, energy, observedFlux, observedError, Ne, Ntimes, times, source, currentF);

		}
		currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, Ntimes, times, source);

		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/
		for (int i = 0; i < my_Nparams; ++i) {
			prevVector[i] = currentVector[i];
			currentVector[i] = vector[i];
		}
		printf("optimization function = %g\n", currentF);
		for (int i = 0; i < my_Nparams; ++i) {
			printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
		}
	}
	printf("finish optimization\n");
}

GridEnumRadiationTimeOptimizer::GridEnumRadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, const int* Npoints) : RadiationTimeOptimizer(evaluator, minParameters, maxParameters, Nparams)
{
	my_Npoints = new int[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_Npoints[i] = Npoints[i];
	}

	my_points = new double* [my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_points[i] = new double[my_Npoints[i]];
		my_points[i][0] = my_minVector[i];
		my_points[i][my_Npoints[i] - 1] = 1.0;
		if (my_Npoints[i] > 1) {
			double delta = (my_points[i][my_Npoints[i] - 1] - my_points[i][0]) / (my_Npoints[i] - 1);
			for (int j = 1; j < my_Npoints[i] - 1; ++j) {
				my_points[i][j] = my_points[i][j - 1] + delta;
			}
		}
	}
}

GridEnumRadiationTimeOptimizer::~GridEnumRadiationTimeOptimizer()
{
	for (int i = 0; i < my_Nparams; ++i) {
		delete[] my_points[i];
	}
	delete[] my_points;
	delete[] my_Npoints;
}

void GridEnumRadiationTimeOptimizer::optimize(double* vector, bool* optPar, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source)
{
	double* tempVector = new double[my_Nparams];

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction(vector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
	int N = 1;
	for (int i = 0; i < my_Nparams; ++i) {
		N *= my_Npoints[i];
	}

	int* index = new int[my_Nparams];
	for (int i = 0; i < N; ++i) {
		int tempN = i;
		for (int j = 0; j < my_Nparams; ++j) {
			index[j] = tempN % my_Npoints[j];
			tempN = tempN / my_Npoints[j];
			if (optPar[j]) {
				tempVector[j] = my_points[j][index[j]];
			}
		}
		double tempF = evaluateOptimizationFunction(tempVector, energy, observedFlux, observedError, Ne, Ntimes, times, source);
		if (tempF < currentF) {
			currentF = tempF;
			for (int j = 0; j < my_Nparams; ++j) {
				if (optPar[j]) {
					vector[j] = tempVector[j];
				}
			}
		}
		printf("optimization function = %g\n", tempF);
	}

	delete[] index;
	delete[] tempVector;
}
