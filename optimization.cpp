#include "stdio.h"
#include "math.h"

#include "util.h"
#include "constants.h"
#include "radiationSource.h"
#include "electronDistribution.h"
#include "photonDistribution.h"
#include "synchrotron.h"
#include "inverseCompton.h"

#include "optimization.h"


SynchrotronOptimizer::SynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale) {
	my_evaluator = evaluator;
	my_Nparams = Nparams;
	my_errorScale = errorScale;
	my_minParameters = new double[my_Nparams];
	my_maxParameters = new double[my_Nparams];
	my_minVector = new double[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_minParameters[i] = minParameters[i];
		my_maxParameters[i] = maxParameters[i];
		my_minVector[i] = my_minParameters[i] / my_maxParameters[i];
	}
}
SynchrotronOptimizer::~SynchrotronOptimizer() {
	delete[] my_minParameters;
	delete[] my_maxParameters;
	delete[] my_minVector;
}

double SynchrotronOptimizer::evaluateOptimizationFunction(const double* vector, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source) {
	double* totalInu = new double[Nnu];
	
	source->resetParameters(vector, my_maxParameters);

	for (int i = 0; i < Nnu; ++i) {
		//1E26 from Jansky
		totalInu[i] = my_evaluator->evaluateSynchrotronFluxFromSource(source, nu[i])*1E26;
	}
	double err = 0;
	for (int j = 0; j < Nnu; ++j) {
		double err1 = 0;
		if (my_errorScale == ErrorScale::LINEAR) {
			err1 = sqr(totalInu[j] - observedInu[j]) / sqr(observedError[j]);
		}
		else {
			//todo
			err1 = sqr(log(totalInu[j]) - log(observedInu[j])) / sqr(observedError[j] / observedInu[j]);
		}
		err = err + err1;
	}

	delete[] totalInu;

	return err;
}

void SynchrotronOptimizer::optimize(double* vector, bool* optPar, double* nu, double* observedInu, int Nnu, RadiationSource* source) {
	double* observedError = new double[Nnu];
	for (int i = 0; i < Nnu; ++i) {
		observedError[i] = 1.0;
	}

	optimize(vector, optPar, nu, observedInu, observedError, Nnu, source);

	delete[] observedError;
}

GradientDescentSynchrotronOptimizer::GradientDescentSynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, ErrorScale errorScale) :SynchrotronOptimizer(evaluator, minParameters, maxParameters, Nparams, errorScale) {
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

GradientDescentSynchrotronOptimizer::~GradientDescentSynchrotronOptimizer() {
	delete[] tempVector;
	delete[] tempVector1;
	delete[] tempVector2;
	delete[] prevVector;
	delete[] currentVector;
	delete[] valley1;
	delete[] valley2;
	delete[] grad;
}

void GradientDescentSynchrotronOptimizer::findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source, const double& currentF) {
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

	double f = evaluateOptimizationFunction(tempVector1, nu, observedInu, observedError, Nnu, source);
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

			f = evaluateOptimizationFunction(tempVector1, nu, observedInu, observedError, Nnu, source);
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

	f = evaluateOptimizationFunction(vector, nu, observedInu, observedError, Nnu, source);

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
	double f3 = evaluateOptimizationFunction(tempVector1, nu, observedInu, observedError, Nnu, source);

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

		double f1 = evaluateOptimizationFunction(tempVector1, nu, observedInu, observedError, Nnu, source);
		double f2 = evaluateOptimizationFunction(tempVector2, nu, observedInu, observedError, Nnu, source);

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

	double tempF = evaluateOptimizationFunction(tempVector1, nu, observedInu, observedError, Nnu, source);

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

void GradientDescentSynchrotronOptimizer::optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source) {
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
	double currentF = evaluateOptimizationFunction(vector, nu, observedInu, observedError, Nnu, source);
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

			double tempF = evaluateOptimizationFunction(tempVector, nu, observedInu, observedError, Nnu, source);
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
					dx = min(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2);
				}
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction(tempVector, nu, observedInu, observedError, Nnu, source);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction(tempVector, nu, observedInu, observedError, Nnu, source);

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


		findMinParametersAtDirection(vector, optPar, grad, nu, observedInu, observedError, Nnu, source, currentF);

		currentF = evaluateOptimizationFunction(vector, nu, observedInu, observedError, Nnu, source);
		//valley second step

		for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				for (int j = 0; j < my_Nparams; ++j) {
					tempVector[j] = vector[j];
				}
				double dx = fabs(vector[i]) / 100;
				if (k > 0) {
					dx = min(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2);
				}
				tempVector[i] = vector[i] + dx;
				double f = evaluateOptimizationFunction(tempVector, nu, observedInu, observedError, Nnu, source);
				tempVector[i] = vector[i] - dx;
				double f1 = evaluateOptimizationFunction(tempVector, nu, observedInu, observedError, Nnu, source);

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

		findMinParametersAtDirection(vector, optPar, grad, nu, observedInu, observedError, Nnu, source, currentF);

		currentF = evaluateOptimizationFunction(vector, nu, observedInu, observedError, Nnu, source);

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

			findMinParametersAtDirection(vector, optPar, grad, nu, observedInu, observedError, Nnu, source, currentF);

		}
		currentF = evaluateOptimizationFunction(vector, nu, observedInu, observedError, Nnu, source);

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

EnumSynchrotronOptimizer::EnumSynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale, const int* Npoints) : SynchrotronOptimizer(evaluator, minParameters, maxParameters, Nparams, errorScale)
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

EnumSynchrotronOptimizer::~EnumSynchrotronOptimizer()
{
	for (int i = 0; i < my_Nparams; ++i) {
		delete[] my_points[i];
	}
	delete[] my_points;
	delete[] my_Npoints;
}

void EnumSynchrotronOptimizer::optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source)
{
	double* tempVector = new double[my_Nparams];

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction(vector, nu, observedInu, observedError, Nnu, source);
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
			tempVector[j] = my_points[j][index[j]];
		}
		double tempF = evaluateOptimizationFunction(tempVector, nu, observedInu, observedError, Nnu, source);
		if (tempF < currentF) {
			currentF = tempF;
			for (int j = 0; j < my_Nparams; ++j) {
				vector[j] = tempVector[j];
			}
		}
		printf("optimization function = %g\n", tempF);
	}

	delete[] index;
	delete[] tempVector;
}
