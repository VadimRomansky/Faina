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


RadiationOptimizer::RadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, KPIevaluator* KPIevaluator) {
	my_evaluator = evaluator;
	my_KPIevaluator = KPIevaluator;
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

double RadiationOptimizer::evaluateOptimizationFunction(const double* vector) {
	return my_KPIevaluator->evaluate(vector, my_maxParameters, my_evaluator);
}


void RadiationOptimizer::outputProfileDiagrams(const double* vector, int Npoints)
{
	if (my_Nparams < 2) {
		printf("cannot output error profilefor Nparams < 2\n");
		printLog("cannot output error profilefor Nparams < 2\n");
	}

	double** paramPoints = new double* [my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		paramPoints[i] = new double[Npoints];
		paramPoints[i][0] = my_minParameters[i];
		paramPoints[i][Npoints-1] = my_maxParameters[i];
		double factor = pow(my_maxParameters[i] / my_minParameters[i], 1.0 / (Npoints - 1));
		for (int j = 1; j < Npoints-1; ++j) {
			paramPoints[i][j] = paramPoints[i][j - 1] * factor;
		}
	}

	for (int i = 0; i < my_Nparams; ++i) {
		std::string fileNumber = convertIntToString(i);
		std::string fileName = "parameter";
		FILE* paramFile = fopen((fileName + fileNumber + ".dat").c_str(), "w");
		for (int j = 0; j < Npoints; ++j) {
			fprintf(paramFile, "%g\n", paramPoints[i][j]);
		}
		fclose(paramFile);
	}

	double* tempVector = new double[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		tempVector[i] = vector[i];
	}

	for (int ip1 = 0; ip1 < my_Nparams-1; ++ip1) {
		for (int ip2 = ip1 + 1; ip2 < my_Nparams; ++ip2) {
			std::string fileNumber1 = convertIntToString(ip1);
			std::string fileNumber2 = convertIntToString(ip2);
			std::string fileName = "error_";
			FILE* errorFile = fopen((fileName + fileNumber1 +"_" + fileNumber2 + ".dat").c_str(), "w");

			for (int i = 0; i < my_Nparams; ++i) {
				tempVector[i] = vector[i];
			}

			for (int j1 = 0; j1 < Npoints; ++j1) {
				for (int j2 = 0; j2 < Npoints; ++j2) {
					tempVector[ip1] = paramPoints[ip1][j1]/my_maxParameters[ip1];
					tempVector[ip2] = paramPoints[ip2][j2]/my_maxParameters[ip2];

					double error = evaluateOptimizationFunction(tempVector);

					fprintf(errorFile, "%g", error);
					if (j2 < Npoints - 1) {
						fprintf(errorFile, " ");
					}
				}
				fprintf(errorFile, "\n");
			}

			fclose(errorFile);
		}
	}

	delete[] tempVector;

	for (int i = 0; i < my_Nparams; ++i) {
		delete[] paramPoints[i];
	}
	delete[] paramPoints;
}

void RadiationOptimizer::outputOptimizedProfileDiagram(const double* vector, bool* optPar, int Npoints, int Nparam1, int Nparam2)
{
	double** paramPoints = new double* [my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		paramPoints[i] = new double[Npoints];
		paramPoints[i][0] = my_minParameters[i];
		paramPoints[i][Npoints - 1] = my_maxParameters[i];
		double factor = pow(my_maxParameters[i] / my_minParameters[i], 1.0 / (Npoints - 1));
		for (int j = 1; j < Npoints - 1; ++j) {
			paramPoints[i][j] = paramPoints[i][j - 1] * factor;
		}
	}

	for (int i = 0; i < my_Nparams; ++i) {
		std::string fileNumber = convertIntToString(i);
		std::string fileName = "parameter";
		FILE* paramFile = fopen((fileName + fileNumber + ".dat").c_str(), "w");
		for (int j = 0; j < Npoints; ++j) {
			fprintf(paramFile, "%g\n", paramPoints[i][j]);
		}
		fclose(paramFile);
	}

	double* tempVector = new double[my_Nparams];
	bool* tempOptPar = new bool[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		tempVector[i] = vector[i];
		tempOptPar[i] = optPar[i];
	}
	tempOptPar[Nparam1] = false;
	tempOptPar[Nparam2] = false;


			std::string fileNumber1 = convertIntToString(Nparam1);
			std::string fileNumber2 = convertIntToString(Nparam2);
			std::string fileName = "error_optimized_";
			FILE* errorFile = fopen((fileName + fileNumber1 + "_" + fileNumber2 + ".dat").c_str(), "w");

			for (int i = 0; i < my_Nparams; ++i) {
				tempVector[i] = vector[i];
			}

			for (int j1 = 0; j1 < Npoints; ++j1) {
				for (int j2 = 0; j2 < Npoints; ++j2) {
					printf("outputOptimizedProfile j1 = %d j2 = %d\n", j1, j2);
					printLog("outputOptimizedProfile j1 = %d j2 = %d\n", j1, j2);
					for (int i = 0; i < my_Nparams; ++i) {
						tempVector[i] = vector[i];
					}
					tempVector[Nparam1] = paramPoints[Nparam1][j1] / my_maxParameters[Nparam1];
					tempVector[Nparam2] = paramPoints[Nparam2][j2] / my_maxParameters[Nparam2];

					optimize(tempVector, tempOptPar);

					double error = evaluateOptimizationFunction(tempVector);

					fprintf(errorFile, "%g", error);
					if (j2 < Npoints - 1) {
						fprintf(errorFile, " ");
					}
				}
				fprintf(errorFile, "\n");
			}

			fclose(errorFile);


	delete[] tempVector;

	for (int i = 0; i < my_Nparams; ++i) {
		delete[] paramPoints[i];
	}
	delete[] paramPoints;
}

GradientDescentRadiationOptimizer::GradientDescentRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, KPIevaluator* KPIevaluator) : RadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, KPIevaluator) {
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

void GradientDescentRadiationOptimizer::findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, const double& currentF) {
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
				printLog("tempvector[i] != tempVector[i]\n");
				exit(0);
			}
		}
	}

	double f = evaluateOptimizationFunction(tempVector1);
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
						printLog("tempvector[i] != tempVector[i]\n");
						exit(0);
					}
				}
			}

			f = evaluateOptimizationFunction(tempVector1);
		}

		if (f < currentF) {
			maxLambda = 2 * step;
			//todo!
			/*for (int i = 0; i < my_Nparams; ++i) {
				vector[i] = tempVector1[i];
			}*/
			//return;
		}
		else {
			return;
		}
	}

	f = evaluateOptimizationFunction(vector);

	/*maxLambda = sqrt(1.0 * Npar);
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - my_minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}*/

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
	double f3 = evaluateOptimizationFunction(tempVector1);

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

		double f1 = evaluateOptimizationFunction(tempVector1);
		double f2 = evaluateOptimizationFunction(tempVector2);

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

	double tempF = evaluateOptimizationFunction(tempVector1);

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

double GradientDescentRadiationOptimizer::getDerivativeByCoordinate(double* vector, int direction)
{
	for (int j = 0; j < my_Nparams; ++j) {
		tempVector[j] = vector[j];
	}
	double dx = fabs(vector[direction]) / 100;

	tempVector[direction] = vector[direction] + dx;
	double f = evaluateOptimizationFunction(tempVector);
	tempVector[direction] = vector[direction] - dx;
	double f1 = evaluateOptimizationFunction(tempVector);

	double der = (f - f1) / (2 * dx);
	return der;
}

void GradientDescentRadiationOptimizer::gradientStep(int iterationNumber, const double& currentF, double* vector, bool* optPar) {
	for (int i = 0; i < my_Nparams; ++i) {
		grad[i] = 0;
	}

	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			for (int j = 0; j < my_Nparams; ++j) {
				tempVector[j] = vector[j];
			}
			double dx = fabs(vector[i]) / 1000000;
			if (iterationNumber > 0) {
				dx = min(dx, max(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2));
			}
			tempVector[i] = vector[i] + dx;
			double f = evaluateOptimizationFunction(tempVector);
			tempVector[i] = vector[i] - dx;
			double f1 = evaluateOptimizationFunction(tempVector);

			grad[i] = (f - f1) / (2 * dx);
			if (grad[i] != grad[i]) {
				printf("grad[i] = NaN\n");
				printLog("grad[i] = NaN\n");
				exit(0);
			}
			if (grad[i] > 0) {
				if (vector[i] - my_minVector[i] < 0.000001* my_minVector[i]) {
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
		return;
	}
	if (gradNorm != gradNorm) {
		printf("gradNorm = NaN\n");
		printLog("gradNorm = NaN\n");
		exit(0);
	}
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			grad[i] = grad[i] / gradNorm;
		}
	}


	findMinParametersAtDirection(vector, optPar, grad, currentF);

}

void GradientDescentRadiationOptimizer::optimize(double* vector, bool* optPar) {
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
	double currentF = evaluateOptimizationFunction(vector);
	printf("optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}
	
	for (int k = 0; k < my_Niterations; ++k) {
		///randomization;
		for (int j = 0; j < 5; ++j) {
			for(int i = 0; i < my_Nparams; ++i) {
				tempVector[i] = vector[i];
				if(optPar[i]){
					tempVector[i] = my_minVector[i] + (1.0 - my_minVector[i])*uniformDistribution();
					//tempVector[i] = vector[i] + 0.2*my_minVector[i]*(0.5 - uniformDistribution());
					if(tempVector[i] > 1.0) {
						tempVector[i] = 1.0;
					}
					if(tempVector[i] < my_minVector[i]) {
						tempVector[i] = my_minVector[i];
					}
				}
			}

			double tempF = evaluateOptimizationFunction(tempVector);
			if(tempF < currentF){
				currentF = tempF;
				for(int i = 0; i < my_Nparams; ++i) {
					vector[i] = tempVector[i];
				}
				printf("optimization function = %g\n", currentF);
				printLog("optimization function = %g\n", currentF);
				printf("random search succeed\n");
				printLog("random search succeed\n");
			}
		}
		double prevF = currentF;
		//
		//valley first step
		for (int i = 0; i < my_Nparams; ++i) {
			valley1[i] = vector[i];
			grad[i] = 0;
		}

		printf("optimization k = %d\n", k);
		printLog("optimization k = %d\n", k);

		gradientStep(k, currentF, vector, optPar);

		currentF = evaluateOptimizationFunction(vector);
		//valley second step
		gradientStep(k, currentF, vector, optPar);

		currentF = evaluateOptimizationFunction(vector);

		for (int i = 0; i < my_Nparams; ++i) {
			valley2[i] = vector[i];
			grad[i] = 0;
		}

		//// valley third step

		/*for (int i = 0; i < my_Nparams; ++i) {
			if (optPar[i]) {
				grad[i] = valley1[i] - valley2[i];
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

		}*/
		currentF = evaluateOptimizationFunction(vector);

		/*if(fabs(currentF - prevF) < 0.00000001){
			break;
		}*/
		for (int i = 0; i < my_Nparams; ++i) {
			prevVector[i] = currentVector[i];
			currentVector[i] = vector[i];
		}
		printf("optimization function = %g\n", currentF);
		printLog("optimization function = %g\n", currentF);
		for (int i = 0; i < my_Nparams; ++i) {
			printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
			printLog("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
		}
	}
	printf("GradientDescentOptimizer: final optimization function = %g\n", currentF);
	printLog("GradientDescentOptimizer: final optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
		printLog("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}
	printf("finish optimization\n");
	printLog("finish optimization\n");
}

GridEnumRadiationOptimizer::GridEnumRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, const int* Npoints, KPIevaluator* KPIevaluator) : RadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, KPIevaluator)
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

void GridEnumRadiationOptimizer::optimize(double* vector, bool* optPar)
{
	double* tempVector = new double[my_Nparams];

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector[i] = vector[i];
	}
	double currentF = evaluateOptimizationFunction(vector);
	int N = 1;
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			N *= my_Npoints[i];
		}
	}

	int* index = new int[my_Nparams];
	for (int i = 0; i < N; ++i) {
		int tempN = i;
		for (int j = 0; j < my_Nparams; ++j) {
			if (optPar[j]) {
				index[j] = tempN % my_Npoints[j];
				tempN = tempN / my_Npoints[j];
				tempVector[j] = my_points[j][index[j]];
			}
		}
		double tempF = evaluateOptimizationFunction(tempVector);
		if (tempF < currentF) {
			currentF = tempF;
			for (int j = 0; j < my_Nparams; ++j) {
				if (optPar[j]) {
					vector[j] = tempVector[j];
				}
			}
			printf("optimization function = %g\n", tempF);
			printLog("optimization function = %g\n", tempF);
		}
	}

	printf("GridENumOptimizer: final optimization function = %g\n", currentF);
	printLog("GridENumOptimizer: final optimization function = %g\n", currentF);
	for (int i = 0; i < my_Nparams; ++i) {
		printf("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
		printLog("parameter[%d] = %g\n", i, vector[i] * my_maxParameters[i]);
	}
	printf("finish optimization\n");
	printLog("finish optimization\n");

	delete[] index;
	delete[] tempVector;
}

void CoordinateRadiationOptimizer::findMinParametersAtCoordinate(double* vector, bool* optPar, const double& der, int ncoord, const double& currentF)
{
	if (!optPar[ncoord]) {
		return;
	}

	double maxLambda = 1.0;
	if (der > 0) {
		maxLambda = min(maxLambda, fabs((vector[ncoord] - my_minVector[ncoord]) / der));
	}
	if (der < 0) {
		maxLambda = min(maxLambda, fabs((1.0 - vector[ncoord]) / der));
	}


	double step = maxLambda / 2;

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}

	tempVector1[ncoord] = vector[ncoord] - der * step;
	if (tempVector1[ncoord] != tempVector1[ncoord]) {
		printf("tempvector[i] != tempVector[i]\n");
		printLog("tempvector[i] != tempVector[i]\n");
		exit(0);
	}


	double f = evaluateOptimizationFunction(tempVector1);
	if (f > currentF) {
		int count = 0;
		while (f > currentF && count < 20) {
			count++;
			step = step / 2;

			for (int i = 0; i < my_Nparams; ++i) {
				tempVector1[i] = vector[i];
			}

			tempVector1[ncoord] = vector[ncoord] - der * step;
			if (tempVector1[ncoord] != tempVector1[ncoord]) {
				printf("tempvector[i] != tempVector[i]\n");
				exit(0);
			}

			f = evaluateOptimizationFunction(tempVector1);
		}

		if (f < currentF) {
			maxLambda = 2 * step;
			//todo!
			/*for (int i = 0; i < my_Nparams; ++i) {
				vector[i] = tempVector1[i];
			}*/
			//return;
		}
		else {
			return;
		}
	}

	f = evaluateOptimizationFunction(vector);

	/*maxLambda = sqrt(1.0 * Npar);
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			if (grad[i] > 0) {
				maxLambda = min(maxLambda, fabs((vector[i] - my_minVector[i]) / grad[i]));
			}
			if (grad[i] < 0) {
				maxLambda = min(maxLambda, fabs((1.0 - vector[i]) / grad[i]));
			}
		}
	}*/

	for (int i = 0; i < my_Nparams; ++i) {
		tempVector1[i] = vector[i];
	}

	tempVector1[ncoord] = vector[ncoord] - maxLambda * der;
	if (tempVector1[ncoord] > 1.0) {
		tempVector1[ncoord] = 1.0;
	}
	if (tempVector1[ncoord] < my_minVector[ncoord]) {
		tempVector1[ncoord] = my_minVector[ncoord];
	}

	double f3 = evaluateOptimizationFunction(tempVector1);

	double leftLambda = 0;
	double rightLambda = maxLambda;


	for (int j = 0; j < 10; ++j) {
		double lambda1 = leftLambda + (rightLambda - leftLambda) / 3.0;
		double lambda2 = leftLambda + (rightLambda - leftLambda) * 2.0 / 3.0;

		for (int i = 0; i < my_Nparams; ++i) {
			tempVector1[i] = vector[i];
			tempVector2[i] = vector[i];
		}

		tempVector1[ncoord] = vector[ncoord] - lambda1 * der;
		tempVector2[ncoord] = vector[ncoord] - lambda2 * der;

		double f1 = evaluateOptimizationFunction(tempVector1);
		double f2 = evaluateOptimizationFunction(tempVector2);

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

	tempVector1[ncoord] = vector[ncoord] - lambda * der;
	if (tempVector1[ncoord] > 1.0) {
		tempVector1[ncoord] = 1.0;
	}
	if (tempVector1[ncoord] < my_minVector[ncoord]) {
		tempVector1[ncoord] = my_minVector[ncoord];
	}

	double tempF = evaluateOptimizationFunction(tempVector1);

	if (f < tempF) {
		if (f < f3) {
			return;
		}
		else {
			vector[ncoord] = vector[ncoord] - maxLambda * der;
		}
	}
	else {
		if (tempF < f3) {
			vector[ncoord] = vector[ncoord] - lambda * der;
		}
		else {
			vector[ncoord] = vector[ncoord] - maxLambda * der;
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

void CoordinateRadiationOptimizer::gradientStep(int iterationNumber, const double& currentF, double* vector, bool* optPar)
{
	for (int i = 0; i < my_Nparams; ++i) {
		if (optPar[i]) {
			for (int j = 0; j < my_Nparams; ++j) {
				tempVector[j] = vector[j];
			}
			double dx = fabs(vector[i]) / 1000000;
			if (iterationNumber > 0) {
				dx = min(dx, max(max(0.000001, fabs(currentVector[i] - prevVector[i]) / 10), fabs(currentVector[i] - my_minVector[i]) / 2));
			}
			tempVector[i] = vector[i] + dx;
			double f = evaluateOptimizationFunction(tempVector);
			tempVector[i] = vector[i] - dx;
			double f1 = evaluateOptimizationFunction(tempVector);

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


			double gradNorm = fabs(grad[i]);

			if (gradNorm <= 0) {
				printf("gradNorm = %g\n", gradNorm);
				continue;
			}
			if (gradNorm != gradNorm) {
				printf("gradNorm = NaN\n");
				printLog("gradNorm = NaN\n");
				exit(0);
			}

			grad[i] = grad[i] / gradNorm;


			findMinParametersAtCoordinate(vector, optPar, grad[i], i, currentF);
		}
	}
}

CoordinateRadiationOptimizer::CoordinateRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, KPIevaluator* KPIevaluator) : GradientDescentRadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator) {

}
CoordinateRadiationOptimizer::~CoordinateRadiationOptimizer(){

}

CombinedRadiationOptimizer::CombinedRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, const int* Npoints, KPIevaluator* KPIevaluator) : RadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, KPIevaluator)
{
	my_EnumOptimizer = new GridEnumRadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	my_GradientOptimzer = new GradientDescentRadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);

}

/*CombinedRadiationOptimizer::CombinedRadiationOptimizer(GridEnumRadiationOptimizer* enumOptimizer, GradientDescentRadiationOptimizer* gradientOptimizer)
{
	my_EnumOptimizer = enumOptimizer;
	my_GradientOptimzer = gradientOptimizer;

	my_evaluator = enumOptimizer->my_evaluator;
	my_Nparams = Nparams;
	my_minParameters = new double[my_Nparams];
	my_maxParameters = new double[my_Nparams];
	my_minVector = new double[my_Nparams];
	for (int i = 0; i < my_Nparams; ++i) {
		my_minParameters[i] = minParameters[i];
		my_maxParameters[i] = maxParameters[i];
		my_minVector[i] = my_minParameters[i] / my_maxParameters[i];
	}
}*/

CombinedRadiationOptimizer::~CombinedRadiationOptimizer()
{
	delete my_EnumOptimizer;
	delete my_GradientOptimzer;
}

void CombinedRadiationOptimizer::optimize(double* vector, bool* optPar)
{
	my_EnumOptimizer->optimize(vector, optPar);
	my_GradientOptimzer->optimize(vector, optPar);
}
