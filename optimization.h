#ifndef optimization_h
#define optimization_h

#include "electronDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "inverseCompton.h"

enum ErrorScale{LINEAR, LOG};

class SynchrotronOptimizer {
protected:
	int my_Nparams;
	int my_Niterations;
	ErrorScale my_errorScale;
	SynchrotronEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
	SynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, ErrorScale errorScale);
	~SynchrotronOptimizer();
	double evaluateOptimizationFunction(const double* vector, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source) = 0;
	void optimize(double* vector, bool* optPar, double* nu, double* observedInu, int Nnu, RadiationSource* source);
};

class GradientDescentSynchrotronOptimizer:public SynchrotronOptimizer{
protected:
	double* tempVector;
	double* tempVector1;
	double* tempVector2;
	double* prevVector;
	double* currentVector;
	double* valley1;
	double* valley2;
	double* grad;
	void findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source, const double& currentF);
public:
	GradientDescentSynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, ErrorScale errorScale);
	~GradientDescentSynchrotronOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
};

#endif
