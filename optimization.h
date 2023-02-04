#ifndef optimization_h
#define optimization_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "inverseCompton.h"

enum ErrorScale{LINEAR, LOG};

class SynchrotronOptimizer {
protected:
	int my_Nparams;
	ErrorScale my_errorScale;
	SynchrotronEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
	SynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale);
	~SynchrotronOptimizer();
	double evaluateOptimizationFunction(const double* vector, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source) = 0;
	void optimize(double* vector, bool* optPar, double* nu, double* observedInu, int Nnu, RadiationSource* source);
};

class GradientDescentSynchrotronOptimizer:public SynchrotronOptimizer{
protected:
	int my_Niterations;

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

class GridEnumSynchrotronOptimizer : public SynchrotronOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
	GridEnumSynchrotronOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale, const int* Npoints);
	~GridEnumSynchrotronOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
};

class SynchrotronTimeOptimizer {
protected:
	int my_Nparams;
	ErrorScale my_errorScale;
	SynchrotronEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
	SynchrotronTimeOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale);
	~SynchrotronTimeOptimizer();
	double evaluateOptimizationFunction(const double* vector, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
	virtual void optimize(double* vector, bool* optPar, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source) = 0;
	void optimize(double* vector, bool* optPar, double** nu, double** observedInu, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

class GradientDescentSynchrotronTimeOptimizer : public SynchrotronTimeOptimizer {
protected:
	int my_Niterations;

	double* tempVector;
	double* tempVector1;
	double* tempVector2;
	double* prevVector;
	double* currentVector;
	double* valley1;
	double* valley2;
	double* grad;
	void findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source, const double& currentF);
public:
	GradientDescentSynchrotronTimeOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, ErrorScale errorScale);
	~GradientDescentSynchrotronTimeOptimizer();
	virtual void optimize(double* vector, bool* optPar, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

class GridEnumSynchrotronTimeOptimizer :public SynchrotronTimeOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
	GridEnumSynchrotronTimeOptimizer(SynchrotronEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale, const int* Npoints);
	~GridEnumSynchrotronTimeOptimizer();
	virtual void optimize(double* vector, bool* optPar, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

#endif
