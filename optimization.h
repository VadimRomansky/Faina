#ifndef optimization_h
#define optimization_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "inverseCompton.h"
#include "radiation.h"

enum ErrorScale{LINEAR, LOG};

class RadiationOptimizer {
protected:
	int my_Nparams;
	ErrorScale my_errorScale;
    RadiationEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
    RadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale);
    virtual ~RadiationOptimizer();
	double evaluateOptimizationFunction(const double* vector, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source) = 0;
	void optimize(double* vector, bool* optPar, double* nu, double* observedInu, int Nnu, RadiationSource* source);
};

class GradientDescentRadiationOptimizer:public RadiationOptimizer{
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
    GradientDescentRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, ErrorScale errorScale);
    virtual ~GradientDescentRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
};

class GridEnumRadiationOptimizer : public RadiationOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
    GridEnumRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale, const int* Npoints);
    virtual ~GridEnumRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* nu, double* observedInu, double* observedError, int Nnu, RadiationSource* source);
};

class RadiationTimeOptimizer {
protected:
	int my_Nparams;
	ErrorScale my_errorScale;
    RadiationEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
    RadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale);
    virtual ~RadiationTimeOptimizer();
	double evaluateOptimizationFunction(const double* vector, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
	virtual void optimize(double* vector, bool* optPar, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source) = 0;
	void optimize(double* vector, bool* optPar, double** nu, double** observedInu, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

class GradientDescentRadiationTimeOptimizer : public RadiationTimeOptimizer {
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
    GradientDescentRadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, ErrorScale errorScale);
    virtual ~GradientDescentRadiationTimeOptimizer();
	virtual void optimize(double* vector, bool* optPar, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

class GridEnumRadiationTimeOptimizer :public RadiationTimeOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
    GridEnumRadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, ErrorScale errorScale, const int* Npoints);
    virtual ~GridEnumRadiationTimeOptimizer();
	virtual void optimize(double* vector, bool* optPar, double** nu, double** observedInu, double** observedError, int* Nnu, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

#endif
