#ifndef optimization_h
#define optimization_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "inverseCompton.h"
#include "radiation.h"

class RadiationOptimizer {
protected:
	int my_Nparams;
    RadiationEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
    RadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams);
    virtual ~RadiationOptimizer();
	double evaluateOptimizationFunction(const double* vector, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
	virtual void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source) = 0;
	void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, int Ne, RadiationSource* source);
	void outputProfileDiagrams(const double* vector, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source, int Npoints);
	void outputOptimizedProfileDiagram(const double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source, int Npoints, int Nparam1, int Nparam2);
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
	void findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source, const double& currentF);
	double getDerivativeByCoordinate(double* vector, int direction, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
	virtual void gradientStep(int iterationNumber, const double& currentF, double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
public:
    GradientDescentRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations);
    virtual ~GradientDescentRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
};

class CoordinateRadiationOptimizer : public GradientDescentRadiationOptimizer {
protected:
	void findMinParametersAtCoordinate(double* vector, bool* optPar, const double& der, int ncoord, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source, const double& currentF);
	virtual void gradientStep(int iterationNumber, const double& currentF, double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
public:
	CoordinateRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations);
	virtual ~CoordinateRadiationOptimizer();
	//virtual void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
};

class GridEnumRadiationOptimizer : public RadiationOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
    GridEnumRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, const int* Npoints);
    virtual ~GridEnumRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
};

class CombinedRadiationOptimizer : public RadiationOptimizer {
protected:
	GradientDescentRadiationOptimizer* my_GradientOptimzer;
	GridEnumRadiationOptimizer* my_EnumOptimizer;
public:
	CombinedRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, const int* Npoints);
	//CombinedRadiationOptimizer(GridEnumRadiationOptimizer* enumOptimizer, GradientDescentRadiationOptimizer* gradientOptimizer);
	virtual ~CombinedRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
};

class RadiationTimeOptimizer {
protected:
	int my_Nparams;
    RadiationEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
public:
    RadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams);
    virtual ~RadiationTimeOptimizer();
	double evaluateOptimizationFunction(const double* vector, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source);
	virtual void optimize(double* vector, bool* optPar, double** energy, double** observedFLux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source) = 0;
	void optimize(double* vector, bool* optPar, double** energy, double** observedFlux, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source);
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
	void findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source, const double& currentF);
	double getDerivativeByCoordinate(double* vector, int direction, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source);
public:
    GradientDescentRadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations);
    virtual ~GradientDescentRadiationTimeOptimizer();
	virtual void optimize(double* vector, bool* optPar, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

class GridEnumRadiationTimeOptimizer :public RadiationTimeOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
    GridEnumRadiationTimeOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, const int* Npoints);
    virtual ~GridEnumRadiationTimeOptimizer();
	virtual void optimize(double* vector, bool* optPar, double** energy, double** observedFlux, double** observedError, int* Ne, int Ntimes, double* times, RadiationTimeDependentSource* source);
};

#endif
