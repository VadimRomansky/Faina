#ifndef optimization_h
#define optimization_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "inverseCompton.h"
#include "radiation.h"
#include "KPIevaluator.h"

class RadiationOptimizer {
protected:
	int my_Nparams;
    RadiationEvaluator* my_evaluator;
	double* my_minParameters;
	double* my_maxParameters;

	double* my_minVector;
	LossEvaluator* my_lossEvaluator;
public:
    RadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, LossEvaluator* KPIevaluator);
    virtual ~RadiationOptimizer();
	virtual double evaluateOptimizationFunction(const double* vector);
	virtual void optimize(double* vector, bool* optPar) = 0;
	//void optimize(double* vector, bool* optPar, RadiationSource* source);
	void outputProfileDiagrams(const double* vector, int Npoints);
	void outputOptimizedProfileDiagram(const double* vector, bool* optPar, int Npoints, int Nparam1, int Nparam2);
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
	void findMinParametersAtDirection(double* vector, bool* optPar, const double* grad, const double& currentF);
	double getDerivativeByCoordinate(double* vector, int direction);
	virtual void gradientStep(int iterationNumber, const double& currentF, double* vector, bool* optPar);
public:
    GradientDescentRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, LossEvaluator* lossEvaluator);
    virtual ~GradientDescentRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar);
};

class CoordinateRadiationOptimizer : public GradientDescentRadiationOptimizer {
protected:
	void findMinParametersAtCoordinate(double* vector, bool* optPar, const double& der, int ncoord, const double& currentF);
	virtual void gradientStep(int iterationNumber, const double& currentF, double* vector, bool* optPar);
public:
	CoordinateRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, LossEvaluator* lossEvaluator);
	virtual ~CoordinateRadiationOptimizer();
	//virtual void optimize(double* vector, bool* optPar, double* energy, double* observedFlux, double* observedError, int Ne, RadiationSource* source);
};

class GridEnumRadiationOptimizer : public RadiationOptimizer {
protected:
	int* my_Npoints;
	double** my_points;
public:
    GridEnumRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, const int* Npoints, LossEvaluator* lossEvaluator);
    virtual ~GridEnumRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar);
};

class CombinedRadiationOptimizer : public RadiationOptimizer {
protected:
	GradientDescentRadiationOptimizer* my_GradientOptimzer;
	GridEnumRadiationOptimizer* my_EnumOptimizer;
public:
	CombinedRadiationOptimizer(RadiationEvaluator* evaluator, const double* minParameters, const double* maxParameters, int Nparams, int Niterations, const int* Npoints, LossEvaluator* lossEvaluator);
	//CombinedRadiationOptimizer(GridEnumRadiationOptimizer* enumOptimizer, GradientDescentRadiationOptimizer* gradientOptimizer);
	virtual ~CombinedRadiationOptimizer();
	virtual void optimize(double* vector, bool* optPar);
};

#endif
