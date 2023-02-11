#ifndef RADIATION_H
#define RADIATION_H

#include "math.h"

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"

class RadiationEvaluator{
protected:
    int my_Ne;
    double my_Emin;
    double my_Emax;

    double* my_Ee;
public:
    RadiationEvaluator(int Ne, const double& Emin, const double& Emax){
        my_Ne = Ne;
        my_Emin = Emin;
        my_Emax = Emax;
        my_Ee = new double[my_Ne];
        double factor = pow(my_Emax / my_Emin, 1.0 / (my_Ne - 1));

        my_Ee[0] = my_Emin;
        for (int i = 1; i < my_Ne; ++i) {
            my_Ee[i] = my_Ee[i - 1] * factor;
        }
    }
    virtual ~RadiationEvaluator(){
        delete[] my_Ee;
    }
    virtual double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance) = 0;
    virtual double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source) = 0;
};

#endif // RADIATION_H
