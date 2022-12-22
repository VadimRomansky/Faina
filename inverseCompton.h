#ifndef inverse_compton_h
#define inverse_compton_h

#include "electronDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"

class InverseComptonEvaluator {
protected:
	int my_Ne;
	int my_Nmu;
	int my_Nphi;
	double my_Emin;
	double my_Emax;

	double* my_cosTheta;
	double* my_cosThetaLeft;
	double* my_dcosTheta;
	double* my_Ee;
public:
	InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax);
	~InverseComptonEvaluator();
	double evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, ElectronDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonLuminocityIsotropicFunction(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, ElectronIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonIsotropicFluxFromSource(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, RadiationSource* source);
};

//correct name?
double evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, ElectronDistribution* electronDistribution, const double& volume, const double& distance, const double& Emin, const double& Emax, const int Ne, const int Nmu, const int Nphi);

#endif