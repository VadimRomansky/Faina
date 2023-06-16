#ifndef inverse_compton_h
#define inverse_compton_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

enum ComptonSolverType {ISOTROPIC_THOMSON, ISOTROPIC_JONES, ISOTROPIC_KLEIN_NISHINA, ISOTROPIC_KLEIN_NISHINA1, ANISOTROPIC_KLEIN_NISHINA};

class InverseComptonEvaluator : public RadiationEvaluator{
protected:
	ComptonSolverType my_solverType;

	int my_Nmu;
	int my_Nphi;

	double my_Ephmin;
	double my_Ephmax;

	double* my_theta;
	double* my_cosTheta;
	double* my_cosThetaLeft;
	double* my_dcosTheta;

    PhotonIsotropicDistribution* my_photonDistribution;

	double evaluateComptonFluxThomsonIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxJonesIsotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaIsotropic(const double& photonFinalEnergyi, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaIsotropic1(const double& photonFinalEnergyi, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaAnisotropic(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
public:
    InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, double Ephmin, double Ephmax, PhotonIsotropicDistribution* photonDistribution, ComptonSolverType solverType);
	~InverseComptonEvaluator();
    void resetParameters(const double *parameters, const double *normalizationUnits);

	void outputDifferentialFlux(const char* fileName);
	void outputDifferentialFluxJones(const char* fileName, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution);
	double evaluateDifferentialFlux(const double& photonFinalEnergy, const double& photonFinalCosTheta, const double& photonFinalPhi, const double& electronInitialEnergy, const double& mu_e, const double& phi_e, const double& mu_ph, const double& phi_ph);
	double evaluateDifferentialFluxJones(const double& photonFinalEnergy, const double& electronInitialEnergy, const double& photonInitialEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution);

	
	double evaluateComptonFluxKleinNishinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);

	//double evaluateComptonFlux(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source);
	double evaluateFluxFromSourceAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, RadiationSource* source);

	virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

#endif
