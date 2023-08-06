#ifndef inverse_compton_h
#define inverse_compton_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

//enum ComptonSolverType {ISOTROPIC_THOMSON, ISOTROPIC_JONES, ISOTROPIC_KLEIN_NISHINA2, ISOTROPIC_KLEIN_NISHINA, ANISOTROPIC_KLEIN_NISHINA};
enum ComptonSolverType {ISOTROPIC_THOMSON, ISOTROPIC_JONES, ISOTROPIC_KLEIN_NISHINA, ANISOTROPIC_KLEIN_NISHINA, ANISOTROPIC_KLEIN_NISHINA2, ANISOTROPIC_KLEIN_NISHINA3};

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

    PhotonDistribution* my_photonDistribution;

	double evaluateComptonFluxThomsonIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxJonesIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	//double evaluateComptonFluxKleinNishinaIsotropic2(const double& photonFinalEnergyi, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaIsotropic(const double& photonFinalEnergyi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaAnisotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	//double evaluateComptonFluxKleinNishinaAnisotropic2(const double& photonFinalEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution, const double& volume, const double& distance);

	virtual PhotonDistribution* getPhotonDistribution(const double& rho, const double& z, const double& phi);
public:
    InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, double Ephmin, double Ephmax, PhotonDistribution* photonDistribution, ComptonSolverType solverType);
	virtual ~InverseComptonEvaluator();
    void resetParameters(const double *parameters, const double *normalizationUnits);

	void outputDifferentialFlux(const char* fileName);
	void outputDifferentialFluxJones(const char* fileName, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution);
	double evaluateDifferentialFlux(const double& photonFinalEnergy, const double& photonFinalCosTheta, const double& photonFinalPhi, const double& electronInitialEnergy, const double& mu_e, const double& phi_e, const double& mu_ph, const double& phi_ph);
	double evaluateDifferentialFluxJones(const double& photonFinalEnergy, const double& electronInitialEnergy, const double& photonInitialEnergy, PhotonIsotropicDistribution* photonDistribution, MassiveParticleIsotropicDistribution* electronDistribution);

	
	double evaluateComptonFluxKleinNishinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaAnisotropic2(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateComptonFluxKleinNishinaAnisotropic3(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);

	//double evaluateComptonFlux(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
    double evaluateFluxFromFunction(const double& photonFinalEnergy, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance, PhotonDistribution* photonDistribution);
	double evaluateFluxFromSourceAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, RadiationSource* source, ComptonSolverType solverType);

	virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

class InverseComptonEvaluatorWithSource : public InverseComptonEvaluator {
protected:
	double my_sourceR;
	double my_sourceZ;
	double my_sourcePhi;
	double my_defaultPhotonConcentration;
	virtual PhotonDistribution* getPhotonDistribution(const double& rho, const double& z, const double& phi);
public:
	InverseComptonEvaluatorWithSource(int Ne, int Nmu, int Nphi, double Emin, double Emax, double Ephmin, double Ephmax, PhotonDistribution* photonDistribution, ComptonSolverType solverType, const double& sourceR, const double& sourceZ, const double& sourcePhi);
	virtual ~InverseComptonEvaluatorWithSource();
};

#endif
