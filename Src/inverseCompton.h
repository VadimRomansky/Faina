#ifndef inverse_compton_h
#define inverse_compton_h

#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "radiationSource.h"
#include "radiation.h"

enum ComptonSolverType {ISOTROPIC_THOMSON, ISOTROPIC_JONES, ISOTROPIC_KLEIN_NISHINA2, ISOTROPIC_KLEIN_NISHINA, ANISOTROPIC_KLEIN_NISHINA};
//enum ComptonSolverType {ISOTROPIC_THOMSON, ISOTROPIC_JONES, ISOTROPIC_KLEIN_NISHINA, ANISOTROPIC_KLEIN_NISHINA, ANISOTROPIC_KLEIN_NISHINA2, ANISOTROPIC_KLEIN_NISHINA3};

class InverseComptonEvaluator : public RadiationEvaluator{
protected:
	ComptonSolverType my_solverType;

	int my_Nmu;
	int my_Nphi;

	int my_Nph;

	double my_Ephmin;
	double my_Ephmax;

	double* my_theta;
	double* my_cosTheta;
	double* my_cosThetaLeft;
	double* my_dcosTheta;

    PhotonDistribution* my_photonDistribution;
	double my_photonConcentration;

	double evaluateComptonEmissivityThomsonIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
	double evaluateComptonEmissivityJonesIsotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
	double evaluateComptonFluxKleinNishinaIsotropic(const double& photonFinalEnergyi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
	double evaluateComptonEmissivityKleinNishinaIsotropic2(const double& photonFinalEnergyi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
	double evaluateComptonEmissivityKleinNishinaAnisotropic(const double& photonFinalEnergy, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);


	virtual PhotonDistribution* getPhotonDistribution(const double& rho, const double& z, const double& phi);
	virtual double getPhotonConcentration(const double& rho, const double& z, const double& phi);

	double evaluateComptonEmissivityKleinNishinaAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
	double evaluateComptonEmissivityKleinNishinaAnisotropic2(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
	double evaluateComptonEmissivityKleinNishinaAnisotropic3(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& photonConcentration, const double& electronConcentration);
public:
    InverseComptonEvaluator(int Ne, int Nmu, int Nphi, double Emin, double Emax, int Nph, double Ephmin, double Ephmax, PhotonDistribution* photonDistribution, const double& photonConcentration, ComptonSolverType solverType, bool absorption = false, bool doppler = false);
	virtual ~InverseComptonEvaluator();
    void resetParameters(const double *parameters, const double *normalizationUnits);

	//double evaluateComptonFlux(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, MassiveParticleDistribution* electronDistribution, const double& volume, const double& distance);
	double evaluateFluxFromSourceAnisotropic(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, const double& photonConcentration, RadiationSourceInCylindrical* source, ComptonSolverType solverType);
	double evaluateTotalFluxInEnergyRangeAnisotropic(const double& Ephmin, const double& Ephmax, const double& photonFinalTheta, const double& photonFinalPhi, int Nph, PhotonDistribution* photonDistribution, const double& photonConcentration, RadiationSourceInCylindrical* source, ComptonSolverType solverType);

	virtual double evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);

	virtual double evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source);

	//virtual double evaluateFluxFromSourceAtPoint(const double& photonFinalEnergy, RadiationSource* source, int irho, int iphi);
};

class InverseComptonEvaluatorWithSource : public InverseComptonEvaluator {
protected:
	double my_sourceR;
	double my_sourceZ;
	double my_sourcePhi;
	double my_defaultPhotonConcentration;
	virtual PhotonDistribution* getPhotonDistribution(const double& rho, const double& z, const double& phi);
	virtual double getPhotonConcentration(const double& rho, const double& z, const double& phi);
public:
	InverseComptonEvaluatorWithSource(int Ne, int Nmu, int Nphi, double Emin, double Emax, int Nph, double Ephmin, double Ephmax, PhotonDistribution* photonDistribution, const double& photonConcentration, ComptonSolverType solverType, const double& sourceR, const double& sourceZ, const double& sourcePhi, bool absorption = false, bool doppler = false);
	virtual ~InverseComptonEvaluatorWithSource();
};

#endif
