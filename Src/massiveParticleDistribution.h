#ifndef electron_distribution_h
#define electron_distribution_h
#include "particleDistribution.h"

/** abstract class, containing electron energy distribution function, normalized to the concentration
* number of particles dN = f(E, mu, phi) dE dmu dphi dV where mu = cos theta
*/

enum DistributionInputType {ENERGY_FE, ENERGY_KIN_FE, GAMMA_FGAMMA, GAMMA_KIN_FGAMMA, MOMENTUM_FP};

class MassiveParticleDistribution : public ParticleDistribution{
protected:
	double my_mass;
public:
    virtual ~MassiveParticleDistribution(){

    };
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi) = 0;
	virtual void resetConcentration(const double& concentration) = 0;
	double getMass() {
		return my_mass;
	}
};

class MassiveParticleIsotropicDistribution : public MassiveParticleDistribution {
public:
    virtual ~MassiveParticleIsotropicDistribution(){

    };
	double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double distribution(const double& energy) {
		return my_concentration * distributionNormalized(energy);
	}
	virtual double distributionNormalized(const double& energy) = 0;
	void writeDistribution(const char* fileName, int Ne, const double& Emin, const double& Emax);

	double distributionNormalizedTransformedByLosses(const double& energy, const double& lossRate, const double& time);

};

class MassiveParticleMonoenergeticDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_E0;
	double my_dE;
public:
	virtual ~MassiveParticleMonoenergeticDistribution() {

	}
	MassiveParticleMonoenergeticDistribution(const double& mass, const double& Energy, const double& halfWidth, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
};

class MassiveParticleMonoenergeticDirectedDistribution : public MassiveParticleDistribution {
protected:
	double my_E0;
	double my_dE;

	double my_theta0;
	double my_phi0;
	double my_deltaTheta;
public:
	virtual ~MassiveParticleMonoenergeticDirectedDistribution() {

	}
	MassiveParticleMonoenergeticDirectedDistribution(const double& mass, const double& Energy, const double& halfWidth, const double& concentration, const double& theta0, const double& phi0, const double& deltaTheta);
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
};

class MassiveParticlePowerLawDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_index;
	double my_E0;
	double my_A;
public:
    virtual ~MassiveParticlePowerLawDistribution(){

    };
	MassiveParticlePowerLawDistribution(const double& mass, const double& index, const double& E0, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	double getIndex();
	double getE0();
};

class MassiveParticlePowerLawCutoffDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_index;
	double my_beta;
	double my_E0;
	double my_Ecut;
	double my_A;
public:
    virtual ~MassiveParticlePowerLawCutoffDistribution(){

    };
	MassiveParticlePowerLawCutoffDistribution(const double& mass, const double& index, const double& E0, const double& beta, const double& Ecut, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	virtual void resetEcut(const double& Ecut);
	double getIndex();
	double getBeta();
	double getE0();
	double getEcutoff();
};

class MassiveParticleBrokenPowerLawDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_index1;
	double my_index2;
	double my_Etran;
	double my_E0;
	double my_A;
	double my_B;
public:
    virtual ~MassiveParticleBrokenPowerLawDistribution(){

    };
	MassiveParticleBrokenPowerLawDistribution(const double& mass, const double& index1, const double& index2, const double& E0, const double& Etran, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	double getIndex1();
	double getIndex2();
	double getE0();
	double getEtran();
};

class MassiveParticleMaxwellDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
    virtual ~MassiveParticleMaxwellDistribution(){

    };
	MassiveParticleMaxwellDistribution(const double& mass, const double& temperature, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	double getTemperature();
};

class MassiveParticleMaxwellJuttnerDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
	MassiveParticleMaxwellJuttnerDistribution(const double& mass, const double& temperature, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	double getTemperature();
};

class MassiveParticleTabulatedIsotropicDistribution : public MassiveParticleIsotropicDistribution {
protected:
	int my_Ne;
	double* my_energy;
	double* my_distribution;
	DistributionInputType my_inputType;
protected:
	void setDistributionAtPoint(int i, const double& energy, const double& distriution);
	void normalizeDistribution();
public:
	MassiveParticleTabulatedIsotropicDistribution(const MassiveParticleTabulatedIsotropicDistribution& distribution);

	MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* fileName, const int N, const double& concentration, DistributionInputType inputType);
	MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* energyFileName, const char* distributionFileName, const int N, const double& concentration, DistributionInputType inputType);
	MassiveParticleTabulatedIsotropicDistribution(const double& mass, const double* energy, const double* distribution, const int N, const double& concentration, DistributionInputType inputType);
    virtual ~MassiveParticleTabulatedIsotropicDistribution();
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	int getN();
	double getEmin();
	double getEmax();
	void rescaleDistribution(const double& k);
	void addPowerLaw(const double& Epower, const double& index);
	void prolongEnergyRange(const double& Emax, int N);
	void addExponentialCutoff(const double& E);
	void setToZeroAboveE(const double& E);
	void transformToLosses(const double& lossRate, const double& time);
};

class MassiveParticleTabulatedPolarDistribution : public MassiveParticleDistribution {
protected:
	int my_Ne;
	int my_Nmu;
	double* my_energy;
	//note that my_mu is supposed to start at 1 and and at -1, grid may be not uniform
	double* my_mu;
	double** my_distribution;
	DistributionInputType my_inputType;
protected:
	//mu is supposed to be already set
	void setDistributionAtPoint(int i, int j, const double& energy, const double& distribution);
	void normalizeDistribution();
public:
	MassiveParticleTabulatedPolarDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const double& concentration, DistributionInputType inputType);
	MassiveParticleTabulatedPolarDistribution(const double& mass, const double* energy, const double* mu, const double** distribution, const int Ne, const int Nmu, const double& concentration, DistributionInputType inputType);
    virtual ~MassiveParticleTabulatedPolarDistribution();
	virtual double getMeanEnergy();
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual void resetConcentration(const double& concentration);
	int getNe();
	double getEmin();
	double getEmax();
	int getNmu();
	void rescaleDistribution(const double& k);
};

class MassiveParticleTabulatedAnisotropicDistribution : public MassiveParticleDistribution {
protected:
	int my_Ne;
	int my_Nmu;
	int my_Nphi;
	double* my_energy;
	//note that my_mu is supposed to start at 1 and and at -1, grid may be not uniform
	//while phi is supposed uniform
	double* my_mu;
	double* my_phi;
	double*** my_distribution;
	DistributionInputType my_inputType;
protected:
	//mu and phi are supposed to be already set
	void setDistributionAtPoint(int i, int j, int k, const double& energy, const double& distribution);
	void normalizeDistribution();
public:
	MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const int Nphi, const double& concentration, DistributionInputType inputType);
	MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const double* energy, const double* mu, const double*** distribution, const int Ne, const int Nmu, const int Nphi, const double& concentration, DistributionInputType inputType);
    virtual ~MassiveParticleTabulatedAnisotropicDistribution();
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
	int getNe();
	double getEmin();
	double getEmax();
	int getNmu();
	int getNphi();
	void rescaleDistribution(const double& k);
};

class CompoundMassiveParticleDistribution : public MassiveParticleDistribution {
private:
	int my_Ndistr;
	MassiveParticleDistribution** my_distributions;
public:
	//todo what if different masses
	CompoundMassiveParticleDistribution(int N, MassiveParticleDistribution** distributions);
	CompoundMassiveParticleDistribution(MassiveParticleDistribution* dist1, MassiveParticleDistribution* dist2);
	CompoundMassiveParticleDistribution(MassiveParticleDistribution* dist1, MassiveParticleDistribution* dist2, MassiveParticleDistribution* dist3);
    virtual ~CompoundMassiveParticleDistribution();

	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
};

class CompoundWeightedMassiveParticleDistribution : public MassiveParticleDistribution {
private:
	int my_Ndistr;
	double* my_weights;
	MassiveParticleDistribution** my_distributions;
public:
	CompoundWeightedMassiveParticleDistribution(int N, const double* weights, MassiveParticleDistribution** distributions);
	CompoundWeightedMassiveParticleDistribution(MassiveParticleDistribution* dist1, const double& w1, MassiveParticleDistribution* dist2, const double& w2);
	CompoundWeightedMassiveParticleDistribution(MassiveParticleDistribution* dist1, const double& w1, MassiveParticleDistribution* dist2, const double& w2, MassiveParticleDistribution* dist3, const double& w3);
    virtual ~CompoundWeightedMassiveParticleDistribution();

	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
	virtual void resetConcentration(const double& concentration);
};

class MassiveParticleDistributionFactory {
public:
	static MassiveParticleDistribution** readTabulatedIsotropicDistributions(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentrations, int Ne);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributions(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentrations, int Ne);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentration, int Ne, const double& Epower, const double& index);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentration, int Ne, const double& Epower, const double& index);
};

#endif
