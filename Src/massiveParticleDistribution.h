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
	virtual double minEnergy() = 0;
	virtual double maxEnergy() = 0;
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi) = 0;
	double getMass() {
		return my_mass;
	}
};

class MassiveParticleIsotropicDistribution : public MassiveParticleDistribution {
public:
    virtual ~MassiveParticleIsotropicDistribution(){

    };
	double distributionNormalized(const double& energy, const double& mu, const double& phi);

	virtual double distributionNormalized(const double& energy) = 0;
	void writeDistribution(const char* fileName, int Ne, const double& Emin, const double& Emax);

	double distributionNormalizedWithLosses(const double& energy, const double& lossRate, const double& time);

};

class MassiveParticleMonoenergeticDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_E0;
	double my_dE;
public:
	virtual ~MassiveParticleMonoenergeticDistribution() {

	}
	MassiveParticleMonoenergeticDistribution(const double& mass, const double& Energy, const double& halfWidth);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
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
	MassiveParticleMonoenergeticDirectedDistribution(const double& mass, const double& Energy, const double& halfWidth, const double& theta0, const double& phi0, const double& deltaTheta);
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
};

class MassiveParticlePowerLawDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_index;
	double my_E0;
	double my_A;
public:
    virtual ~MassiveParticlePowerLawDistribution(){

    };
	MassiveParticlePowerLawDistribution(const double& mass, const double& index, const double& E0);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
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
	MassiveParticlePowerLawCutoffDistribution(const double& mass, const double& index, const double& E0, const double& beta, const double& Ecut);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
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
	MassiveParticleBrokenPowerLawDistribution(const double& mass, const double& index1, const double& index2, const double& E0, const double& Etran);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
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
	MassiveParticleMaxwellDistribution(const double& mass, const double& temperature);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
	double getTemperature();
};

class MassiveParticleMaxwellJuttnerDistribution : public MassiveParticleIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
	MassiveParticleMaxwellJuttnerDistribution(const double& mass, const double& temperature);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
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

	MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* fileName, DistributionInputType inputType);
	MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* energyFileName, const char* distributionFileName, DistributionInputType inputType);
	MassiveParticleTabulatedIsotropicDistribution(const double& mass, const double* energy, const double* distribution, const int N, DistributionInputType inputType);
	MassiveParticleTabulatedIsotropicDistribution(MassiveParticleIsotropicDistribution* distribution, const double& minE, const double& maxE, const int N);
    virtual ~MassiveParticleTabulatedIsotropicDistribution();
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
	int getN();
	void rescaleDistribution(const double& k);
	void addPowerLaw(const double& Epower, const double& index);
	void prolongEnergyRange(const double& Emax, int N);
	void addExponentialCutoff(const double& E);
	void setToZeroAboveE(const double& E);
	void transformToLosses(const double& lossRate, const double& time);
	void transformToLosses2(const double& k, const double& l1, const double& l2);
	void transformToThickRegime(const double& Uph);

	double* getEnergyArray();
	double* getDistributionArray();
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
	MassiveParticleTabulatedPolarDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, DistributionInputType inputType);
	MassiveParticleTabulatedPolarDistribution(const double& mass, const double* energy, const double* mu, const double** distribution, const int Ne, const int Nmu, DistributionInputType inputType);
    virtual ~MassiveParticleTabulatedPolarDistribution();
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	int getNe();
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
	MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const int Nphi, DistributionInputType inputType);
	MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const double* energy, const double* mu, const double*** distribution, const int Ne, const int Nmu, const int Nphi, DistributionInputType inputType);
    virtual ~MassiveParticleTabulatedAnisotropicDistribution();
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
	int getNe();
	int getNmu();
	int getNphi();
	void rescaleDistribution(const double& k);
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
	virtual double minEnergy();
	virtual double maxEnergy();
};

class CompoundWeightedMassiveParticleIsotropicDistribution : public MassiveParticleIsotropicDistribution {
private:
	int my_Ndistr;
	double* my_weights;
	MassiveParticleIsotropicDistribution** my_distributions;
public:
	CompoundWeightedMassiveParticleIsotropicDistribution(int N, const double* weights, MassiveParticleIsotropicDistribution** distributions);
	CompoundWeightedMassiveParticleIsotropicDistribution(MassiveParticleIsotropicDistribution* dist1, const double& w1, MassiveParticleIsotropicDistribution* dist2, const double& w2);
	CompoundWeightedMassiveParticleIsotropicDistribution(MassiveParticleIsotropicDistribution* dist1, const double& w1, MassiveParticleIsotropicDistribution* dist2, const double& w2, MassiveParticleIsotropicDistribution* dist3, const double& w3);
	virtual ~CompoundWeightedMassiveParticleIsotropicDistribution();

	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
};

class MassiveParticleDistributionFactory {
public:
	static double evaluateNorm(double* energy, double* distribution, int Ne);
	static void readTabulatedIsotropicDistributionFromMonteCarlo(const double& mass, const char* fileName, MassiveParticleTabulatedIsotropicDistribution*& outputDistribution, double& outputConcentration);
	static void readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(const double& mass, const char* xfileName, const char* pfileName, const char* distributionFileName, double* & xgrid, double* & energy, double** & distributions, double* & concentration, int& Nenergy, int& Nx);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributions(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributions(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne, const double& Epower, const double& index);
	static MassiveParticleDistribution** readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne, const double& Epower, const double& index);
	static void evaluateDistributionAdvectionWithLosses(const double& mass, double* energy, double* distribution, double** outEnergy, double** outDistribution, const int Ne, const int Nx, double* x, double advectionV, double* B, double Uph1, double Eph1, double Uph2, double Eph2);
	static void evaluateDistributionAdvectionWithLossesBetweenTwoPoints(const double& mass, const double* inEnergy, const double* inDistribution, double* outEnergy, double* outDistribution, const int Ne, const double& deltaX, const double& advectionV, const double& B, const double& Uph1, const double& Eph1, const double& Uph2, const double& Eph2, const int Nstep);
	static void redistributeArray(const double* inEnergy, const double* inDistribution, double* outEnergy, double* outDistribution, int Ne);
	static void redistributeArrayWithFixedEnergy(const double* inEnergy, const double* inDistribution, const double* outEnergy, double* outDistribution, int Ne);
	static double getDistribution(const double& E, const double* energy, const double* distribution, int Ne);
};

class MassiveParticleMovingDistribution : public MassiveParticleDistribution{
private:
	MassiveParticleDistribution* my_distribution;
	double my_velocity;
	double my_gamma;
public:
	MassiveParticleMovingDistribution(MassiveParticleDistribution* distribution, const double& velocity);
	virtual ~MassiveParticleMovingDistribution() {

	};
	virtual double getMeanEnergy();
	virtual double minEnergy();
	virtual double maxEnergy();
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
};

#endif
