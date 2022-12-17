#ifndef electron_distribution_h
#define electron_distribution_h
#include "particleDistribution.h"

enum ElectronInputType {ENERGY_FE, GAMMA_FGAMMA, MOMENTUM_FP};

class ElectronDistribution : public ParticleDistribution{
	virtual double distribution(const double& energy, const double& mu, const double& phi) = 0;
};

class ElectronIsotropicDistribution : public ElectronDistribution {
	double distribution(const double& energy, const double& mu, const double& phi);
	virtual double distribution(const double& energy) = 0;
};

class ElectronPowerLawDistribution : public ElectronIsotropicDistribution {
protected:
	double my_index;
	double my_E0;
	double my_A;
public:
	ElectronPowerLawDistribution(const double& index, const double& E0, const double& concentration);
	double distribution(const double& energy);
	double getIndex();
	double getE0();
};

class ElectronMaxwellDistribution : public ElectronIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
	ElectronMaxwellDistribution(const double& temperature, const double& concentration);
	double distribution(const double& energy);
	double getTemperature();
};

class ElectronMaxwellJuttnerDistribution : public ElectronIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
	ElectronMaxwellJuttnerDistribution(const double& temperature, const double& concentration);
	double distribution(const double& energy);
	double getTemperature();
};

class ElectronTabulatedIsotropicDistribution : public ElectronIsotropicDistribution {
protected:
	int my_N;
	double* my_energy;
	double* my_distribution;
	ElectronInputType my_inputType;
protected:
	void setDistributionAtPoint(int i, const double& energy, const double& distriution);
	void normalizeDistribution();
public:
	ElectronTabulatedIsotropicDistribution(const char* fileName, const int N, const double& concentration, ElectronInputType inputType);
	ElectronTabulatedIsotropicDistribution(const char* energyFileName, const char* distributionFileName, const int N, const double& concentration, ElectronInputType inputType);
	ElectronTabulatedIsotropicDistribution(const double* energy, const double* distribution, const int N, const double& concentration, ElectronInputType inputType);
	~ElectronTabulatedIsotropicDistribution();
	double distribution(const double& energy);
	int getN();
};

class ElectronTabulatedAzimutalDistribution : public ElectronDistribution {
protected:
	int my_Ne;
	int my_Nmu;
	double* my_energy;
	//note that my_mu is supposed to start at 1 and and at -1, grid may be not uniform
	double* my_mu;
	double** my_distribution;
	ElectronInputType my_inputType;
protected:
	//mu is supposed to be already set
	void setDistributionAtPoint(int i, int j, const double& energy, const double& distribution);
	void normalizeDistribution();
public:
	ElectronTabulatedAzimutalDistribution(const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const double& concentration, ElectronInputType inputType);
	ElectronTabulatedAzimutalDistribution(const double* energy, const double* mu, const double** distribution, const int Ne, const int Nmu, const double& concentration, ElectronInputType inputType);
	~ElectronTabulatedAzimutalDistribution();
	double distribution(const double& energy, const double& mu, const double& phi);
	int getNe();
	int getNmu();
};

class ElectronTabulatedAnisotropicDistribution : public ElectronDistribution {
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
	ElectronInputType my_inputType;
protected:
	//mu and phi are supposed to be already set
	void setDistributionAtPoint(int i, int j, int k, const double& energy, const double& distribution);
	void normalizeDistribution();
public:
	ElectronTabulatedAnisotropicDistribution(const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const int Nphi, const double& concentration, ElectronInputType inputType);
	ElectronTabulatedAnisotropicDistribution(const double* energy, const double* mu, const double*** distribution, const int Ne, const int Nmu, const int Nphi, const double& concentration, ElectronInputType inputType);
	~ElectronTabulatedAnisotropicDistribution();
	double distribution(const double& energy, const double& mu, const double& phi);
	int getNe();
	int getNmu();
	int getNphi();
};

#endif