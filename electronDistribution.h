#ifndef electron_distribution_h
#define electron_distribution_h
#include "particleDistribution.h"

/** abstract class, containing electron energy distribution function, normalized to the concentration
* number of particles dN = f(E, mu, phi) dE dmu dphi dV where mu = cos theta
*/

enum ElectronInputType {ENERGY_FE, ENERGY_KIN_FE, GAMMA_KIN_FGAMMA, GAMMA_FGAMMA, MOMENTUM_FP};

class ElectronDistribution : public ParticleDistribution{
public:
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi) = 0;
	virtual void resetConcentration(const double& concentration) = 0;
};

class ElectronIsotropicDistribution : public ElectronDistribution {
public:
	double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double distribution(const double& energy) {
		return my_concentration * distributionNormalized(energy);
	}
	virtual double distributionNormalized(const double& energy) = 0;
};

class ElectronPowerLawDistribution : public ElectronIsotropicDistribution {
protected:
	double my_index;
	double my_E0;
	double my_A;
public:
	ElectronPowerLawDistribution(const double& index, const double& E0, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual void resetConcentration(const double& concentration);
	double getIndex();
	double getE0();
};

class ElectronMaxwellDistribution : public ElectronIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
	ElectronMaxwellDistribution(const double& temperature, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual void resetConcentration(const double& concentration);
	double getTemperature();
};

class ElectronMaxwellJuttnerDistribution : public ElectronIsotropicDistribution {
protected:
	double my_temperature;
	double my_A;
public:
	ElectronMaxwellJuttnerDistribution(const double& temperature, const double& concentration);
	virtual double distributionNormalized(const double& energy);
	virtual void resetConcentration(const double& concentration);
	double getTemperature();
};

class ElectronTabulatedIsotropicDistribution : public ElectronIsotropicDistribution {
protected:
	int my_Ne;
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
	virtual double distributionNormalized(const double& energy);
	virtual void resetConcentration(const double& concentration);
	int getN();
	void rescaleDistribution(const double& k);
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
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual void resetConcentration(const double& concentration);
	int getNe();
	int getNmu();
	void rescaleDistribution(const double& k);
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
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual void resetConcentration(const double& concentration);
	int getNe();
	int getNmu();
	int getNphi();
	void rescaleDistribution(const double& k);
};

class CompoundElectronDistribution : public ElectronDistribution {
private:
	int my_Ndistr;
	ElectronDistribution** my_distributions;
public:
	CompoundElectronDistribution(int N, ElectronDistribution** distributions);
	CompoundElectronDistribution(ElectronDistribution* dist1, ElectronDistribution* dist2);
	CompoundElectronDistribution(ElectronDistribution* dist1, ElectronDistribution* dist2, ElectronDistribution* dist3);
	~CompoundElectronDistribution();

	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual void resetConcentration(const double& concentration);
};

class CompoundWeightedElectronDistribution : public ElectronDistribution {
private:
	int my_Ndistr;
	double* my_weights;
	ElectronDistribution** my_distributions;
public:
	CompoundWeightedElectronDistribution(int N, const double* weights, ElectronDistribution** distributions);
	CompoundWeightedElectronDistribution(ElectronDistribution* dist1, const double& w1, ElectronDistribution* dist2, const double& w2);
	CompoundWeightedElectronDistribution(ElectronDistribution* dist1, const double& w1, ElectronDistribution* dist2, const double& w2, ElectronDistribution* dist3, const double& w3);
	~CompoundWeightedElectronDistribution();

	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual void resetConcentration(const double& concentration);
};

class ElectronDistributionFactory {
public:
	static ElectronIsotropicDistribution** readTabulatedIsotropicDistributions(const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, ElectronInputType inputType, const double& electronConcentrations, int Ne);
	static ElectronIsotropicDistribution** readTabulatedIsotropicDistributions(const char* fileName, const char* fileExtension, int Nfiles, ElectronInputType inputType, const double& electronConcentrations, int Ne);
};

#endif