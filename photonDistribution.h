#ifndef photon_distribution
#define photon_distribution
#include "particleDistribution.h"

/** abstract class, containing photon energy distribution function, normalized to the concentration
* number of particles dN = f(E, mu, phi) dE dmu dphi dV where mu = cos theta
*/

class PhotonDistribution : public ParticleDistribution{
public:
	virtual double distribution(const double& energy, const double& mu, const double& phi) = 0;
};

class PhotonIsotropicDistribution : public PhotonDistribution {
public:
	double distribution(const double& energy, const double& mu, const double& phi);
	virtual double distribution(const double& energy) = 0;
};

class PhotonPowerLawDistribution : public PhotonIsotropicDistribution {
private:
	double my_index;
	double my_E0;
	double my_A;
public:
	PhotonPowerLawDistribution(const double& index, const double& E0, const double& concentration);
	virtual double distribution(const double& energy);

	double getIndex();
	double getE0();
};

class PhotonPlankDistribution : public PhotonIsotropicDistribution {
private:
	double my_temperature;
	double my_A;

	static PhotonPlankDistribution* my_CMBradiation;
public:

	PhotonPlankDistribution(const double& temperature, const double& amplitude);
	virtual double distribution(const double& energy);

	double getTemperature();

	static PhotonPlankDistribution* getCMBradiation();
};

class PhotonMultiPlankDistribution : public PhotonIsotropicDistribution {
private:
	double my_Nplank;
	double* my_temperatures;
	double* my_concentrations;
	double* my_A;

	static PhotonMultiPlankDistribution* my_GalacticField;
public:
	PhotonMultiPlankDistribution(int Nplank, const double* const temperatures, const double* const amplitudes);
	~PhotonMultiPlankDistribution();
	virtual double distribution(const double& energy);

	//Mathis 1983?
	static PhotonMultiPlankDistribution* getGalacticField();
};

class CompoundPhotonDistribution : public PhotonDistribution {
private:
	int my_Ndistr;
	PhotonDistribution** my_distributions;
public:
	CompoundPhotonDistribution(int N, PhotonDistribution** distributions);
	CompoundPhotonDistribution(PhotonDistribution* dist1, PhotonDistribution* dist2);
	CompoundPhotonDistribution(PhotonDistribution* dist1, PhotonDistribution* dist2, PhotonDistribution* dist3);
	~CompoundPhotonDistribution();

	virtual double distribution(const double& energy, const double& mu, const double& phi);
};

#endif
