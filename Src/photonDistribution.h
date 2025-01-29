#ifndef photon_distribution
#define photon_distribution
#include "particleDistribution.h"

/** abstract class, containing photon energy distribution function, normalized to the concentration
* number of particles dN = f(E, mu, phi) dE dmu dphi dV where mu = cos theta
*/

class PhotonDistribution : public ParticleDistribution{
public:
    virtual ~PhotonDistribution(){

    };
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi) = 0;
};

class PhotonIsotropicDistribution : public PhotonDistribution {
public:
    virtual ~PhotonIsotropicDistribution(){

    };
	double distributionNormalized(const double& energy, const double& mu, const double& phi);

	virtual double distributionNormalized(const double& energy) = 0;
	void writeDistribution(const char* fileName, int Ne, const double& Emin, const double& Emax);
};

class PhotonPowerLawDistribution : public PhotonIsotropicDistribution {
private:
	double my_index;
	double my_E0;
	double my_A;
public:
	PhotonPowerLawDistribution(const double& index, const double& E0);
    virtual ~PhotonPowerLawDistribution(){

    };
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();

	double getIndex();
	double getE0();
};

class PhotonPlankDistribution : public PhotonIsotropicDistribution {
private:
	double my_temperature;
	double my_A;
	double my_concentration;

	static PhotonPlankDistribution* my_CMBradiation;
public:

	PhotonPlankDistribution(const double& temperature, const double& amplitude);
    virtual ~PhotonPlankDistribution(){

    }
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();

	double getConcentration();

	double getTemperature();

	static PhotonPlankDistribution* getCMBradiation();
};

class PhotonMultiPlankDistribution : public PhotonIsotropicDistribution {
private:
	int my_Nplank;
	double* my_temperatures;
	double* my_A;
	double* my_concentrations;
	double my_concentration;

	static PhotonMultiPlankDistribution* my_GalacticField;
public:
	PhotonMultiPlankDistribution(const double& temperature1, const double& amplitude1, const double& temperature2, const double& amplitude2);
	PhotonMultiPlankDistribution(int Nplank, const double* const temperatures, const double* const amplitudes);
    virtual ~PhotonMultiPlankDistribution();
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
	virtual double getConcentration();

	//Mathis 1983?
	static PhotonMultiPlankDistribution* getGalacticField();
};

class CompoundWeightedPhotonDistribution : public PhotonDistribution {
private:
	int my_Ndistr;
	double* my_weights;
	PhotonDistribution** my_distributions;
public:
	CompoundWeightedPhotonDistribution(int N, const double* weights, PhotonDistribution** distributions);
	CompoundWeightedPhotonDistribution(PhotonDistribution* dist1, const double& w1, PhotonDistribution* dist2, const double& w2);
	CompoundWeightedPhotonDistribution(PhotonDistribution* dist1, const double& w1, PhotonDistribution* dist2, const double& w2, PhotonDistribution* dist3,const double& w3);
    virtual ~CompoundWeightedPhotonDistribution();

	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();
};

class PhotonMonoenergeticDistribution : public PhotonIsotropicDistribution {
protected:
	double my_E0;
	double my_dE;
public:
	virtual ~PhotonMonoenergeticDistribution() {

	}
	PhotonMonoenergeticDistribution(const double& Energy, const double& halfWidth);
	virtual double distributionNormalized(const double& energy);
	virtual double getMeanEnergy();
};

class PhotonPlankDirectedDistribution : public PhotonDistribution {
private:
	double my_temperature;
	double my_A;
	double my_concentration;

	double my_theta0;
	double my_phi0;
	double my_deltaTheta;
public:

	PhotonPlankDirectedDistribution(const double& temperature, const double& theta0, const double& phi0, const double& deltaTheta, const double& amplitude);
	virtual ~PhotonPlankDirectedDistribution();
	double distributionNormalized(const double& energy, const double& mu, const double& phi);
	virtual double getMeanEnergy();

	double getConcentration();
	double getTemperature();
};

#endif
