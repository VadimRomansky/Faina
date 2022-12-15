#ifndef photon_distribution
#define photon_distribution

/** abstract class, containing photon energy distribution function, normalized to the concentration
* number of particles dN = f(E, mu, phi) dE dmu dphi dV where mu = cos theta
*/

class PhotonDistribution {
protected:
	double my_concentration;
public:
	virtual double distribution(const double& energy, const double& mu, const double& phi) = 0;
	double getConcentration();
};

class PhotonPowerLawDistribution : public PhotonDistribution {
private:
	double my_index;
	double my_E0;
	double my_A;
public:
	PhotonPowerLawDistribution(const double& index, const double& E0, const double& concentration);
	double distribution(const double& energy, const double& mu, const double& phi);

	double getIndex();
	double getE0();
};

class PhotonPlankDistribution : public PhotonDistribution {
private:
	double my_temperature;
public:
	PhotonPlankDistribution(const double& temperature, const double& concentration);
	double distribution(const double& energy, const double& mu, const double& phi);

	double getTemperature();
};

class PhotonMultiPlankDistribution : public PhotonDistribution {
private:
	double my_Nplank;
	double* my_temperatures;
	double* my_concentrations;
	double* my_A;
public:
	PhotonMultiPlankDistribution(int Nplank, double* temperatures, double* concentrations);
	double distribution(const double& energy, const double& mu, const double& phi);
};

#endif
