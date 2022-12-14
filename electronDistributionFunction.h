#ifndef electron_distribution_h
#define electron_distribution_h

class ElectronDistribution {
public:
	virtual double distribution(const double& energy, const double& mu, const double& phi);
};

class ElectronPowerLawDistribution : public ElectronDistribution {
private:
	double my_index;
	double my_E0;
	double my_concentration;
	double my_A;
public:
	ElectronPowerLawDistribution(const double& index, const double& E0, const double& concentration);
	double distribution(const double& energy, const double& mu, const double& phi);
};

#endif
