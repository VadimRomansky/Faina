#ifndef partcle_distribution_h
#define particle_distribution_h

class ParticleDistribution {
protected:
	double my_concentration;
public:
	virtual double distribution(const double& energy, const double& mu, const double& phi) = 0;
	double getConcentration() {
		return my_concentration;
	};
};

#endif