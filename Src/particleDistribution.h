#ifndef particle_distribution_h
#define particle_distribution_h

class ParticleDistribution {
protected:
	//double my_concentration;
public:
    virtual ~ParticleDistribution(){

    }
	/*virtual double distribution(const double& energy, const double& mu, const double& phi) {
		return my_concentration * distributionNormalized(energy, mu, phi);
	};*/
	virtual double distributionNormalized(const double& energy, const double& mu, const double& phi) = 0;
	virtual double getMeanEnergy() = 0;
	/*double getConcentration() {
		return my_concentration;
	};*/
	/*virtual void resetConcentration(const double& concentration) {
		my_concentration = concentration;
	}*/
};

#endif
