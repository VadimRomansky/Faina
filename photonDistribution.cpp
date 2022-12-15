#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "photonDistribution.h"

double PhotonDistribution::getConcentration()
{
    return my_concentration;
}

PhotonPowerLawDistribution::PhotonPowerLawDistribution(const double& index, const double& E0, const double& concentration)
{
	if (index <= 1.0) {
		printf("photon spectrum index <= 1.0, contains infinit energy\n");
		printLog("photon spectrum index <= 1.0, contains infinit energy\n");
		exit(0);
	}
	my_index = index;
	my_E0 = E0;
	if (concentration <= 0) {
		printf("photons concentration <= 0\n");
		printLog("photons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_A = my_concentration * (my_index - 1) / (my_E0 * 4 * pi);
}

double PhotonPowerLawDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	if (energy < 0) {
		printf("photon energy < 0\n");
		printLog("photon energy < 0\n");
		exit(0);
	}
    return my_A/pow(energy/my_E0, my_index);
}

double PhotonPowerLawDistribution::getIndex()
{
    return my_index;
}

double PhotonPowerLawDistribution::getE0()
{
    return my_E0;
}

PhotonPlankDistribution::PhotonPlankDistribution(const double& temperature, const double& concentration)
{
}
