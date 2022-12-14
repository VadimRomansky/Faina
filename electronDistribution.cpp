#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "electronDistributionFunction.h"

ElectronPowerLawDistribution::ElectronPowerLawDistribution(const double& index, const double& E0, const double& concentration) {
	if (index <= 1.0) {
		printf("electron spectrum index <= 1.0, contains infinit energy\n");
		printLog("electron spectrum index <= 1.0, contains infinit energy\n");
		exit(0);
	}
	my_index = index;
	if (E0 < me_c2) {
		printf("electrons minimum energy is less than m c^2\n");
		printLog("electrons minimum energy is less than m c^2\n");
		exit(0);
	}
	my_E0 = E0;
	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_A = my_concentration * (my_index - 1) / my_E0;
}

double ElectronPowerLawDistribution::distribution(const double& energy, const double& mu, const double& phi) {
	if (energy < me_c2) {
		printf("warning: energy is less than m c^2\n");
		printLog("warning: energy is less than m c^2\n");
		//exit(0);
		return 0;
	}
	return my_A / pow(energy / my_E0, my_index);
}
