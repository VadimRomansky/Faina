#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "photonDistribution.h"

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

PhotonPlankDistribution::PhotonPlankDistribution(const double& temperature, const double& amplitude)
{
	my_temperature = temperature;
	my_A = amplitude;

	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;

	my_concentration = my_A*intPlank2*(8 * pi / cube(hplank * speed_of_light)) * cube(kBoltzman * my_temperature);
}

double PhotonPlankDistribution::distribution(const double& energy, const double& mu, const double& phi) {
	double theta = energy / (kBoltzman * my_temperature);
	return my_A*(2 * energy * energy / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
}

double PhotonPlankDistribution::getTemperature() {
	return my_temperature;
}

PhotonMultiPlankDistribution::PhotonMultiPlankDistribution(int Nplank, double* temperatures, double* amplitudes)
{
	my_Nplank = Nplank;
	my_temperatures = new double[my_Nplank];
	my_A = new double[my_Nplank];
	my_concentrations = new double[my_Nplank];

	my_concentration = 0;
	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;

	for (int i = 0; i < my_Nplank; ++i) {
		my_temperatures[i] = temperatures[i];
		my_A[i] = amplitudes[i];
		my_concentrations[i] = my_A[i] * intPlank2 * (8 * pi / cube(hplank * speed_of_light)) * cube(kBoltzman * my_temperatures[i]);
		my_concentration += my_concentrations[i];
	}
}

PhotonMultiPlankDistribution::~PhotonMultiPlankDistribution()
{
	delete[] my_concentrations;
	delete[] my_A;
	delete[] my_temperatures;
}

double PhotonMultiPlankDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	double result = 0;
	for (int i = 0; i < my_Nplank; ++i) {
		double theta = energy / (kBoltzman * my_temperatures[i]);
		result += my_A[i] * (2 * energy * energy / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
	}

	return result;
}
