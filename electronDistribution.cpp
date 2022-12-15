#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "electronDistribution.h"

double ElectronDistribution::getConcentration()
{
	return my_concentration;
}

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

	my_A = my_concentration * (my_index - 1) / (my_E0*4*pi);
}

double ElectronPowerLawDistribution::distribution(const double& energy, const double& mu, const double& phi) {
	if (energy < 0) {
		printf("electron energy < 0\n");
		printLog("electron energy < 0\n");
		exit(0);
	}
	if (energy < me_c2) {
		printf("warning: energy is less than m c^2\n");
		printLog("warning: energy is less than m c^2\n");
		//exit(0);
		return 0;
	}
	return my_A / pow(energy / my_E0, my_index);
}

double ElectronPowerLawDistribution::getIndex()
{
	return my_index;
}

double ElectronPowerLawDistribution::getE0()
{
	return my_E0;
}

ElectronMaxwellDistribution::ElectronMaxwellDistribution(const double& temperature, const double& concentration)
{
	if (temperature <= 0) {
		printf("electrons temperature <= 0\n");
		printLog("electrons temperature <= 0\n");
		exit(0);
	}
	my_temperature = temperature;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_A = my_concentration / (2 * sqrt(cube(pi * kBoltzman * my_temperature)));
}

double ElectronMaxwellDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	return my_A*sqrt(energy)*exp(-energy/(kBoltzman*my_temperature));
}

double ElectronMaxwellDistribution::getTemperature()
{
	return my_temperature;
}

ElectronMaxwellJuttnerDistribution::ElectronMaxwellJuttnerDistribution(const double& temperature, const double& concentration)
{
	if (temperature <= 0) {
		printf("electrons temperature <= 0\n");
		printLog("electrons temperature <= 0\n");
		exit(0);
	}
	my_temperature = temperature;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	double theta = kBoltzman * my_temperature / me_c2;
	my_A = my_concentration / (4*pi*cube(me_c2)*theta*McDonaldFunction(2, 1/theta));
}

double ElectronMaxwellJuttnerDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	if (energy < me_c2) {
		printf("electron energy is less than m c^2 in maxwell-juttner distribution\n");
		printLog("electron energy is less than m c^2 in maxwell-juttner distribution\n");
		exit(0);
	}
	return my_A * sqrt(energy) * sqrt(energy*energy - me_c2*me_c2)*energy*exp(-energy / (kBoltzman * my_temperature));
}

double ElectronMaxwellJuttnerDistribution::getTemperature()
{
	return my_temperature;
}

ElectronTabulatedSphericalDistribution::ElectronTabulatedSphericalDistribution(const char* fileName, const int N, const double& concentration, ElectronInputType inputType) {
	if (N <= 0) {
		printf("grid number <= 0 in tabulated spherical distribution\n");
		printLog("grid number <= 0 in tabulated spherical distribution\n");
		exit(0);
	}
	my_N = N;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_N];
	my_distribution = new double[my_N];

	FILE* file = fopen(fileName, "r");
	for (int i = 0; i < my_N; ++i) {
		double x;
		double y;
		fscanf(file, "%lf", &x);
		fscanf(file, "%lf", &y);
		if (my_inputType == ElectronInputType::ENERGY_FE) {
			my_energy[i] = x;
			my_distribution[i] = y;
		} else if (my_inputType == ElectronInputType::GAMMA_FGAMMA) {
			my_energy[i] = x * me_c2;
			my_distribution[i] = y / me_c2;
		} else if (my_inputType == ElectronInputType::MOMENTUM_FP) {
			my_energy[i] = sqrt(x * x * speed_of_light2 + me_c2 * me_c2);
			my_distribution[i] = y * x * my_energy[i] / speed_of_light2;
		} else {
			printf("unknown electron input type\n");
			printLog("unknown electron input type\n");
			exit(0);
		}
	}

	double norm = my_distribution[0] * (my_energy[1] - my_energy[0]);
	for (int i = 1; i < my_N; ++i) {
		norm += my_distribution[i] * (my_energy[i] - my_energy[i - 1]);
	}
	for (int i = 0; i < my_N; ++i) {
		my_distribution[i] *= my_concentration / norm;
	}

	fclose(file);
}

ElectronTabulatedSphericalDistribution::~ElectronTabulatedSphericalDistribution() {
	delete[] my_distribution;
	delete[] my_energy;
}

double ElectronTabulatedSphericalDistribution::distribution(const double& energy, const double& mu, const double& phi) {
	if (energy <= my_energy[0]) {
		printf("warning: energy is less than minimum energy\n");
		printLog("warning: energy is less than minimum energy\n");
		//exit(0);
		return 0;
	}
	else if (energy >= my_energy[my_N-1]) {
		printf("warning: energy is greater than maximum energy\n");
		printLog("warning: energy isgreater than maximum energy\n");
		return 0;
	}
	else {
		int currentI = 0;
		int nextI = 1;
		for (int i = 1; i < my_N; ++i) {
			currentI = i - 1;
			nextI = i;
			if (my_energy[i] > energy) {
				break;
			}
		}
		double result;
		if (my_distribution[currentI] <= 0 || my_distribution[nextI] <= 0) {
			result = (my_distribution[currentI] * (my_energy[nextI] - energy) + my_distribution[nextI] * (energy - my_energy[currentI])) / (my_energy[nextI] - my_energy[currentI]);
		}
		else {
			result = my_distribution[currentI] * exp(log(my_distribution[nextI] / my_distribution[currentI]) * ((energy - my_energy[currentI]) / (my_energy[nextI] - my_energy[currentI])));
		}
		if (result != result) {
			printf("result = NaN in electron tabulated distribution\n");
			printLog("result = NaN in electron tabulated distribution\n");
			exit(0);
		}
		return result;
	}
}

int ElectronTabulatedSphericalDistribution::getN() {
	return my_N;
}
