#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "electronDistribution.h"

double ElectronIsotropicDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	return distribution(energy);
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

double ElectronPowerLawDistribution::distribution(const double& energy) {
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

double ElectronMaxwellDistribution::distribution(const double& energy)
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

double ElectronMaxwellJuttnerDistribution::distribution(const double& energy)
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

void ElectronTabulatedIsotropicDistribution::setDistributionAtPoint(int i, const double& energy, const double& distribution)
{
	if (my_inputType == ElectronInputType::ENERGY_FE) {
		my_energy[i] = energy;
		my_distribution[i] = distribution;
	}
	else if (my_inputType == ElectronInputType::GAMMA_FGAMMA) {
		my_energy[i] = energy * me_c2;
		my_distribution[i] = distribution / me_c2;
	}
	else if (my_inputType == ElectronInputType::MOMENTUM_FP) {
		my_energy[i] = sqrt(energy * energy * speed_of_light2 + me_c2 * me_c2);
		my_distribution[i] = distribution * energy * my_energy[i] / speed_of_light2;
	}
	else {
		printf("unknown electron input type\n");
		printLog("unknown electron input type\n");
		exit(0);
	}
}

void ElectronTabulatedIsotropicDistribution::normalizeDistribution()
{
	double norm = my_distribution[0] * (my_energy[1] - my_energy[0]);
	for (int i = 1; i < my_N; ++i) {
		norm += my_distribution[i] * (my_energy[i] - my_energy[i - 1]);
	}
	norm *= 4 * pi;
	for (int i = 0; i < my_N; ++i) {
		my_distribution[i] *= my_concentration / norm;
	}
}

ElectronTabulatedIsotropicDistribution::ElectronTabulatedIsotropicDistribution(const char* fileName, const int N, const double& concentration, ElectronInputType inputType) {
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
		if (y < 0) {
			printf("input distribution < 0\n");
			printLog("input distribution < 0\n");
			exit(0);
		}
		setDistributionAtPoint(i, x, y);
	}

	normalizeDistribution();

	fclose(file);
}

ElectronTabulatedIsotropicDistribution::ElectronTabulatedIsotropicDistribution(const char* energyFileName, const char* distributionFileName, const int N, const double& concentration, ElectronInputType inputType) {
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

	FILE* energyFile = fopen(energyFileName, "r");
	FILE* distributionFile = fopen(distributionFileName, "r");
	for (int i = 0; i < my_N; ++i) {
		double x;
		double y;
		fscanf(energyFile, "%lf", &x);
		fscanf(distributionFile, "%lf", &y);
		if (y < 0) {
			printf("input distribution < 0\n");
			printLog("input distribution < 0\n");
			exit(0);
		}
		setDistributionAtPoint(i, x, y);
	}

	normalizeDistribution();

	fclose(energyFile);
	fclose(distributionFile);
}

ElectronTabulatedIsotropicDistribution::ElectronTabulatedIsotropicDistribution(const double* energy, const double* distribution, const int N, const double& concentration, ElectronInputType inputType) {
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

	for (int i = 0; i < my_N; ++i) {
		double x = energy[i];
		double y = distribution[i];
		if (y < 0) {
			printf("input distribution < 0\n");
			printLog("input distribution < 0\n");
			exit(0);
		}
		setDistributionAtPoint(i, x, y);
	}

	normalizeDistribution();
}

ElectronTabulatedIsotropicDistribution::~ElectronTabulatedIsotropicDistribution() {
	delete[] my_distribution;
	delete[] my_energy;
}

double ElectronTabulatedIsotropicDistribution::distribution(const double& energy) {
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
		// linear
		/*for (int i = 1; i < my_N; ++i) {
			currentI = i - 1;
			nextI = i;
			if (my_energy[i] > energy) {
				break;
			}
		}*/

		//log
		int minI = 0;
		int maxI = my_N-1;
		while (maxI - minI > 1) {
			int tempI = (minI + maxI) / 2;
			if (my_energy[tempI] > energy) {
				maxI = tempI;
			}
			else {
				minI = tempI;
			}
		}
		currentI = minI;
		nextI = maxI;


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

int ElectronTabulatedIsotropicDistribution::getN() {
	return my_N;
}

void ElectronTabulatedAzimutalDistribution::setDistributionAtPoint(int i, int j, const double& energy, const double& distribution)
{
	if (my_inputType == ElectronInputType::ENERGY_FE) {
		my_energy[i] = energy;
		my_distribution[i][j] = distribution;
	}
	else if (my_inputType == ElectronInputType::GAMMA_FGAMMA) {
		my_energy[i] = energy * me_c2;
		my_distribution[i][j] = distribution / me_c2;
	}
	else if (my_inputType == ElectronInputType::MOMENTUM_FP) {
		my_energy[i] = sqrt(energy * energy * speed_of_light2 + me_c2 * me_c2);
		my_distribution[i][j] = distribution * energy * my_energy[i] / speed_of_light2;
	}
	else {
		printf("unknown electron input type\n");
		printLog("unknown electron input type\n");
		exit(0);
	}
}

void ElectronTabulatedAzimutalDistribution::normalizeDistribution()
{
	double norm = 0;
	for (int imu = 0; imu < my_Nmu; ++imu) {
		double dmu;
		if (imu == 0) {
			dmu = 1 - (my_mu[0] + my_mu[1]) / 2;
		}
		else if (imu == my_Nmu - 1) {
			dmu = (my_mu[my_Nmu - 2] + my_mu[my_Nmu - 1]) / 2 + 1;
		}
		else {
			dmu = (my_mu[imu - 1] - my_mu[imu + 1]) / 2;
		}
		norm += my_distribution[0][imu] * (my_energy[1] - my_energy[0]) * dmu;
		for (int i = 1; i < my_Ne; ++i) {
			norm += my_distribution[i][imu] * (my_energy[i] - my_energy[i - 1]) * dmu;
		}
	}

	for (int imu = 0; imu < my_Nmu; ++imu) {
		for (int i = 0; i < my_Ne; ++i) {
			my_distribution[i][imu] *= my_concentration / norm;
		}
	}
}

ElectronTabulatedAzimutalDistribution::ElectronTabulatedAzimutalDistribution(const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const double& concentration, ElectronInputType inputType)
{
	if (Ne <= 0) {
		printf("energy grid number <= 0 in tabulated azimutal distribution\n");
		printLog("energy grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Ne = Ne;

	if (Nmu <= 0) {
		printf("mu grid number <= 0 in tabulated azimutal distribution\n");
		printLog("mu grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Nmu = Nmu;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_mu = new double[my_Nmu];
	my_distribution = new double*[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = new double[my_Nmu];
	}

	FILE* energyFile = fopen(energyFileName, "r");
	FILE* muFile = fopen(muFileName, "r");
	FILE* distributionFile = fopen(distributionFileName, "r");
	for (int i = 0; i < my_Nmu; ++i) {
		fscanf(muFile, "%lf", &my_mu[i]);
	}

	for (int i = 0; i < my_Ne; ++i) {
		double x;
		double y;
		fscanf(energyFile, "%lf", &x);
		for (int j = 0; j < my_Nmu; ++j) {
			fscanf(distributionFile, "%lf", &y);
			if (y < 0) {
				printf("input distribution < 0\n");
				printLog("input distribution < 0\n");
				exit(0);
			}
			setDistributionAtPoint(i, j, x, y);
		}
	}

	normalizeDistribution();

	fclose(energyFile);
	fclose(muFile);
	fclose(distributionFile);
}

ElectronTabulatedAzimutalDistribution::ElectronTabulatedAzimutalDistribution(const double* energy, const double* mu, const double** distribution, const int Ne, const int Nmu, const double& concentration, ElectronInputType inputType)
{
	if (Ne <= 0) {
		printf("energy grid number <= 0 in tabulated azimutal distribution\n");
		printLog("energy grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Ne = Ne;

	if (Nmu <= 0) {
		printf("mu grid number <= 0 in tabulated azimutal distribution\n");
		printLog("mu grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Nmu = Nmu;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_mu = new double[my_Nmu];
	my_distribution = new double* [my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = new double[my_Nmu];
	}

	for (int i = 0; i < my_Nmu; ++i) {
		my_mu[i] = mu[i];
	}

	for (int i = 0; i < my_Ne; ++i) {
		double x;
		double y;
		x = energy[i];
		for (int j = 0; j < my_Nmu; ++j) {
			y = distribution[i][j];
			if (y < 0) {
				printf("input distribution < 0\n");
				printLog("input distribution < 0\n");
				exit(0);
			}
			setDistributionAtPoint(i, j, x, y);
		}
	}

	normalizeDistribution();
}

ElectronTabulatedAzimutalDistribution::~ElectronTabulatedAzimutalDistribution()
{
	delete[] my_energy;
	delete[] my_mu;
	for (int i = 0; i < my_Ne; ++i) {
		delete[] my_distribution[i];
	}
	delete[] my_distribution;
}

double ElectronTabulatedAzimutalDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	if (mu > 1.0) {
		printf("mu = %lf > 1.0 in azimutal distribution\n", mu);
		printLog("mu> 1.0 in azimutal distribution\n");
		exit(0);
	}
	if (mu < -1.0) {
		printf("mu = %lf < -1.0 in azimutal distribution\n", mu);
		printLog("mu < -1.0 in azimutal distribution\n");
		exit(0);
	}
	if (energy <= my_energy[0]) {
		printf("warning: energy is less than minimum energy\n");
		printLog("warning: energy is less than minimum energy\n");
		//exit(0);
		return 0;
	}
	else if (energy >= my_energy[my_Ne - 1]) {
		printf("warning: energy is greater than maximum energy\n");
		printLog("warning: energy isgreater than maximum energy\n");
		return 0;
	}
	else {
		int currentI = 0;
		int nextI = 1;
		// linear
		/*for (int i = 1; i < my_N; ++i) {
			currentI = i - 1;
			nextI = i;
			if (my_energy[i] > energy) {
				break;
			}
		}*/

		//log
		int minI = 0;
		int maxI = my_Ne - 1;
		while (maxI - minI > 1) {
			int tempI = (minI + maxI) / 2;
			if (my_energy[tempI] > energy) {
				maxI = tempI;
			}
			else {
				minI = tempI;
			}
		}
		currentI = minI;
		nextI = maxI;

		int currentJ = 0;
		int nextJ = 1;
		if (mu > my_mu[0]) {
			currentJ = 0;
			nextJ = 1;
		}
		else if (mu < my_mu[my_Nmu]) {
			currentJ = my_Nmu - 2;
			nextJ = my_Nmu - 1;
		}
		else {
			int minJ = 0;
			int maxJ = my_Nmu - 1;
			while (maxJ - minJ > 1) {
				int tempJ = (minJ + maxJ) / 2;
				if (my_mu[tempJ] < mu) {
					maxJ = tempJ;
				}
				else {
					minJ = tempJ;
				}
			}
			currentJ = minJ;
			nextJ = maxJ;
		}


		double result;
		double totalArea = (my_energy[nextI] - my_energy[currentI]) * (my_mu[currentJ] - my_mu[nextJ]);
		double weight00 = (my_energy[nextI] - energy) * (mu - my_mu[nextJ])/totalArea;
		double weight01 = (my_energy[nextI] - energy) * (my_mu[currentJ] - mu)/totalArea;
		double weight10 = (energy - my_energy[currentI]) * (mu - my_mu[nextJ])/totalArea;
		double weight11 = (energy - my_energy[currentI]) * (my_mu[currentJ] - mu)/totalArea;
		if (my_distribution[currentI][currentJ] <= 0 || my_distribution[nextI][currentJ] <= 0 || my_distribution[currentI][nextJ] <= 0 || my_distribution[nextI][nextJ] <= 0) {
			result = weight00 * my_distribution[currentI][currentJ] +
				weight01 * my_distribution[currentI][nextJ] +
				weight10 * my_distribution[nextI][currentJ] +
				weight11 * my_distribution[nextI][nextJ];
			if (result < 0) {
				printf("warning: linear aproximated distribution < 0\n");
				printLog("warning: linear aproximated distribution < 0\n");
				result = 0;
			}
		}
		else {
			double logResult = weight00 * 0 +
				weight01 * log(my_distribution[currentI][nextJ] / my_distribution[currentI][currentJ]) +
				weight10 * log(my_distribution[nextI][currentJ] / my_distribution[currentI][currentJ]) +
				weight11 * log(my_distribution[nextI][nextJ] / my_distribution[currentI][currentJ]);
			result = my_distribution[currentI][currentJ] * exp(logResult);
		}
		if (result != result) {
			printf("result = NaN in electron tabulated distribution\n");
			printLog("result = NaN in electron tabulated distribution\n");
			exit(0);
		}
		return result;
	}
}

int ElectronTabulatedAzimutalDistribution::getNe()
{
	return my_Ne;
}

int ElectronTabulatedAzimutalDistribution::getNmu()
{
	return my_Nmu;
}

void ElectronTabulatedAnisotropicDistribution::setDistributionAtPoint(int i, int j, int k, const double& energy, const double& distribution)
{
	if (my_inputType == ElectronInputType::ENERGY_FE) {
		my_energy[i] = energy;
		my_distribution[i][j][k] = distribution;
	}
	else if (my_inputType == ElectronInputType::GAMMA_FGAMMA) {
		my_energy[i] = energy * me_c2;
		my_distribution[i][j][k] = distribution / me_c2;
	}
	else if (my_inputType == ElectronInputType::MOMENTUM_FP) {
		my_energy[i] = sqrt(energy * energy * speed_of_light2 + me_c2 * me_c2);
		my_distribution[i][j][k] = distribution * energy * my_energy[i] / speed_of_light2;
	}
	else {
		printf("unknown electron input type\n");
		printLog("unknown electron input type\n");
		exit(0);
	}
}

void ElectronTabulatedAnisotropicDistribution::normalizeDistribution()
{
	double norm = 0;
	double dphi = 2 * pi / my_Nphi;
	for (int imu = 0; imu < my_Nmu; ++imu) {
		double dmu;
		if (imu == 0) {
			dmu = 1 - (my_mu[0] + my_mu[1]) / 2;
		}
		else if (imu == my_Nmu - 1) {
			dmu = (my_mu[my_Nmu - 2] + my_mu[my_Nmu - 1]) / 2 + 1;
		}
		else {
			dmu = (my_mu[imu - 1] - my_mu[imu + 1]) / 2;
		}
		for (int k = 0; k < my_Nphi; ++k) {
			norm += my_distribution[0][imu][k] * (my_energy[1] - my_energy[0]) * dmu * dphi;
			for (int i = 1; i < my_Ne; ++i) {
				norm += my_distribution[i][imu][k] * (my_energy[i] - my_energy[i - 1]) * dmu * dphi;
			}
		}
	}

	for (int imu = 0; imu < my_Nmu; ++imu) {
		for (int i = 0; i < my_Ne; ++i) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_distribution[i][imu][iphi] *= my_concentration / norm;
			}
		}
	}
}

ElectronTabulatedAnisotropicDistribution::ElectronTabulatedAnisotropicDistribution(const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const int Nphi, const double& concentration, ElectronInputType inputType)
{
	if (Ne <= 0) {
		printf("energy grid number <= 0 in tabulated azimutal distribution\n");
		printLog("energy grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Ne = Ne;

	if (Nmu <= 0) {
		printf("mu grid number <= 0 in tabulated azimutal distribution\n");
		printLog("mu grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Nmu = Nmu;

	if (Nphi <= 0) {
		printf("phi grid number <= 0 in tabulated azimutal distribution\n");
		printLog("phi grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Nphi = Nphi;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_mu = new double[my_Nmu];
	my_phi = new double[my_Nphi];
	my_distribution = new double** [my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = new double*[my_Nmu];
		for (int j = 0; j < my_Nmu; ++j) {
			my_distribution[i][j] = new double[my_Nphi];
		}
	}

	FILE* energyFile = fopen(energyFileName, "r");
	FILE* muFile = fopen(muFileName, "r");
	FILE* distributionFile = fopen(distributionFileName, "r");
	for (int i = 0; i < my_Nmu; ++i) {
		fscanf(muFile, "%lf", &my_mu[i]);
	}
	for (int i = 0; i < my_Nphi; ++i) {
		my_phi[i] = 2 * pi * (i + 0.5) / my_Nphi;
	}

	for (int i = 0; i < my_Ne; ++i) {
		double x;
		double y;
		fscanf(energyFile, "%lf", &x);
		for (int j = 0; j < my_Nmu; ++j) {
			for (int k = 0; k < my_Nphi; ++k) {
				fscanf(distributionFile, "%lf", &y);
				if (y < 0) {
					printf("input distribution < 0\n");
					printLog("input distribution < 0\n");
					exit(0);
				}
				setDistributionAtPoint(i, j, k, x, y);
			}
		}
	}

	normalizeDistribution();

	fclose(energyFile);
	fclose(muFile);
	fclose(distributionFile);
}

ElectronTabulatedAnisotropicDistribution::ElectronTabulatedAnisotropicDistribution(const double* energy, const double* mu, const double*** distribution, const int Ne, const int Nmu, const int Nphi, const double& concentration, ElectronInputType inputType)
{
	if (Ne <= 0) {
		printf("energy grid number <= 0 in tabulated azimutal distribution\n");
		printLog("energy grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Ne = Ne;

	if (Nmu <= 0) {
		printf("mu grid number <= 0 in tabulated azimutal distribution\n");
		printLog("mu grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Nmu = Nmu;

	if (Nphi <= 0) {
		printf("phi grid number <= 0 in tabulated azimutal distribution\n");
		printLog("phi grid number <= 0 in tabulated azimutal distribution\n");
		exit(0);
	}
	my_Nphi = Nphi;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_mu = new double[my_Nmu];
	my_phi = new double[my_Nphi];
	my_distribution = new double** [my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = new double* [my_Nmu];
		for (int j = 0; j < my_Nmu; ++j) {
			my_distribution[i][j] = new double[my_Nphi];
		}
	}

	for (int i = 0; i < my_Nmu; ++i) {
		my_mu[i] = mu[i];
	}
	for (int i = 0; i < my_Nphi; ++i) {
		my_phi[i] = 2 * pi * (i + 0.5) / my_Nphi;
	}

	for (int i = 0; i < my_Ne; ++i) {
		double x;
		double y;
		x = energy[i];
		for (int j = 0; j < my_Nmu; ++j) {
			for (int k = 0; k < my_Nphi; ++k) {
				y = distribution[i][j][k];
				if (y < 0) {
					printf("input distribution < 0\n");
					printLog("input distribution < 0\n");
					exit(0);
				}
				setDistributionAtPoint(i, j, k, x, y);
			}
		}
	}

	normalizeDistribution();
}

ElectronTabulatedAnisotropicDistribution::~ElectronTabulatedAnisotropicDistribution()
{
	for (int i = 0; i < my_Ne; ++i) {
		for (int j = 0; j < my_Nmu; ++j) {
			delete[] my_distribution[i][j];
		}
		delete[] my_distribution[i];
	}
	delete[] my_distribution;
	delete[] my_mu;
	delete[] my_phi;
}

double ElectronTabulatedAnisotropicDistribution::distribution(const double& energy, const double& mu, const double& phi)
{
	if (mu > 1.0) {
		printf("mu = %lf > 1.0 in azimutal distribution\n", mu);
		printLog("mu> 1.0 in azimutal distribution\n");
		exit(0);
	}
	if (mu < -1.0) {
		printf("mu = %lf < -1.0 in azimutal distribution\n", mu);
		printLog("mu < -1.0 in azimutal distribution\n");
		exit(0);
	}
	if (energy <= my_energy[0]) {
		printf("warning: energy is less than minimum energy\n");
		printLog("warning: energy is less than minimum energy\n");
		//exit(0);
		return 0;
	}
	else if (energy >= my_energy[my_Ne - 1]) {
		printf("warning: energy is greater than maximum energy\n");
		printLog("warning: energy isgreater than maximum energy\n");
		return 0;
	}
	double tempPhi = phi;
	if (phi < 0) {
		printf("warning: phi = %lf < 0\n", phi);
		printLog("warning: phi < 0\n");
		while (tempPhi < 0) {
			tempPhi += 2 * pi;
		}
	}
	if (phi >= 2 * pi) {
		printf("warning: phi = %lf > 2 pi\n", phi);
		printLog("warning: phi > 2 pi\n");
		while (tempPhi >= 2*pi) {
			tempPhi -= 2 * pi;
		}
	}
	else {
		int currentI = 0;
		int nextI = 1;
		// linear
		/*for (int i = 1; i < my_N; ++i) {
			currentI = i - 1;
			nextI = i;
			if (my_energy[i] > energy) {
				break;
			}
		}*/

		//log
		int minI = 0;
		int maxI = my_Ne - 1;
		while (maxI - minI > 1) {
			int tempI = (minI + maxI) / 2;
			if (my_energy[tempI] > energy) {
				maxI = tempI;
			}
			else {
				minI = tempI;
			}
		}
		currentI = minI;
		nextI = maxI;

		int currentJ = 0;
		int nextJ = 1;
		if (mu > my_mu[0]) {
			currentJ = 0;
			nextJ = 1;
		}
		else if (mu < my_mu[my_Nmu-1]) {
			currentJ = my_Nmu - 2;
			nextJ = my_Nmu - 1;
		}
		else {
			int minJ = 0;
			int maxJ = my_Nmu - 1;
			while (maxJ - minJ > 1) {
				int tempJ = (minJ + maxJ) / 2;
				if (my_mu[tempJ] < mu) {
					maxJ = tempJ;
				}
				else {
					minJ = tempJ;
				}
			}
			currentJ = minJ;
			nextJ = maxJ;
		}

		int currentK = 0;
		int nextK = 1;
		if (tempPhi < my_phi[0]) {
			currentK = 0;
			nextK = 1;
		}
		else if (tempPhi > my_phi[my_Nphi-1]) {
			currentK = my_Nphi - 2;
			nextK = my_Nphi - 1;
		}
		else {
			int minK = 0;
			int maxK = my_Nphi - 1;
			while (maxK - minK > 1) {
				int tempK = (minK + maxK) / 2;
				if (my_phi[tempK] > tempPhi) {
					maxK = tempK;
				}
				else {
					minK = tempK;
				}
			}
			currentK = minK;
			nextK = maxK;
		}

		double dphi = 2 * pi / my_Nphi;
		double result;
		double totalVolume = (my_energy[nextI] - my_energy[currentI]) * (my_mu[currentJ] - my_mu[nextJ])*dphi;
		double weight000 = (my_energy[nextI] - energy) * (mu - my_mu[nextJ]) * (my_phi[nextK] - tempPhi) / totalVolume;
		double weight010 = (my_energy[nextI] - energy) * (my_mu[currentJ] - mu) * (my_phi[nextK] - tempPhi) / totalVolume;
		double weight100 = (energy - my_energy[currentI]) * (mu - my_mu[nextJ]) * (my_phi[nextK] - tempPhi) / totalVolume;
		double weight110 = (energy - my_energy[currentI]) * (my_mu[currentJ] - mu) * (my_phi[nextK] - tempPhi) / totalVolume;
		double weight001 = (my_energy[nextI] - energy) * (mu - my_mu[nextJ]) * (tempPhi - my_phi[currentK]) / totalVolume;
		double weight011 = (my_energy[nextI] - energy) * (my_mu[currentJ] - mu) * (tempPhi - my_phi[currentK]) / totalVolume;
		double weight101 = (energy - my_energy[currentI]) * (mu - my_mu[nextJ]) * (tempPhi - my_phi[currentK]) / totalVolume;
		double weight111 = (energy - my_energy[currentI]) * (my_mu[currentJ] - mu) * (tempPhi - my_phi[currentK]) / totalVolume;
		if (my_distribution[currentI][currentJ][currentK] <= 0 || my_distribution[nextI][currentJ][currentK] <= 0 || my_distribution[currentI][nextJ][currentK] <= 0 || my_distribution[nextI][nextJ][currentK] <= 0 ||
			my_distribution[currentI][currentJ][nextK] <= 0 || my_distribution[nextI][currentJ][nextK] <= 0 || my_distribution[currentI][nextJ][nextK] <= 0 || my_distribution[nextI][nextJ][nextK] <= 0) {
			result = weight000 * my_distribution[currentI][currentJ][currentK] +
				weight010 * my_distribution[currentI][nextJ][currentK] +
				weight100 * my_distribution[nextI][currentJ][currentK] +
				weight110 * my_distribution[nextI][nextJ][currentK] +
				weight001 * my_distribution[currentI][currentJ][nextK] +
				weight011 * my_distribution[currentI][nextJ][nextK] +
				weight101 * my_distribution[nextI][currentJ][nextK] +
				weight111 * my_distribution[nextI][nextJ][nextK];
			if (result < 0) {
				printf("warning: linear aproximated distribution < 0\n");
				printLog("warning: linear aproximated distribution < 0\n");
				result = 0;
			}
		}
		else {
			double logResult = weight000 * 0 +
				weight010 * log(my_distribution[currentI][nextJ][currentK]) +
				weight100 * log(my_distribution[nextI][currentJ][currentK]) +
				weight110 * log(my_distribution[nextI][nextJ][currentK]) +
				weight001 * log(my_distribution[currentI][currentJ][nextK]) +
				weight011 * log(my_distribution[currentI][nextJ][nextK]) +
				weight101 * log(my_distribution[nextI][currentJ][nextK]) +
				weight111 * log(my_distribution[nextI][nextJ][nextK]);
			result = my_distribution[currentI][currentJ][currentK] * exp(logResult);
		}
		if (result != result) {
			printf("result = NaN in electron tabulated distribution\n");
			printLog("result = NaN in electron tabulated distribution\n");
			exit(0);
		}
		return result;
	}
}

int ElectronTabulatedAnisotropicDistribution::getNe()
{
	return my_Ne;
}

int ElectronTabulatedAnisotropicDistribution::getNmu()
{
	return my_Nmu;
}

int ElectronTabulatedAnisotropicDistribution::getNphi()
{
	return my_Nphi;
}
