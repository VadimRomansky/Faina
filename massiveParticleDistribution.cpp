#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "massiveParticleDistribution.h"

double MassiveParticleIsotropicDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	return distributionNormalized(energy);
}

void MassiveParticleIsotropicDistribution::writeDistribution(const char* fileName, int Ne, const double& Emin, const double& Emax) {
	double* energy = new double[Ne];
	double factor = pow(Emax / Emin, 1.0 / (Ne - 1));
	energy[0] = Emin;
	for (int i = 1; i < Ne; ++i) {
		energy[i] = energy[i - 1] * factor;
	}
	FILE* outFile = fopen(fileName, "w");
	for (int i = 0; i < Ne; ++i) {
		fprintf(outFile, "%g %g\n", energy[i], distribution(energy[i]));
	}
	fclose(outFile);
	delete[] energy;
}

MassiveParticlePowerLawDistribution::MassiveParticlePowerLawDistribution(const double& mass, const double& index, const double& E0, const double& concentration) {
	my_mass = mass;
	if (index <= 1.0) {
		printf("electron spectrum index <= 1.0, contains infinit energy\n");
		printLog("electron spectrum index <= 1.0, contains infinit energy\n");
		exit(0);
	}
	my_index = index;
	if (E0 < my_mass*speed_of_light2) {
		printf("particle minimum energy is less than m c^2\n");
		printLog("particle minimum energy is less than m c^2\n");
		exit(0);
	}
	my_E0 = E0;
	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_A = (my_index - 1) / (my_E0*4*pi);
}

double MassiveParticlePowerLawDistribution::distributionNormalized(const double& energy) {
	//for debug
	/*if (energy > 1000 * getMass() * speed_of_light2) {
		return 0;
	}*/
	if (energy < 0) {
		printf("electron energy < 0\n");
		printLog("electron energy < 0\n");
		exit(0);
	}
	if (energy < my_mass*speed_of_light2) {
		printf("warning: energy is less than m c^2\n");
		printLog("warning: energy is less than m c^2\n");
		//exit(0);
		return 0;
	}
	return my_A / pow(energy / my_E0, my_index);
}

double MassiveParticlePowerLawDistribution::getMeanEnergy()
{
	return my_A*4*pi*my_E0*my_E0/(my_index - 2);
}

void MassiveParticlePowerLawDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
}

double MassiveParticlePowerLawDistribution::getIndex()
{
	return my_index;
}

double MassiveParticlePowerLawDistribution::getE0()
{
	return my_E0;
}

MassiveParticlePowerLawCutoffDistribution::MassiveParticlePowerLawCutoffDistribution(const double& mass, const double& index, const double& E0, const double& beta, const double& Ecut, const double& concentration) {
	my_mass = mass;
	if (index <= 1.0) {
		printf("electron spectrum index <= 1.0, contains infinit energy\n");
		printLog("electron spectrum index <= 1.0, contains infinit energy\n");
		exit(0);
	}
	my_index = index;
	if (E0 < my_mass * speed_of_light2) {
		printf("particle minimum energy is less than m c^2\n");
		printLog("particle minimum energy is less than m c^2\n");
		exit(0);
	}
	my_E0 = E0;
	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_beta = beta;
	my_Ecut = Ecut;
	if (my_Ecut < my_E0) {
		printf("Ecutoff < E0\n");
		printLog("Ecutoff < E0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_A = (my_index - 1) / (my_E0 * 4 * pi);

	double Emin = my_E0;
	double Emax = 10 * my_Ecut;
	int N = 100000;
	double factor = pow(Emax / Emin, 1.0 / (N - 1));
	double norm = 0;
	double E = Emin;
	for (int i = 0; i < N; ++i) {
		double dE = E * (factor - 1.0);
		norm += 4 * pi * distributionNormalized(E) * dE;
		E = E * factor;
	}
	my_A = my_A / norm;

}

double MassiveParticlePowerLawCutoffDistribution::distributionNormalized(const double& energy) {
	if (energy < 0) {
		printf("electron energy < 0\n");
		printLog("electron energy < 0\n");
		exit(0);
	}
	if (energy < my_mass * speed_of_light2) {
		printf("warning: energy is less than m c^2\n");
		printLog("warning: energy is less than m c^2\n");
		//exit(0);
		return 0;
	}
	return (my_A / pow(energy / my_E0, my_index))*exp(-pow(energy/my_Ecut, my_beta));
}

double MassiveParticlePowerLawCutoffDistribution::getMeanEnergy()
{
	double Emin = my_E0;
	double Emax = 10 * my_Ecut;
	double Emean = 0;
	int N = 100000;
	double factor = pow(Emax / Emin, 1.0 / (N - 1));
	double norm = 0;
	double E = Emin;
	for (int i = 0; i < N; ++i) {
		double dE = E * (factor - 1.0);
		Emean += 4 * pi * distributionNormalized(E) * E * dE;
		E = E * factor;
	}
	return Emean;
}

void MassiveParticlePowerLawCutoffDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
}

double MassiveParticlePowerLawCutoffDistribution::getIndex()
{
	return my_index;
}

double MassiveParticlePowerLawCutoffDistribution::getBeta()
{
	return my_beta;
}

double MassiveParticlePowerLawCutoffDistribution::getE0()
{
	return my_E0;
}

double MassiveParticlePowerLawCutoffDistribution::getEcutoff()
{
	return my_Ecut;
}

MassiveParticleBrokenPowerLawDistribution::MassiveParticleBrokenPowerLawDistribution(const double& mass, const double& index1, const double& index2, const double& E0, const double& Etran, const double& concentration) {
	my_mass = mass;
	if (index2 <= 1.0) {
		printf("electron spectrum index2 <= 1.0, contains infinit energy\n");
		printLog("electron spectrum index2 <= 1.0, contains infinit energy\n");
		exit(0);
	}
	my_index1 = index1;
	my_index2 = index2;
	if (E0 < my_mass * speed_of_light2) {
		printf("particle minimum energy is less than m c^2\n");
		printLog("particle minimum energy is less than m c^2\n");
		exit(0);
	}
	if (Etran <= E0) {
		printf("Etran < E0\n");
		printLog("Etran < E0\n");
		exit(0);
	}
	my_E0 = E0;
	my_Etran = Etran;
	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	double normCoef;
	if (my_index1 == 1.0) {
		normCoef = log(my_Etran / my_E0) + 1.0 / (my_index2 - 1.0);
	} else {
		normCoef = ((1.0 - pow(my_E0 / my_Etran, my_index1 - 1)) / (my_index1 - 1.0)) + pow(my_E0 / my_Etran, my_index1 - 1) / (my_index2 - 1);
	}
	my_A = 1.0 / (my_E0 * 4 * pi * normCoef);
	my_B = my_A * pow(my_E0 / my_Etran, my_index1 - my_index2);
}

double MassiveParticleBrokenPowerLawDistribution::distributionNormalized(const double& energy) {
	if (energy < 0) {
		printf("electron energy < 0\n");
		printLog("electron energy < 0\n");
		exit(0);
	}
	if (energy < my_mass * speed_of_light2) {
		printf("warning: energy is less than m c^2\n");
		printLog("warning: energy is less than m c^2\n");
		//exit(0);
		return 0;
	}
	if (energy < my_Etran) {
		return my_A / pow(energy / my_E0, my_index1);
	} else {
		return my_B / pow(energy / my_E0, my_index2);
	}
}

double MassiveParticleBrokenPowerLawDistribution::getMeanEnergy()
{
	double first = 0;
	if(my_index1 == 2.0){
		first = my_E0 * my_E0 * my_A * log(my_Etran / my_E0);
	}
	else {
		first = my_E0 * my_E0 * my_A * (1.0 - pow(my_E0 / my_Etran, my_index1 - 2)) / (my_index1 - 2.0);
	}
	double second = 0;
	if (my_index2 <= 2.0) {
		printf("MassiveParticleBrokenPowerLawDistribution has infinite energy\n");
		printLog("MassiveParticleBrokenPowerLawDistribution has infinite energy\n");
		exit(0);
	}
	second = my_E0 * my_E0 * my_B * pow(my_E0/my_Etran, my_index2 - 2.0) / (my_index2 - 2.0);
	return first + second;
}

void MassiveParticleBrokenPowerLawDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
}

double MassiveParticleBrokenPowerLawDistribution::getIndex1()
{
	return my_index1;
}

double MassiveParticleBrokenPowerLawDistribution::getIndex2()
{
	return my_index2;
}

double MassiveParticleBrokenPowerLawDistribution::getE0()
{
	return my_E0;
}

double MassiveParticleBrokenPowerLawDistribution::getEtran()
{
	return my_Etran;
}

MassiveParticleMaxwellDistribution::MassiveParticleMaxwellDistribution(const double& mass, const double& temperature, const double& concentration)
{
	my_mass = mass;
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

	my_A = 1.0 / (2 * sqrt(cube(pi * kBoltzman * my_temperature)));
}

double MassiveParticleMaxwellDistribution::distributionNormalized(const double& energy)
{
	return my_A*sqrt(energy)*exp(-energy/(kBoltzman*my_temperature));
}

double MassiveParticleMaxwellDistribution::getMeanEnergy()
{
	return 3.0 * kBoltzman * my_temperature / 2.0;
}

void MassiveParticleMaxwellDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
}

double MassiveParticleMaxwellDistribution::getTemperature()
{
	return my_temperature;
}

MassiveParticleMaxwellJuttnerDistribution::MassiveParticleMaxwellJuttnerDistribution(const double& mass, const double& temperature, const double& concentration)
{
	my_mass = mass;
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

	double m_c2 = my_mass * speed_of_light2;

	double theta = kBoltzman * my_temperature / m_c2;
	my_A = 1.0 / (4*pi*cube(m_c2)*theta*McDonaldFunction(2, 1/theta));
}

double MassiveParticleMaxwellJuttnerDistribution::distributionNormalized(const double& energy)
{
	double m_c2 = my_mass * speed_of_light2;
	if (energy < m_c2) {
		printf("particle energy is less than m c ^ 2 in maxwell - juttner distribution\n");
		printLog("particle energy is less than m c^2 in maxwell-juttner distribution\n");
		exit(0);
	}
	return my_A * sqrt(energy*energy - m_c2*m_c2)*energy*exp(-energy / (kBoltzman * my_temperature));
}

double MassiveParticleMaxwellJuttnerDistribution::getMeanEnergy()
{
	double m_c2 = my_mass * speed_of_light2;
	double Emin = m_c2;
	double Emax = 100 * kBoltzman*my_temperature;
	double Emean = 0;
	int N = 100000;
	double factor = pow(Emax / Emin, 1.0 / (N - 1));
	double norm = 0;
	double E = Emin;
	for (int i = 0; i < N; ++i) {
		double dE = E * (factor - 1.0);
		Emean += 4 * pi * distributionNormalized(E) * E * dE;
		E = E * factor;
	}
	return Emean;
}

void MassiveParticleMaxwellJuttnerDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
}

double MassiveParticleMaxwellJuttnerDistribution::getTemperature()
{
	return my_temperature;
}

void MassiveParticleTabulatedIsotropicDistribution::setDistributionAtPoint(int i, const double& energy, const double& distribution)
{
	double m_c2 = my_mass * speed_of_light2;
	if (my_inputType == DistributionInputType::ENERGY_FE) {
		my_energy[i] = energy;
		my_distribution[i] = distribution;
	} else if (my_inputType == DistributionInputType::ENERGY_KIN_FE) {
		my_energy[i] = energy + m_c2;
		my_distribution[i] = distribution;
	}
	else if (my_inputType == DistributionInputType::GAMMA_FGAMMA) {
		my_energy[i] = energy * m_c2;
		my_distribution[i] = distribution / m_c2;
	}
	else if (my_inputType == DistributionInputType::GAMMA_KIN_FGAMMA) {
		my_energy[i] = (energy + 1) * m_c2;
		my_distribution[i] = distribution / m_c2;
	}
	else if (my_inputType == DistributionInputType::MOMENTUM_FP) {
		my_energy[i] = sqrt(energy * energy * speed_of_light2 + m_c2 * m_c2);
		my_distribution[i] = distribution * energy * my_energy[i] / speed_of_light2;
	}
	else {
		printf("unknown electron input type\n");
		printLog("unknown electron input type\n");
		exit(0);
	}
}

void MassiveParticleTabulatedIsotropicDistribution::normalizeDistribution()
{
	double norm = my_distribution[0] * (my_energy[1] - my_energy[0]);
	for (int i = 1; i < my_Ne; ++i) {
		norm += my_distribution[i] * (my_energy[i] - my_energy[i - 1]);
	}
	norm *= 4 * pi;
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] *= 1.0 / norm;
	}
}

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* fileName, const int N, const double& concentration, DistributionInputType inputType) {
	my_mass = mass;
	
	if (N <= 0) {
		printf("grid number <= 0 in tabulated spherical distribution\n");
		printLog("grid number <= 0 in tabulated spherical distribution\n");
		exit(0);
	}
	my_Ne = N;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_distribution = new double[my_Ne];

	FILE* file = fopen(fileName, "r");
	for (int i = 0; i < my_Ne; ++i) {
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

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* energyFileName, const char* distributionFileName, const int N, const double& concentration, DistributionInputType inputType) {
	my_mass = mass;
	if (N <= 0) {
		printf("grid number <= 0 in tabulated spherical distribution\n");
		printLog("grid number <= 0 in tabulated spherical distribution\n");
		exit(0);
	}
	my_Ne = N;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_distribution = new double[my_Ne];

	FILE* energyFile = fopen(energyFileName, "r");
	FILE* distributionFile = fopen(distributionFileName, "r");
	for (int i = 0; i < my_Ne; ++i) {
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

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const double& mass, const double* energy, const double* distribution, const int N, const double& concentration, DistributionInputType inputType) {
	my_mass = mass;
	if (N <= 0) {
		printf("grid number <= 0 in tabulated spherical distribution\n");
		printLog("grid number <= 0 in tabulated spherical distribution\n");
		exit(0);
	}
	my_Ne = N;

	if (concentration <= 0) {
		printf("electrons concentration <= 0\n");
		printLog("electrons concentration <= 0\n");
		exit(0);
	}
	my_concentration = concentration;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_distribution = new double[my_Ne];

	for (int i = 0; i < my_Ne; ++i) {
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

MassiveParticleTabulatedIsotropicDistribution::~MassiveParticleTabulatedIsotropicDistribution() {
	delete[] my_distribution;
	delete[] my_energy;
}

double MassiveParticleTabulatedIsotropicDistribution::distributionNormalized(const double& energy) {
	if (energy <= my_energy[0]) {
		//printf("warning: energy is less than minimum energy\n");
		//printLog("warning: energy is less than minimum energy\n");
		return 0;
	}
	else if (energy >= my_energy[my_Ne-1]) {
		//printf("warning: energy is greater than maximum energy\n");
		//printLog("warning: energy is greater than maximum energy\n");
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
		int maxI = my_Ne-1;
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

double MassiveParticleTabulatedIsotropicDistribution::getMeanEnergy()
{
	double Emean = 0;
	Emean += 4*pi*my_distribution[0] * (my_energy[1] - my_energy[0]) * (my_energy[1] + my_energy[0])/2.0;
	for (int i = 1; i < my_Ne; ++i) {
		Emean += 4*pi*my_distribution[i] * (my_energy[i] - my_energy[i - 1]) * (my_energy[i] + my_energy[i - 1])/2.0;
	}
	return Emean;
}

void MassiveParticleTabulatedIsotropicDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
	//normalizeDistribution();
}

int MassiveParticleTabulatedIsotropicDistribution::getN() {
	return my_Ne;
}

double MassiveParticleTabulatedIsotropicDistribution::getEmin() {
	return my_energy[0];
}

double MassiveParticleTabulatedIsotropicDistribution::getEmax() {
	return my_energy[my_Ne - 1];
}

void MassiveParticleTabulatedIsotropicDistribution::rescaleDistribution(const double& k)
{
	double m_c2 = my_mass * speed_of_light2;
	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = m_c2 + (my_energy[i] - m_c2) * k;
	}
	normalizeDistribution();
}

void MassiveParticleTabulatedIsotropicDistribution::addPowerLaw(const double& Epower, const double& index)
{
	int ie = my_Ne-1;
	for (int i = 0; i < my_Ne; ++i) {
		if (my_energy[i] > Epower) {
			ie = i;
			break;
		}
	}

	for (int i = ie; i < my_Ne; ++i) {
		my_distribution[i] = my_distribution[ie] * pow(my_energy[ie] / my_energy[i], index);
	}

	normalizeDistribution();
}

void MassiveParticleTabulatedPolarDistribution::setDistributionAtPoint(int i, int j, const double& energy, const double& distribution)
{
	double m_c2 = my_mass * speed_of_light2;
	if (my_inputType == DistributionInputType::ENERGY_FE) {
		my_energy[i] = energy;
		my_distribution[i][j] = distribution;
	}
	else if (my_inputType == DistributionInputType::ENERGY_KIN_FE) {
		my_energy[i] = energy + m_c2;
		my_distribution[i][j] = distribution;
	}
	else if (my_inputType == DistributionInputType::GAMMA_FGAMMA) {
		my_energy[i] = energy * m_c2;
		my_distribution[i][j] = distribution / m_c2;
	}
	else if (my_inputType == DistributionInputType::GAMMA_KIN_FGAMMA) {
		my_energy[i] = (energy + 1) * m_c2;
		my_distribution[i][j] = distribution / m_c2;
	}
	else if (my_inputType == DistributionInputType::MOMENTUM_FP) {
		my_energy[i] = sqrt(energy * energy * speed_of_light2 + m_c2 * m_c2);
		my_distribution[i][j] = distribution * energy * my_energy[i] / speed_of_light2;
	}
	else {
		printf("unknown electron input type\n");
		printLog("unknown electron input type\n");
		exit(0);
	}
}

void MassiveParticleTabulatedPolarDistribution::normalizeDistribution()
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
			my_distribution[i][imu] *= 1.0 / norm;
		}
	}
}

MassiveParticleTabulatedPolarDistribution::MassiveParticleTabulatedPolarDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const double& concentration, DistributionInputType inputType)
{
	my_mass = mass;
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

MassiveParticleTabulatedPolarDistribution::MassiveParticleTabulatedPolarDistribution(const double& mass, const double* energy, const double* mu, const double** distribution, const int Ne, const int Nmu, const double& concentration, DistributionInputType inputType)
{
	my_mass = mass;
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

MassiveParticleTabulatedPolarDistribution::~MassiveParticleTabulatedPolarDistribution()
{
	delete[] my_energy;
	delete[] my_mu;
	for (int i = 0; i < my_Ne; ++i) {
		delete[] my_distribution[i];
	}
	delete[] my_distribution;
}

double MassiveParticleTabulatedPolarDistribution::getMeanEnergy()
{
	double Emean = 0;
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
		Emean += my_distribution[0][imu] * (my_energy[1] - my_energy[0]) * dmu * (my_energy[1] + my_energy[0])/2.0;
		for (int i = 1; i < my_Ne; ++i) {
			Emean += my_distribution[i][imu] * (my_energy[i] - my_energy[i - 1]) * dmu * (my_energy[i] + my_energy[i - 1])/2.0;
		}
	}
	return Emean;
}

double MassiveParticleTabulatedPolarDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
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
		//printf("warning: energy is less than minimum energy\n");
		//printLog("warning: energy is less than minimum energy\n");
		return 0;
	}
	else if (energy >= my_energy[my_Ne - 1]) {
		//printf("warning: energy is greater than maximum energy\n");
		//printLog("warning: energy is greater than maximum energy\n");
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

void MassiveParticleTabulatedPolarDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
	//normalizeDistribution();
}

int MassiveParticleTabulatedPolarDistribution::getNe()
{
	return my_Ne;
}

double MassiveParticleTabulatedPolarDistribution::getEmin() {
	return my_energy[0];
}

double MassiveParticleTabulatedPolarDistribution::getEmax() {
	return my_energy[my_Ne - 1];
}

int MassiveParticleTabulatedPolarDistribution::getNmu()
{
	return my_Nmu;
}

void MassiveParticleTabulatedPolarDistribution::rescaleDistribution(const double& k)
{
	double m_c2 = my_mass * speed_of_light2;
	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = m_c2 + (my_energy[i] - m_c2);
	}
	normalizeDistribution();
}

void MassiveParticleTabulatedAnisotropicDistribution::setDistributionAtPoint(int i, int j, int k, const double& energy, const double& distribution)
{
	double m_c2 = my_mass * speed_of_light2;
	if (my_inputType == DistributionInputType::ENERGY_FE) {
		my_energy[i] = energy;
		my_distribution[i][j][k] = distribution;
	} else if (my_inputType == DistributionInputType::ENERGY_KIN_FE) {
		my_energy[i] = energy + m_c2;
		my_distribution[i][j][k] = distribution;
	}
	else if (my_inputType == DistributionInputType::GAMMA_FGAMMA) {
		my_energy[i] = energy * m_c2;
		my_distribution[i][j][k] = distribution / m_c2;
	}
	else if (my_inputType == DistributionInputType::GAMMA_KIN_FGAMMA) {
		my_energy[i] = (energy + 1) * m_c2;
		my_distribution[i][j][k] = distribution / m_c2;
	}
	else if (my_inputType == DistributionInputType::MOMENTUM_FP) {
		my_energy[i] = sqrt(energy * energy * speed_of_light2 + m_c2 * m_c2);
		my_distribution[i][j][k] = distribution * energy * my_energy[i] / speed_of_light2;
	}
	else {
		printf("unknown electron input type\n");
		printLog("unknown electron input type\n");
		exit(0);
	}
}

void MassiveParticleTabulatedAnisotropicDistribution::normalizeDistribution()
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
				my_distribution[i][imu][iphi] *= 1.0 / norm;
			}
		}
	}
}

MassiveParticleTabulatedAnisotropicDistribution::MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const int Nphi, const double& concentration, DistributionInputType inputType)
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

MassiveParticleTabulatedAnisotropicDistribution::MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const double* energy, const double* mu, const double*** distribution, const int Ne, const int Nmu, const int Nphi, const double& concentration, DistributionInputType inputType)
{
	my_mass = mass;
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

MassiveParticleTabulatedAnisotropicDistribution::~MassiveParticleTabulatedAnisotropicDistribution()
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

double MassiveParticleTabulatedAnisotropicDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	if (mu > 1.0) {
		printf("mu = %lf > 1.0 in azimutal distribution\n", mu);
        printLog("mu = %lf > 1.0 in azimutal distribution\n", mu);
		exit(0);
	}
	if (mu < -1.0) {
		printf("mu = %lf < -1.0 in azimutal distribution\n", mu);
        printLog("mu = %lf < -1.0 in azimutal distribution\n", mu);
		exit(0);
	}
	if (energy <= my_energy[0]) {
		//printf("warning: energy is less than minimum energy\n");
		//printLog("warning: energy is less than minimum energy\n");
		return 0;
	}
	else if (energy >= my_energy[my_Ne - 1]) {
		//printf("warning: energy is greater than maximum energy\n");
		//printLog("warning: energy is greater than maximum energy\n");
		return 0;
	}
	double tempPhi = phi;
	if (phi < 0) {
		printf("warning: phi = %lf < 0\n", phi);
        printLog("warning: phi = %lf < 0\n", phi);
		while (tempPhi < 0) {
			tempPhi += 2 * pi;
		}
	}
	if (phi >= 2 * pi) {
		printf("warning: phi = %lf > 2 pi\n", phi);
        printLog("warning: phi = %lf > 2 pi\n", phi);
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

double MassiveParticleTabulatedAnisotropicDistribution::getMeanEnergy()
{
	double Emean = 0;
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
			Emean += my_distribution[0][imu][k] * (my_energy[1] - my_energy[0]) * dmu * dphi * (my_energy[1] + my_energy[0])/2.0;
			for (int i = 1; i < my_Ne; ++i) {
				Emean += my_distribution[i][imu][k] * (my_energy[i] - my_energy[i - 1]) * dmu * dphi * (my_energy[i] + my_energy[i - 1])/2.0;
			}
		}
	}
	return Emean;
}

void MassiveParticleTabulatedAnisotropicDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
	//normalizeDistribution();
}

int MassiveParticleTabulatedAnisotropicDistribution::getNe()
{
	return my_Ne;
}

double MassiveParticleTabulatedAnisotropicDistribution::getEmin() {
	return my_energy[0];
}

double MassiveParticleTabulatedAnisotropicDistribution::getEmax() {
	return my_energy[my_Ne - 1];
}

int MassiveParticleTabulatedAnisotropicDistribution::getNmu()
{
	return my_Nmu;
}

int MassiveParticleTabulatedAnisotropicDistribution::getNphi()
{
	return my_Nphi;
}

void MassiveParticleTabulatedAnisotropicDistribution::rescaleDistribution(const double& k)
{
	double m_c2 = my_mass * speed_of_light2;
	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = m_c2 + (my_energy[i] - m_c2) * k;
	}
	normalizeDistribution();
}

CompoundMassiveParticleDistribution::CompoundMassiveParticleDistribution(int N, MassiveParticleDistribution** distributions)
{
	my_mass = distributions[0]->getMass();
	my_Ndistr = N;

	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	my_concentration = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		if (distributions[0]->getMass() != my_mass) {
			printf("masses in compound distribution are not equal\n");
			printLog("masses in compound distribution are not equal\n");
			exit(0);
		}
		my_distributions[i] = distributions[i];
		my_concentration += my_distributions[i]->getConcentration();
	}
}

CompoundMassiveParticleDistribution::CompoundMassiveParticleDistribution(MassiveParticleDistribution* dist1, MassiveParticleDistribution* dist2)
{
	my_mass = dist1->getMass();
	if (my_mass != dist2->getMass()) {
		printf("masses in compound distribution are not equal\n");
		printLog("masses in compound distribution are not equal\n");
		exit(0);
	}
	my_Ndistr = 2;
	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	my_concentration = dist1->getConcentration() + dist2->getConcentration();
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
}

CompoundMassiveParticleDistribution::CompoundMassiveParticleDistribution(MassiveParticleDistribution* dist1, MassiveParticleDistribution* dist2, MassiveParticleDistribution* dist3)
{
	my_mass = dist1->getMass();
	if (my_mass != dist2->getMass()) {
		printf("masses in compound distribution are not equal\n");
		printLog("masses in compound distribution are not equal\n");
		exit(0);
	}
	if (my_mass != dist3->getMass()) {
		printf("masses in compound distribution are not equal\n");
		printLog("masses in compound distribution are not equal\n");
		exit(0);
	}
	my_Ndistr = 3;
	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	my_concentration = dist1->getConcentration() + dist2->getConcentration() + dist3->getConcentration();
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_distributions[2] = dist3;
}

CompoundMassiveParticleDistribution::~CompoundMassiveParticleDistribution()
{
	delete[] my_distributions;
}



double CompoundMassiveParticleDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	double result = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		result += my_distributions[i]->distribution(energy, mu, phi);
	}
	return result/my_concentration;
}

double CompoundMassiveParticleDistribution::getMeanEnergy()
{
	double Emean = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		Emean += my_distributions[i]->getConcentration()*my_distributions[i]->getMeanEnergy();
	}
	return Emean / my_concentration;
}

void CompoundMassiveParticleDistribution::resetConcentration(const double& concentration)
{
	double ratio = concentration / my_concentration;
	for (int i = 0; i < my_Ndistr; ++i) {
		my_distributions[i]->resetConcentration(my_distributions[i]->getConcentration() * ratio);
	}
	my_concentration = concentration;
}

CompoundWeightedMassiveParticleDistribution::CompoundWeightedMassiveParticleDistribution(int N, const double* weights, MassiveParticleDistribution** distributions)
{
	my_mass = distributions[0]->getMass();
	my_Ndistr = N;

	my_weights = new double[my_Ndistr];
	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	my_concentration = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		if (distributions[0]->getMass() != my_mass) {
			printf("masses in compound weighted distribution are not equal\n");
			printLog("masses in compound weighted distribution are not equal\n");
			exit(0);
		}
		my_weights[i] = weights[i];
		my_distributions[i] = distributions[i];
		my_concentration += my_weights[i]*my_distributions[i]->getConcentration();
	}
}

CompoundWeightedMassiveParticleDistribution::CompoundWeightedMassiveParticleDistribution(MassiveParticleDistribution* dist1, const double& w1, MassiveParticleDistribution* dist2, const double& w2)
{
	my_mass = dist1->getMass();
	if (my_mass != dist2->getMass()) {
		printf("masses in compound weighted distribution are not equal\n");
		printLog("masses in compound weighted distribution are not equal\n");
		exit(0);
	}
	my_Ndistr = 2;
	my_weights = new double[my_Ndistr];
	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	my_concentration = w1*dist1->getConcentration() + w2*dist2->getConcentration();
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_weights[0] = w1;
	my_weights[1] = w2;
}

CompoundWeightedMassiveParticleDistribution::CompoundWeightedMassiveParticleDistribution(MassiveParticleDistribution* dist1, const double& w1, MassiveParticleDistribution* dist2, const double& w2, MassiveParticleDistribution* dist3, const double& w3)
{
	my_mass = dist1->getMass();
	if (my_mass != dist2->getMass()) {
		printf("masses in compound weighted distribution are not equal\n");
		printLog("masses in compound weighted distribution are not equal\n");
		exit(0);
	}
	if (my_mass != dist3->getMass()) {
		printf("masses in compound weighted distribution are not equal\n");
		printLog("masses in compound weighted distribution are not equal\n");
		exit(0);
	}
	my_Ndistr = 3;
	my_weights = new double[my_Ndistr];
	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	my_concentration = dist1->getConcentration() + dist2->getConcentration() + dist3->getConcentration();
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_distributions[2] = dist3;
	my_weights[0] = w1;
	my_weights[1] = w2;
	my_weights[2] = w3;
}

CompoundWeightedMassiveParticleDistribution::~CompoundWeightedMassiveParticleDistribution()
{
	delete[] my_weights;
	delete[] my_distributions;
}

double CompoundWeightedMassiveParticleDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	double result = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		result += my_weights[i]*my_distributions[i]->distributionNormalized(energy, mu, phi);
	}
	return result;
}

double CompoundWeightedMassiveParticleDistribution::getMeanEnergy()
{
	double Emean = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		Emean += my_weights[i] * my_distributions[i]->getMeanEnergy();
	}
	return Emean;
}

void CompoundWeightedMassiveParticleDistribution::resetConcentration(const double& concentration)
{
	//double ratio = concentration / my_concentration;
	for (int i = 0; i < my_Ndistr; ++i) {
		my_distributions[i]->resetConcentration(concentration);
	}
	my_concentration = concentration;
}

MassiveParticleIsotropicDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentration, int Ne)
{
	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileNameE = energyFileName;
		const std::string fileNameF = distributionFileName;
		std::string a = fileNameE + fileNumber + fileExtension;
		std::string b = fileNameF + fileNumber + fileExtension;
		const char* energyFileName = a.c_str();
		const char* distributionFileName = b.c_str();
		distributions[i] = new MassiveParticleTabulatedIsotropicDistribution(mass, energyFileName, distributionFileName, Ne, electronConcentration, inputType);
	}
	return distributions;
}

MassiveParticleIsotropicDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentration, int Ne)
{
	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileName = fileName;
		std::string a = fileName + fileNumber + fileExtension;
		distributions[i] = new MassiveParticleTabulatedIsotropicDistribution(mass, a.c_str(), Ne, electronConcentration, inputType);
	}
	return distributions;
}

MassiveParticleIsotropicDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentration, int Ne, const double& Epower, const double& index)
{
	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileName = fileName;
		std::string a = fileName + fileNumber + fileExtension;
		MassiveParticleTabulatedIsotropicDistribution* distributionTabulated = new MassiveParticleTabulatedIsotropicDistribution(mass, a.c_str(), Ne, electronConcentration, inputType);
		distributionTabulated->addPowerLaw(Epower, index);
		distributions[i] = distributionTabulated;
	}
	return distributions;
}

MassiveParticleIsotropicDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, const double& electronConcentration, int Ne, const double& Epower, const double& index)
{
	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileNameE = energyFileName;
		const std::string fileNameF = distributionFileName;
		std::string a = fileNameE + fileNumber + fileExtension;
		std::string b = fileNameF + fileNumber + fileExtension;
		const char* energyFileName = a.c_str();
		const char* distributionFileName = b.c_str();
		MassiveParticleTabulatedIsotropicDistribution* distributionTabulated = new MassiveParticleTabulatedIsotropicDistribution(mass, energyFileName, distributionFileName, Ne, electronConcentration, inputType);
		distributionTabulated->addPowerLaw(Epower, index);
		distributions[i] = distributionTabulated;
	}
	return distributions;
}
