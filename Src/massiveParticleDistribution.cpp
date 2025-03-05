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
		fprintf(outFile, "%g %g\n", energy[i], distributionNormalized(energy[i]));
	}
	fclose(outFile);
	delete[] energy;
}

double MassiveParticleIsotropicDistribution::distributionNormalizedWithLosses(const double& energy, const double& lossRate, const double& time)
{
	double factor = energy * lossRate * time;
	if (factor >= 1) {
		return 0;
	}
	return distributionNormalized(energy / (1 - factor)) / sqr(1 - factor);
}

MassiveParticlePowerLawDistribution::MassiveParticlePowerLawDistribution(const double& mass, const double& index, const double& E0) {
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

	if (energy < my_E0) {
		return 0;
	}
	return my_A / pow(energy / my_E0, my_index);
}

double MassiveParticlePowerLawDistribution::getMeanEnergy()
{
	return my_A*4*pi*my_E0*my_E0/(my_index - 2);
}

double MassiveParticlePowerLawDistribution::minEnergy()
{
	return my_E0;
}

double MassiveParticlePowerLawDistribution::maxEnergy()
{
	return -1.0;
}

double MassiveParticlePowerLawDistribution::getIndex()
{
	return my_index;
}

double MassiveParticlePowerLawDistribution::getE0()
{
	return my_E0;
}

MassiveParticlePowerLawCutoffDistribution::MassiveParticlePowerLawCutoffDistribution(const double& mass, const double& index, const double& E0, const double& beta, const double& Ecut) {
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
	my_beta = beta;
	my_Ecut = Ecut;
	if (my_Ecut < my_E0) {
		printf("Ecutoff < E0\n");
		printLog("Ecutoff < E0\n");
		exit(0);
	}

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
	if (energy < my_E0) {
		return 0;
	}
	if (pow(energy / my_Ecut, my_beta) > 300) {
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

double MassiveParticlePowerLawCutoffDistribution::minEnergy()
{
	return my_E0;
}

double MassiveParticlePowerLawCutoffDistribution::maxEnergy()
{
	return -1.0;
}

void MassiveParticlePowerLawCutoffDistribution::resetEcut(const double& Ecut)
{
	//todo need to renorm?
	my_Ecut = Ecut;
	if (my_Ecut < my_E0) {
		printf("Ecutoff < E0\n");
		printLog("Ecutoff < E0\n");
		exit(0);
	}
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

MassiveParticleBrokenPowerLawDistribution::MassiveParticleBrokenPowerLawDistribution(const double& mass, const double& index1, const double& index2, const double& E0, const double& Etran) {
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
	if (energy < my_E0) {
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
		first = my_E0 * my_E0 * my_A * (1.0 - pow(my_E0 / my_Etran, my_index1 - 2.0)) / (my_index1 - 2.0);
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

double MassiveParticleBrokenPowerLawDistribution::minEnergy()
{
	return my_E0;
}

double MassiveParticleBrokenPowerLawDistribution::maxEnergy()
{
	return -1.0;
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

MassiveParticleMaxwellDistribution::MassiveParticleMaxwellDistribution(const double& mass, const double& temperature)
{
	my_mass = mass;
	if (temperature <= 0) {
		printf("electrons temperature <= 0\n");
		printLog("electrons temperature <= 0\n");
		exit(0);
	}
	my_temperature = temperature;

	my_A = 1.0 / (2 * sqrt(cube(pi * kBoltzman * my_temperature)));
}

double MassiveParticleMaxwellDistribution::distributionNormalized(const double& energy)
{
	return my_A*sqrt((energy- my_mass*speed_of_light2))*exp(-(energy-my_mass*speed_of_light2)/(kBoltzman*my_temperature));
}

double MassiveParticleMaxwellDistribution::getMeanEnergy()
{
	return my_mass*speed_of_light2*3.0 * kBoltzman * my_temperature / 2.0;
}

double MassiveParticleMaxwellDistribution::minEnergy()
{
	return my_mass*speed_of_light2;
}

double MassiveParticleMaxwellDistribution::maxEnergy()
{
	return -1.0;
}

double MassiveParticleMaxwellDistribution::getTemperature()
{
	return my_temperature;
}

MassiveParticleMaxwellJuttnerDistribution::MassiveParticleMaxwellJuttnerDistribution(const double& mass, const double& temperature)
{
	my_mass = mass;
	if (temperature <= 0) {
		printf("electrons temperature <= 0\n");
		printLog("electrons temperature <= 0\n");
		exit(0);
	}
	my_temperature = temperature;

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

double MassiveParticleMaxwellJuttnerDistribution::minEnergy()
{
	return my_mass*speed_of_light2;
}

double MassiveParticleMaxwellJuttnerDistribution::maxEnergy()
{
	return -1.0;
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
	my_kineticEnergy[i] = my_energy[i] - my_mass * speed_of_light2;
}

void MassiveParticleTabulatedIsotropicDistribution::normalizeDistribution()
{
	//double norm = my_distribution[0] * (my_energy[1] - my_energy[0]);
	double norm = my_distribution[0] * (my_kineticEnergy[1] - my_kineticEnergy[0]);
	for (int i = 1; i < my_Ne; ++i) {
		norm += my_distribution[i] * (my_kineticEnergy[i] - my_kineticEnergy[i - 1]);
	}
	norm *= 4 * pi;
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] *= 1.0 / norm;
	}
}

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const MassiveParticleTabulatedIsotropicDistribution& distribution) {
	my_mass = distribution.my_mass;
	my_Ne = distribution.my_Ne;

	my_energy = new double[my_Ne];
	my_kineticEnergy = new double[my_Ne];
	my_distribution = new double[my_Ne];
	my_inputType = distribution.my_inputType;

	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = distribution.my_energy[i];
		my_kineticEnergy[i] = my_energy[i] - my_mass * speed_of_light2;
		my_distribution[i] = distribution.my_distribution[i];
	}

	normalizeDistribution();
}

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* fileName, DistributionInputType inputType) {
	my_mass = mass;


	int n = 0;
	FILE* file = fopen(fileName, "r");
	while (!feof(file)) {
		double a;
		double b;
		fscanf(file, "%lf %lf", &a, &b);
		n = n + 1;
	}
	fclose(file);
	my_Ne = n - 1;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_kineticEnergy = new double[my_Ne];
	my_distribution = new double[my_Ne];

	file = fopen(fileName, "r");
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

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const double& mass, const char* energyFileName, const char* distributionFileName, DistributionInputType inputType) {
	my_mass = mass;

	int n = 0;
	FILE* energyFile = fopen(energyFileName, "r");
	while (!feof(energyFile)) {
		double a;
		fscanf(energyFile, "%lf", &a);
		n = n + 1;
	}
	fclose(energyFile);
	my_Ne = n-1;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_kineticEnergy = new double[my_Ne];
	my_distribution = new double[my_Ne];

	energyFile = fopen(energyFileName, "r");
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

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(const double& mass, const double* energy, const double* distribution, const int N, DistributionInputType inputType) {
	my_mass = mass;
	if (N <= 0) {
		printf("grid number <= 0 in tabulated spherical distribution\n");
		printLog("grid number <= 0 in tabulated spherical distribution\n");
		exit(0);
	}
	my_Ne = N;

	my_inputType = inputType;

	my_energy = new double[my_Ne];
	my_kineticEnergy = new double[my_Ne];
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

MassiveParticleTabulatedIsotropicDistribution::MassiveParticleTabulatedIsotropicDistribution(MassiveParticleIsotropicDistribution* distribution, const double& minE, const double& maxE, const int N)
{
	my_mass = distribution->getMass();
	if (N <= 0) {
		printf("grid number <= 0 in tabulated spherical distribution\n");
		printLog("grid number <= 0 in tabulated spherical distribution\n");
		exit(0);
	}
	my_Ne = N;

	my_inputType = DistributionInputType::ENERGY_FE;

	my_energy = new double[my_Ne];
	my_kineticEnergy = new double[my_Ne];
	my_distribution = new double[my_Ne];
	double factor = pow(maxE / minE, 1.0 / (my_Ne - 1.0));
	double x = minE;
	double y = distribution->distributionNormalized(x);
	setDistributionAtPoint(0, x, y);
	for (int i = 1; i < my_Ne; ++i) {
		x = x * factor;
		y = distribution->distributionNormalized(x);
		setDistributionAtPoint(i, x, y);
	}

	normalizeDistribution();
}

MassiveParticleTabulatedIsotropicDistribution::~MassiveParticleTabulatedIsotropicDistribution() {
	delete[] my_distribution;
	delete[] my_energy;
	delete[] my_kineticEnergy;
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
			result = my_distribution[currentI] * exp(log(my_distribution[nextI] / my_distribution[currentI]) * ((log(energy/my_energy[currentI])) / (log(my_energy[nextI]/my_energy[currentI]))));
			//result = (my_distribution[currentI] * (my_energy[nextI] - energy) + my_distribution[nextI] * (energy - my_energy[currentI])) / (my_energy[nextI] - my_energy[currentI]);
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

double MassiveParticleTabulatedIsotropicDistribution::minEnergy()
{
	return my_energy[0];
}

double MassiveParticleTabulatedIsotropicDistribution::maxEnergy()
{
	return my_energy[my_Ne-1];
}

int MassiveParticleTabulatedIsotropicDistribution::getN() {
	return my_Ne;
}

void MassiveParticleTabulatedIsotropicDistribution::rescaleDistribution(const double& k)
{
	double m_c2 = my_mass * speed_of_light2;
	for (int i = 0; i < my_Ne; ++i) {
		my_kineticEnergy[i] = my_kineticEnergy[i] * k;
		my_energy[i] = m_c2 + my_kineticEnergy[i];
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

void MassiveParticleTabulatedIsotropicDistribution::prolongEnergyRange(const double& Emax, int N)
{
	double* tempEnergy = my_energy;
	double* tempKineticEnergy = my_kineticEnergy;
	double* tempDistribution = my_distribution;
	int tempNe = my_Ne;

	my_Ne = my_Ne + N;
	my_energy = new double[my_Ne];
	my_kineticEnergy = new double[my_Ne];
	my_distribution = new double[my_Ne];
	for (int i = 0; i < tempNe; ++i) {
		my_energy[i] = tempEnergy[i];
		my_kineticEnergy[i] = tempKineticEnergy[i];
		my_distribution[i] = tempDistribution[i];
	}

	double factor = pow(Emax / my_energy[tempNe - 1], 1.0 / N);
	for (int i = tempNe; i < my_Ne; ++i) {
		my_energy[i] = my_energy[i - 1] * factor;
		my_kineticEnergy[i] = my_energy[i] - my_mass * speed_of_light2;
		my_distribution[i] = 0;
	}
	
	delete[] tempEnergy;
	delete[] tempKineticEnergy;
	delete[] tempDistribution;
}

void MassiveParticleTabulatedIsotropicDistribution::addExponentialCutoff(const double& E)
{
	if (E < 0) {
		printf("MassiveParticleTabulatedIsotropicDistribution::addExponentialCutoff E < 0\n");
		printLog("MassiveParticleTabulatedIsotropicDistribution::addExponentialCutoff E < 0\n");
		exit(0);
	}

	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = my_distribution[i] * exp(-my_energy[i] / E);
	}

	normalizeDistribution();
}

void MassiveParticleTabulatedIsotropicDistribution::setToZeroAboveE(const double& E) {
	if (E <= my_energy[0]) {
		printf("MassiveParticleTabulatedIsotropicDistribution: canot set whole distribution to zero\n");
		printLog("MassiveParticleTabulatedIsotropicDistribution: canot set whole distribution to zero\n");
		exit(0);
	}

	int index = 1;
	for (int i = 1; i < my_Ne; ++i) {
		if (my_energy[i] > E) {
			index = i;
			break;
		}
	}

	for (int i = index; i < my_Ne; ++i) {
		my_distribution[i] = 0;
	}

	normalizeDistribution();
}

void MassiveParticleTabulatedIsotropicDistribution::transformToLosses(const double& lossRate, const double& time)
{
	for (int i = 0; i < my_Ne; ++i) {
		//double factor = (my_energy[i] - my_mass * speed_of_light2) * lossRate * time;
		double factor = my_kineticEnergy[i] * lossRate * time;
		my_kineticEnergy[i] = my_kineticEnergy[i] / (1 + factor);
		my_energy[i] = my_kineticEnergy[i] + my_mass * speed_of_light2;
		my_distribution[i] = my_distribution[i] * sqr(1 + factor);
		if (my_distribution[i] != my_distribution[i]) {
			printf("distribution = NaN in transformToLosses\n");
			printLog("distribution = NaN in transformToLosses\n");
			exit(0);
		}
	}
	normalizeDistribution();
}

//todo what to do with concentration?
void MassiveParticleTabulatedIsotropicDistribution::transformToLosses2(const double& k, const double& l1, const double& l2)
{
	//l1 > l2!!!
	double* tempDistribution = new double[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		tempDistribution[i] = 0;
	}
	for (int i = 0; i < my_Ne; ++i) {
		double E = my_energy[i];
		if (k * l2 * E > 1.0) {
			tempDistribution[i] = 0;
		}
		else {
			if(l1 < 1.5*l2){
				double tempE = E / (1 - k * l2 * E);
				tempDistribution[i] = distributionNormalized(tempE) / sqr(1 - k * l2 * E);
			}
			else {
				double Emin = E / (1 - k * l2 * E);
				double Emax;
				if (k * l1 * E > 1.0) {
					Emax = my_energy[my_Ne - 1];
				}
				else {
					Emax = E / (1 - k * l1 * E);
				}
				int N = 100;
				double dE = (Emax - Emin) / N;
				double factor = pow(Emax / Emin, 1.0 / (N - 1));
				double tempE = Emin;
				for (int j = 0; j < N; ++j) {
					tempDistribution[i] += distributionNormalized(tempE) * tempE*(factor-1) / ((l1 - l2)*k * tempE * tempE);
					tempE = tempE * factor;
				}
			}
			//double factor = (my_energy[i] - my_mass * speed_of_light2) * lossRate * time;
			//my_energy[i] = my_mass * speed_of_light2 + (my_energy[i] - my_mass * speed_of_light2) / (1 + factor);
			//my_distribution[i] = my_distribution[i] * sqr(1 + factor);
		}
	}
	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = tempDistribution[i];
	}
	delete[] tempDistribution;
	normalizeDistribution();
}

void MassiveParticleTabulatedIsotropicDistribution::transformToThickRegime(const double& Uph)
{
	double sigmaT = (8.0 * pi / 3.0) * sqr(electron_charge*electron_charge / (my_mass * speed_of_light2));
	double coef = (4.0 / 3.0) * sigmaT / (my_mass * my_mass * speed_of_light * speed_of_light * speed_of_light);
	double* tempDistribution = new double[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		double E = my_energy[i];
		double intJdE = 0;
		for (int j = i; j < my_Ne; ++j) {
			double dE = 0;
			if (j == 0) {
				//dE = my_energy[1] - my_energy[0];
				dE = my_kineticEnergy[1] - my_kineticEnergy[0];
			}
			else {
				//dE = my_energy[j] - my_energy[j - 1];
				dE = my_kineticEnergy[j] - my_kineticEnergy[j - 1];
			}

			intJdE = intJdE + my_distribution[j] * dE;
			if ((intJdE != intJdE) || (0*intJdE != 0*intJdE)) {
				printf("intJdE = NaN\n");
				printLog("intJdE = NaN\n");
				exit(0);
			}
		}
		tempDistribution[i] = intJdE / (coef * E * E);
		if ((tempDistribution[i] != tempDistribution[i]) || (0 * tempDistribution[i] != 0 * tempDistribution[i])) {
			printf("temp distribution = NaN in transform to thick regime\n");
			printLog("temp distribution = NaN in transform to thick regime\n");
			exit(0);
		}
	}

	for (int i = 0; i < my_Ne; ++i) {
		my_distribution[i] = tempDistribution[i];
	}
	delete[] tempDistribution;
	normalizeDistribution();
}

double* MassiveParticleTabulatedIsotropicDistribution::getEnergyArray()
{
	double* energy = new double[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		energy[i] = my_energy[i];
	}

	return energy;
}

double* MassiveParticleTabulatedIsotropicDistribution::getDistributionArray()
{
	double* distribution = new double[my_Ne];
	for (int i = 0; i < my_Ne; ++i) {
		distribution[i] = my_distribution[i];
	}

	return distribution;
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

MassiveParticleTabulatedPolarDistribution::MassiveParticleTabulatedPolarDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, DistributionInputType inputType)
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

MassiveParticleTabulatedPolarDistribution::MassiveParticleTabulatedPolarDistribution(const double& mass, const double* energy, const double* mu, const double** distribution, const int Ne, const int Nmu, DistributionInputType inputType)
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

double MassiveParticleTabulatedPolarDistribution::minEnergy()
{
	return my_energy[0];
}

double MassiveParticleTabulatedPolarDistribution::maxEnergy()
{
	return my_energy[my_Ne-1];
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

int MassiveParticleTabulatedPolarDistribution::getNe()
{
	return my_Ne;
}

int MassiveParticleTabulatedPolarDistribution::getNmu()
{
	return my_Nmu;
}

void MassiveParticleTabulatedPolarDistribution::rescaleDistribution(const double& k)
{
	double m_c2 = my_mass * speed_of_light2;
	for (int i = 0; i < my_Ne; ++i) {
		my_energy[i] = m_c2 + (my_energy[i] - m_c2)*k;
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

MassiveParticleTabulatedAnisotropicDistribution::MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const char* energyFileName, const char* muFileName, const char* distributionFileName, const int Ne, const int Nmu, const int Nphi, DistributionInputType inputType)
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

MassiveParticleTabulatedAnisotropicDistribution::MassiveParticleTabulatedAnisotropicDistribution(const double& mass, const double* energy, const double* mu, const double*** distribution, const int Ne, const int Nmu, const int Nphi, DistributionInputType inputType)
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

double MassiveParticleTabulatedAnisotropicDistribution::minEnergy()
{
	return my_energy[0];
}

double MassiveParticleTabulatedAnisotropicDistribution::maxEnergy()
{
	return my_energy[my_Ne-1];
}

int MassiveParticleTabulatedAnisotropicDistribution::getNe()
{
	return my_Ne;
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

CompoundWeightedMassiveParticleDistribution::CompoundWeightedMassiveParticleDistribution(int N, const double* weights, MassiveParticleDistribution** distributions)
{
	my_mass = distributions[0]->getMass();
	my_Ndistr = N;

	my_weights = new double[my_Ndistr];
	my_distributions = new MassiveParticleDistribution * [my_Ndistr];
	for (int i = 0; i < my_Ndistr; ++i) {
		if (distributions[0]->getMass() != my_mass) {
			printf("masses in compound weighted distribution are not equal\n");
			printLog("masses in compound weighted distribution are not equal\n");
			exit(0);
		}
		my_weights[i] = weights[i];
		my_distributions[i] = distributions[i];
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
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_distributions[2] = dist3;
	my_weights[0] = w1;
	my_weights[1] = w2;
	my_weights[2] = w3;
}

CompoundWeightedMassiveParticleDistribution::~CompoundWeightedMassiveParticleDistribution()
{
	//todo delete or not?
	delete[] my_weights;
	//delete[] my_distributions;
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

double CompoundWeightedMassiveParticleDistribution::minEnergy()
{
	double tempEmin = my_distributions[0]->minEnergy();
	for (int i = 0; i < my_Ndistr; ++i) {
		double emin = my_distributions[i]->minEnergy();
		if (tempEmin > emin) {
			tempEmin = emin;
		}
	}
	return tempEmin;
}

double CompoundWeightedMassiveParticleDistribution::maxEnergy()
{
	double tempEmax = my_distributions[0]->maxEnergy();
	for (int i = 0; i < my_Ndistr; ++i) {
		double emax = my_distributions[i]->maxEnergy();
		if (emax < 0) {
			return -1.0;
		}
		if (tempEmax < emax) {
			tempEmax = emax;
		}
	}
	return tempEmax;
}

CompoundWeightedMassiveParticleIsotropicDistribution::CompoundWeightedMassiveParticleIsotropicDistribution(int N, const double* weights, MassiveParticleIsotropicDistribution** distributions)
{
	my_mass = distributions[0]->getMass();
	my_Ndistr = N;

	my_weights = new double[my_Ndistr];
	my_distributions = new MassiveParticleIsotropicDistribution * [my_Ndistr];
	for (int i = 0; i < my_Ndistr; ++i) {
		if (distributions[0]->getMass() != my_mass) {
			printf("masses in compound weighted distribution are not equal\n");
			printLog("masses in compound weighted distribution are not equal\n");
			exit(0);
		}
		my_weights[i] = weights[i];
		my_distributions[i] = distributions[i];
	}
}

CompoundWeightedMassiveParticleIsotropicDistribution::CompoundWeightedMassiveParticleIsotropicDistribution(MassiveParticleIsotropicDistribution* dist1, const double& w1, MassiveParticleIsotropicDistribution* dist2, const double& w2)
{
	my_mass = dist1->getMass();
	if (my_mass != dist2->getMass()) {
		printf("masses in compound weighted distribution are not equal\n");
		printLog("masses in compound weighted distribution are not equal\n");
		exit(0);
	}
	my_Ndistr = 2;
	my_weights = new double[my_Ndistr];
	my_distributions = new MassiveParticleIsotropicDistribution * [my_Ndistr];
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_weights[0] = w1;
	my_weights[1] = w2;
}

CompoundWeightedMassiveParticleIsotropicDistribution::CompoundWeightedMassiveParticleIsotropicDistribution(MassiveParticleIsotropicDistribution* dist1, const double& w1, MassiveParticleIsotropicDistribution* dist2, const double& w2, MassiveParticleIsotropicDistribution* dist3, const double& w3)
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
	my_distributions = new MassiveParticleIsotropicDistribution * [my_Ndistr];
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_distributions[2] = dist3;
	my_weights[0] = w1;
	my_weights[1] = w2;
	my_weights[2] = w3;
}

CompoundWeightedMassiveParticleIsotropicDistribution::~CompoundWeightedMassiveParticleIsotropicDistribution()
{
	delete[] my_weights;
	delete[] my_distributions;
}

double CompoundWeightedMassiveParticleIsotropicDistribution::distributionNormalized(const double& energy)
{
	double result = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		result += my_weights[i] * my_distributions[i]->distributionNormalized(energy);
	}
	return result;
}

double CompoundWeightedMassiveParticleIsotropicDistribution::getMeanEnergy()
{
	double Emean = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		Emean += my_weights[i] * my_distributions[i]->getMeanEnergy();
	}
	return Emean;
}

double CompoundWeightedMassiveParticleIsotropicDistribution::minEnergy()
{
	double tempEmin = my_distributions[0]->minEnergy();
	for (int i = 0; i < my_Ndistr; ++i) {
		double emin = my_distributions[i]->minEnergy();
		if (tempEmin > emin) {
			tempEmin = emin;
		}
	}
	return tempEmin;
}

double CompoundWeightedMassiveParticleIsotropicDistribution::maxEnergy()
{
	double tempEmax = my_distributions[0]->maxEnergy();
	for (int i = 0; i < my_Ndistr; ++i) {
		double emax = my_distributions[i]->maxEnergy();
		if (emax < 0) {
			return -1.0;
		}
		if (tempEmax < emax) {
			tempEmax = emax;
		}
	}
	return tempEmax;
}

double MassiveParticleDistributionFactory::evaluateNorm(double* energy, double* distribution, int Ne)
{
	double norm = distribution[0] * (energy[1] - energy[0]);
	for (int i = 1; i < Ne; ++i) {
		norm += distribution[i] * (energy[i] - energy[i - 1]);
	}
	return norm;
}

void MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(const double& mass, const char* fileName, MassiveParticleTabulatedIsotropicDistribution*& outputDistribution, double& outputConcentration)
{
	int n = 0;
	FILE* file = fopen(fileName, "r");
	while (!feof(file)) {
		double a;
		double b;
		fscanf(file, "%lf %lf", &a, &b);
		n = n + 1;
	}
	fclose(file);
	n = n - 1;


	double* energy = new double[n];
	double* distribution = new double[n];

	file = fopen(fileName, "r");
	for (int i = 0; i < n; ++i) {
		double x;
		double y;
		fscanf(file, "%lf", &x);
		fscanf(file, "%lf", &y);
		if (y < 0) {
			printf("input distribution < 0\n");
			printLog("input distribution < 0\n");
			exit(0);
		}
		double m_c2 = mass * speed_of_light2;


			//energy[i] = sqrt(x * x * speed_of_light2 + m_c2 * m_c2);
			double p = x * massProton * speed_of_light;
			energy[i] = sqrt(p * p * speed_of_light2 + m_c2 * m_c2);
			//distribution[i] = y * x * energy[i] / speed_of_light2;
			distribution[i] = (y * p * energy[i] / speed_of_light2)/(x*x*x*x*cube(massProton*speed_of_light));

	}
	fclose(file);

	double norm = distribution[0] * (energy[1] - energy[0]);
	for (int i = 1; i < n; ++i) {
		norm += distribution[i] * (energy[i] - energy[i - 1]);
	}
	norm *= 4 * pi;
	for (int i = 0; i < n; ++i) {
		distribution[i] *= 1.0 / norm;
	}

	outputConcentration = norm;

	outputDistribution = new MassiveParticleTabulatedIsotropicDistribution(mass, energy, distribution, n, DistributionInputType::ENERGY_FE);
}

void MassiveParticleDistributionFactory::readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(const double& mass, const char* xfileName, const char* pfileName, const char* distributionFileName, double* & xgrid, double* & energy, double** & distributions, double* & concentration, int& Nenergy, int& Nx)
{
	Nx = 0;
	FILE* xfile = fopen(xfileName, "r");
	while (!feof(xfile)) {
		double a;
		fscanf(xfile, "%lf", &a);
		Nx = Nx + 1;
	}
	fclose(xfile);
	Nx = Nx - 1;

	Nenergy = 0;
	FILE* pfile = fopen(pfileName, "r");
	while (!feof(xfile)) {
		double a;
		fscanf(pfile, "%lf", &a);
		Nenergy = Nenergy + 1;
	}
	fclose(pfile);
	Nenergy = Nenergy - 1;

	energy = new double[Nenergy];
	xgrid = new double[Nx];
	concentration = new double[Nx];
	distributions = new double* [Nx];
	for (int i = 0; i < Nx; ++i) {
		distributions[i] = new double[Nenergy];
	}

	xfile = fopen(xfileName, "r");
	pfile = fopen(pfileName, "r");
	FILE* distributionFile = fopen(distributionFileName, "r");
	for (int j = 0; j < Nenergy; ++j) {
		double p;
		fscanf(pfile, "%lf", &p);
		double m_c2 = mass * speed_of_light2;
		p = p * massProton * speed_of_light;
		energy[j] = sqrt(p * p * speed_of_light2 + m_c2 * m_c2);
	}
	for (int i = 0; i < Nx; ++i) {
		fscanf(xfile, "%lf", &xgrid[i]);


		//energy[i] = sqrt(x * x * speed_of_light2 + m_c2 * m_c2);
		//distribution[i] = y * x * energy[i] / speed_of_light2;
		for (int j = 0; j < Nenergy; ++j) {
			double m_c2 = mass * speed_of_light2;
			double p = sqrt(energy[j] * energy[j] - m_c2 * m_c2) / speed_of_light;
			double x = p / (massProton * speed_of_light);
			double a;
			fscanf(distributionFile, "%lf %lf", &a, &distributions[i][j]);
			distributions[i][j] = (distributions[i][j] * p * energy[j] / speed_of_light2) / (x * x * x * x * cube(massProton * speed_of_light));
		}

		double norm = distributions[i][0] * (energy[1] - energy[0]);
		for (int j = 1; j < Nenergy; ++j) {
			norm += distributions[i][j] * (energy[j] - energy[j - 1]);
		}
		norm *= 4 * pi;
		if (norm > 0) {
			for (int j = 0; j < Nenergy; ++j) {
				distributions[i][j] *= 1.0 / norm;
			}
			concentration[i] = norm;
		}
		else {
			concentration[i] = 0;
			distributions[i][0] = 1.0 / (energy[2] - energy[0]);
			distributions[i][1] = 1.0 / (energy[2] - energy[0]);
		}
	}
	fclose(xfile);
	fclose(pfile);
	fclose(distributionFile);
}

MassiveParticleDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne)
{
	MassiveParticleDistribution** distributions = new MassiveParticleDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileNameE = energyFileName;
		const std::string fileNameF = distributionFileName;
		std::string a = fileNameE + fileNumber + fileExtension;
		std::string b = fileNameF + fileNumber + fileExtension;
		const char* energyFileName = a.c_str();
		const char* distributionFileName = b.c_str();
		distributions[i] = new MassiveParticleTabulatedIsotropicDistribution(mass, energyFileName, distributionFileName, inputType);
	}
	return distributions;
}

MassiveParticleDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne)
{
	MassiveParticleDistribution** distributions = new MassiveParticleDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileName = fileName;
		std::string a = fileName + fileNumber + fileExtension;
		distributions[i] = new MassiveParticleTabulatedIsotropicDistribution(mass, a.c_str(), inputType);
	}
	return distributions;
}

MassiveParticleDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* fileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne, const double& Epower, const double& index)
{
	MassiveParticleDistribution** distributions = new MassiveParticleDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileName = fileName;
		std::string a = fileName + fileNumber + fileExtension;
		MassiveParticleTabulatedIsotropicDistribution* distributionTabulated = new MassiveParticleTabulatedIsotropicDistribution(mass, a.c_str(), inputType);
		distributionTabulated->addPowerLaw(Epower, index);
		distributions[i] = distributionTabulated;
	}
	return distributions;
}

MassiveParticleDistribution** MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(const double& mass, const char* energyFileName, const char* distributionFileName, const char* fileExtension, int Nfiles, DistributionInputType inputType, int Ne, const double& Epower, const double& index)
{
	MassiveParticleDistribution** distributions = new MassiveParticleDistribution * [Nfiles];
	for (int i = 0; i < Nfiles; ++i) {
		std::string fileNumber = convertIntToString(i);
		const std::string fileNameE = energyFileName;
		const std::string fileNameF = distributionFileName;
		std::string a = fileNameE + fileNumber + fileExtension;
		std::string b = fileNameF + fileNumber + fileExtension;
		const char* energyFileName = a.c_str();
		const char* distributionFileName = b.c_str();
		MassiveParticleTabulatedIsotropicDistribution* distributionTabulated = new MassiveParticleTabulatedIsotropicDistribution(mass, energyFileName, distributionFileName, inputType);
		distributionTabulated->addPowerLaw(Epower, index);
		distributions[i] = distributionTabulated;
	}
	return distributions;
}

void MassiveParticleDistributionFactory::evaluateDistributionAdvectionWithLosses(const double& mass, double* energy, double* distribution, double** outEnergy, double** outDistribution, const int Ne, const int Nx, double* x, double advectionV, double* B, double Uph1, double Eph1, double Uph2, double Eph2)
{
	double dx;
	dx = x[0];

	//double* tempEnergy = new double[Ne];
	//double* tempDistribution = new double[Ne];
	//double* dEold = new double[Ne];
	//double* dE = new double[Ne];

	double mc2 = mass * speed_of_light * speed_of_light;

	double r2 = sqr(electron_charge * electron_charge / mc2);
	double sigmaT = 8 * pi * r2 / 3.0;
	
	double maxLossRate = ((4.0 / 3.0) * sigmaT * energy[Ne-1] * energy[Ne-1] / (mass * mass * cube(speed_of_light))) * (B[0] * B[0] / (8 * pi) + Uph1 / pow(1 + 4 * energy[Ne-1] * Eph1 / sqr(mc2), 3.0 / 2.0) + Uph2 / pow(1 + 4 * energy[Ne-1] * Eph2 / sqr(mc2), 3.0 / 2.0));
	double maxLossRate2 = ((4.0 / 3.0) * sigmaT * energy[Ne-2] * energy[Ne-2] / (mass * mass * cube(speed_of_light))) * (B[0] * B[0] / (8 * pi) + Uph1 / pow(1 + 4 * energy[Ne-2] * Eph1 / sqr(mc2), 3.0 / 2.0) + Uph2 / pow(1 + 4 * energy[Ne-2] * Eph2 / sqr(mc2), 3.0 / 2.0));
	double maxdE = energy[Ne - 1] - energy[Ne - 2];
	double dt = dx / advectionV;
	int Nstep = 10 * floor((maxLossRate - maxLossRate2) * dt / maxdE) + 3;
	evaluateDistributionAdvectionWithLossesBetweenTwoPoints(mass, energy, distribution, outEnergy[0], outDistribution[0], Ne, dx, advectionV, B[0], Uph1, Eph1, Uph2, Eph2, Nstep);

	for (int j = 1; j < Nx; ++j) {
		maxLossRate = ((4.0 / 3.0) * sigmaT * outEnergy[j-1][Ne - 1] * outEnergy[j-1][Ne - 1] / (mass * mass * cube(speed_of_light))) * (B[j] * B[j] / (8 * pi) + Uph1 / pow(1 + 4 * outEnergy[j-1][Ne - 1] * Eph1 / sqr(mc2), 3.0 / 2.0) + Uph2 / pow(1 + 4 * outEnergy[j-1][Ne - 1] * Eph2 / sqr(mc2), 3.0 / 2.0));
		maxLossRate2 = ((4.0 / 3.0) * sigmaT * outEnergy[j-1][Ne - 2] * outEnergy[j-1][Ne - 2] / (mass * mass * cube(speed_of_light))) * (B[j] * B[j] / (8 * pi) + Uph1 / pow(1 + 4 * outEnergy[j-1][Ne - 2] * Eph1 / sqr(mc2), 3.0 / 2.0) + Uph2 / pow(1 + 4 * outEnergy[j-1][Ne - 2] * Eph2 / sqr(mc2), 3.0 / 2.0));
		maxdE = outEnergy[j-1][Ne - 1] - outEnergy[j-1][Ne - 2];
		dt = dx / advectionV;
		int Nstep = 10 * floor(maxLossRate * dt / maxdE) + 3;
		printf("j = %d\n", j);
		printLog("j = %d\n", j);
		printf("Nstep = %d\n", Nstep);
		printLog("Nstep = %d\n", Nstep);
		dx = x[j] - x[j - 1];
		evaluateDistributionAdvectionWithLossesBetweenTwoPoints(mass, outEnergy[j - 1], outDistribution[j - 1], outEnergy[j], outDistribution[j], Ne, dx, advectionV, B[j], Uph1, Eph1, Uph2, Eph2, Nstep);
	}
}

void MassiveParticleDistributionFactory::evaluateDistributionAdvectionWithLossesBetweenTwoPoints(const double& mass, const double* inEnergy, const double* inDistribution, double* outEnergy, double* outDistribution, const int Ne, const double& deltaX, const double& advectionV, const double& B, const double& Uph1, const double& Eph1, const double& Uph2, const double& Eph2, const int Nstep)
{
	double dx = deltaX / Nstep;
	double dt = dx / advectionV;

	double mc2 = mass * speed_of_light * speed_of_light;

	double r2 = sqr(electron_charge * electron_charge / mc2);
	double sigmaT = 8 * pi * r2 / 3.0;

	double* tempEnergy = new double[Ne];
	double* tempDistribution = new double[Ne];
	for (int i = 0; i < Ne; ++i) {
		tempEnergy[i] = inEnergy[i];
		tempDistribution[i] = inDistribution[i];
	}

	double* dEold = new double[Ne];
	double* dE = new double[Ne];
	int i;
	for (int j = 0; j < Nstep; ++j) {
#pragma omp parallel for private(i) shared(tempEnergy, tempDistribution, dEold, outEnergy, outDistribution, dE, dt, sigmaT, B, Uph1, Eph1, Uph2, Eph2)
		for (i = 0; i < Ne; ++i) {
			if (i == 0) {
				dEold[i] = (tempEnergy[1] - tempEnergy[0]) / 2;
				outEnergy[i] = tempEnergy[i];
			}
			else if (i < Ne - 1) {
				dEold[i] = (tempEnergy[i + 1] - tempEnergy[i - 1]) / 2;

				double E = tempEnergy[i];
				double gamma = E / mc2;
				double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
				double lossRate = ((4.0 / 3.0) * sigmaT * beta * beta * E * E / (mass * mass * cube(speed_of_light))) * (B * B / (8 * pi) + Uph1 / pow(1 + 4 * E * Eph1 / sqr(mc2), 3.0 / 2.0) + Uph2 / pow(1 + 4 * E * Eph2 / sqr(mc2), 3.0 / 2.0));

				outEnergy[i] = E - lossRate * dt;
			}
			else {
				dEold[i] = (tempEnergy[i] - tempEnergy[i - 1]) / 2;

				double E = tempEnergy[i];
				double gamma = E / mc2;
				double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
				double lossRate = ((4.0 / 3.0) * sigmaT * beta * beta * E * E / (mass * mass * cube(speed_of_light))) * (B * B / (8 * pi) + Uph1 / pow(1 + 4 * E * Eph1 / sqr(mc2), 3.0 / 2.0) + Uph2 / pow(1 + 4 * E * Eph2 / sqr(mc2), 3.0 / 2.0));

				outEnergy[i] = E - lossRate * dt;
			}
		}

		for (i = 0; i < Ne; ++i) {
			if (i == 0) {
				dE[i] = (outEnergy[1] - outEnergy[0]) / 2;
			}
			else if (i < Ne - 1) {
				dE[i] = (outEnergy[i + 1] - outEnergy[i - 1]) / 2;
			}
			else {
				dE[i] = (outEnergy[i] - outEnergy[i - 1]) / 2;
			}

			if (dE[i] < 0) {
				printf("dE < 0 in evaluate distribution advection with losses\n");
				printLog("dE < 0 in evaluate distribution advection with losses\n");
				exit(0);
			}

			outDistribution[i] = tempDistribution[i] * dEold[i] / dE[i];
		}

		int minIndex = 0;
		for (int i = 0; i < Ne; ++i) {
			if (outDistribution[i] > 0) {
				minIndex = i;
				break;
			}
		}
		int maxIndex = Ne - 1;
		for (int i = Ne - 1; i >= 0; --i) {
			if (outDistribution[i] > 0) {
				maxIndex = i;
				break;
			}
		}

		tempEnergy[0] = outEnergy[0];
		tempEnergy[Ne - 1] = outEnergy[Ne - 1];
		if (minIndex > 0) {
			tempEnergy[minIndex] = outEnergy[minIndex];
			double factor1 = pow(tempEnergy[minIndex] / tempEnergy[0], 1.0 / minIndex);
			for (int i = 1; i < minIndex; ++i) {
				//tempEnergy[i] = tempEnergy[0] + i * (tempEnergy[minIndex] - tempEnergy[0]) / minIndex;
				tempEnergy[i] = tempEnergy[i - 1] * factor1;
			}
		}
		if (maxIndex < Ne - 1) {
			tempEnergy[maxIndex] = outEnergy[maxIndex];
			double factor1 = pow(tempEnergy[Ne - 1] / tempEnergy[maxIndex], 1.0 / (Ne - 1 - maxIndex));
			for (int i = maxIndex + 1; i < Ne; ++i) {
				//tempEnergy[i] = tempEnergy[maxIndex] + (i - maxIndex) * (tempEnergy[Ne - 1] - tempEnergy[maxIndex]) / (Ne - 1 - maxIndex);
				tempEnergy[i] = tempEnergy[i - 1] * factor1;
			}
		}

		double factor = pow(tempEnergy[maxIndex] / tempEnergy[minIndex], 1.0 / (maxIndex - minIndex));
		for (int i = minIndex + 1; i < maxIndex; ++i) {
			tempEnergy[i] = tempEnergy[i - 1] * factor;
		}

		redistributeArrayWithFixedEnergy(outEnergy, outDistribution, tempEnergy, tempDistribution, Ne);
	}

	for (int i = 0; i < Ne; ++i) {
		outEnergy[i] = tempEnergy[i];
		outDistribution[i] = tempDistribution[i];
	}

	delete[] dEold;
	delete[] dE;
	delete[] tempEnergy;
	delete[] tempDistribution;
}

void MassiveParticleDistributionFactory::redistributeArray(const double* inEnergy, const double* inDistribution, double* outEnergy, double* outDistribution, int Ne)
{
	double minE = inEnergy[0];
	double maxE = inEnergy[Ne - 1];
	double factor = pow(maxE / minE, 1.0 / (Ne - 1));

	outEnergy[0] = minE;
	for (int i = 1; i < Ne; ++i) {
		outEnergy[i] = outEnergy[i - 1] * factor;
	}
	outEnergy[Ne - 1] = maxE;

	outDistribution[0] = inDistribution[0];
	for (int i = 1; i < Ne - 1; ++i) {
		outDistribution[i] = MassiveParticleDistributionFactory::getDistribution(outEnergy[i], inEnergy, inDistribution, Ne);
	}
	outDistribution[Ne - 1] = inDistribution[Ne - 1];
}

void MassiveParticleDistributionFactory::redistributeArrayWithFixedEnergy(const double* inEnergy, const double* inDistribution, const double* outEnergy, double* outDistribution, int Ne)
{
	for (int i = 0; i < Ne; ++i) {
		outDistribution[i] = MassiveParticleDistributionFactory::getDistribution(outEnergy[i], inEnergy, inDistribution, Ne);
	}
}

double MassiveParticleDistributionFactory::getDistribution(const double& E, const double* energy, const double* distribution, int Ne)
{
	if (E < energy[0]) {
		return 0;
	}

	if (E > energy[Ne - 1]) {
		return 0;
	}

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
	int maxI = Ne - 1;
	while (maxI - minI > 1) {
		int tempI = (minI + maxI) / 2;
		if (energy[tempI] == E) {
			return distribution[tempI];
		} else if (energy[tempI] > E) {
			maxI = tempI;
		}
		else {
			minI = tempI;
		}
	}
	currentI = minI;
	nextI = maxI;


	double result;
	if (distribution[currentI] <= 0 || distribution[nextI] <= 0) {
		result = (distribution[currentI] * (energy[nextI] - E) + distribution[nextI] * (E - energy[currentI])) / (energy[nextI] - energy[currentI]);
		result = 0;
	}
	else {
		result = distribution[currentI] * exp(log(distribution[nextI] / distribution[currentI]) * ((log(E / energy[currentI])) / (log(energy[nextI] / energy[currentI]))));
		//result = (my_distribution[currentI] * (my_energy[nextI] - energy) + my_distribution[nextI] * (energy - my_energy[currentI])) / (my_energy[nextI] - my_energy[currentI]);
	}
	if (result != result) {
		printf("result = NaN in get distribution\n");
		printLog("result = NaN in get distribution\n");
		exit(0);
	}
	return result;
}

MassiveParticleMonoenergeticDistribution::MassiveParticleMonoenergeticDistribution(const double& mass, const double& Energy, const double& halfWidth)
{
	my_mass = mass;
	my_E0 = Energy;
	my_dE = halfWidth;
	if (halfWidth > my_E0) {
		printf("dE > E0 in monoenergetic distribution\n");
		printLog("dE > E0 in monoenergetic distribution\n");
		exit(0);
	}
}

double MassiveParticleMonoenergeticDistribution::distributionNormalized(const double& energy)
{
	if (energy > my_E0 + my_dE || energy < my_E0 - my_dE) {
		return 0;
	}
	return 1.0/(4*pi*2*my_dE);
}

double MassiveParticleMonoenergeticDistribution::getMeanEnergy()
{
	//only for narrow beam
	return my_E0;
}

double MassiveParticleMonoenergeticDistribution::minEnergy()
{
	return my_E0 - my_dE;
}

double MassiveParticleMonoenergeticDistribution::maxEnergy()
{
	return my_E0 + my_dE;
}

MassiveParticleMonoenergeticDirectedDistribution::MassiveParticleMonoenergeticDirectedDistribution(const double& mass, const double& Energy, const double& halfWidth, const double& theta0, const double& phi0, const double& deltaTheta)
{
	my_mass = mass;

	my_E0 = Energy;
	my_dE = halfWidth;

	my_theta0 = theta0;
	my_phi0 = phi0;
	my_deltaTheta = deltaTheta;
}

double MassiveParticleMonoenergeticDirectedDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	double sinTheta = sin(acos(mu));

	double cosDelta = mu * cos(my_theta0) + sinTheta * sin(my_theta0) * cos(phi - my_phi0);

	if (cosDelta < cos(my_deltaTheta)) {
		return 0;
	}

	if (energy > my_E0 + my_dE || energy < my_E0 - my_dE) {
		return 0;
	}

	double sphericalFraction = 0.5 * (1.0 - cos(my_deltaTheta));

	return 1.0 / (4 * pi * 2 * my_dE * sphericalFraction);
}

double MassiveParticleMonoenergeticDirectedDistribution::getMeanEnergy()
{
	return my_E0;
}

double MassiveParticleMonoenergeticDirectedDistribution::minEnergy()
{
	return my_E0 - my_dE;
}

double MassiveParticleMonoenergeticDirectedDistribution::maxEnergy()
{
	return my_E0 + my_dE;
}

//todo what to do with concentration?
MassiveParticleMovingDistribution::MassiveParticleMovingDistribution(MassiveParticleDistribution* distribution, const double& velocity)
{
	my_distribution = distribution;
	my_mass = distribution->getMass();
	if (velocity > speed_of_light || velocity < -speed_of_light) {
		printf("velocity > c in MassiveParticleMovingDistribution\n");
		printLog("velocity > c in MassiveParticleMovingDistribution\n");
		exit(0);
	}
	my_velocity = velocity;
	my_gamma = 1.0 / sqrt(1.0 - my_velocity * my_velocity / speed_of_light2);
}

double MassiveParticleMovingDistribution::getMeanEnergy()
{
	double tempEmin = my_distribution->getMeanEnergy();
	double Emean = 0;
	int Nphi = 20;
	int Nmu = 20;
	int Ne = 100;
	double dphi = 2 * pi / Nphi;
	double dmu = 2 / Nmu;
	double Emin = my_mass * speed_of_light2;
	double Emax = 10 * tempEmin * my_gamma;
	double factor = pow(Emax / Emin, 1.0 / (Ne - 1));
	for (int imu = 0; imu < Nmu; ++imu) {
		double mu = -1 + (imu + 0.5) * dmu;
		for (int k = 0; k < Nphi; ++k) {
			double currentEnergy = Emin*1.00000001;
			for (int i = 0; i < Ne; ++i) {
				Emean += distributionNormalized(currentEnergy, mu, (k+0.5)*dphi) * currentEnergy*(factor-1.0) * dmu * dphi * currentEnergy;
			}
		}
	}
	return Emean;
}

double MassiveParticleMovingDistribution::minEnergy()
{
	return my_mass*speed_of_light2;
}

double MassiveParticleMovingDistribution::maxEnergy()
{
	double emax = my_distribution->maxEnergy();
	if (emax < 0) {
		return -1.0;
	}

	double pmax = sqrt(emax * emax - my_mass * my_mass * speed_of_light2) / speed_of_light;
	double tempemax = my_gamma * emax + my_velocity * my_gamma * pmax;

	return tempemax;
}

double MassiveParticleMovingDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	if (energy < my_mass * speed_of_light2) {
		printf("energy < m*c2 in MassiveOarticleMovingDistribution\n");
		printLog("energy < m*c2 in MassiveOarticleMovingDistribution\n");
		exit(0);
	}

	double p = sqrt(energy * energy - my_mass * my_mass * speed_of_light4) / speed_of_light;
	double pparallel = p * mu;
	double pnorm = p * sqrt(1.0 - mu * mu);

	double beta = my_velocity / speed_of_light;

	double energy1 = my_gamma * energy - my_gamma * beta * mu * p * speed_of_light;
	double p1 = sqrt(energy1 * energy1 - my_mass * my_mass * speed_of_light4) / speed_of_light;
	double mu1 = (-my_gamma * beta * energy + my_gamma * mu * p * speed_of_light) / (p1 * speed_of_light);

	double jacobian = (p/p1)*(beta*beta + 1)/(my_gamma*my_gamma);//todo

	return my_distribution->distributionNormalized(energy1, mu1, phi)*jacobian;
}
