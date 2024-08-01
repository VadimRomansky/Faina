#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "photonDistribution.h"

PhotonPlankDistribution* PhotonPlankDistribution::my_CMBradiation = 0;
PhotonMultiPlankDistribution* PhotonMultiPlankDistribution::my_GalacticField = 0;

double PhotonIsotropicDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	return distributionNormalized(energy);
}

void PhotonIsotropicDistribution::writeDistribution(const char* fileName, int Ne, const double& Emin, const double& Emax) {
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

	my_A = (my_index - 1) / (my_E0 * 4 * pi);
}

double PhotonPowerLawDistribution::distributionNormalized(const double& energy)
{
	if (energy < 0) {
		printf("photon energy < 0\n");
		printLog("photon energy < 0\n");
		exit(0);
	}
	if (energy < my_E0) {
		return 0;
	}
    return my_A/pow(energy/my_E0, my_index);
}

double PhotonPowerLawDistribution::getMeanEnergy() {
	return my_A * my_E0 * my_E0 / (my_index - 2);
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

double PhotonPlankDistribution::distributionNormalized(const double& energy) {
	//for debug
	/*if (energy > 10 * kBoltzman * my_temperature) {
		return 0;
	}*/
	double theta = energy / (kBoltzman * my_temperature);
	if (theta < 1E-12) {
		return my_A * (2 * kBoltzman * my_temperature * energy / cube(hplank * speed_of_light)) / my_concentration;
	}
	return my_A*(2 * energy * energy / cube(hplank * speed_of_light)) / (exp(theta) - 1.0)/my_concentration;
}

double PhotonPlankDistribution::getMeanEnergy()
{
	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;
	double intPlank3 = pi*pi*pi*pi/15;
	return kBoltzman*my_temperature*intPlank3/intPlank2;
}

double PhotonPlankDistribution::getTemperature() {
	return my_temperature;
}

PhotonPlankDistribution* PhotonPlankDistribution::getCMBradiation()
{
	if (!my_CMBradiation) {
		my_CMBradiation = new PhotonPlankDistribution(2.725, 1.0);
	}
	return my_CMBradiation;
}

PhotonMultiPlankDistribution::PhotonMultiPlankDistribution(int Nplank, const double* const temperatures, const double* const amplitudes)
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

double PhotonMultiPlankDistribution::distributionNormalized(const double& energy)
{
	double result = 0;
	for (int i = 0; i < my_Nplank; ++i) {
		double theta = energy / (kBoltzman * my_temperatures[i]);
		result += my_A[i] * (2 * energy * energy / cube(hplank * speed_of_light)) / (exp(theta) - 1.0);
	}

	return result/my_concentration;
}

double PhotonMultiPlankDistribution::getMeanEnergy()
{
	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;
	double intPlank3 = pi * pi * pi * pi / 15;
	double result = 0;
	for (int i = 0; i < my_Nplank; ++i) {
		result += my_concentrations[i] * my_temperatures[i] * kBoltzman * intPlank3 / intPlank2;
	}
	return result/getConcentration();
}

PhotonMultiPlankDistribution* PhotonMultiPlankDistribution::getGalacticField()
{
	if (!my_GalacticField) {
		double temperatures[5] = { 2.725, 20, 3000, 4000, 7500 };
		double amplitudes[5] = { 1.0, 4E-4, 4E-13, 1.65E-13, 1E-14 };
		my_GalacticField = new PhotonMultiPlankDistribution(5, temperatures, amplitudes);
	}

	return my_GalacticField;
}

CompoundPhotonDistribution::CompoundPhotonDistribution(int N, PhotonDistribution** distributions)
{
	my_Ndistr = N;

	my_distributions = new PhotonDistribution * [my_Ndistr];
	my_concentration = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		my_distributions[i] = distributions[i];
		my_concentration += my_distributions[i]->getConcentration();
	}
}

CompoundPhotonDistribution::CompoundPhotonDistribution(PhotonDistribution* dist1, PhotonDistribution* dist2)
{
	my_Ndistr = 2;
	my_distributions = new PhotonDistribution * [my_Ndistr];
	my_concentration = dist1->getConcentration() + dist2->getConcentration();
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
}

CompoundPhotonDistribution::CompoundPhotonDistribution(PhotonDistribution* dist1, PhotonDistribution* dist2, PhotonDistribution* dist3)
{
	my_Ndistr = 3;
	my_distributions = new PhotonDistribution * [my_Ndistr];
	my_concentration = dist1->getConcentration() + dist2->getConcentration() + dist3->getConcentration();
	my_distributions[0] = dist1;
	my_distributions[1] = dist2;
	my_distributions[2] = dist3;
}

CompoundPhotonDistribution::~CompoundPhotonDistribution()
{
	delete[] my_distributions;
}

double CompoundPhotonDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	double result = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		result += my_distributions[i]->distribution(energy, mu, phi);
	}
	return result/my_concentration;
}

double CompoundPhotonDistribution::getMeanEnergy()
{
	double result = 0;
	for (int i = 0; i < my_Ndistr; ++i) {
		result += my_distributions[i]->getConcentration() * my_distributions[i]->getMeanEnergy();
	}
	return result / my_concentration;
}

PhotonPlankDirectedDistribution::PhotonPlankDirectedDistribution(const double& temperature, const double& amplitude, const double& theta0, const double& phi0, const double& deltaTheta)
{
	my_temperature = temperature;

	my_theta0 = theta0;
	my_phi0 = phi0;
	my_deltaTheta = deltaTheta;

	double sphericalFraction = 0.5 * (1.0 - cos(my_deltaTheta));
	my_A = amplitude/sphericalFraction;

	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;

	my_concentration = my_A * intPlank2 * (8 * pi / cube(hplank * speed_of_light)) * cube(kBoltzman * my_temperature)*sphericalFraction;
}

double PhotonPlankDirectedDistribution::distributionNormalized(const double& energy, const double& mu, const double& phi)
{
	//for debug
	/*if (energy > 10 * kBoltzman * my_temperature) {
		return 0;
	}*/

	double sinTheta = sin(acos(mu));

	double cosDelta = mu * cos(my_theta0) + sinTheta * sin(my_theta0) * cos(phi - my_phi0);

	if (cosDelta < cos(my_deltaTheta)) {
		return 0;
	}
	//debug
	//double factor1 = (1.0/ (fabs(cos(my_deltaTheta) - 1)*sqrt(pi)))*exp(-sqr((cosDelta - 1) / (cos(my_deltaTheta) - 1)));
	double factor1 = 1;

	double theta = energy / (kBoltzman * my_temperature);
	if (theta < 1E-12) {
		return factor1*my_A * (2 * kBoltzman * my_temperature * energy / cube(hplank * speed_of_light)) / my_concentration;
	}
	return factor1*my_A * (2 * energy * energy / cube(hplank * speed_of_light)) / (exp(theta) - 1.0) / my_concentration;
}

double PhotonPlankDirectedDistribution::getMeanEnergy()
{
	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;
	double intPlank3 = pi*pi*pi*pi/15;
	return kBoltzman*my_temperature*intPlank3/intPlank2;
}

double PhotonPlankDirectedDistribution::getTemperature()
{
	return my_temperature;
}

PhotonPlankDirectedDistribution::~PhotonPlankDirectedDistribution() {

}

PhotonMonoenergeticDistribution::PhotonMonoenergeticDistribution(const double& Energy, const double& halfWidth, const double& concentration)
{
	my_E0 = Energy;
	my_dE = halfWidth;
	if (halfWidth > my_E0) {
		printf("dE > E0 in monoenergetic distribution\n");
		printLog("dE > E0 in monoenergetic distribution\n");
		exit(0);
	}
	my_concentration = concentration;
}

double PhotonMonoenergeticDistribution::distributionNormalized(const double& energy)
{
	if (energy > my_E0 + my_dE || energy < my_E0 - my_dE) {
		return 0;
	}
	return 1.0 / (4 * pi * 2 * my_dE);
}

double PhotonMonoenergeticDistribution::getMeanEnergy()
{
	return my_E0;
}

void PhotonMonoenergeticDistribution::resetConcentration(const double& concentration)
{
	my_concentration = concentration;
}
