#include "stdio.h"
#include "math.h"
#include <memory>

#include "constants.h"
#include "util.h"
#include "massiveParticleDistribution.h"

#include "radiationSource.h"

RadiationSource::RadiationSource(int Nrho, int Nz, int Nphi, double distance, double redShift) {
	my_Nrho = Nrho;
	my_Nz = Nz;
	my_Nphi = Nphi;
	my_distance = distance;
	my_redShift = redShift;
}

double RadiationSource::evaluateAverageVelocity()
{
	double v = 0;
	double m = 0;
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for(int iphi = 0; iphi < my_Nphi; ++iphi){
				double tempv, theta, phi;
				getVelocity(irho, iz, iphi, tempv, theta, phi);
				v = v + tempv * getConcentration(irho, iz, iphi) * getVolume(irho, iz, iphi);
				m = m + getConcentration(irho, iz, iphi) * getVolume(irho, iz, iphi);
			}
		}
	}
	return v/m;
}

int RadiationSource::getNrho() {
	return my_Nrho;
}
int RadiationSource::getNz() {
	return my_Nz;
}
int RadiationSource::getNphi() {
	return my_Nphi;
}
double RadiationSource::getDistance() {
	return my_distance;
}

double RadiationSource::getRedShift()
{
	return my_redShift;
}

double RadiationSource::getVolume(int irho, int iz, int iphi)
{
	return getArea(irho, iz, iphi)*getLength(irho, iz, iphi);
}

DiskSource::DiskSource(int Nrho, int Nz, int Nphi, const double& rho, const double& z, const double& distance, const double& redShift) : RadiationSource(Nrho, Nz, Nphi, distance, redShift)
{
	my_rho = rho;
	my_z = z;
}

double DiskSource::getArea(int irho, int iz, int iphi) {
	double rho0 = irho * getMaxRho() / my_Nrho;
	double rho1 = (irho + 1) * getMaxRho() / my_Nrho;
	return pi * (rho1 * rho1 - rho0 * rho0) / my_Nphi;
}

double DiskSource::getCrossSectionArea(int irho, int iphi)
{
	double rho1 = (irho + 1) * my_rho / my_Nrho;
	double rho0 = irho*my_rho / my_Nrho;
	return (pi/my_Nphi)*(rho1*rho1 - rho0*rho0);
}

double DiskSource::getTotalVolume()
{
	return pi * my_z * my_rho * my_rho;
}
double DiskSource::getMaxRho() {
	return my_rho;
}
double DiskSource::getMinRho() {
	return 0;
}
double DiskSource::getMinZ() {
	return 0;
}
double DiskSource::getMaxZ() {
	return my_z;
}

SimpleFlatSource::SimpleFlatSource(MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& concentration, const double& phi, const double& rho, const double& z, const double& distance, const double& velocity, const double& redShift) : DiskSource(1,1,1, rho, z,distance, redShift) {
	my_distribution = electronDistribution;
	my_B = B;
	my_theta = theta;
	my_phi = phi;
	my_concentration = concentration;
	my_velocity = velocity;
}

double SimpleFlatSource::getRho(int irho) {
	return 0;
}
double SimpleFlatSource::getZ(int iz) {
	return 0;
}
double SimpleFlatSource::getPhi(int iphi) {
	return 0;
}

int SimpleFlatSource::getRhoIndex(const double& rho) {
	if (rho < 0) {
		printf("rho < 0 in get rhoIndex\n");
		printLog("rho < 0 in get rhoIndex\n");
		exit(0);
	}
	if (rho > my_rho) {
		printf("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		printLog("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		exit(0);
	}

	return 0;
}

bool SimpleFlatSource::isSource(int irho, int iphi) {
	return true;
}

double SimpleFlatSource::getB(int irho, int iz, int iphi)
{
	return my_B;
}
double SimpleFlatSource::getMaxB()
{
	return my_B;
}
double SimpleFlatSource::getMaxOuterB()
{
	return my_B;
}
double SimpleFlatSource::getAverageSigma()
{
	return my_B*my_B/(4*pi*my_concentration*massProton*speed_of_light2);
}
double SimpleFlatSource::getAverageConcentration()
{
	return my_concentration;
}
double SimpleFlatSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration;
}
void SimpleFlatSource::getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi)
{
	velocity = my_velocity;
	if (velocity > speed_of_light) {
		printf("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		printLog("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		exit(0);
	}
	theta = 0;
	phi = 0;
}
double SimpleFlatSource::getSinTheta(int irho, int iz, int iphi)
{
	return sin(my_theta);
}
double SimpleFlatSource::getBTheta(int irho, int iz, int iphi)
{
	return my_theta;
}
double SimpleFlatSource::getBPhi(int irho, int iz, int iphi)
{
	return my_phi;
}
/*void SimpleFlatSource::resetConcentration(const double& concentration)
{
	my_distribution->resetConcentration(concentration);
}*/
void SimpleFlatSource::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
	* R = parameters[0]
	* sigma = parameters[1]
	* n = parameters[2]
	* z = R*parameters[3]
	* v = parameters[4]
	*/
	double sigma = parameters[1] * normalizationUnits[1];
	my_rho = parameters[0] * normalizationUnits[0];
	my_concentration = parameters[2] * normalizationUnits[2];
	my_z = my_rho * parameters[3] * normalizationUnits[3];
	my_B = sqrt(sigma * 4 * pi * massProton * my_concentration * speed_of_light2);
	my_velocity = parameters[4] * normalizationUnits[4];
}
double SimpleFlatSource::getLength(int irho, int iz, int iphi) {
	return my_z;
}
MassiveParticleDistribution* SimpleFlatSource::getParticleDistribution(int irho, int iz, int iphi) {
	my_distribution->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distribution;
}

SimpleFlatSource2::SimpleFlatSource2(int Ndistributions, double* velocities, MassiveParticleIsotropicDistribution** electronDistributions, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& velocity, const double& redShift) : SimpleFlatSource(NULL, B, theta, phi, concentration, rho, z, distance, velocity, redShift){
	my_Ndistributions = Ndistributions;
	my_maxThreads = omp_get_max_threads();
	my_velocities = new double[my_Ndistributions];
	my_outputDistributions = new MassiveParticleDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_outputDistributions[i] = NULL;
	}
	my_distributions = new MassiveParticleDistribution * [my_Ndistributions];
	for (int i = 0; i < my_Ndistributions; ++i) {
		my_velocities[i] = velocities[i];
		my_distributions[i] = electronDistributions[i];
	}
}

MassiveParticleDistribution* SimpleFlatSource2::getParticleDistribution(int irho, int iz, int iphi)
{
	int num_threads = omp_get_num_threads();
	//todo not thread safe!!!
	if (my_outputDistributions[num_threads] != NULL) {
		return my_outputDistributions[num_threads];
	}

	/*for (int i = 0; i < my_Ndistributions; ++i) {
		my_distributions[i]->resetConcentration(getConcentration(irho, iz, iphi));
	}*/
	if (my_velocity <= my_velocities[0]) {
		MassiveParticleIsotropicDistribution* iso1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[0]);
		if (iso1 != NULL) {
			my_outputDistributions[num_threads] = new CompoundWeightedMassiveParticleIsotropicDistribution(iso1, 0.5, iso1, 0.5);
		}
		else {
			my_outputDistributions[num_threads] = new CompoundWeightedMassiveParticleDistribution(my_distributions[0], 0.5, my_distributions[0], 0.5);
		}
	}
	else if (my_velocity >= my_velocities[my_Ndistributions - 1]) {
		MassiveParticleIsotropicDistribution* iso1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[my_Ndistributions - 1]);
		if (iso1 != NULL) {
			my_outputDistributions[num_threads] = new CompoundWeightedMassiveParticleIsotropicDistribution(iso1, 0.5, iso1, 0.5);
		}
		else {
			my_outputDistributions[num_threads] = new CompoundWeightedMassiveParticleDistribution(my_distributions[my_Ndistributions - 1], 0.5, my_distributions[my_Ndistributions - 1], 0.5);
		}
	}
	else {
		for (int i = 1; i < my_Ndistributions; ++i) {
			if (my_velocities[i] > my_velocity) {
				double left = (my_velocities[i] - my_velocity) / (my_velocities[i] - my_velocities[i - 1]);
				double right = (my_velocity - my_velocities[i - 1]) / (my_velocities[i] - my_velocities[i - 1]);
				MassiveParticleIsotropicDistribution* iso1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[i - 1]);
				MassiveParticleIsotropicDistribution* iso2 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[i]);
				if ((iso1 != NULL) && (iso2 != NULL)) {
					my_outputDistributions[num_threads] = new CompoundWeightedMassiveParticleIsotropicDistribution(iso1, left, iso2, right);
				}
				else {
					my_outputDistributions[num_threads] = new CompoundWeightedMassiveParticleDistribution(my_distributions[i - 1], left, my_distributions[i], right);
				}
				return my_outputDistributions[num_threads];
			}
		}
	}
	return my_outputDistributions[num_threads];
}

void SimpleFlatSource2::resetParameters(const double* parameters, const double* normalizationUnits)
{
	if (my_distribution != NULL) {
		delete my_distribution;
		my_distribution = NULL;
	}

	for (int i = 0; i < my_maxThreads; ++i) {
		if (my_outputDistributions[i] != NULL) {
			delete my_outputDistributions[i];
			my_outputDistributions[i] = NULL;
		}
	}

	SimpleFlatSource::resetParameters(parameters, normalizationUnits);
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, const double& velocity, const double& redShift) : DiskSource(Nrho, Nz, Nphi, rho, z, distance, redShift) {
	my_distribution = electronDistribution;

	my_velocity = velocity;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		my_v[irho] = new double* [my_Nz];
		my_vtheta[irho] = new double* [my_Nz];
		my_vphi[irho] = new double* [my_Nphi];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			my_v[irho][iz] = new double[my_Nphi];
			my_vtheta[irho][iz] = new double[my_Nphi];
			my_vphi[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_theta[irho][iz][iphi] = theta[irho][iz][iphi];
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				my_v[irho][iz][iphi] = my_velocity;
				my_vtheta[irho][iz][iphi] = 0;
				my_vphi[irho][iz][iphi] = 0;
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& velocity, const double& redShift) : DiskSource(Nrho, Nz, Nphi, rho, z, distance, redShift)
{
	my_distribution = electronDistribution;
	my_velocity = velocity;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		my_v[irho] = new double* [my_Nz];
		my_vtheta[irho] = new double* [my_Nz];
		my_vphi[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			my_v[irho][iz] = new double[my_Nphi];
			my_vtheta[irho][iz] = new double[my_Nphi];
			my_vphi[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_theta[irho][iz][iphi] = theta;
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
				my_v[irho][iz][iphi] = my_velocity;
				my_vtheta[irho][iz][iphi] = 0;
				my_vphi[irho][iz][iphi] = 0;
			}
		}
	}

	my_isSource = new bool*[my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}
TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : DiskSource(Nrho, Nz, Nphi, rho, z, distance, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		my_v[irho] = new double* [my_Nz];
		my_vtheta[irho] = new double* [my_Nz];
		my_vphi[irho] = new double* [my_Nphi];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			my_v[irho][iz] = new double[my_Nphi];
			my_vtheta[irho][iz] = new double[my_Nphi];
			my_vphi[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_theta[irho][iz][iphi] = theta[irho][iz][iphi];
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				my_v[irho][iz][iphi] = velocity[irho][iz][iphi];
				my_vtheta[irho][iz][iphi] = vtheta[irho][iz][iphi];
				my_vphi[irho][iz][iphi] = vphi[irho][iz][iphi];
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}

	my_velocity = evaluateAverageVelocity();
}
TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : DiskSource(Nrho, Nz, Nphi, rho, z, distance, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		my_v[irho] = new double* [my_Nz];
		my_vtheta[irho] = new double* [my_Nz];
		my_vphi[irho] = new double* [my_Nphi];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			my_v[irho][iz] = new double[my_Nphi];
			my_vtheta[irho][iz] = new double[my_Nphi];
			my_vphi[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_theta[irho][iz][iphi] = theta;
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
				my_v[irho][iz][iphi] = velocity[irho][iz][iphi];
				my_vtheta[irho][iz][iphi] = vtheta[irho][iz][iphi];
				my_vphi[irho][iz][iphi] = vphi[irho][iz][iphi];
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}

	my_velocity = evaluateAverageVelocity();
}
TabulatedDiskSource::~TabulatedDiskSource()
{
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_B[irho][iz];
			delete[] my_theta[irho][iz];
			delete[] my_phi[irho][iz];
			delete[] my_concentration[irho][iz];
			delete[] my_v[irho][iz];
			delete[] my_vtheta[irho][iz];
			delete[] my_vphi[irho][iz];
		}
		delete[] my_B[irho];
		delete[] my_theta[irho];
		delete[] my_phi[irho];
		delete[] my_concentration[irho];
		delete[] my_isSource[irho];
		delete[] my_v[irho];
		delete[] my_vtheta[irho];
		delete[] my_vphi[irho];
	}
	delete[] my_B;
	delete[] my_theta;
	delete[] my_phi;
	delete[] my_concentration;
	delete[] my_isSource;
	delete[] my_v;
	delete[] my_vtheta;
	delete[] my_vphi;
}

double TabulatedDiskSource::getRho(int irho) {
	return (irho + 0.5) * my_rho / my_Nrho;
}
double TabulatedDiskSource::getZ(int iz) {
	return (iz + 0.5) * my_z / my_Nz;
}
double TabulatedDiskSource::getPhi(int iphi) {
	return (iphi + 0.5) * 2 * pi / my_Nphi;
}

int TabulatedDiskSource::getRhoIndex(const double& rho) {
	if (rho < 0) {
		printf("rho < 0 in get rhoIndex\n");
		printLog("rho < 0 in get rhoIndex\n");
		exit(0);
	}
	if (rho > my_rho) {
		printf("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		printLog("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		exit(0);
	}
	if (rho == my_rho) {
		return my_Nrho - 1;
	}
	double drho = my_rho / my_Nrho;
	return floor(rho / drho);
}

void TabulatedDiskSource::setMask(bool** mask) {
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			my_isSource[irho][iphi] = mask[irho][iphi];
		}
	}
}

bool TabulatedDiskSource::isSource(int irho, int iphi) {
	return my_isSource[irho][iphi];
}

double TabulatedDiskSource::getAverageBsquared()
{
	double magneticEnergy = 0;
	double volume = 0;


	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				magneticEnergy += my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * getVolume(irho, iz, iphi);
				volume += getVolume(irho, iz, iphi);
			}
		}
	}

	return sqrt(magneticEnergy / volume);
}
double TabulatedDiskSource::getB(int irho, int iz, int iphi)
{
	return my_B[irho][iz][iphi];
}
double TabulatedDiskSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration[irho][iz][iphi];
}
void TabulatedDiskSource::getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi)
{
	velocity = my_v[irho][iz][iphi];
	if (velocity > speed_of_light) {
		printf("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		printLog("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		exit(0);
	}
	theta = my_vtheta[irho][iz][iphi];
	phi = my_vphi[irho][iz][iphi];
}
double TabulatedDiskSource::getSinTheta(int irho, int iz, int iphi)
{
	return sin(my_theta[irho][iz][iphi]);
}
double TabulatedDiskSource::getBTheta(int irho, int iz, int iphi)
{
	return my_theta[irho][iz][iphi];
}
double TabulatedDiskSource::getBPhi(int irho, int iz, int iphi)
{
	return my_phi[irho][iz][iphi];
}
/*void TabulatedDiskSource::resetConcentration(const double& concentration)
{
	my_distribution->resetConcentration(concentration);
}*/
void TabulatedDiskSource::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* sigma[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* z = R*parameters[3]
*/
	my_rho = parameters[0] * normalizationUnits[0];
	double sigma = parameters[1] * normalizationUnits[1];
	//double B0 = my_B[my_Nrho - 1][0][0];
	//double n0 = my_concentration[my_Nrho - 1][0][0];
	double n0 = getAverageConcentration();
	//double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	double sigma0 = getAverageSigma();

	my_z = my_rho * parameters[3] * normalizationUnits[3];
	double old_velocity = my_velocity;
	my_velocity = parameters[4] * normalizationUnits[4];

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double sigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi]*speed_of_light2);
				sigma *= parameters[1] * normalizationUnits[1]/sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2]/n0;
				my_B[irho][iz][iphi] = sqrt(sigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				if(old_velocity > 0){
					my_v[irho][iz][iphi] *= my_velocity / old_velocity;
				}
				else {
					my_v[irho][iz][iphi] = my_velocity;
				}
			}
		}
	}
	
}
double TabulatedDiskSource::getLength(int irho, int iz, int iphi) {
	return my_z/my_Nz;
}
double TabulatedDiskSource::getMaxB()
{
	double Bmax = 0;
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (my_B[irho][iz][iphi] > Bmax) {
					Bmax = my_B[irho][iz][iphi];
				}
			}
		}
	}
	return Bmax;
}
double TabulatedDiskSource::getMaxOuterB()
{
	double Bmax = 0;
	int irho = my_Nrho - 1;
	int iz = 0;
	for (iz = 0; iz < my_Nz; ++iz) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}
		}
	}

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			iz = 0;
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}

			iz = my_Nz - 1;
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}
		}
	}

	return Bmax;
}
double TabulatedDiskSource::getAverageSigma()
{
	double magneticEnergy = 0;
	double restEnergy = 0;


	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				magneticEnergy += my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * getVolume(irho, iz, iphi) / (4 * pi);
				restEnergy += my_concentration[irho][iz][iphi] * massProton * speed_of_light2 * getVolume(irho, iz, iphi);
			}
		}
	}

	return magneticEnergy / restEnergy;
}
double TabulatedDiskSource::getAverageConcentration()
{
	double concentration = 0;
	double volume = 0;

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				concentration += my_concentration[irho][iz][iphi] * getVolume(irho, iz, iphi);
				volume += getVolume(irho, iz, iphi);
			}
		}
	}

	double result = concentration / volume;
	if (result != result) {
		printf("averageConcentration = NaN\n");
		printLog("averageConcentration = NaN\n");
		exit(0);
	}
	return result;
}
MassiveParticleDistribution* TabulatedDiskSource::getParticleDistribution(int irho, int iz, int iphi) {
	my_distribution->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distribution;
}

double SphericalLayerSource::evaluateLength(int irho, int iz, int iphi) {
	double dr = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;
	double length = 0;
	double r = (irho + 0.5) * dr;
	double z1 = -sqrt(my_rho * my_rho - r * r);
	double z2 = 0;
	if (my_rhoin > r) {
		z2 = -sqrt(my_rhoin * my_rhoin - r * r);
	}
	double z3 = -z2;
	double z4 = -z1;
	double zmin = -my_rho + iz * dz;
	double zmax = zmin + dz;
	if (zmin >= 0) {
		if (z3 > zmax) {
			length = 0;
		}
		else
			if (z4 < zmin) {
				length = 0;
			}
			else {
				double lowz = max(zmin, z3);
				double topz = min(zmax, z4);
				length = topz - lowz;
			}
	}
	else {
		if (z1 > zmax) {
			length = 0;
		}
		else if (z2 < zmin) {
			length = 0;
		}
		else {
			double lowz = max(zmin, z1);
			double topz = min(zmax, z2);
			length = topz - lowz;
		}
	}
	return length;
}

double SphericalLayerSource::evaluateArea(int irho, int iz, int iphi) {
	double dr = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;
	double area = 0;
	double rho0 = irho * getMaxRho() / my_Nrho;
	double rho1 = (irho + 1) * getMaxRho() / my_Nrho;
	double rmin = rho0;
	double rmax = rho1;
	if (iz >= my_Nz / 2) {
		//upper hemisphere
		double zmax = (iz + 1 - my_Nz / 2) * (2 * my_rho / my_Nz);
		double zmin = (iz - my_Nz / 2) * (2 * my_rho / my_Nz);
		rmin = rho0;
		if (zmax < my_rhoin) {
			rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
		}
		rmax = rho1;
		rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
	}
	else {
		//lower hemisphere, z inversed
		double zmax = fabs((iz - my_Nz / 2) * (2 * my_rho / my_Nz));
		double zmin = fabs((iz + 1 - my_Nz / 2) * (2 * my_rho / my_Nz));

		rmin = rho0;
		if (zmax < my_rhoin) {
			rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
		}
		rmax = rho1;
		rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
	}
	if (rmax < rho0 || rmin > rho1 || rmax < rmin) {
		area = 0;
	}
	else {
		area = pi * (rmax * rmax - rmin * rmin) / my_Nphi;
	}

	return area;
}

void SphericalLayerSource::evaluateLengthAndArea()
{
	//length
	//FILE* lengthFile = fopen("length.dat", "w");
	//FILE* areaFile = fopen("area.dat", "w");
	//FILE* volumeFile = fopen("volume.dat", "w");
	double dr = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;
	for (int irho = 0; irho < my_Nrho; ++irho) {
		double r = (irho + 0.5) * dr;
		double z1 = -sqrt(my_rho * my_rho - r * r);
		double z2 = 0;
		if (my_rhoin > r) {
			z2 = -sqrt(my_rhoin * my_rhoin - r * r);
		}
		double z3 = -z2;
		double z4 = -z1;
		for (int iz = 0; iz < my_Nz; ++iz) {
			double zmin = -my_rho + iz * dz;
			double zmax = zmin + dz;
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (zmin >= 0) {
					if (z3 > zmax) {
						my_length[irho][iz][iphi] = 0;
					} else
					if (z4 < zmin) {
						my_length[irho][iz][iphi] = 0;
					}
					else {
						double lowz = max(zmin, z3);
						double topz = min(zmax, z4);
						my_length[irho][iz][iphi] = topz - lowz;
						if (my_length[irho][iz][iphi] != my_length[irho][iz][iphi]) {
							printf("my_length = NaN\n");
							printLog("my_length = NaN\n");
							exit(0);
						}
					}
				}
				else {
					if (z1 > zmax) {
						my_length[irho][iz][iphi] = 0;
					} else if (z2 < zmin) {
						my_length[irho][iz][iphi] = 0;
					}
					else {
						double lowz = max(zmin, z1);
						double topz = min(zmax, z2);
						my_length[irho][iz][iphi] = topz - lowz;
						if (my_length[irho][iz][iphi] != my_length[irho][iz][iphi]) {
							printf("my_length = NaN\n");
							printLog("my_length = NaN\n");
							exit(0);
						}
					}
				}
				//fprintf(lengthFile, "%g\n", my_length[irho][iz][iphi]);
			}
		}
	}

	//area
	for (int irho = 0; irho < my_Nrho; ++irho) {
		double rho0 = irho * my_rho / my_Nrho;
		double rho1 = (irho + 1) * my_rho / my_Nrho;
		double rmin = rho0;
		double rmax = rho1;
		for (int iz = 0; iz < my_Nz; ++iz) {
			if (iz >= my_Nz / 2) {
				//upper hemisphere
				double zmin = -my_rho + iz * dz;
				double zmax = zmin + dz;
				rmin = rho0;
				if (zmax < my_rhoin) {
					rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
				}
				rmax = rho1;
				rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
			}
			else {
				//lower hemisphere, z inversed
				double zmax = fabs( - my_rho + iz * dz);
				double zmin = zmax - dz;

				rmin = rho0;
				if (zmax < my_rhoin) {
					rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
				}
				rmax = rho1;
				rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
			}
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (rmax < rho0 || rmin > rho1 || rmax < rmin) {
					my_area[irho][iz][iphi] = 0;
				}
				else {
					my_area[irho][iz][iphi] = pi * (rmax * rmax - rmin * rmin) / my_Nphi;
					if (my_area[irho][iz][iphi] != my_area[irho][iz][iphi]) {
						printf("my_area = NaN\n");
						printLog("my_area = NaN\n");
						exit(0);
					}
				}
				//fprintf(areaFile, "%g\n", my_area[irho][iz][iphi]);
			}
		}
	}

	/*for (int i = 0; i < my_Nrho; ++i) {
		for (int j = 0; j < my_Nz; ++j) {
			for (int k = 0; k < my_Nphi; ++k) {
				fprintf(volumeFile, "%g\n", my_length[i][j][k] * my_area[i][j][k]);
			}
		}
	}*/

	//fclose(lengthFile);
	//fclose(areaFile);
	//fclose(volumeFile);

	my_geometryCashed = true;
}

SphericalLayerSource::SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance, const double& velocity, const double& redShift) : RadiationSource(Nrho, Nz, Nphi, distance, redShift) {
	if (rhoin > rho) {
		printf("rhoin > rho in spherical layer source\n");
		printLog("rhoin > rho in spherical layer source\n");
		exit(0);
	}
	if (rho - rhoin < rho / Nrho) {
		printf("warning: rho - rhoin < dr\n");
		printLog("warning: rho - rhoin < dr\n");
	}
	if (Nz % 2 != 0) {
		printf("Nz in spherical source must be even\n");
		printLog("Nz in spherical source must be even\n");
		exit(0);
	}
	my_rho = rho;
	my_rhoin = rhoin;

	double dr = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;

	my_velocity = velocity;

	my_geometryCashed = false;
	my_area = new double** [my_Nrho];
	my_length = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_area[i] = new double* [my_Nz];
		my_length[i] = new double* [my_Nz];
		my_v[i] = new double* [my_Nz];
		my_vtheta[i] = new double* [my_Nz];
		my_vphi[i] = new double* [my_Nz];
		for (int j = 0; j < my_Nz; ++j) {
			my_area[i][j] = new double [my_Nphi];
			my_length[i][j] = new double [my_Nphi];
			my_v[i][j] = new double[my_Nphi];
			my_vtheta[i][j] = new double[my_Nphi];
			my_vphi[i][j] = new double[my_Nphi];
			for (int k = 0; k < my_Nphi; ++k) {
				my_area[i][j][k] = 0;
				my_length[i][j][k] = 0;
			}
		}
	}

	evaluateLengthAndArea();


	for (int i = 0; i < my_Nrho; ++i) {
		double r = (i + 0.5) * dr;
		for (int j = 0; j < my_Nz; ++j) {
			double z = -my_rho + j * dz;
			for (int k = 0; k < my_Nphi; ++k) {
				double phi = (k + 0.5) * 2 * pi / my_Nphi;
				my_v[i][j][k] = my_velocity;
				my_vtheta[i][j][k] = acos(z / sqr(r * r + z * z));
				my_vphi[i][j][k] = phi;
			}
		}
	}

}

SphericalLayerSource::SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : RadiationSource(Nrho, Nz, Nphi, distance, redShift)
{
	if (rhoin > rho) {
		printf("rhoin > rho in spherical layer source\n");
		printLog("rhoin > rho in spherical layer source\n");
		exit(0);
	}
	if (rho - rhoin < rho / Nrho) {
		printf("warning: rho - rhoin < dr\n");
		printLog("warning: rho - rhoin < dr\n");
	}
	if (Nz % 2 != 0) {
		printf("Nz in spherical source must be even\n");
		printLog("Nz in spherical source must be even\n");
		exit(0);
	}
	my_rho = rho;
	my_rhoin = rhoin;

	double dr = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;

	my_geometryCashed = false;
	my_area = new double** [my_Nrho];
	my_length = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_area[i] = new double* [my_Nz];
		my_length[i] = new double* [my_Nz];
		my_v[i] = new double* [my_Nz];
		my_vtheta[i] = new double* [my_Nz];
		my_vphi[i] = new double* [my_Nz];
		for (int j = 0; j < my_Nz; ++j) {
			my_area[i][j] = new double[my_Nphi];
			my_length[i][j] = new double[my_Nphi];
			my_v[i][j] = new double[my_Nphi];
			my_vtheta[i][j] = new double[my_Nphi];
			my_vphi[i][j] = new double[my_Nphi];
			for (int k = 0; k < my_Nphi; ++k) {
				my_area[i][j][k] = 0;
				my_length[i][j][k] = 0;
			}
		}
	}

	evaluateLengthAndArea();


	for (int i = 0; i < my_Nrho; ++i) {
		double r = (i + 0.5) * dr;
		for (int j = 0; j < my_Nz; ++j) {
			double z = -my_rho + j * dz;
			for (int k = 0; k < my_Nphi; ++k) {
				double phi = (k + 0.5) * 2 * pi / my_Nphi;
				my_v[i][j][k] = velocity[i][j][k];
				my_vtheta[i][j][k] = vtheta[i][j][k];
				my_vphi[i][j][k] = vphi[i][j][k];
			}
		}
	}

	my_velocity = evaluateAverageVelocity();
}

SphericalLayerSource::~SphericalLayerSource()
{
	for (int i = 0; i < my_Nrho; ++i) {
		for (int j = 0; j < my_Nz; ++j) {
			delete[] my_area[i][j];
			delete[] my_length[i][j];
			delete[] my_v[i][j];
			delete[] my_vtheta[i][j];
			delete[] my_vphi[i][j];
		}
		delete[] my_area[i];
		delete[] my_length[i];
		delete[] my_v[i];
		delete[] my_vtheta[i];
		delete[] my_vphi[i];
	}
	delete[] my_area;
	delete[] my_length;
	delete[] my_v;
	delete[] my_vtheta;
	delete[] my_vphi;
}

/*double SphericalLayerSource::getLength(int irho, int iz, int iphi) {
	double dr = my_rho / my_Nrho;
	double r = (irho + 0.5) * dr;
	double z1 = -sqrt(my_rho * my_rho - r * r);
	double z2 = 0;
	if (my_rhoin > r) {
		z2 = -sqrt(my_rhoin * my_rhoin - r * r);
	}
	double z3 = -z2;
	double z4 = -z1;
	double dz = 2 * my_rho / my_Nz;
	double zmin = -my_rho + iz * dz;
	double zmax = zmin + dz;
	if (zmin >= 0) {
		if (z3 > zmax) {
			return 0;
		}
		if (z4 < zmin) {
			return 0;
		}
		double lowz = max(zmin, z3);
		double topz = min(zmax, z4);
		return topz - lowz;
	}
	else {
		if (z1 > zmax) {
			return 0;
		}
		if (z2 < zmin) {
			return 0;
		}
		double lowz = max(zmin, z1);
		double topz = min(zmax, z2);
		return topz - lowz;
	}
}

double SphericalLayerSource::getArea(int irho, int iz, int iphi)
{
	double rho0 = irho * getMaxRho() / my_Nrho;
	double rho1 = (irho + 1) * getMaxRho() / my_Nrho;
	double rmin = rho0;
	double rmax = rho1;
	if (iz >= my_Nz / 2) {
		//upper hemisphere
		double zmax = (iz + 1 - my_Nz / 2) * (2 * my_rho / my_Nz);
		double zmin = (iz - my_Nz / 2) * (2 * my_rho / my_Nz);
		rmin = rho0;
		if (zmax < my_rhoin) {
			rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
		}
		rmax = rho1;
		rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
	}
	else {
		//lower hemisphere, z inversed
		double zmax = fabs((iz - my_Nz / 2) * (2 * my_rho / my_Nz));
		double zmin = fabs((iz + 1 - my_Nz / 2) * (2 * my_rho / my_Nz));

		rmin = rho0;
		if (zmax < my_rhoin) {
			rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
		}
		rmax = rho1;
		rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
	}
	if (rmax < rho0 || rmin > rho1 || rmax < rmin) {
		return 0;
	}

	return 2 * pi * (rmax * rmax - rmin * rmin) / my_Nphi;
}*/

double SphericalLayerSource::getRho(int irho) {
	return (irho + 0.5) * my_rho / my_Nrho;
}
double SphericalLayerSource::getZ(int iz) {
	return -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
}
double SphericalLayerSource::getPhi(int iphi) {
	return iphi * 2 * pi / my_Nphi;
}

int SphericalLayerSource::getRhoIndex(const double& rho) {
	if (rho < 0) {
		printf("rho < 0 in get rhoIndex\n");
		printLog("rho < 0 in get rhoIndex\n");
		exit(0);
	}
	if (rho > my_rho) {
		printf("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		printLog("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		exit(0);
	}
	if (rho == my_rho) {
		return my_Nrho - 1;
	}
	double drho = my_rho / my_Nrho;
	return floor(rho / drho);
}

double SphericalLayerSource::getLength(int irho, int iz, int iphi) {
	if (!my_geometryCashed) {
		evaluateLengthAndArea();
	}
	return my_length[irho][iz][iphi];
	//return evaluateLength(irho, iz, iphi);
}

double SphericalLayerSource::getArea(int irho, int iz, int iphi)
{
	if (!my_geometryCashed) {
		evaluateLengthAndArea();
	}
	return my_area[irho][iz][iphi];
	//return evaluateArea(irho, iz, iphi);
}

double SphericalLayerSource::getCrossSectionArea(int irho, int iphi)
{
	double rho1 = (irho + 1) * my_rho / my_Nrho;
	double rho0 = irho * my_rho / my_Nrho;
	return (pi/my_Nphi)*(rho1*rho1 - rho0*rho0);
}

double SphericalLayerSource::getMaxRho() {
	return my_rho;
}
double SphericalLayerSource::getMinRho() {
	return 0;
}
double SphericalLayerSource::getMinZ() {
	return -my_rho;
}
double SphericalLayerSource::getMaxZ() {
	return my_rho;
}
double SphericalLayerSource::getTotalVolume() {
	return 4 * pi * (my_rho * my_rho * my_rho - my_rhoin * my_rhoin * my_rhoin) / 3;
}
double SphericalLayerSource::getInnerRho() {
	return my_rhoin;
}

TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity, const double& redShift) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance, velocity, redShift) {
	my_distribution = electronDistribution;
	
	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			double rho1 = (irho + 0.5) * my_rho / my_Nrho;
			double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;

			double r = sqrt(z * z + rho1 * rho1);

			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				//my_B[irho][iz][iphi] = B[irho][iz][iphi]*(my_rho/r);
				my_theta[irho][iz][iphi] = theta[irho][iz][iphi];
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				//my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi] * sqr(my_rho / r);
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}
TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity, const double& redShift) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance, velocity, redShift) {
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_theta[irho][iz][iphi] = theta;
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance, velocity, vtheta, vphi, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			double rho1 = (irho + 0.5) * my_rho / my_Nrho;
			double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;

			double r = sqrt(z * z + rho1 * rho1);

			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				//my_B[irho][iz][iphi] = B[irho][iz][iphi]*(my_rho/r);
				my_theta[irho][iz][iphi] = theta[irho][iz][iphi];
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				//my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi] * sqr(my_rho / r);
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance, velocity, vtheta, vphi, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_theta[irho][iz][iphi] = theta;
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSphericalLayerSource::~TabulatedSphericalLayerSource() {
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_B[irho][iz];
			delete[] my_theta[irho][iz];
			delete[] my_phi[irho][iz];
		}
		delete[] my_B[irho];
		delete[] my_theta[irho];
		delete[] my_phi[irho];
		delete[] my_isSource[irho];
	}
	delete[] my_B;
	delete[] my_theta;
	delete[] my_phi;
	delete[] my_isSource;
}

bool TabulatedSphericalLayerSource::rayTraceToNextCell(const double& rho0, const double& z0, int iphi, const double& theta, double& rho1, double& z1, double& lB2) {
	double drho = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;
	int irho = floor(rho0 / drho);

	double r = sqrt(rho0 * rho0 + z0 * z0);
	if (r >= my_rho) {
		rho1 = rho0;
		z1 = z0;
		lB2 = 0;
		return true;
	}

	double cosTheta = cos(theta);
	double sinTheta = sin(theta);

	double nextRho = drho * (irho + 1);
	double deltaRho = nextRho - rho0;
	if (deltaRho <= 0) {
		irho = irho + 1;
		nextRho = drho * (irho + 1);
		deltaRho = nextRho - rho0;
	}

	if (cosTheta > 0) {
		int iz = floor((z0 + my_rho) / dz);
		double B2 = sqr(getB(irho, iz, iphi));
		double nextZ = dz * (iz + 1) - my_rho;
		double deltaZ = nextZ - z0;
		if (deltaZ <= 0) {
			nextZ = dz * (iz + 2) - my_rho;
			deltaZ = dz;
			B2 = sqr(getB(irho, iz + 1, iphi));
		}

		double rterm = cosTheta * deltaRho;
		double zterm = sinTheta * deltaZ;

		double l = 0;

		if (rterm > zterm) {
			//sin can be = 0;
			
			l = deltaZ / cosTheta;
			z1 = nextZ;
			rho1 = rho0 + sqrt(l * l - deltaZ * deltaZ);
		}
		else if (zterm > rterm) {
			//cos can be = 0;
			l = deltaRho / sinTheta;
			rho1 = nextRho;
			z1 = z0 + sqrt(l * l - deltaRho * deltaRho);
		}
		else {
			z1 = nextZ;
			rho1 = nextRho;
			l = sqrt(deltaRho * deltaRho + deltaZ * deltaZ);
		}
		double nextR = sqrt(rho1 * rho1 + z1 * z1);
		if (nextR > my_rho) {
			l = l - (nextR - my_rho);
			if (l < 0) {
				printf("l < 0 in rayTraceToNextCell\n");
				printLog("l < 0 in rayTraceToNextCell\n");
				exit(0);
			}
			//todo change rgo1 and z1
			lB2 = l * B2;
			return true;
		}
		lB2 = l * B2;
		return false;
	} if (cosTheta < 0) {
		int iz = my_Nz/2 - floor(fabs((z0) / dz));
		double B2 = sqr(getB(irho, iz, iphi));
		double nextZ = dz * (iz - 1) - my_rho;
		double deltaZ = -(nextZ - z0);
		if (deltaZ <= 0) {
			nextZ = dz * (iz - 2) - my_rho;
			deltaZ = dz;
			B2 = sqr(getB(irho, iz-1, iphi));
		}

		double rterm = -cosTheta * deltaRho;
		double zterm = sinTheta * deltaZ;

		double l = 0;

		if (rterm > zterm) {
			//sin can be = 0;

			l = - deltaZ / cosTheta;
			z1 = nextZ;
			rho1 = rho0 + sqrt(l * l - deltaZ * deltaZ);
		}
		else if (zterm > rterm) {
			//cos can be = 0;
			l = deltaRho / sinTheta;
			rho1 = nextRho;
			z1 = z0 + sqrt(l * l - deltaRho * deltaRho);
		}
		else {
			z1 = nextZ;
			rho1 = nextRho;
			l = sqrt(deltaRho*deltaRho + deltaZ*deltaZ);
		}
		double nextR = sqrt(rho1*rho1 + z1*z1);
		if (nextR > my_rho) {
			l = l - (nextR - my_rho);
			if (l < 0) {
				printf("l < 0 in rayTraceToNextCell\n");
				printLog("l < 0 in rayTraceToNextCell\n");
				exit(0);
			}
			//todo change rgo1 and z1
			lB2 = l * B2;
			return true;
		}
		lB2 = l * B2;
		return false;
	}
	else {
		z1 = z0;
		rho1 = rho0;
		double l = deltaRho;
		//todo
		double B2 = sqr(getB(irho, my_Nz, iphi));
		double nextR = sqrt(rho1 * rho1 + z1 * z1);
		if (nextR > my_rho) {
			l = l - (nextR - my_rho);
			if (l < 0) {
				printf("l < 0 in rayTraceToNextCell\n");
				printLog("l < 0 in rayTraceToNextCell\n");
				exit(0);
			}
			//todo change rgo1 and z1
			lB2 = l * B2;
			return true;
		}
		lB2 = l * B2;
		return false;
	}
}

double TabulatedSphericalLayerSource::evaluateTotalLB2fromPoint(const double& rho0, const double& z0, int iphi, const double& theta) {
	double r = sqrt(rho0 * rho0 + z0 * z0);
	if (r > my_rho) {
		return 0.0;
	}
	double rho1;
	double z1;
	bool reachedSurface = false;
	double totalLB2 = 0;
	double lB2 = 0;
	double tempRho = rho0;
	double tempZ = z0;
	int numberOfIteration = 0;
	while (!reachedSurface) {
		if (numberOfIteration > 1.5 * (my_Nrho + my_Nz)) {
			printf("aaa\n");
		}
		reachedSurface = rayTraceToNextCell(tempRho, tempZ, iphi, theta, rho1, z1, lB2);
		totalLB2 += lB2;
		tempRho = rho1;
		tempZ = z1;
		numberOfIteration++;
		if (numberOfIteration > 2 * (my_Nrho + my_Nz)) {
			printf("infinite number of iterations in evaluateTotalLB2fromPoint\n");
			printLog("infinite number of iterations in evaluateTotalLB2fromPoint\n");
			exit(0);
		}
	}
	return totalLB2;
}

void TabulatedSphericalLayerSource::setMask(bool** mask) {
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			my_isSource[irho][iphi] = mask[irho][iphi];
		}
	}
}

bool TabulatedSphericalLayerSource::isSource(int irho, int iphi) {
	return my_isSource[irho][iphi];
}

double TabulatedSphericalLayerSource::getB(int irho, int iz, int iphi) {
	return my_B[irho][iz][iphi];
}
double TabulatedSphericalLayerSource::getMaxB()
{
	double Bmax = 0;
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (my_B[irho][iz][iphi] > Bmax) {
					Bmax = my_B[irho][iz][iphi];
				}
			}
		}
	}
	return Bmax;
}
double TabulatedSphericalLayerSource::getMaxOuterB()
{
	double Bmax = 0;
	int irho = my_Nrho - 1;
	int iz = 0;
	for (iz = 0; iz < my_Nz; ++iz) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}
		}
	}

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			iz = 0;
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}

			iz = my_Nz - 1;
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}
		}
	}

	return Bmax;
}
double TabulatedSphericalLayerSource::getAverageSigma()
{
	double magneticEnergy = 0;
	double restEnergy = 0;


	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				magneticEnergy += my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * getVolume(irho, iz, iphi) / (4 * pi);
				restEnergy += my_concentration[irho][iz][iphi] * massProton * speed_of_light2 * getVolume(irho, iz, iphi);
			}
		}
	}

	return magneticEnergy/restEnergy;
}
double TabulatedSphericalLayerSource::getAverageBsquared()
{
	double magneticEnergy = 0;
	double volume = 0;


	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				magneticEnergy += my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * getVolume(irho, iz, iphi);
				volume += getVolume(irho, iz, iphi);
			}
		}
	}

	return sqrt(magneticEnergy / volume);
}
double TabulatedSphericalLayerSource::getAverageConcentration()
{
	double concentration = 0;
	double volume = 0;

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				concentration += my_concentration[irho][iz][iphi] * getVolume(irho, iz, iphi);
				volume += getVolume(irho, iz, iphi);
			}
		}
	}

	double result =  concentration/volume;
	if (result != result) {
		printf("averageConcentration = NaN\n");
		printLog("averageConcentration = NaN\n");
		exit(0);
	}
	return result;
}
double TabulatedSphericalLayerSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration[irho][iz][iphi];
}
void TabulatedSphericalLayerSource::getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi)
{
	velocity = my_v[irho][iz][iphi];
	if (velocity > speed_of_light) {
		printf("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		printLog("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		exit(0);
	}
	theta = my_vtheta[irho][iz][iphi];
	phi = my_vphi[irho][iz][iphi];
}
double TabulatedSphericalLayerSource::getSinTheta(int irho, int iz, int iphi) {
	return sin(my_theta[irho][iz][iphi]);
}
double TabulatedSphericalLayerSource::getBTheta(int irho, int iz, int iphi)
{
	return my_theta[irho][iz][iphi];
}
double TabulatedSphericalLayerSource::getBPhi(int irho, int iz, int iphi)
{
	return my_phi[irho][iz][iphi];
}
/*void TabulatedSphericalLayerSource::resetConcentration(const double& concentration)
{
	my_distribution->resetConcentration(concentration);
}*/
void TabulatedSphericalLayerSource::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* B[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* Rin = R*(1 - parameters[3])
*/
	my_geometryCashed = false;

	my_rho = parameters[0] * normalizationUnits[0];
	my_rhoin = my_rho * (1.0 - parameters[3] * normalizationUnits[3]);
	evaluateLengthAndArea();
	double sigma = parameters[1] * normalizationUnits[1];
	//double B0 = my_B[my_Nrho - 1][0][0];
	//double n0 = my_concentration[my_Nrho - 1][0][0];
	double n0 = getAverageConcentration();
	//double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	double sigma0 = getAverageSigma();
	double old_velocity = my_velocity;
	my_velocity = parameters[4] * normalizationUnits[4];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double localSigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				localSigma *= sigma / sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2] / n0;
				if (my_concentration[irho][iz][iphi] != my_concentration[irho][iz][iphi]) {
					printf("my_concentration = NaN\n");
					printLog("my_concentration = NaN\n");
					exit(0);
				}
				my_B[irho][iz][iphi] = sqrt(localSigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				if (old_velocity > 0) {
					my_v[irho][iz][iphi] *= my_velocity / old_velocity;
				}
				else {
					my_v[irho][iz][iphi] = my_velocity;
				}
			}
		}
	}
}
MassiveParticleDistribution* TabulatedSphericalLayerSource::getParticleDistribution(int irho, int iz, int iphi) {
	my_distribution->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distribution;
}

TabulatedSphericalLayerSource2::TabulatedSphericalLayerSource2(int Ndistributions, double* velocities, MassiveParticleIsotropicDistribution** distributions, int Nrho, int Nz, int Nphi, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, theta, phi, concentration, rho, rhoin, distance, velocity, redShift) {
	my_Ndistributions = Ndistributions;
	my_maxThreads = omp_get_max_threads();
	my_velocities = new double[my_Ndistributions];
	my_outputDistributions = new MassiveParticleDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_outputDistributions[i] = NULL;
	}
	my_distributions = new MassiveParticleDistribution * [my_Ndistributions];
	for (int i = 0; i < my_Ndistributions; ++i) {
		my_velocities[i] = velocities[i];
		my_distributions[i] = distributions[i];
	}
}

MassiveParticleDistribution* TabulatedSphericalLayerSource2::getParticleDistribution(int irho, int iz, int iphi)
{
	int num_thread = omp_get_thread_num();
	//todo not thread safe!!!
	if (my_outputDistributions[num_thread] != NULL) {
		return my_outputDistributions[num_thread];
	}

	/*for (int i = 0; i < my_Ndistributions; ++i) {
		my_distributions[i]->resetConcentration(getConcentration(irho, iz, iphi));
	}*/
	if (my_velocity <= my_velocities[0]) {
		MassiveParticleIsotropicDistribution* iso1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[0]);
		if (iso1 != NULL) {
			my_outputDistributions[num_thread] = new CompoundWeightedMassiveParticleIsotropicDistribution(iso1, 0.5, iso1, 0.5);
		}
		else {
			my_outputDistributions[num_thread] = new CompoundWeightedMassiveParticleDistribution(my_distributions[0], 0.5, my_distributions[0], 0.5);
		}
	}
	else if (my_velocity >= my_velocities[my_Ndistributions - 1]) {
		MassiveParticleIsotropicDistribution* iso1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[my_Ndistributions - 1]);
		if (iso1 != NULL) {
			my_outputDistributions[num_thread] = new CompoundWeightedMassiveParticleIsotropicDistribution(iso1, 0.5, iso1, 0.5);
		}
		else {
			my_outputDistributions[num_thread] = new CompoundWeightedMassiveParticleDistribution(my_distributions[my_Ndistributions - 1], 0.5, my_distributions[my_Ndistributions - 1], 0.5);
		}
	}
	else {
		for (int i = 1; i < my_Ndistributions; ++i) {
			if (my_velocities[i] > my_velocity) {
				double left = (my_velocities[i] - my_velocity) / (my_velocities[i] - my_velocities[i - 1]);
				double right = (my_velocity - my_velocities[i - 1]) / (my_velocities[i] - my_velocities[i - 1]);
				MassiveParticleIsotropicDistribution* iso1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[i - 1]);
				MassiveParticleIsotropicDistribution* iso2 = dynamic_cast<MassiveParticleIsotropicDistribution*>(my_distributions[i]);
				if ((iso1 != NULL) && (iso2 != NULL)) {
					my_outputDistributions[num_thread] = new CompoundWeightedMassiveParticleIsotropicDistribution(iso1, left, iso2, right);
				}
				else {
					my_outputDistributions[num_thread] = new CompoundWeightedMassiveParticleDistribution(my_distributions[i - 1], left, my_distributions[i], right);
				}
				return my_outputDistributions[num_thread];
			}
		}
	}
	return my_outputDistributions[num_thread];
}

void TabulatedSphericalLayerSource2::resetParameters(const double* parameters, const double* normalizationUnits)
{
	if (my_distribution != NULL) {
		delete my_distribution;
		my_distribution = NULL;
	}
	for (int i = 0; i < my_maxThreads; ++i) {
		if (my_outputDistributions[i] != NULL) {
			delete my_outputDistributions[i];
			my_outputDistributions[i] = NULL;
		}
	}

	TabulatedSphericalLayerSource::resetParameters(parameters, normalizationUnits);
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, phi, concentration, rho, rhoin, distance, velocity, redShift)
{
	my_Ntheta = Ntheta;
	my_distributions = new MassiveParticleDistribution*[my_Ntheta];

	for (int i = 0; i < my_Ntheta; ++i) {
		my_distributions[i] = electronDistributions[i];
	}

	my_shockWaveAngle = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_shockWaveAngle[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_shockWaveAngle[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				double rho = (irho + 0.5) * my_rho / my_Nrho;
				double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
				double alpha = (iphi + 0.5) * 2 * pi / my_Nphi;

				double r = sqrt(z * z + rho * rho);
				double nz = z / r;
				double nx = rho * cos(alpha) / r;
				double ny = rho * sin(alpha) / r;

				double Bz = my_B[irho][iz][iphi] * cos(my_theta[irho][iz][iphi]);
				double Bx = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * cos(my_phi[irho][iz][iphi]);
				double By = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * sin(my_phi[irho][iz][iphi]);

				double cosTheta = (nz * Bz + nx * Bx + ny * By) / my_B[irho][iz][iphi];

				if (cosTheta > 1.0 || cosTheta < -1.0) {
					printf("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
                    printLog("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
				exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
				//my_shockWaveAngle[irho][iz][iphi] = pi*4/9;
			}
		}
	}
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, phi, concentration, rho, rhoin, distance, velocity, redShift)
{
	my_Ntheta = Ntheta;
	my_distributions = new MassiveParticleDistribution * [my_Ntheta];

	for (int i = 0; i < my_Ntheta; ++i) {
		my_distributions[i] = electronDistributions[i];
	}

	my_shockWaveAngle = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_shockWaveAngle[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_shockWaveAngle[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_concentration[irho][iz][iphi] = concentration;
				double rho = (irho + 0.5) * my_rho / my_Nrho;
				double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
				double alpha = (iphi + 0.5) * 2 * pi / my_Nphi;

				double r = sqrt(z * z + rho * rho);
				double nz = z / r;
				double nx = rho * cos(alpha) / r;
				double ny = rho * sin(alpha) / r;

				double Bz = my_B[irho][iz][iphi] * cos(my_theta[irho][iz][iphi]);
				double Bx = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * cos(my_phi[irho][iz][iphi]);
				double By = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * sin(my_phi[irho][iz][iphi]);

				double cosTheta = (nz * Bz + nx * Bx + ny * By) / my_B[irho][iz][iphi];

				if (cosTheta > 1.0 || cosTheta < -1.0) {
					printf("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
                    printLog("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
				//my_shockWaveAngle[irho][iz][iphi] = pi*4/9;
			}
		}
	}
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, phi, concentration, rho, rhoin, distance, velocity, vtheta, vphi, redShift)
{
	my_Ntheta = Ntheta;
	my_distributions = new MassiveParticleDistribution * [my_Ntheta];

	for (int i = 0; i < my_Ntheta; ++i) {
		my_distributions[i] = electronDistributions[i];
	}

	my_shockWaveAngle = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_shockWaveAngle[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_shockWaveAngle[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				double rho = (irho + 0.5) * my_rho / my_Nrho;
				double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
				double alpha = (iphi + 0.5) * 2 * pi / my_Nphi;

				double r = sqrt(z * z + rho * rho);
				double nz = z / r;
				double nx = rho * cos(alpha) / r;
				double ny = rho * sin(alpha) / r;

				double Bz = my_B[irho][iz][iphi] * cos(my_theta[irho][iz][iphi]);
				double Bx = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * cos(my_phi[irho][iz][iphi]);
				double By = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * sin(my_phi[irho][iz][iphi]);

				double cosTheta = (nz * Bz + nx * Bx + ny * By) / my_B[irho][iz][iphi];

				if (cosTheta > 1.0 || cosTheta < -1.0) {
					printf("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					printLog("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
				//my_shockWaveAngle[irho][iz][iphi] = pi*4/9;
			}
		}
	}
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, phi, concentration, rho, rhoin, distance, velocity, vtheta, vphi, redShift)
{
	my_Ntheta = Ntheta;
	my_distributions = new MassiveParticleDistribution * [my_Ntheta];

	for (int i = 0; i < my_Ntheta; ++i) {
		my_distributions[i] = electronDistributions[i];
	}

	my_shockWaveAngle = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_shockWaveAngle[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_shockWaveAngle[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_concentration[irho][iz][iphi] = concentration;
				double rho = (irho + 0.5) * my_rho / my_Nrho;
				double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
				double alpha = (iphi + 0.5) * 2 * pi / my_Nphi;

				double r = sqrt(z * z + rho * rho);
				double nz = z / r;
				double nx = rho * cos(alpha) / r;
				double ny = rho * sin(alpha) / r;

				double Bz = my_B[irho][iz][iphi] * cos(my_theta[irho][iz][iphi]);
				double Bx = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * cos(my_phi[irho][iz][iphi]);
				double By = my_B[irho][iz][iphi] * sin(my_theta[irho][iz][iphi]) * sin(my_phi[irho][iz][iphi]);

				double cosTheta = (nz * Bz + nx * Bx + ny * By) / my_B[irho][iz][iphi];

				if (cosTheta > 1.0 || cosTheta < -1.0) {
					printf("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					printLog("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
				//my_shockWaveAngle[irho][iz][iphi] = pi*4/9;
			}
		}
	}
}

AngleDependentElectronsSphericalSource::~AngleDependentElectronsSphericalSource()
{
	delete[] my_distributions;

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_shockWaveAngle[irho][iz];
			delete[] my_concentration[irho][iz];
		}
		delete[] my_shockWaveAngle[irho];
		delete[] my_concentration[irho];
	}
	delete[] my_shockWaveAngle;
	delete[] my_concentration;
}

/*void AngleDependentElectronsSphericalSource::resetConcentration(const double& concentration)
{
	for(int i = 0; i < my_Ntheta; ++i){
		my_distributions[i]->resetConcentration(concentration);
}
}*/

double AngleDependentElectronsSphericalSource::getShockWaveAngle(int irho, int iz, int iphi)
{
	return my_shockWaveAngle[irho][iz][iphi];
}

MassiveParticleDistribution* AngleDependentElectronsSphericalSource::getParticleDistribution(int irho, int iz, int iphi)
{
	double dtheta = 0.5*pi / (my_Ntheta-1);
	double theta = my_shockWaveAngle[irho][iz][iphi];
	if (theta > pi / 2) {
		theta = pi - theta;
	}

	int angleIndex = floor((theta + 0.5*dtheta) / dtheta);
	if (angleIndex == my_Ntheta) {
		angleIndex = my_Ntheta - 1;
	}

	//for debug
	//angleIndex = 3;

	my_distributions[angleIndex]->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distributions[angleIndex];

}

ExpandingRemnantSource::ExpandingRemnantSource(const double& R0, const double& B0, const double& concentration0, const double& v, const double& widthFraction, RadiationSource* source, const double& t0, const double& velocityPower, const double& Bpower, const double& concentrationPower, const double& widthPower) : RadiationTimeDependentSource(source, t0) {
	my_R0 = R0;
	my_B0 = B0;
	my_concentration0 = concentration0;
	my_widthFraction = widthFraction;
	my_v = v;
	my_velocityPower = velocityPower;
	my_Bpower = Bpower;
	my_concentrationPower = concentrationPower;
	my_widthPower = widthPower;
}

void ExpandingRemnantSource::resetParameters(const double* parameters, const double* normalizationUnits) {
	my_R0 = parameters[0] * normalizationUnits[0];
	double sigma = parameters[1] * normalizationUnits[1];
	my_concentration0 = parameters[2] * normalizationUnits[2];
	my_widthFraction = parameters[3] * normalizationUnits[3];
	my_v = parameters[4] * normalizationUnits[4];
	my_B0 = sqrt(sigma * 4 * pi * massProton * my_concentration0 * speed_of_light2);
	my_velocityPower = parameters[5] * normalizationUnits[5] - 1.0;
	my_concentrationPower = parameters[6] * normalizationUnits[6] - 1.0;
	my_Bpower = parameters[7] * normalizationUnits[7] - 1.0;
	my_widthPower = parameters[8] * normalizationUnits[8] - 1.0;
}

//just one possible example
RadiationSource* ExpandingRemnantSource::getRadiationSource(double& time, const double* normalizationUnits) {
	//double R = my_R0 + my_v * (time - my_t0);
	double R;
	if (my_velocityPower == 1.0) {
		R = my_R0 + my_v * my_t0 * log(time / my_t0);
	}
	else {
		R = my_R0 + my_v * my_t0 * (1.0 - pow(time / my_t0, 1.0 - my_velocityPower))/(my_velocityPower - 1.0);
	}
	//double sigma = sqr(my_B0) / (4 * pi * massProton * my_concentration0 * speed_of_light2)/(my_R0/R);
	double sigma = sqr(my_B0) / (4 * pi * massProton * my_concentration0 * speed_of_light2)*pow(my_R0/R, 2*my_Bpower- my_concentrationPower);
	//double B = my_B0;
	//double n = my_concentration0 * sqr(my_R0 / R);
	//double n = my_concentration0;
	//double n = my_concentration0*cube(my_R0/R);
	double n = my_concentration0*pow(my_R0/R, my_concentrationPower);
	double fracton = my_widthFraction * pow(my_R0/R, my_widthPower);
	//double fracton = my_widthFraction;

	double localV = my_v * pow(time / my_t0, -my_velocityPower);

	double parameters[5];
	parameters[0] = R / normalizationUnits[0];
	parameters[1] = sigma / normalizationUnits[1];
	parameters[2] = n / normalizationUnits[2];
	parameters[3] = fracton / normalizationUnits[3];
	parameters[4] = localV / normalizationUnits[4];
	my_radiationSource->resetParameters(parameters, normalizationUnits);
	return my_radiationSource;
}

TabulatedSLSourceWithSynchCutoff::TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, const double& velocity, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, distance, velocity, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution*[my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
}

TabulatedSLSourceWithSynchCutoff::TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, const double& velocity, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, distance, velocity, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
}

TabulatedSLSourceWithSynchCutoff::TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, distance, velocity, vtheta, vphi, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
}

TabulatedSLSourceWithSynchCutoff::TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, distance, velocity, vtheta, vphi, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
}

TabulatedSLSourceWithSynchCutoff::~TabulatedSLSourceWithSynchCutoff()
{
	for (int i = 0; i < my_maxThreads; ++i) {
		if (my_localDistribution[i] != NULL) {
			delete my_localDistribution[i];
		}
	}
	delete[] my_localDistribution;
	delete3dArray(my_LB2, my_Nrho, my_Nz, my_Nphi);
}

void TabulatedSLSourceWithSynchCutoff::updateLB2() {
	double drho = my_rho / my_Nrho;
	double dz = 2 * my_rho / my_Nz;
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double rho = (irho + 0.5) * drho;
				double z = (iz + 0.5) * dz - my_rho;
				double theta = acos(z / sqrt(rho * rho + z * z));
				//my_LB2[irho][iz][iphi] = evaluateTotalLB2fromPoint(rho, z, iphi, theta);
				double l = my_rho - sqrt(rho*rho + z*z);
				if (l < 0) {
					l = 0;
				}
				my_LB2[irho][iz][iphi] = my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * l;
			}
		}
	}
}

void TabulatedSLSourceWithSynchCutoff::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* B[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* Rin = R*(1 - parameters[3])
*/
	my_geometryCashed = false;

	my_rho = parameters[0] * normalizationUnits[0];
	my_rhoin = my_rho * (1.0 - parameters[3] * normalizationUnits[3]);
	evaluateLengthAndArea();
	double sigma = parameters[1] * normalizationUnits[1];
	//double B0 = my_B[my_Nrho - 1][0][0];
	//double n0 = my_concentration[my_Nrho - 1][0][0];
	double n0 = getAverageConcentration();
	//double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	double sigma0 = getAverageSigma();
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double localSigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				localSigma *= sigma / sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2] / n0;
				if (my_concentration[irho][iz][iphi] != my_concentration[irho][iz][iphi]) {
					printf("my_concentration = NaN\n");
					printLog("my_concentration = NaN\n");
					exit(0);
				}
				my_B[irho][iz][iphi] = sqrt(localSigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
			}
		}
	}
	my_velocity = parameters[4] * normalizationUnits[4];
	my_meanB = getAverageBsquared();
	updateLB2();
	for (int i = 0; i < my_maxThreads; ++i) {
		if (my_localDistribution[i] != NULL) {
			delete my_localDistribution[i];
			my_localDistribution[i] = NULL;
		}
	}
}

MassiveParticleDistribution* TabulatedSLSourceWithSynchCutoff::getParticleDistribution(int irho, int iz, int iphi)
{
	int numthreads = omp_get_thread_num();
	if (my_localDistribution[numthreads] != NULL) {
		delete my_localDistribution[numthreads];
		my_localDistribution[numthreads] = NULL;
	}
	my_localDistribution[numthreads] = new MassiveParticleTabulatedIsotropicDistribution(*my_cutoffDistribution);
	double LB2 = my_LB2[irho][iz][iphi];
	if (LB2 <= 0) {
		/*printf("l <= 0 in TabulatedSLSourceWithSynchCutoff::getParticleDistribution irho = %d iz = %d\n", irho, iz);
		printf("rho = %g z = %g r = %g R = %g\n", rho, z, r, my_rho);
		printLog("l <= 0 in TabulatedSLSourceWithSynchCutoff::getParticleDistribution irho = %d iz = %d\n", irho, iz);
		printLog("rho = %g z = %g r = %g R = %g\n", rho, z, r, my_rho);*/
		//exit(0);
	}
	if (LB2 > 0) {
		double mass = my_cutoffDistribution->getMass();
		double Ecut = mass * mass * mass * mass * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * LB2);
		if (Ecut < mass * speed_of_light2) {
			//todo
			Ecut = mass * speed_of_light2 * 2.0;
		}
		//my_cutoffDistribution->resetEcut(Ecut);
		//my_localDistribution->setToZeroAboveE(Ecut);
		//my_localDistribution->addExponentialCutoff(Ecut);
		double rho = getRho(irho);
		double z = getZ(iz);
		double l = my_rho - sqrt(rho * rho + z * z);
		if (l < 0) {
			l = 0;
		}
		double lossRate = (4.0 / 9.0) * electron_charge * electron_charge * electron_charge * electron_charge * my_meanB * my_meanB / (mass * mass * mass * mass * pow(speed_of_light, 7.0));
		double time = l / my_downstreamVelocity;
		my_localDistribution[numthreads]->transformToLosses(lossRate, time);
	}
	my_localDistribution[numthreads]->resetConcentration(getConcentration(irho, iz, iphi));
	return my_localDistribution[numthreads];
	//return std::unique_ptr<MassiveParticleDistribution>(new MassiveParticleMaxwellDistribution(massElectron, 1000, 1.0));
}

TabulatedDiskSourceWithSynchCutoff::TabulatedDiskSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, const double& velocity, const double& redShift) : TabulatedDiskSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, z, distance, velocity, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution*[my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
}

TabulatedDiskSourceWithSynchCutoff::TabulatedDiskSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, const double& velocity, const double& redShift) : TabulatedDiskSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, z, distance, velocity, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
}

TabulatedDiskSourceWithSynchCutoff::TabulatedDiskSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedDiskSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, z, distance, velocity, vtheta, vphi, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
}

TabulatedDiskSourceWithSynchCutoff::TabulatedDiskSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& concentration, const double& theta, const double& phi, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedDiskSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, z, distance, velocity, vtheta, vphi, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticleTabulatedIsotropicDistribution\n");
		exit(0);
	}
	my_maxThreads = omp_get_max_threads();
	my_localDistribution = new MassiveParticleTabulatedIsotropicDistribution * [my_maxThreads];
	for (int i = 0; i < my_maxThreads; ++i) {
		my_localDistribution[i] = NULL;
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
}

TabulatedDiskSourceWithSynchCutoff::~TabulatedDiskSourceWithSynchCutoff()
{
	for (int i = 0; i < my_maxThreads; ++i) {
		if (my_localDistribution[i] != NULL) {
			delete my_localDistribution[i];
		}
	}
	delete[] my_localDistribution;
}

void TabulatedDiskSourceWithSynchCutoff::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* sigma[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* z = R*parameters[3]
*/
	my_rho = parameters[0] * normalizationUnits[0];
	double sigma = parameters[1] * normalizationUnits[1];
	//double B0 = my_B[my_Nrho - 1][0][0];
	//double n0 = my_concentration[my_Nrho - 1][0][0];
	double n0 = getAverageConcentration();
	//double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	double sigma0 = getAverageSigma();
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double sigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				sigma *= parameters[1] * normalizationUnits[1] / sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2] / n0;
				my_B[irho][iz][iphi] = sqrt(sigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
			}
		}
	}
	my_z = my_rho * parameters[3] * normalizationUnits[3];
	my_velocity = parameters[4] * normalizationUnits[4];
	my_meanB = getAverageBsquared();

	for (int i = 0; i < my_maxThreads; ++i) {
		if (my_localDistribution[i] != NULL) {
			delete my_localDistribution[i];
			my_localDistribution[i] = NULL;
		}
	}
}

MassiveParticleDistribution* TabulatedDiskSourceWithSynchCutoff::getParticleDistribution(int irho, int iz, int iphi)
{
	int numthread = omp_get_thread_num();
	if (my_localDistribution[numthread] != NULL) {
		delete my_localDistribution[numthread];
		my_localDistribution[numthread] = NULL;
	}
	my_localDistribution[numthread] = new MassiveParticleTabulatedIsotropicDistribution(*my_cutoffDistribution);
	//my_cutoffDistribution->resetEcut(my_defaultCutoff);
	double rho = (irho + 0.5) * my_rho / my_Nrho;
	double z = (iz + 0.5) * my_z / my_Nz;
	//double l1 = my_rho - rho;
	//double l2 = my_rho - sqrt(rho*rho + z*z);
	double l2 = my_z - z;
	//printf("l2 = %g\n", l2);
	//double l3 = z;
	//double l = min3(l1, l2, l3);
	double l = l2;
	if (l <= 0) {
		/*printf("l <= 0 in TabulatedSLSourceWithSynchCutoff::getParticleDistribution irho = %d iz = %d\n", irho, iz);
		printf("rho = %g z = %g r = %g R = %g\n", rho, z, r, my_rho);
		printLog("l <= 0 in TabulatedSLSourceWithSynchCutoff::getParticleDistribution irho = %d iz = %d\n", irho, iz);
		printLog("rho = %g z = %g r = %g R = %g\n", rho, z, r, my_rho);*/
		//exit(0);
	}
	if (l > 0) {
		double mass = my_cutoffDistribution->getMass();
		double Ecut = mass * mass * mass * mass * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * my_meanB * my_meanB * l);
		if (Ecut < mass * speed_of_light2) {
			//todo
			Ecut = mass * speed_of_light2 * 2.0;
		}
		double gammaCut = Ecut / (me_c2);
		//printf("gammaCut = %g\n", gammaCut);
		//my_cutoffDistribution->resetEcut(Ecut);
		//my_localDistribution->setToZeroAboveE(Ecut);
		//my_localDistribution->addExponentialCutoff(Ecut);
		double lossRate = (4.0/9.0)*electron_charge * electron_charge * electron_charge * electron_charge * my_meanB * my_meanB / (mass * mass * mass * mass * pow(speed_of_light, 7.0));
		double time = l / my_downstreamVelocity;
		my_localDistribution[numthread]->transformToLosses(lossRate, time);
	}
	my_localDistribution[numthread]->resetConcentration(getConcentration(irho, iz, iphi));
	return my_localDistribution[numthread];
}

SectoralSphericalLayerSource::SectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& velocity, const double& redShift) : RadiationSource(Nrho, Nz, Nphi, distance, redShift)
{
	my_rho = rho;
	my_rhoin = rhoin;
	my_minrho = minrho;
	my_phi_sectoral = phi_sectoral;
	my_drho = (my_rho - my_minrho) / my_Nrho;
	my_z = sqrt(my_rho * my_rho - my_minrho * my_minrho);
	my_dz = 2*my_z / my_Nz;
	my_dphi = phi_sectoral / my_Nphi;

	my_velocity = velocity;

	my_geometryCashed = false;
	my_area = new double** [my_Nrho];
	my_length = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_area[i] = new double* [my_Nz];
		my_length[i] = new double* [my_Nz];
		my_v[i] = new double* [my_Nz];
		my_vtheta[i] = new double* [my_Nz];
		my_vphi[i] = new double* [my_Nz];
		for (int j = 0; j < my_Nz; ++j) {
			my_area[i][j] = new double[my_Nphi];
			my_length[i][j] = new double[my_Nphi];
			my_v[i][j] = new double [my_Nphi];
			my_vtheta[i][j] = new double [my_Nphi];
			my_vphi[i][j] = new double [my_Nphi];
			for (int k = 0; k < my_Nphi; ++k) {
				my_area[i][j][k] = 0;
				my_length[i][j][k] = 0;
			}
		}
	}

	evaluateLengthAndArea();


	for (int i = 0; i < my_Nrho; ++i) {
		double r = my_minrho + (i + 0.5) * my_drho;
		for (int j = 0; j < my_Nz; ++j) {
			double z = -my_rho + j * my_dz;
			for (int k = 0; k < my_Nphi; ++k) {
				double phi = (k + 0.5) * my_dphi;
				my_v[i][j][k] = my_velocity;
				my_vtheta[i][j][k] = acos(z / sqr(r * r + z * z));
				my_vphi[i][j][k] = phi;
			}
		}
	}
}

SectoralSphericalLayerSource::SectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : RadiationSource(Nrho, Nz, Nphi, distance, redShift)
{
	my_rho = rho;
	my_rhoin = rhoin;
	my_minrho = minrho;
	my_phi_sectoral = phi_sectoral;
	my_drho = (my_rho - my_minrho) / my_Nrho;
	my_z = sqrt(my_rho * my_rho - my_minrho * my_minrho);
	my_dz = 2 * my_z / my_Nz;
	my_dphi = phi_sectoral / my_Nphi;

	my_geometryCashed = false;
	my_area = new double** [my_Nrho];
	my_length = new double** [my_Nrho];
	my_v = new double** [my_Nrho];
	my_vtheta = new double** [my_Nrho];
	my_vphi = new double** [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_area[i] = new double* [my_Nz];
		my_length[i] = new double* [my_Nz];
		my_v[i] = new double* [my_Nz];
		my_vtheta[i] = new double* [my_Nz];
		my_vphi[i] = new double* [my_Nz];
		for (int j = 0; j < my_Nz; ++j) {
			my_area[i][j] = new double[my_Nphi];
			my_length[i][j] = new double[my_Nphi];
			my_v[i][j] = new double[my_Nphi];
			my_vtheta[i][j] = new double[my_Nphi];
			my_vphi[i][j] = new double[my_Nphi];
			for (int k = 0; k < my_Nphi; ++k) {
				my_area[i][j][k] = 0;
				my_length[i][j][k] = 0;
			}
		}
	}

	evaluateLengthAndArea();


	for (int i = 0; i < my_Nrho; ++i) {
		double r = my_minrho + (i + 0.5) * my_drho;
		for (int j = 0; j < my_Nz; ++j) {
			double z = -my_rho + j * my_dz;
			for (int k = 0; k < my_Nphi; ++k) {
				double phi = (k + 0.5) * my_dphi;
				my_v[i][j][k] = velocity[i][j][k];
				my_vtheta[i][j][k] = vtheta[i][j][k];
				my_vphi[i][j][k] = vphi[i][j][k];
			}
		}
	}
}

SectoralSphericalLayerSource::~SectoralSphericalLayerSource()
{
	for (int i = 0; i < my_Nrho; ++i) {
		for (int j = 0; j < my_Nz; ++j) {
			delete[] my_area[i][j];
			delete[] my_length[i][j];
			delete[] my_v[i][j];
			delete[] my_vtheta[i][j];
			delete[] my_vphi[i][j];
		}
		delete[] my_area[i];
		delete[] my_length[i];
		delete[] my_v[i];
		delete[] my_vtheta[i];
		delete[] my_vphi[i];
	}
	delete[] my_area;
	delete[] my_length;
	delete[] my_v;
	delete[] my_vtheta;
	delete[] my_vphi;
}

double SectoralSphericalLayerSource::getRho(int irho) {
	return my_minrho + (irho + 0.5) * my_drho;
}
double SectoralSphericalLayerSource::getZ(int iz) {
	return -my_z + (iz + 0.5) * my_dz;
}
double SectoralSphericalLayerSource::getPhi(int iphi) {
	return iphi * my_dphi;
}

int SectoralSphericalLayerSource::getRhoIndex(const double& rho) {
	if (rho < my_minrho) {
		printf("rho = %g < my_minrho = %g in get rhoIndex\n", rho, my_minrho);
		printLog("rho = %g < my_minrho = %g in get rhoIndex\n", rho, my_minrho);
		exit(0);
	}
	if (rho > my_rho) {
		printf("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		printLog("rho = %g > my_rho = %g in getRhoIndex\n", rho, my_rho);
		exit(0);
	}
	if (rho == my_rho) {
		return my_Nrho - 1;
	}

	return floor((rho-my_minrho) / my_drho);
}

double SectoralSphericalLayerSource::getMaxRho()
{
	return my_rho;
}

double SectoralSphericalLayerSource::getMinRho() {
	return my_minrho;
}

double SectoralSphericalLayerSource::getRhoin()
{
	return my_rhoin;
}

double SectoralSphericalLayerSource::getMinZ()
{
	return -my_z;
}

double SectoralSphericalLayerSource::getMaxZ()
{
	return my_z;
}

double SectoralSphericalLayerSource::getPhiSectoral()
{
	return my_phi_sectoral;
}

double SectoralSphericalLayerSource::getTotalVolume()
{
	if (my_minrho >= my_rhoin) {
		return (2*my_phi_sectoral / 3) * (my_rho - my_minrho) * (my_rho - my_minrho) * (my_rho - my_minrho);
	}
	else {
		double a = (2 * my_phi_sectoral / 3) * (my_rho - my_minrho) * (my_rho - my_minrho) * (my_rho - my_minrho);
		double b = (2 * my_phi_sectoral / 3) * (my_rhoin - my_minrho) * (my_rhoin - my_minrho) * (my_rhoin - my_minrho);
		return a - b;
	}
}

void SectoralSphericalLayerSource::getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi)
{
	velocity = my_v[irho][iz][iphi];
	if (velocity > speed_of_light) {
		printf("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		printLog("v > c in get velocity %g, irho = %d, iz = %d, iphi = %d\n", velocity / speed_of_light, irho, iz, iphi);
		exit(0);
	}
	theta = my_vtheta[irho][iz][iphi];
	phi = my_vphi[irho][iz][iphi];
}

double SectoralSphericalLayerSource::getArea(int irho, int iz, int iphi)
{
	if (!my_geometryCashed) {
		evaluateLengthAndArea();
	}
	return my_area[irho][iz][iphi];
	//return evaluateArea(irho, iz, iphi);
}

double SectoralSphericalLayerSource::getLength(int irho, int iz, int iphi)
{
	if (!my_geometryCashed) {
		evaluateLengthAndArea();
	}
	return my_length[irho][iz][iphi];
	//return evaluateLength(irho, iz, iphi);
}

double SectoralSphericalLayerSource::getCrossSectionArea(int irho, int iphi)
{
	double rho1 = my_minrho + (irho + 1) * my_drho;
	double rho0 = my_minrho + irho * my_drho;
	return (0.5*my_dphi)*(rho1*rho1 - rho0*rho0);
}

double SectoralSphericalLayerSource::evaluateLength(int irho, int iz, int iphi) {
	double length = 0;
	double r = my_minrho + (irho + 0.5) * my_drho;
	double z1 = -sqrt(my_rho * my_rho - r * r);
	double z2 = 0;
	if (my_rhoin > r) {
		z2 = -sqrt(my_rhoin * my_rhoin - r * r);
	}
	double z3 = -z2;
	double z4 = -z1;
	double zmin = -my_z + iz * my_dz;
	double zmax = zmin + my_dz;
	if (zmin >= 0) {
		if (z3 > zmax) {
			length = 0;
		}
		else
			if (z4 < zmin) {
				length = 0;
			}
			else {
				double lowz = max(zmin, z3);
				double topz = min(zmax, z4);
				length = topz - lowz;
			}
	}
	else {
		if (z1 > zmax) {
			length = 0;
		}
		else if (z2 < zmin) {
			length = 0;
		}
		else {
			double lowz = max(zmin, z1);
			double topz = min(zmax, z2);
			length = topz - lowz;
		}
	}
	return length;
}

double SectoralSphericalLayerSource::evaluateArea(int irho, int iz, int iphi) {
	double area = 0;
	double rho0 = my_minrho + irho * my_drho;
	double rho1 = my_minrho + (irho+1) * my_drho;
	double rmin = rho0;
	double rmax = rho1;
	if (iz >= my_Nz / 2) {
		//upper hemisphere
		double zmax = (iz + 1 - my_Nz / 2) * my_dz;
		double zmin = (iz - my_Nz / 2) * my_dz;
		rmin = rho0;
		if (zmax < my_rhoin) {
			rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
		}
		rmax = rho1;
		rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
	}
	else {
		//lower hemisphere, z inversed
		double zmax = fabs((iz - my_Nz / 2) * (2 * my_z / my_Nz));
		double zmin = fabs((iz + 1 - my_Nz / 2) * (2 * my_z / my_Nz));

		rmin = rho0;
		if (zmax < my_rhoin) {
			rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
		}
		rmax = rho1;
		rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
	}
	if (rmax < rho0 || rmin > rho1 || rmax < rmin) {
		area = 0;
	}
	else {
		area = (rmax * rmax - rmin * rmin) * my_dphi;
	}

	return area;
}

void SectoralSphericalLayerSource::evaluateLengthAndArea()
{
	//length
	for (int irho = 0; irho < my_Nrho; ++irho) {
		double r = my_minrho + (irho + 0.5) * my_drho;
		double z1 = -sqrt(my_rho * my_rho - r * r);
		double z2 = 0;
		if (my_rhoin > r) {
			z2 = -sqrt(my_rhoin * my_rhoin - r * r);
		}
		double z3 = -z2;
		double z4 = -z1;
		for (int iz = 0; iz < my_Nz; ++iz) {
			double zmin = -my_z + iz * my_dz;
			double zmax = zmin + my_dz;
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (zmin >= 0) {
					if (z3 > zmax) {
						my_length[irho][iz][iphi] = 0;
					}
					else
						if (z4 < zmin) {
							my_length[irho][iz][iphi] = 0;
						}
						else {
							double lowz = max(zmin, z3);
							double topz = min(zmax, z4);
							my_length[irho][iz][iphi] = topz - lowz;
							if (my_length[irho][iz][iphi] != my_length[irho][iz][iphi]) {
								printf("my_length = NaN\n");
								printLog("my_length = NaN\n");
								exit(0);
							}
						}
				}
				else {
					if (z1 > zmax) {
						my_length[irho][iz][iphi] = 0;
					}
					else if (z2 < zmin) {
						my_length[irho][iz][iphi] = 0;
					}
					else {
						double lowz = max(zmin, z1);
						double topz = min(zmax, z2);
						my_length[irho][iz][iphi] = topz - lowz;
						if (my_length[irho][iz][iphi] != my_length[irho][iz][iphi]) {
							printf("my_length = NaN\n");
							printLog("my_length = NaN\n");
							exit(0);
						}
					}
				}
			}
		}
	}

	//area
	for (int irho = 0; irho < my_Nrho; ++irho) {
		double rho0 = my_minrho + irho * my_drho;
		double rho1 = my_minrho + (irho + 1) * my_drho;
		double rmin = rho0;
		double rmax = rho1;
		for (int iz = 0; iz < my_Nz; ++iz) {
			if (iz >= my_Nz / 2) {
				//upper hemisphere
				double zmin = -my_z + iz * my_dz;
				double zmax = zmin + my_dz;
				rmin = rho0;
				if (zmax < my_rhoin) {
					rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
				}
				rmax = rho1;
				rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
			}
			else {
				//lower hemisphere, z inversed
				double zmax = fabs(-my_z + iz * my_dz);
				double zmin = zmax - my_dz;

				rmin = rho0;
				if (zmax < my_rhoin) {
					rmin = max(rho0, sqrt(my_rhoin * my_rhoin - zmax * zmax));
				}
				rmax = rho1;
				rmax = min(rho1, sqrt(my_rho * my_rho - zmin * zmin));
			}
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (rmax < rho0 || rmin > rho1 || rmax < rmin) {
					my_area[irho][iz][iphi] = 0;
				}
				else {
					my_area[irho][iz][iphi] = (rmax * rmax - rmin * rmin) * my_dphi;
					if (my_area[irho][iz][iphi] != my_area[irho][iz][iphi]) {
						printf("my_area = NaN\n");
						printLog("my_area = NaN\n");
						exit(0);
					}
				}
			}
		}
	}

	my_geometryCashed = true;
}

bool TabulatedSectoralSphericalLayerSource::rayTraceToNextCell(const double& rho0, const double& z0, int iphi, const double& theta, double& rho1, double& z1, double& lB2)
{
	int irho = floor((rho0 - my_minrho) / my_drho);

	double r = sqrt(rho0 * rho0 + z0 * z0);
	if (r >= my_rho) {
		rho1 = rho0;
		z1 = z0;
		lB2 = 0;
		return true;
	}

	double cosTheta = cos(theta);
	double sinTheta = sin(theta);

	double nextRho = my_minrho + my_drho * (irho + 1);
	double deltaRho = nextRho - rho0;
	if (deltaRho <= 0) {
		irho = irho + 1;
		nextRho = my_minrho + my_drho * (irho + 1);
		deltaRho = nextRho - rho0;
	}

	if (cosTheta > 0) {
		int iz = floor((z0 + my_z) / my_dz);
		double B2 = sqr(getB(irho, iz, iphi));
		double nextZ = my_dz * (iz + 1) - my_z;
		double deltaZ = nextZ - z0;
		if (deltaZ <= 0) {
			nextZ = my_dz * (iz + 2) - my_z;
			deltaZ = my_dz;
			B2 = sqr(getB(irho, iz + 1, iphi));
		}

		double rterm = cosTheta * deltaRho;
		double zterm = sinTheta * deltaZ;

		double l = 0;

		if (rterm > zterm) {
			//sin can be = 0;

			l = deltaZ / cosTheta;
			z1 = nextZ;
			rho1 = rho0 + sqrt(l * l - deltaZ * deltaZ);
		}
		else if (zterm > rterm) {
			//cos can be = 0;
			l = deltaRho / sinTheta;
			rho1 = nextRho;
			z1 = z0 + sqrt(l * l - deltaRho * deltaRho);
		}
		else {
			z1 = nextZ;
			rho1 = nextRho;
			l = sqrt(deltaRho * deltaRho + deltaZ * deltaZ);
		}
		double nextR = sqrt(rho1 * rho1 + z1 * z1);
		if (nextR > my_rho) {
			l = l - (nextR - my_rho);
			if (l < 0) {
				printf("l < 0 in rayTraceToNextCell\n");
				printLog("l < 0 in rayTraceToNextCell\n");
				exit(0);
			}
			//todo change rgo1 and z1
			lB2 = l * B2;
			return true;
		}
		lB2 = l * B2;
		return false;
	} if (cosTheta < 0) {
		int iz = my_Nz / 2 - floor(fabs((z0) / my_dz));
		double B2 = sqr(getB(irho, iz, iphi));
		double nextZ = my_dz * (iz - 1) - my_z;
		double deltaZ = -(nextZ - z0);
		if (deltaZ <= 0) {
			nextZ = my_dz * (iz - 2) - my_z;
			deltaZ = my_dz;
			B2 = sqr(getB(irho, iz - 1, iphi));
		}

		double rterm = -cosTheta * deltaRho;
		double zterm = sinTheta * deltaZ;

		double l = 0;

		if (rterm > zterm) {
			//sin can be = 0;

			l = -deltaZ / cosTheta;
			z1 = nextZ;
			rho1 = rho0 + sqrt(l * l - deltaZ * deltaZ);
		}
		else if (zterm > rterm) {
			//cos can be = 0;
			l = deltaRho / sinTheta;
			rho1 = nextRho;
			z1 = z0 + sqrt(l * l - deltaRho * deltaRho);
		}
		else {
			z1 = nextZ;
			rho1 = nextRho;
			l = sqrt(deltaRho * deltaRho + deltaZ * deltaZ);
		}
		double nextR = sqrt(rho1 * rho1 + z1 * z1);
		if (nextR > my_rho) {
			l = l - (nextR - my_rho);
			if (l < 0) {
				printf("l < 0 in rayTraceToNextCell\n");
				printLog("l < 0 in rayTraceToNextCell\n");
				exit(0);
			}
			//todo change rgo1 and z1
			lB2 = l * B2;
			return true;
		}
		lB2 = l * B2;
		return false;
	}
	else {
		z1 = z0;
		rho1 = rho0;
		double l = deltaRho;
		//todo
		double B2 = sqr(getB(irho, my_Nz, iphi));
		double nextR = sqrt(rho1 * rho1 + z1 * z1);
		if (nextR > my_rho) {
			l = l - (nextR - my_rho);
			if (l < 0) {
				printf("l < 0 in rayTraceToNextCell\n");
				printLog("l < 0 in rayTraceToNextCell\n");
				exit(0);
			}
			//todo change rgo1 and z1
			lB2 = l * B2;
			return true;
		}
		lB2 = l * B2;
		return false;
	}
}

double TabulatedSectoralSphericalLayerSource::evaluateTotalLB2fromPoint(const double& rho0, const double& z0, int iphi, const double& theta)
{
	double r = sqrt(rho0 * rho0 + z0 * z0);
	if (r > my_rho) {
		return 0.0;
	}
	double rho1;
	double z1;
	bool reachedSurface = false;
	double totalLB2 = 0;
	double lB2 = 0;
	double tempRho = rho0;
	double tempZ = z0;
	int numberOfIteration = 0;
	while (!reachedSurface) {
		if (numberOfIteration > 1.5 * (my_Nrho + my_Nz)) {
			printf("aaa\n");
		}
		reachedSurface = rayTraceToNextCell(tempRho, tempZ, iphi, theta, rho1, z1, lB2);
		totalLB2 += lB2;
		tempRho = rho1;
		tempZ = z1;
		numberOfIteration++;
		if (numberOfIteration > 2 * (my_Nrho + my_Nz)) {
			printf("infinite number of iterations in evaluateTotalLB2fromPoint\n");
			printLog("infinite number of iterations in evaluateTotalLB2fromPoint\n");
			exit(0);
		}
	}
	return totalLB2;
}

TabulatedSectoralSphericalLayerSource::TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& velocity, const double& redShift) : SectoralSphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, minrho, phi_sectoral, distance, velocity, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_theta[irho][iz][iphi] = theta;
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSectoralSphericalLayerSource::TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& velocity, const double& redShift) : SectoralSphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, minrho, phi_sectoral, distance, velocity, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			double rho1 = (irho + 0.5) * my_rho / my_Nrho;
			double z = -my_z + (iz + 0.5) * 2 * my_rho / my_Nz;

			double r = sqrt(z * z + rho1 * rho1);

			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				//my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_B[irho][iz][iphi] = B[irho][iz][iphi] * (my_rho / r);
				my_theta[irho][iz][iphi] = theta[irho][iz][iphi];
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				//my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi] * sqr(my_rho / r);
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSectoralSphericalLayerSource::TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : SectoralSphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, minrho, phi_sectoral, distance, velocity, vtheta, vphi, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_theta[irho][iz][iphi] = theta;
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSectoralSphericalLayerSource::TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : SectoralSphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, minrho, phi_sectoral, distance, velocity, vtheta, vphi, redShift)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_theta = new double** [my_Nrho];
	my_phi = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_theta[irho] = new double* [my_Nz];
		my_phi[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_theta[irho][iz] = new double[my_Nphi];
			my_phi[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_theta[irho][iz][iphi] = theta[irho][iz][iphi];
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
			}
		}
	}

	my_isSource = new bool* [my_Nrho];
	for (int i = 0; i < my_Nrho; ++i) {
		my_isSource[i] = new bool[my_Nphi];
		for (int j = 0; j < my_Nphi; ++j) {
			my_isSource[i][j] = true;
		}
	}
}

TabulatedSectoralSphericalLayerSource::~TabulatedSectoralSphericalLayerSource()
{
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_B[irho][iz];
			delete[] my_theta[irho][iz];
		}
		delete[] my_B[irho];
		delete[] my_theta[irho];
		delete[] my_isSource[irho];
	}
	delete[] my_B;
	delete[] my_theta;
	delete[] my_isSource;
}

void TabulatedSectoralSphericalLayerSource::setMask(bool** mask)
{
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			my_isSource[irho][iphi] = mask[irho][iphi];
		}
	}
}

bool TabulatedSectoralSphericalLayerSource::isSource(int irho, int iphi)
{
	return my_isSource[irho][iphi];
}

double TabulatedSectoralSphericalLayerSource::getB(int irho, int iz, int iphi)
{
	return my_B[irho][iz][iphi];
}

double TabulatedSectoralSphericalLayerSource::getMaxB()
{
	double Bmax = 0;
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				if (my_B[irho][iz][iphi] > Bmax) {
					Bmax = my_B[irho][iz][iphi];
				}
			}
		}
	}
	return Bmax;
}

double TabulatedSectoralSphericalLayerSource::getMaxOuterB()
{
	double Bmax = 0;
	int irho = my_Nrho - 1;
	int iz = 0;
	for (iz = 0; iz < my_Nz; ++iz) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}
		}
	}

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iphi = 0; iphi < my_Nphi; ++iphi) {
			iz = 0;
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}

			iz = my_Nz - 1;
			if (my_B[irho][iz][iphi] > Bmax) {
				Bmax = my_B[irho][iz][iphi];
			}
		}
	}

	return Bmax;
}
double TabulatedSectoralSphericalLayerSource::getAverageSigma()
{
	double magneticEnergy = 0;
	double restEnergy = 0;


	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				magneticEnergy += my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * getVolume(irho, iz, iphi) / (4 * pi);
				restEnergy += my_concentration[irho][iz][iphi] * massProton * speed_of_light2 * getVolume(irho, iz, iphi);
			}
		}
	}

	return magneticEnergy / restEnergy;
}
double TabulatedSectoralSphericalLayerSource::getAverageBsquared()
{
	double magneticEnergy = 0;
	double volume = 0;


	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				magneticEnergy += my_B[irho][iz][iphi] * my_B[irho][iz][iphi] * getVolume(irho, iz, iphi);
				volume += getVolume(irho, iz, iphi);
			}
		}
	}

	return sqrt(magneticEnergy / volume);
}
double TabulatedSectoralSphericalLayerSource::getAverageConcentration()
{
	double concentration = 0;
	double volume = 0;

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				concentration += my_concentration[irho][iz][iphi] * getVolume(irho, iz, iphi);
				volume += getVolume(irho, iz, iphi);
			}
		}
	}

	double result = concentration / volume;
	if (result != result) {
		printf("averageConcentration = NaN\n");
		printLog("averageConcentration = NaN\n");
		exit(0);
	}
	return result;
}
double TabulatedSectoralSphericalLayerSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration[irho][iz][iphi];
}

double TabulatedSectoralSphericalLayerSource::getSinTheta(int irho, int iz, int iphi) {
	return sin(my_theta[irho][iz][iphi]);
}

double TabulatedSectoralSphericalLayerSource::getBTheta(int irho, int iz, int iphi)
{
	return my_theta[irho][iz][iphi];
}

double TabulatedSectoralSphericalLayerSource::getBPhi(int irho, int iz, int iphi)
{
	return my_phi[irho][iz][iphi];
}

void TabulatedSectoralSphericalLayerSource::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* B[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* Rin = R*(1 - parameters[3])
*/
	my_geometryCashed = false;

	my_rho = parameters[0] * normalizationUnits[0];
	my_rhoin = my_rho * (1.0 - parameters[3] * normalizationUnits[3]);
	//todo
	my_minrho = my_rho * (1.0 - 2*parameters[3] * normalizationUnits[3]);
	my_drho = (my_rho - my_minrho) / my_Nrho;
	my_z = sqrt(my_rho * my_rho - my_minrho * my_minrho);
	my_dz = 2 * my_z / my_Nz;

	evaluateLengthAndArea();
	double sigma = parameters[1] * normalizationUnits[1];
	//double B0 = my_B[my_Nrho - 1][0][0];
	//double n0 = my_concentration[my_Nrho - 1][0][0];
	double n0 = getAverageConcentration();
	//double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	double sigma0 = getAverageSigma();
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double localSigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				localSigma *= sigma / sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2] / n0;
				if (my_concentration[irho][iz][iphi] != my_concentration[irho][iz][iphi]) {
					printf("my_concentration = NaN\n");
					printLog("my_concentration = NaN\n");
					exit(0);
				}
				my_B[irho][iz][iphi] = sqrt(localSigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
			}
		}
	}
	my_velocity = parameters[4] * normalizationUnits[4];
}
MassiveParticleDistribution* TabulatedSectoralSphericalLayerSource::getParticleDistribution(int irho, int iz, int iphi) {
	my_distribution->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distribution;
}

TabulatedSectoralSLSourceWithSynchCutoff::TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, const double& velocity, const double& redShift) : TabulatedSectoralSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, minrho, phi_sectoral, distance, velocity, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticlePowerLawCutoffDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		exit(0);
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
	my_defaultCutoff = my_cutoffDistribution->getEcutoff();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	double*** Ecutoff = create3dArray(my_Nrho, my_Nz, my_Nphi, 2);
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				Ecutoff[irho][iz][iphi] = (9.0 * massElectron * massElectron * massElectron * massElectron * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * my_LB2[irho][iz][iphi])) / me_c2;
			}
		}
	}
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
	write3dArrayToFile(Ecutoff, my_Nrho, my_Nz, my_Nphi, "Ecutoff.dat");
	delete3dArray(Ecutoff, my_Nrho, my_Nz, my_Nphi);
}

TabulatedSectoralSLSourceWithSynchCutoff::TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, const double& velocity, const double& redShift) : TabulatedSectoralSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, minrho, phi_sectoral, distance, velocity, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticlePowerLawCutoffDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		exit(0);
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
	my_defaultCutoff = my_cutoffDistribution->getEcutoff();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	double*** Ecutoff = create3dArray(my_Nrho, my_Nz, my_Nphi, 2);
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				Ecutoff[irho][iz][iphi] = (9.0 * massElectron * massElectron * massElectron * massElectron * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * my_LB2[irho][iz][iphi]))/me_c2;
			}
		}
	}
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
	write3dArrayToFile(Ecutoff, my_Nrho, my_Nz, my_Nphi, "Ecutoff.dat");
	delete3dArray(Ecutoff, my_Nrho, my_Nz, my_Nphi);
}

TabulatedSectoralSLSourceWithSynchCutoff::TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedSectoralSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, minrho, phi_sectoral, distance, velocity, vtheta, vphi, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticlePowerLawCutoffDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		exit(0);
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
	my_defaultCutoff = my_cutoffDistribution->getEcutoff();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	double*** Ecutoff = create3dArray(my_Nrho, my_Nz, my_Nphi, 2);
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				Ecutoff[irho][iz][iphi] = (9.0 * massElectron * massElectron * massElectron * massElectron * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * my_LB2[irho][iz][iphi])) / me_c2;
			}
		}
	}
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
	write3dArrayToFile(Ecutoff, my_Nrho, my_Nz, my_Nphi, "Ecutoff.dat");
	delete3dArray(Ecutoff, my_Nrho, my_Nz, my_Nphi);
}

TabulatedSectoralSLSourceWithSynchCutoff::TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift) : TabulatedSectoralSphericalLayerSource(Nrho, Nz, Nphi, electronDistribution, B, theta, phi, concentration, rho, rhoin, minrho, phi_sectoral, distance, velocity, vtheta, vphi, redShift)
{
	my_cutoffDistribution = dynamic_cast<MassiveParticlePowerLawCutoffDistribution*>(electronDistribution);
	if (my_cutoffDistribution == NULL) {
		printf("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		printLog("distribution in TabulatedSLSourceWithSynchCutoff must be only MassiveParticlePowerLawCutoffDistribution\n");
		exit(0);
	}
	my_downstreamVelocity = downstreamVelocity;
	my_meanB = getAverageBsquared();
	my_defaultCutoff = my_cutoffDistribution->getEcutoff();

	my_LB2 = create3dArray(my_Nrho, my_Nz, my_Nphi);
	updateLB2();
	double*** Ecutoff = create3dArray(my_Nrho, my_Nz, my_Nphi, 2);
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				Ecutoff[irho][iz][iphi] = (9.0 * massElectron * massElectron * massElectron * massElectron * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * my_LB2[irho][iz][iphi])) / me_c2;
			}
		}
	}
	write3dArrayToFile(my_LB2, my_Nrho, my_Nz, my_Nphi, "LB2.dat");
	write3dArrayToFile(Ecutoff, my_Nrho, my_Nz, my_Nphi, "Ecutoff.dat");
	delete3dArray(Ecutoff, my_Nrho, my_Nz, my_Nphi);
}

TabulatedSectoralSLSourceWithSynchCutoff::~TabulatedSectoralSLSourceWithSynchCutoff()
{
	delete3dArray(my_LB2, my_Nrho, my_Nz, my_Nphi);
}

void TabulatedSectoralSLSourceWithSynchCutoff::updateLB2() {
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double rho = my_minrho + (irho + 0.5) * my_drho;
				double z = (iz + 0.5) * my_dz - my_z;
				double theta = acos(z / sqrt(rho * rho + z * z));
				my_LB2[irho][iz][iphi] = evaluateTotalLB2fromPoint(rho, z, iphi, theta);
			}
		}
	}
}

void TabulatedSectoralSLSourceWithSynchCutoff::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* B[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* Rin = R*(1 - parameters[3])
*/
	my_geometryCashed = false;

	my_rho = parameters[0] * normalizationUnits[0];
	my_rhoin = my_rho * (1.0 - parameters[3] * normalizationUnits[3]);
	//todo
	my_minrho = my_rho * (1.0 - 2 * parameters[3] * normalizationUnits[3]);
	my_drho = (my_rho - my_minrho) / my_Nrho;
	my_z = sqrt(my_rho * my_rho - my_minrho * my_minrho);
	my_dz = 2 * my_z / my_Nz;
	evaluateLengthAndArea();
	double sigma = parameters[1] * normalizationUnits[1];
	//double B0 = my_B[my_Nrho - 1][0][0];
	//double n0 = my_concentration[my_Nrho - 1][0][0];
	double n0 = getAverageConcentration();
	//double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	double sigma0 = getAverageSigma();
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double localSigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
				localSigma *= sigma / sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2] / n0;
				if (my_concentration[irho][iz][iphi] != my_concentration[irho][iz][iphi]) {
					printf("my_concentration = NaN\n");
					printLog("my_concentration = NaN\n");
					exit(0);
				}
				my_B[irho][iz][iphi] = sqrt(localSigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
			}
		}
	}
	my_velocity = parameters[4] * normalizationUnits[4];
	my_meanB = getAverageBsquared();
	updateLB2();
}

MassiveParticleDistribution* TabulatedSectoralSLSourceWithSynchCutoff::getParticleDistribution(int irho, int iz, int iphi)
{
	my_cutoffDistribution->resetEcut(my_defaultCutoff);
	double LB2 = my_LB2[irho][iz][iphi];
	if (LB2 <= 0) {
		/*printf("l <= 0 in TabulatedSLSourceWithSynchCutoff::getParticleDistribution irho = %d iz = %d\n", irho, iz);
		printf("rho = %g z = %g r = %g R = %g\n", rho, z, r, my_rho);
		printLog("l <= 0 in TabulatedSLSourceWithSynchCutoff::getParticleDistribution irho = %d iz = %d\n", irho, iz);
		printLog("rho = %g z = %g r = %g R = %g\n", rho, z, r, my_rho);*/
		//exit(0);
	}
	if (LB2 > 0) {
		double mass = my_cutoffDistribution->getMass();
		double Ecut = 9.0 * mass * mass * mass * mass * pow(speed_of_light, 7) * my_downstreamVelocity / (electron_charge * electron_charge * electron_charge * electron_charge * LB2);
		if (Ecut < mass * speed_of_light2) {
			//todo
			Ecut = mass * speed_of_light2 * 2.0;
		}
		my_cutoffDistribution->resetEcut(Ecut);
	}
	my_cutoffDistribution->resetConcentration(getConcentration(irho, iz, iphi));
	return my_cutoffDistribution;
}
