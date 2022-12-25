#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "electronDistribution.h"

#include "radiationSource.h"

RadiationSource::RadiationSource(int Nrho, int Nz, int Nphi, double distance) {
	my_Nrho = Nrho;
	my_Nz = Nz;
	my_Nphi = Nphi;
	my_distance = distance;
}

double RadiationSource::getNrho() {
	return my_Nrho;
}
double RadiationSource::getNz() {
	return my_Nz;
}
double RadiationSource::getNphi() {
	return my_Nphi;
}
double RadiationSource::getDistance() {
	return my_distance;
}
double RadiationSource::getArea(int irho) {
	double rho0 = irho*getMaxRho()/my_Nrho;
	double rho1 = (irho+1)*getMaxRho()/my_Nrho;
	return 2 * pi * (rho1 * rho1 - rho0 * rho0) / my_Nphi;
}

double RadiationSource::getVolume(int irho, int iz, int iphi)
{
	return getArea(irho)*getLength(irho, iz, iphi);
}

DiskSource::DiskSource(int Nrho, int Nz, int Nphi, const double& rho, const double& z, const double& distance) : RadiationSource(Nrho, Nz, Nphi, distance)
{
	my_rho = rho;
	my_z = z;
}

double DiskSource::getTotalVolume()
{
	return 2 * pi * my_z * my_rho * my_rho;;
}
double DiskSource::getMaxRho() {
	return my_rho;
}
double DiskSource::getMinZ() {
	return 0;
}
double DiskSource::getMaxZ() {
	return my_z;
}

SimpleFlatSource::SimpleFlatSource(ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& concentration, const double& rho, const double& z, const double& distance) : DiskSource(1,1,1, rho, z,distance) {
	my_distribution = electronDistribution;
	my_B = B;
	my_sinTheta = sinTheta;
	my_concentration = concentration;
}
double SimpleFlatSource::getB(int irho, int iz, int iphi)
{
	return my_B;
}
double SimpleFlatSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration;
}
double SimpleFlatSource::getSinTheta(int irho, int iz, int iphi)
{
	return my_sinTheta;
}
/*void SimpleFlatSource::resetConcentration(const double& concentration)
{
	my_distribution->resetConcentration(concentration);
}*/
void SimpleFlatSource::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
	* R = parameters[0]
	* B = parameters[1]
	* n = parameters[2]
	* z = R*parameters[3]
	*/
	my_rho = parameters[0] * normalizationUnits[0];
	my_B = parameters[1] * normalizationUnits[1];
	my_concentration = parameters[2] * normalizationUnits[2];
	my_z = my_rho * parameters[3] * normalizationUnits[3];
}
double SimpleFlatSource::getLength(int irho, int iz, int iphi) {
	return my_z;
}
ElectronIsotropicDistribution* SimpleFlatSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, double*** concentration, const double& rho, const double& z, const double& distance) : DiskSource(Nrho, Nz, Nphi, rho, z, distance) {
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_sinTheta[irho][iz][iphi] = sinTheta[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
			}
		}
	}
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& concentration, const double& rho, const double& z, const double& distance) : DiskSource(Nrho, Nz, Nphi, rho, z, distance)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_sinTheta[irho][iz][iphi] = sinTheta;
				my_concentration[irho][iz][iphi] = concentration;
			}
		}
	}
}
TabulatedDiskSource::~TabulatedDiskSource()
{
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_B[irho][iz];
			delete[] my_sinTheta[irho][iz];
			delete[] my_concentration[irho][iz];
		}
		delete[] my_B[irho];
		delete[] my_sinTheta[irho];
		delete[] my_concentration[irho];
	}
	delete[] my_B;
	delete[] my_sinTheta;
	delete[] my_concentration;
}
double TabulatedDiskSource::getB(int irho, int iz, int iphi)
{
	return my_B[irho][iz][iphi];
}
double TabulatedDiskSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration[irho][iz][iphi];
}
double TabulatedDiskSource::getSinTheta(int irho, int iz, int iphi)
{
	return my_sinTheta[irho][iz][iphi];
}
/*void TabulatedDiskSource::resetConcentration(const double& concentration)
{
	my_distribution->resetConcentration(concentration);
}*/
void TabulatedDiskSource::resetParameters(const double* parameters, const double* normalizationUnits)
{
	/*parameters must be
* R = parameters[0]
* B[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* z = R*parameters[3]
*/
	my_rho = parameters[0] * normalizationUnits[0];
	double B = my_B[my_Nrho - 1][0][0];
	double n = my_concentration[my_Nrho - 1][0][0];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] *= parameters[1] * normalizationUnits[1]/B;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2]/n;
			}
		}
	}
	my_z = my_rho * parameters[3] * normalizationUnits[3];
}
double TabulatedDiskSource::getLength(int irho, int iz, int iphi) {
	return my_z/my_Nz;
}
ElectronIsotropicDistribution* TabulatedDiskSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}

SphericalLayerSource::SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance) : RadiationSource(Nrho, Nz, Nphi, distance) {
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
}
double SphericalLayerSource::getMaxRho() {
	return my_rho;
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

TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, double*** concentration, const double& rho, const double& rhoin, const double& distance) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance) {
	my_distribution = electronDistribution;
	
	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_sinTheta[irho][iz][iphi] = sinTheta[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
			}
		}
	}
}
TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& concentration, const double& rho, const double& rhoin, const double& distance) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance) {
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_sinTheta[irho][iz][iphi] = sinTheta;
				my_concentration[irho][iz][iphi] = concentration;
			}
		}
	}
}

TabulatedSphericalLayerSource::~TabulatedSphericalLayerSource() {
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_B[irho][iz];
			delete[] my_sinTheta[irho][iz];
		}
		delete[] my_B[irho];
		delete[] my_sinTheta[irho];
	}
	delete[] my_B;
	delete[] my_sinTheta;
}

double TabulatedSphericalLayerSource::getLength(int irho, int iz, int iphi) {
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
double TabulatedSphericalLayerSource::getB(int irho, int iz, int iphi) {
	return my_B[irho][iz][iphi];
}
double TabulatedSphericalLayerSource::getConcentration(int irho, int iz, int iphi)
{
	return my_concentration[irho][iz][iphi];
}
double TabulatedSphericalLayerSource::getSinTheta(int irho, int iz, int iphi) {
	return my_sinTheta[irho][iz][iphi];
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
	my_rho = parameters[0] * normalizationUnits[0];
	double B = my_B[my_Nrho - 1][0][0];
	double n = my_concentration[my_Nrho - 1][0][0];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] *= parameters[1] * normalizationUnits[1] / B;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2] / n;
			}
		}
	}
	my_rhoin = my_rho * (1.0 - parameters[3] * normalizationUnits[3]);
}
ElectronIsotropicDistribution* TabulatedSphericalLayerSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, ElectronIsotropicDistribution** electronDistributions, double*** B, double*** sinTheta, double*** concentration, double*** phi, const double& rho, const double& rhoin, const double& distance) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, concentration, rho, rhoin, distance)
{
	my_Ntheta = Ntheta;
	my_distributions = new ElectronIsotropicDistribution*[my_Ntheta];

	for (int i = 0; i < my_Ntheta; ++i) {
		my_distributions[i] = electronDistributions[i];
	}

	my_phi = new double** [my_Nrho];
	my_shockWaveAngle = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_phi[irho] = new double* [my_Nz];
		my_shockWaveAngle[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_phi[irho][iz] = new double[my_Nphi];
			my_shockWaveAngle[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_phi[irho][iz][iphi] = phi[irho][iz][iphi];
				my_concentration[irho][iz][iphi] = concentration[irho][iz][iphi];
				double rho = (irho + 0.5) * my_rho / my_Nrho;
				double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
				double alpha = (iphi + 0.5) * 2 * pi / my_Nphi;

				double r = sqrt(z * z + rho * rho);
				double nz = z / r;
				double nx = rho * cos(alpha) / r;
				double ny = rho * sin(alpha) / r;

				double Bz = my_B[irho][iz][iphi] * sqrt(1.0 - sqr(my_sinTheta[irho][iz][iphi]));
				double Bx = my_B[irho][iz][iphi] * my_sinTheta[irho][iz][iphi] * cos(my_phi[irho][iz][iphi]);
				double By = my_B[irho][iz][iphi] * my_sinTheta[irho][iz][iphi] * sin(my_phi[irho][iz][iphi]);

				double cosTheta = (nz * Bz + nx * Bx + ny * By) / my_B[irho][iz][iphi];

				if (cosTheta > 1.0 || cosTheta < -1.0) {
					printf("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					printLog("cos theta in angle depended distribution = > 1.0\n");
				exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
			}
		}
	}
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, ElectronIsotropicDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, concentration, rho, rhoin, distance)
{
	my_Ntheta = Ntheta;
	my_distributions = new ElectronIsotropicDistribution * [my_Ntheta];

	for (int i = 0; i < my_Ntheta; ++i) {
		my_distributions[i] = electronDistributions[i];
	}

	my_phi = new double** [my_Nrho];
	my_shockWaveAngle = new double** [my_Nrho];
	my_concentration = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_phi[irho] = new double* [my_Nz];
		my_shockWaveAngle[irho] = new double* [my_Nz];
		my_concentration[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_phi[irho][iz] = new double[my_Nphi];
			my_shockWaveAngle[irho][iz] = new double[my_Nphi];
			my_concentration[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_phi[irho][iz][iphi] = phi;
				my_concentration[irho][iz][iphi] = concentration;
				double rho = (irho + 0.5) * my_rho / my_Nrho;
				double z = -my_rho + (iz + 0.5) * 2 * my_rho / my_Nz;
				double alpha = (iphi + 0.5) * 2 * pi / my_Nphi;

				double r = sqrt(z * z + rho * rho);
				double nz = z / r;
				double nx = rho * cos(alpha) / r;
				double ny = rho * sin(alpha) / r;

				double Bz = my_B[irho][iz][iphi] * sqrt(1.0 - sqr(my_sinTheta[irho][iz][iphi]));
				double Bx = my_B[irho][iz][iphi] * my_sinTheta[irho][iz][iphi] * cos(my_phi[irho][iz][iphi]);
				double By = my_B[irho][iz][iphi] * my_sinTheta[irho][iz][iphi] * sin(my_phi[irho][iz][iphi]);

				double cosTheta = (nz * Bz + nx * Bx + ny * By) / my_B[irho][iz][iphi];

				if (cosTheta > 1.0 || cosTheta < -1.0) {
					printf("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
					printLog("cos theta in angle depended distribution = > 1.0\n");
					exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
			}
		}
	}
}

AngleDependentElectronsSphericalSource::~AngleDependentElectronsSphericalSource()
{
	delete[] my_distributions;

	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			delete[] my_phi[irho][iz];
			delete[] my_shockWaveAngle[irho][iz];
			delete[] my_concentration[irho][iz];
		}
		delete[] my_phi[irho];
		delete[] my_shockWaveAngle[irho];
		delete[] my_concentration[irho];
	}
	delete[] my_phi;
	delete[] my_shockWaveAngle;
	delete[] my_concentration;
}

/*void AngleDependentElectronsSphericalSource::resetConcentration(const double& concentration)
{
	for(int i = 0; i < my_Ntheta; ++i){
		my_distributions[i]->resetConcentration(concentration);
}
}*/

ElectronIsotropicDistribution* AngleDependentElectronsSphericalSource::getElectronDistribution(int irho, int iz, int iphi)
{
	double dtheta = pi / my_Ntheta;
	double theta = my_shockWaveAngle[irho][iz][iphi];
	if (theta > pi / 2) {
		theta = pi - theta;
	}

	int angleIndex = floor(pi / dtheta);
	if (angleIndex == my_Ntheta) {
		angleIndex = my_Ntheta - 1;
	}

	return my_distributions[angleIndex];
}
