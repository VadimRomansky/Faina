#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "massiveParticleDistribution.h"

#include "radiationSource.h"

RadiationSource::RadiationSource(int Nrho, int Nz, int Nphi, double distance) {
	my_Nrho = Nrho;
	my_Nz = Nz;
	my_Nphi = Nphi;
	my_distance = distance;
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

double RadiationSource::getVolume(int irho, int iz, int iphi)
{
	return getArea(irho, iz, iphi)*getLength(irho, iz, iphi);
}

DiskSource::DiskSource(int Nrho, int Nz, int Nphi, const double& rho, const double& z, const double& distance) : RadiationSource(Nrho, Nz, Nphi, distance)
{
	my_rho = rho;
	my_z = z;
}

double DiskSource::getArea(int irho, int iz, int iphi) {
	double rho0 = irho * getMaxRho() / my_Nrho;
	double rho1 = (irho + 1) * getMaxRho() / my_Nrho;
	return 2 * pi * (rho1 * rho1 - rho0 * rho0) / my_Nphi;
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

SimpleFlatSource::SimpleFlatSource(MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& z, const double& distance) : DiskSource(1,1,1, rho, z,distance) {
	my_distribution = electronDistribution;
	my_B = B;
	my_sinTheta = sinTheta;
	my_concentration = my_distribution->getConcentration();
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
	* sigma = parameters[1]
	* n = parameters[2]
	* z = R*parameters[3]
	*/
	double sigma = parameters[1] * normalizationUnits[1];
	my_rho = parameters[0] * normalizationUnits[0];
	my_concentration = parameters[2] * normalizationUnits[2];
	my_z = my_rho * parameters[3] * normalizationUnits[3];
	my_B = sqrt(sigma * 4 * pi * massProton * my_concentration * speed_of_light2);
}
double SimpleFlatSource::getLength(int irho, int iz, int iphi) {
	return my_z;
}
MassiveParticleIsotropicDistribution* SimpleFlatSource::getParticleDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, double*** concentration, const double& rho, const double& z, const double& distance) : DiskSource(Nrho, Nz, Nphi, rho, z, distance) {
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

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& concentration, const double& rho, const double& z, const double& distance) : DiskSource(Nrho, Nz, Nphi, rho, z, distance)
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
* sigma[i][j][k] ~ parameters[1]
* n[i][j][k] ~ parameters[2]
* where B[Nrho-1][0][0] = parameters[1]
* z = R*parameters[3]
*/
	my_rho = parameters[0] * normalizationUnits[0];
	double sigma = parameters[1] * normalizationUnits[1];
	double B0 = my_B[my_Nrho - 1][0][0];
	double n0 = my_concentration[my_Nrho - 1][0][0];
	double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
	for (int irho = 0; irho < my_Nrho; ++irho) {
		for (int iz = 0; iz < my_Nz; ++iz) {
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				double sigma = sqr(my_B[irho][iz][iphi]) / (4 * pi * massProton * my_concentration[irho][iz][iphi]*speed_of_light2);
				sigma *= parameters[1] * normalizationUnits[1]/sigma0;
				my_concentration[irho][iz][iphi] *= parameters[2] * normalizationUnits[2]/n0;
				my_B[irho][iz][iphi] = sqrt(sigma * 4 * pi * massProton * my_concentration[irho][iz][iphi] * speed_of_light2);
			}
		}
	}
	my_z = my_rho * parameters[3] * normalizationUnits[3];
}
double TabulatedDiskSource::getLength(int irho, int iz, int iphi) {
	return my_z/my_Nz;
}
MassiveParticleIsotropicDistribution* TabulatedDiskSource::getParticleDistribution(int irho, int iz, int iphi) {
	my_distribution->resetConcentration(getConcentration(irho, iz, iphi));
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

double SphericalLayerSource::getLength(int irho, int iz, int iphi) {
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
		rmax = min(rho1, sqrt(my_rho*my_rho - zmin*zmin));
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

	return 2 * pi * (rmax * rmax - rmin*rmin) / my_Nphi;
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

TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, double*** concentration, const double& rho, const double& rhoin, const double& distance) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance) {
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
TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& concentration, const double& rho, const double& rhoin, const double& distance) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance) {
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
	double sigma = parameters[1] * normalizationUnits[1];
	double B0 = my_B[my_Nrho - 1][0][0];
	double n0 = my_concentration[my_Nrho - 1][0][0];
	double sigma0 = sqr(my_B[my_Nrho - 1][0][0]) / (4 * pi * massProton * my_concentration[my_Nrho - 1][0][0] * speed_of_light2);
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
	my_rhoin = my_rho * (1.0 - parameters[3] * normalizationUnits[3]);
}
MassiveParticleIsotropicDistribution* TabulatedSphericalLayerSource::getParticleDistribution(int irho, int iz, int iphi) {
	my_distribution->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distribution;
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleIsotropicDistribution** electronDistributions, double*** B, double*** sinTheta, double*** concentration, double*** phi, const double& rho, const double& rhoin, const double& distance) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, concentration, rho, rhoin, distance)
{
	my_Ntheta = Ntheta;
	my_distributions = new MassiveParticleIsotropicDistribution*[my_Ntheta];

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
                    printLog("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
				exit(0);
				}

				my_shockWaveAngle[irho][iz][iphi] = acos(cosTheta);
			}
		}
	}
}

AngleDependentElectronsSphericalSource::AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleIsotropicDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance) : TabulatedSphericalLayerSource(Nrho, Nz, Nphi, NULL, B, sinTheta, concentration, rho, rhoin, distance)
{
	my_Ntheta = Ntheta;
	my_distributions = new MassiveParticleIsotropicDistribution * [my_Ntheta];

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
                    printLog("cos theta in angle depended distribution = %g > 1.0\n", cosTheta);
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

MassiveParticleIsotropicDistribution* AngleDependentElectronsSphericalSource::getParticleDistribution(int irho, int iz, int iphi)
{
	double dtheta = pi / (my_Ntheta-1);
	double theta = my_shockWaveAngle[irho][iz][iphi];
	if (theta > pi / 2) {
		theta = pi - theta;
	}

	int angleIndex = floor((theta + 0.5*dtheta) / dtheta);
	if (angleIndex == my_Ntheta) {
		angleIndex = my_Ntheta - 1;
	}

	my_distributions[angleIndex]->resetConcentration(getConcentration(irho, iz, iphi));
	return my_distributions[angleIndex];
}

ExpandingRemnantSource::ExpandingRemnantSource(const double& R0, const double& B0, const double& concentration0, const double& v, const double& widthFraction, RadiationSource* source, const double& t0) : RadiationTimeDependentSource(source, t0) {
	my_R0 = R0;
	my_B0 = B0;
	my_concentration0 = concentration0;
	my_widthFraction = widthFraction;
	my_v = v;
}

void ExpandingRemnantSource::resetParameters(const double* parameters, const double* normalizationUnits) {
	my_R0 = parameters[0] * normalizationUnits[0];
	double sigma = parameters[1] * normalizationUnits[1];
	my_concentration0 = parameters[2] * normalizationUnits[2];
	my_widthFraction = parameters[3] * normalizationUnits[3];
	my_v = parameters[4] * normalizationUnits[4];
	my_B0 = sqrt(sigma * 4 * pi * massProton * my_concentration0 * speed_of_light2);
}

//just one possible example
RadiationSource* ExpandingRemnantSource::getRadiationSource(double& time, const double* normalizationUnits) {
	double R = my_R0 + my_v * (time - my_t0);
	//double R = my_R0 + (5.0/2.0) * my_v * my_t0 * (pow(time/my_t0, 2.0/5.0) - 1.0);
	double sigma = sqr(my_B0) / (4 * pi * massProton * my_concentration0 * speed_of_light2)/(my_R0/R);
	//double B = my_B0;
	//double n = my_concentration0 * sqr(my_R0 / R);
	//double n = my_concentration0;
	double n = my_concentration0*cube(my_R0/R);
	double fracton = my_widthFraction * my_R0/R;
	//double fracton = my_widthFraction;

	double parameters[4];
	parameters[0] = R / normalizationUnits[0];
	parameters[1] = sigma / normalizationUnits[1];
	parameters[2] = n / normalizationUnits[2];
	parameters[3] = fracton / normalizationUnits[3];
	my_radiationSource->resetParameters(parameters, normalizationUnits);
	return my_radiationSource;
}
