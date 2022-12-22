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

SimpleFlatSource::SimpleFlatSource(ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& z, const double& distance) : DiskSource(1,1,1, rho, z,distance) {
	my_distribution = electronDistribution;
	my_B = B;
	my_sinTheta = sinTheta;
}
double SimpleFlatSource::getB(int irho, int iz, int iphi)
{
	return my_B;
}
double SimpleFlatSource::getSinTheta(int irho, int iz, int iphi)
{
	return my_sinTheta;
}
double SimpleFlatSource::getLength(int irho, int iz, int iphi) {
	return my_z;
}
ElectronIsotropicDistribution* SimpleFlatSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, const double& rho, const double& z, const double& distance) : DiskSource(Nrho, Nz, Nphi, rho, z, distance) {
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_sinTheta[irho][iz][iphi] = sinTheta[irho][iz][iphi];
			}
		}
	}
}

TabulatedDiskSource::TabulatedDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& z, const double& distance) : DiskSource(Nrho, Nz, Nphi, rho, z, distance)
{
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_sinTheta[irho][iz][iphi] = sinTheta;
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
		}
		delete[] my_B[irho];
		delete[] my_sinTheta[irho];
	}
	delete[] my_B;
	delete[] my_sinTheta;
}
double TabulatedDiskSource::getB(int irho, int iz, int iphi)
{
	return my_B[irho][iz][iphi];
}
double TabulatedDiskSource::getSinTheta(int irho, int iz, int iphi)
{
	return my_sinTheta[irho][iz][iphi];
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

TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, const double& rho, const double& rhoin, const double& distance) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance) {
	my_distribution = electronDistribution;
	
	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B[irho][iz][iphi];
				my_sinTheta[irho][iz][iphi] = sinTheta[irho][iz][iphi];
			}
		}
	}
}
TabulatedSphericalLayerSource::TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& rhoin, const double& distance) : SphericalLayerSource(Nrho, Nz, Nphi, rho, rhoin, distance) {
	my_distribution = electronDistribution;

	my_B = new double** [my_Nrho];
	my_sinTheta = new double** [my_Nrho];
	for (int irho = 0; irho < my_Nrho; ++irho) {
		my_B[irho] = new double* [my_Nz];
		my_sinTheta[irho] = new double* [my_Nz];
		for (int iz = 0; iz < my_Nz; ++iz) {
			my_B[irho][iz] = new double[my_Nphi];
			my_sinTheta[irho][iz] = new double[my_Nphi];
			for (int iphi = 0; iphi < my_Nphi; ++iphi) {
				my_B[irho][iz][iphi] = B;
				my_sinTheta[irho][iz][iphi] = sinTheta;
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

}
double TabulatedSphericalLayerSource::getB(int irho, int iz, int iphi) {
	return my_B[irho][iz][iphi];
}
double TabulatedSphericalLayerSource::getSinTheta(int irho, int iz, int iphi) {
	return my_sinTheta[irho][iz][iphi];
}
ElectronIsotropicDistribution* TabulatedSphericalLayerSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}