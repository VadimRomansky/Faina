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
	return 2 * pi * (rho0 * rho0 - rho1 * rho1) / my_Nphi;
}

SimpleFlatSource::SimpleFlatSource(ElectronIsotropicDistribution* electronDistribution, const double& rho, const double& z, const double& distance) : RadiationSource(1,1,1,distance) {
	my_distribution = electronDistribution;
	my_rho = rho;
	my_z = z;
}

double SimpleFlatSource::getMaxRho() {
	return my_rho;
}
double SimpleFlatSource::getMinZ() {
	return 0;
}
double SimpleFlatSource::getMaxZ() {
	return my_z;
}
double SimpleFlatSource::getLength(int irho, int iz, int iphi) {
	return my_z;
}
double SimpleFlatSource::getTotalVolume() {
	return 2 * pi * my_z * my_rho * my_rho;
}
ElectronIsotropicDistribution* SimpleFlatSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}

FlatDiskSource::FlatDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& rho, const double& z, const double& distance) : RadiationSource(Nrho, Nz, Nphi, distance) {
	my_distribution = electronDistribution;
	my_rho = rho;
	my_z = z;
}

double FlatDiskSource::getMaxRho() {
	return my_rho;
}
double FlatDiskSource::getMinZ() {
	return 0;
}
double FlatDiskSource::getMaxZ() {
	return my_z;
}
double FlatDiskSource::getLength(int irho, int iz, int iphi) {
	return my_z/my_Nz;
}
double FlatDiskSource::getTotalVolume() {
	return 2 * pi * my_z * my_rho * my_rho;
}
ElectronIsotropicDistribution* FlatDiskSource::getElectronDistribution(int irho, int iz, int iphi) {
	return my_distribution;
}