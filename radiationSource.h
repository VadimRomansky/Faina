#ifndef radiation_source_h
#define radiation_source_h

#include "electronDistribution.h"

//sources have different shapes in cylindrycal geometry. Z axis to observer

class RadiationSource {
protected:
	int my_Nrho;
	int my_Nz;
	int my_Nphi;
	double my_distance;
public:
	RadiationSource(int Nrho, int Nz, int Nphi, double distance);

	virtual double getMaxRho()=0;
	virtual double getMinZ()=0;
	virtual double getMaxZ()=0;
	double getNrho();
	double getNz();
	double getNphi();
	double getDistance();

	virtual double getTotalVolume()=0;
	virtual double getLength(int irho, int iz, int iphi) = 0;
	virtual double getArea(int irho);
	virtual ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi)=0;
};

class SimpleFlatSource : public RadiationSource {
protected:
	double my_rho;
	double my_z;
	ElectronIsotropicDistribution* my_distribution;
public:
	SimpleFlatSource(ElectronIsotropicDistribution* electronDistribution, const double& rho, const double& z, const double& distance);
	double getMaxRho();
	double getMinZ();
	double getMaxZ();
	double getTotalVolume();
	double getLength(int irho, int iz, int iphi);
	ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi);
};

class FlatDiskSource : public RadiationSource {
protected:
	double my_rho;
	double my_z;
	ElectronIsotropicDistribution* my_distribution;
public:
	FlatDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& rho, const double& z, const double& distance);
	double getMaxRho();
	double getMinZ();
	double getMaxZ();
	double getTotalVolume();
	double getLength(int irho, int iz, int iphi);
	ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi);
};


#endif
