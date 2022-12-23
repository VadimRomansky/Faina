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
	virtual double getArea(int irho);
	virtual double getVolume(int irho, int iz, int iphi);

	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual double getTotalVolume()=0;
	virtual double getLength(int irho, int iz, int iphi) = 0;
	virtual ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi)=0;
};

class DiskSource : public RadiationSource {
protected:
	double my_rho;
	double my_z;
public:
	DiskSource(int Nrho, int Nz, int Nphi, const double& rho, const double& z, const double& distance);
	double getMaxRho();
	double getMinZ();
	double getMaxZ();
	double getTotalVolume();

	virtual double getLength(int irho, int iz, int iphi) = 0;
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi) = 0;
};

class SimpleFlatSource : public DiskSource {
protected:
	double my_B;
	double my_sinTheta;
	ElectronIsotropicDistribution* my_distribution;
public:
	SimpleFlatSource(ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& z, const double& distance);
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getSinTheta(int irho, int iz, int iphi);
	ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi);
};

class TabulatedDiskSource : public DiskSource {
protected:
	double*** my_B;
	double*** my_sinTheta;
	ElectronIsotropicDistribution* my_distribution;
public:
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta , const double& rho, const double& z, const double& distance);
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta , const double& rho, const double& z, const double& distance);
	~TabulatedDiskSource();
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getSinTheta(int irho, int iz, int iphi);
	ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi);
};

class SphericalLayerSource : public RadiationSource {
protected:
	double my_rho;
	double my_rhoin;
public:
	SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance);
	double getMaxRho();
	double getMinZ();
	double getMaxZ();
	double getTotalVolume();
	double getInnerRho();

	virtual double getLength(int irho, int iz, int iphi) = 0;
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi) = 0;
};

class TabulatedSphericalLayerSource : public SphericalLayerSource {
protected:
	double*** my_B;
	double*** my_sinTheta;
	ElectronIsotropicDistribution* my_distribution;
public:
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, const double& rho, const double& rhoin, const double& distance);
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, ElectronIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& rhoin, const double& distance);
	~TabulatedSphericalLayerSource();

	virtual double getLength(int irho, int iz, int iphi);
	virtual double getB(int irho, int iz, int iphi);
	virtual double getSinTheta(int irho, int iz, int iphi);
	virtual ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi);
};

class AngleDependedElectronsSphericalSource : public TabulatedSphericalLayerSource {
protected:
	int my_Ntheta;
	ElectronIsotropicDistribution** my_distributions;
	double*** my_phi;
	double*** my_shockWaveAngle;
public:
	AngleDependedElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, ElectronIsotropicDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, const double& rho, const double& rhoin, const double& distance);
	AngleDependedElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, ElectronIsotropicDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& rho, const double& rhoin, const double& distance);
	~AngleDependedElectronsSphericalSource();

	virtual ElectronIsotropicDistribution* getElectronDistribution(int irho, int iz, int iphi);
};


#endif
