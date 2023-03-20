#ifndef radiation_source_h
#define radiation_source_h

#include "massiveParticleDistribution.h"

//sources have different shapes in cylindrycal geometry. Z axis to observer

class RadiationSource {
protected:
	int my_Nrho;
	int my_Nz;
	int my_Nphi;
	double my_distance;
public:
	RadiationSource(int Nrho, int Nz, int Nphi, double distance);
    virtual ~RadiationSource(){

    }

	virtual double getMaxRho()=0;
	virtual double getMinZ()=0;
	virtual double getMaxZ()=0;
	int getNrho();
	int getNz();
	int getNphi();
	double getDistance();
	virtual double getArea(int irho, int iz, int iphi)=0;
	virtual double getVolume(int irho, int iz, int iphi);

	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getConcentration(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual double getTotalVolume()=0;
	virtual double getLength(int irho, int iz, int iphi) = 0;
	//virtual void resetConcentration(const double& concentration) = 0;
	virtual void resetParameters(const double* parameters, const double* normalizationUnits)=0;
	virtual MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi)=0;
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

	virtual double getArea(int irho, int iz, int iphi);
	virtual double getLength(int irho, int iz, int iphi) = 0;
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi) = 0;
};

class SimpleFlatSource : public DiskSource {
protected:
	double my_B;
	double my_sinTheta;
	double my_concentration;
	MassiveParticleIsotropicDistribution* my_distribution;
public:
	SimpleFlatSource(MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& rho, const double& z, const double& distance);
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getConcentration(int irho, int iz, int iphi);
	double getSinTheta(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class TabulatedDiskSource : public DiskSource {
protected:
	double*** my_B;
	double*** my_sinTheta;
	double*** my_concentration;
	MassiveParticleIsotropicDistribution* my_distribution;
public:
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, double*** concentration, const double& rho, const double& z, const double& distance);
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& sinTheta, const double& concentration , const double& rho, const double& z, const double& distance);
    virtual ~TabulatedDiskSource();
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getConcentration(int irho, int iz, int iphi);
	double getSinTheta(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi);
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

	virtual double getLength(int irho, int iz, int iphi);
	virtual double getArea(int irho, int iz, int iphi);
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi) = 0;
};

class TabulatedSphericalLayerSource : public SphericalLayerSource {
protected:
	double*** my_B;
	double*** my_sinTheta;
	double*** my_concentration;
	MassiveParticleIsotropicDistribution* my_distribution;
public:
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, double*** B, double*** sinTheta, double*** concentration, const double& rho, const double& rhoin, const double& distance);
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& concentration, const double& sinTheta, const double& rho, const double& rhoin, const double& distance);
    virtual ~TabulatedSphericalLayerSource();

	//virtual double getLength(int irho, int iz, int iphi);
	virtual double getB(int irho, int iz, int iphi);
	double getConcentration(int irho, int iz, int iphi);
	virtual double getSinTheta(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class AngleDependentElectronsSphericalSource : public TabulatedSphericalLayerSource {
protected:
	int my_Ntheta;
	MassiveParticleIsotropicDistribution** my_distributions;
	double*** my_phi;
	double*** my_shockWaveAngle;
public:
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleIsotropicDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance);
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleIsotropicDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance);
    virtual ~AngleDependentElectronsSphericalSource();

	//void resetConcentration(const double& concentration);
	virtual MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class RadiationTimeDependentSource {
protected:
	double my_t0;
	RadiationSource* my_radiationSource;
public:
	RadiationTimeDependentSource(RadiationSource* source, const double& t0) {
		my_radiationSource = source;
		my_t0 = t0;
	}
    virtual ~RadiationTimeDependentSource(){

    }
	//note that number of parameters and they sence are on your responsibility
	virtual void resetParameters(const double* parameters, const double* normalizationUnits) = 0;
	virtual RadiationSource* getRadiationSource(double& time, const double* normalizationUnits) = 0;
};

class ExpandingRemnantSource : public RadiationTimeDependentSource {
protected:
	double my_R0;
	double my_B0;
	double my_concentration0;
	double my_v;
	double my_widthFraction;
public:
	ExpandingRemnantSource(const double& R0, const double& B0, const double& concentration0, const double& v, const double& widthFraction, RadiationSource* source, const double& t0);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual RadiationSource* getRadiationSource(double& time, const double* normalizationUnits);
};


#endif
