#ifndef radiation_source_h
#define radiation_source_h

#include "massiveParticleDistribution.h"

//sources have different shapes in cylindrycal geometry. Z axis to observer

class RadiationSource {
protected:
	int my_Nx1;
	int my_Nx2;
	int my_Nz;

	double my_distance;
	double my_redShift;
public:

	RadiationSource(int Nx1, int Nz, int Nx2, double distance, double redShift);
	virtual ~RadiationSource() {

	}

	virtual double getMinZ() = 0;
	virtual double getMaxZ() = 0;
	virtual double getMaxB() = 0;
	virtual double getMaxOuterB() = 0;
	virtual double getAverageSigma() = 0;
	virtual double getAverageConcentration() = 0;
	virtual double getAverageBsquared() = 0;

	int getNz();
	int getNx1();
	int getNx2();

	virtual double getX1(int ix1) = 0;
	virtual double getZ(int iz) = 0;
	virtual double getX2(int ix2) = 0;

	virtual bool isSource(int irho, int iphi) = 0;
	virtual double getArea(int irho, int iz, int iphi) = 0;
	virtual double getVolume(int irho, int iz, int iphi);
	virtual double getCrossSectionArea(int irhi, int iphi) = 0;
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi) = 0;

	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getConcentration(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual double getBTheta(int irho, int iz, int iphi) = 0;
	virtual double getBPhi(int irho, int iz, int iphi) = 0;
	virtual double getTotalVolume() = 0;
	virtual double getLength(int irho, int iz, int iphi) = 0;
	//virtual void resetConcentration(const double& concentration) = 0;
	virtual void resetParameters(const double* parameters, const double* normalizationUnits) = 0;
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi) = 0;

	double getDistance();
	double getRedShift();
};

class RadiationSourceInCartesian : public RadiationSource {
protected:
	double my_minX;
	double my_maxX;
	double my_minY;
	double my_maxY;
	double my_minZ;
	double my_maxZ;

	int my_Nx;
	int my_Ny;
public:
	RadiationSourceInCartesian(int Nx, int Ny, int Nz, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, double distance, double redShift);
	virtual ~RadiationSourceInCartesian() {

	}

	virtual double getMinX();
	virtual double getMaxX();
	virtual double getMinY();
	virtual double getMaxY();
	virtual double getMinZ();
	virtual double getMaxZ();
	int getNx();
	int getNy();
	virtual double getX(int ix) = 0;
	virtual double getZ(int iz) = 0;
	virtual double getY(int iy) = 0;
	virtual double getX1(int ix1);
	virtual double getX2(int ix2);
	virtual int gerXindex(double x) = 0;
	virtual int getYindex(double y) = 0;
	virtual int getZindex(double z) = 0;
};

class RectangularSource : public RadiationSourceInCartesian {
protected:
	double* my_xgrid;
	double* my_ygrid;
	double* my_zgrid;

	double*** my_B;
	double*** my_theta;
	double*** my_phi;
	double*** my_concentration;
	double my_velocity;

	double*** my_v;
	double*** my_vtheta;
	double*** my_vphi;

	bool** my_isSource;
	MassiveParticleDistribution* my_distribution;
public:
	RectangularSource(int Nx, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double B, double theta, double phi, double concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	RectangularSource(int Nx, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	RectangularSource(int Nx, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double B, double theta, double phi, double concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	RectangularSource(int Nx, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	RectangularSource(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double B, double theta, double phi, double concentration, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	RectangularSource(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	RectangularSource(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double B, double theta, double phi, double concentration, double minY, double maxY, double minZ, double maxZ, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	RectangularSource(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, double minY, double maxY, double minZ, double maxZ, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	virtual double getX(int ix);
	virtual double getZ(int iz);
	virtual double getY(int iy);
	virtual int gerXindex(double x);
	virtual int getYindex(double y);
	virtual int getZindex(double z);

	virtual double getMaxB();
	virtual double getMaxOuterB();
	virtual double getAverageSigma();
	virtual double getAverageConcentration();
	virtual double getAverageBsquared();

	virtual bool isSource(int irho, int iphi);
	virtual double getArea(int irho, int iz, int iphi) ;
	virtual double getCrossSectionArea(int irhi, int iphi);
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);

	virtual double getB(int irho, int iz, int iphi);
	virtual double getConcentration(int irho, int iz, int iphi);
	virtual double getSinTheta(int irho, int iz, int iphi);
	virtual double getBTheta(int irho, int iz, int iphi);
	virtual double getBPhi(int irho, int iz, int iphi);
	virtual double getTotalVolume();
	virtual double getLength(int irho, int iz, int iphi);
	//virtual void resetConcentration(const double& concentration) = 0;
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

//be careful uses a lot of memory
class RectangularSourceInhomogenousDistribution : public RectangularSource {
protected:
	MassiveParticleDistribution**** my_distributions;
public:
	RectangularSourceInhomogenousDistribution(int Nx, int Ny, int Nz, MassiveParticleDistribution**** electronDistributions, double B, double theta, double phi, double concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	RectangularSourceInhomogenousDistribution(int Nx, int Ny, int Nz, MassiveParticleDistribution**** electronDistributions, double B, double theta, double phi, double*** concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	RectangularSourceInhomogenousDistribution(int Nx, int Ny, int Nz, MassiveParticleDistribution**** electronDistributions, double*** B, double*** theta, double*** phi, double*** concentration, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);

	RectangularSourceInhomogenousDistribution(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution**** electronDistributions, double*** B, double*** theta, double*** phi, double*** concentration, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);

	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class ThermalRectangularSource : public RectangularSource {
protected:
	double my_mass;
	double*** my_temperature;
	MassiveParticleDistribution** my_localDistribution;
	int my_maxThreads;
public:
	ThermalRectangularSource(int Nx, int Ny, int Nz, double mass, double B, double theta, double phi, double concentration, double temperature, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);
	ThermalRectangularSource(int Nx, int Ny, int Nz, double mass, double*** B, double*** theta, double*** phi, double*** concentration, double*** temperature, double minX, double maxX, double minY, double maxY, double minZ, double maxZ, const double& distance, const double& velocity = 0, const double& redShift = 0);

	virtual double getParticleMass();
	virtual double getTemperature(int ix, int iz, int iy);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class RectangularSourceWithSynchAndComptCutoffFromRight : public RectangularSource {
protected:
	int my_maxThreads;
	double*** my_downstreamVelocity;
	double my_meanB;
	double my_photonEnergyDensity;
	//MassiveParticlePowerLawCutoffDistribution* my_cutoffDistribution;
	MassiveParticleTabulatedIsotropicDistribution* my_cutoffDistribution;
	MassiveParticleTabulatedIsotropicDistribution** my_localDistribution;
	double*** my_LB2;
	void updateLB2();
public:
	RectangularSourceWithSynchAndComptCutoffFromRight(int Nx, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& minX, const double& maxX, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& distance, const double& downstreamVelocity, const double& photonEnergyDensity = 0, const double& velocity = 0, const double& redShift = 0);
	RectangularSourceWithSynchAndComptCutoffFromRight(int Nx, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& minX, const double& maxX, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& distance, const double& downstreamVelocity, const double& photonEnergyDensity = 0, const double& velocity = 0, const double& redShift = 0);
	//sRectangularSourceWithSynchAndComptCutoffFromRight(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& distance, const double& downstreamVelocity, const double& photonEnergyDensity = 0, const double& velocity = 0, const double& redShift = 0);
	RectangularSourceWithSynchAndComptCutoffFromRight(int Nx, double* xgrid, int Ny, int Nz, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& distance, const double& downstreamVelocity1, const double& downstreamVelocity2, const double& photonEnergyDensity = 0, const double& velocity = 0, const double& redShift = 0);
	virtual ~RectangularSourceWithSynchAndComptCutoffFromRight();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};



class RadiationSourceInCylindrical : public RadiationSource{
protected:
	int my_Nrho;
	int my_Nphi;

	double evaluateAverageVelocity();
public:
	RadiationSourceInCylindrical(int Nrho, int Nz, int Nphi, double distance, double redShift);
    virtual ~RadiationSourceInCylindrical(){

    }

	virtual double getMaxRho()=0;
	virtual double getMinRho()=0;
	
	int getNrho();
	int getNphi();
	virtual double getRho(int irho) = 0;
	virtual double getZ(int iz) = 0;
	virtual double getPhi(int iphi) = 0;
	virtual double getX1(int ix1);
	virtual double getX2(int ix2);
	virtual int getRhoIndex(const double& rho) = 0;
};

class DiskSource : public RadiationSourceInCylindrical {
protected:
	double my_rho;
	double my_z;
public:
	DiskSource(int Nrho, int Nz, int Nphi, const double& rho, const double& z, const double& distance, const double& redShift);
	double getMaxRho();
	double getMinRho();
	double getMinZ();
	double getMaxZ();
	double getTotalVolume();

	virtual double getArea(int irho, int iz, int iphi);
	virtual double getLength(int irho, int iz, int iphi) = 0;
	virtual double getCrossSectionArea(int irho, int iphi);
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual double getBTheta(int irho, int iz, int iphi) = 0;
	virtual double getBPhi(int irho, int iz, int iphi) = 0;
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi) = 0;
};

class SimpleFlatSource : public DiskSource {
protected:
	double my_B;
	double my_theta;
	double my_phi;
	double my_concentration;
	double my_velocity;
	MassiveParticleDistribution* my_distribution;
public:
	SimpleFlatSource(MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& velocity = 0, const double& redShift = 0);
	virtual double getRho(int irho);
	virtual double getZ(int iz);
	virtual double getPhi(int iphi);

	virtual int getRhoIndex(const double& rho);
	
	virtual bool isSource(int irho, int iphi);
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getMaxB();
	double getMaxOuterB();
	double getAverageSigma();
	double getAverageConcentration();
	double getAverageBsquared();
	double getConcentration(int irho, int iz, int iphi);
	void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
	double getSinTheta(int irho, int iz, int iphi);
	double getBTheta(int irho, int iz, int iphi);
	double getBPhi(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class SimpleFlatSource2 : public SimpleFlatSource {
	protected:
		int my_maxThreads;
		int my_Ndistributions;
		double* my_velocities;
		MassiveParticleDistribution** my_distributions;
		MassiveParticleDistribution** my_outputDistributions;
	public:
		SimpleFlatSource2(int Ndistributions, double* velocities, MassiveParticleIsotropicDistribution** electronDistributions, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& velocity = 0, const double& redShift = 0);
		virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
		virtual void resetParameters(const double* parameters, const double* normalizationUnits);
};

class TabulatedDiskSource : public DiskSource {
protected:
	double*** my_B;
	double*** my_theta;
	double*** my_phi;
	double*** my_concentration;
	double my_velocity;

	double*** my_v;
	double*** my_vtheta;
	double*** my_vphi;

	bool** my_isSource;
	MassiveParticleDistribution* my_distribution;
public:
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, const double& velocity = 0, const double& redShift = 0);
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration , const double& rho, const double& z, const double& distance, const double& velocity = 0, const double& redShift = 0);
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
    virtual ~TabulatedDiskSource();

	virtual double getRho(int irho);
	virtual double getZ(int iz);
	virtual double getPhi(int iphi);
	virtual int getRhoIndex(const double& rho);
	virtual void setMask(bool** mask);
	virtual bool isSource(int irho, int iphi);
	virtual double getMaxB();
	virtual double getMaxOuterB();
	virtual double getAverageSigma();
	virtual double getAverageConcentration();
	virtual double getAverageBsquared();
	virtual double getLength(int irho, int iz, int iphi);
	virtual double getB(int irho, int iz, int iphi);
	virtual double getConcentration(int irho, int iz, int iphi);
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
	virtual double getSinTheta(int irho, int iz, int iphi);
	virtual double getBTheta(int irho, int iz, int iphi);
	virtual double getBPhi(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class TabulatedDiskSourceWithSynchAndComptCutoff : public TabulatedDiskSource {
protected:
	int my_maxThreads;
	double my_downstreamVelocity;
	double my_meanB;
	double my_photonEnergyDensity;
	//MassiveParticlePowerLawCutoffDistribution* my_cutoffDistribution;
	MassiveParticleTabulatedIsotropicDistribution* my_cutoffDistribution;
	MassiveParticleTabulatedIsotropicDistribution** my_localDistribution;
	double*** my_LB2;
	void updateLB2();
public:
	TabulatedDiskSourceWithSynchAndComptCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, const double& photonEnergyDensity = 0, const double& velocity = 0, const double& redShift = 0);
	TabulatedDiskSourceWithSynchAndComptCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, const double& photonEnergyDensity = 0, const double& velocity = 0, const double& redShift = 0);
	TabulatedDiskSourceWithSynchAndComptCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& photonEnergyDensity = 0, const double& redShift = 0);
	TabulatedDiskSourceWithSynchAndComptCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& z, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& photonEnergyDensity = 0, const double& redShift = 0);
	virtual ~TabulatedDiskSourceWithSynchAndComptCutoff();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class SphericalLayerSource : public RadiationSourceInCylindrical {
protected:
	double my_rho;
	double my_rhoin;

	double*** my_area;
	double*** my_length;
	bool my_geometryCashed;
	double my_velocity;

	double*** my_v;
	double*** my_vtheta;
	double*** my_vphi;

	virtual double evaluateLength(int irho, int iz, int iphi);
	virtual double evaluateArea(int irho, int iz, int iphi);
	virtual void evaluateLengthAndArea();
public:
	SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0, const double& redShift = 0);
	SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	virtual ~SphericalLayerSource();

	virtual double getRho(int irho);
	virtual double getZ(int iz);
	virtual double getPhi(int iphi);
	virtual int getRhoIndex(const double& rho);
	double getMaxRho();
	double getMinRho();
	double getMinZ();
	double getMaxZ();
	double getTotalVolume();
	double getInnerRho();

	virtual double getLength(int irho, int iz, int iphi);
	virtual double getArea(int irho, int iz, int iphi);
	virtual double getCrossSectionArea(int irho, int iphi);
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual double getBTheta(int irho, int iz, int iphi) = 0;
	virtual double getBPhi(int irho, int iz, int iphi) = 0;
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi) = 0;
};

class TabulatedSphericalLayerSource : public SphericalLayerSource {
protected:
	double*** my_B;
	double*** my_theta;
	double*** my_phi;
	double*** my_concentration;
	bool** my_isSource;
	MassiveParticleDistribution* my_distribution;

	bool rayTraceToNextCell(const double& rho0, const double& z0, int iphi, const double& theta, double& rho1, double& z1, double& lB2);
	double evaluateTotalLB2fromPoint(const double& rho0, const double& z0, int iphi, const double& theta);
public:
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0, const double& redShift = 0);
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0, const double& redShift = 0);
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
    virtual ~TabulatedSphericalLayerSource();

	//virtual double getLength(int irho, int iz, int iphi);
	virtual void setMask(bool** mask);
	virtual bool isSource(int irho, int iphi);
	virtual double getB(int irho, int iz, int iphi);
	virtual double getMaxB();
	virtual double getMaxOuterB();
	virtual double getAverageSigma();
	virtual double getAverageBsquared();
	virtual double getAverageConcentration();
	double getConcentration(int irho, int iz, int iphi);
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
	virtual double getSinTheta(int irho, int iz, int iphi);
	virtual double getBTheta(int irho, int iz, int iphi);
	virtual double getBPhi(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class TabulatedSphericalLayerSource2 : public TabulatedSphericalLayerSource {
protected:
	int my_Ndistributions;
	int my_maxThreads;
	double* my_velocities;
	MassiveParticleDistribution** my_distributions;
	MassiveParticleDistribution** my_outputDistributions;

public:
	TabulatedSphericalLayerSource2(int Ndistributions, double* velocities, MassiveParticleIsotropicDistribution** distributions, int Nrho, int Nz, int Nphi, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0, const double& redShift = 0);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
};

class AngleDependentElectronsSphericalSource : public TabulatedSphericalLayerSource {
protected:
	int my_Ntheta;
	MassiveParticleDistribution** my_distributions;
	double*** my_shockWaveAngle;
public:
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0, const double& redShift = 0);
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0, const double& redShift = 0);
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
    virtual ~AngleDependentElectronsSphericalSource();

	//void resetConcentration(const double& concentration);
	double getShockWaveAngle(int irho, int iz, int iphi);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class TabulatedSLSourceWithSynchCutoff : public TabulatedSphericalLayerSource {
protected:
	int my_maxThreads;
	double my_downstreamVelocity;
	double my_meanB;
	double my_defaultCutoff;
	MassiveParticleTabulatedIsotropicDistribution* my_cutoffDistribution;
	MassiveParticleTabulatedIsotropicDistribution** my_localDistribution;
	double*** my_LB2;
	void updateLB2();
public:
	TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, const double& velocity = 0, const double& redShift = 0);
	TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, const double& velocity = 0, const double& redShift = 0);
	TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	TabulatedSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	virtual ~TabulatedSLSourceWithSynchCutoff();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class SectoralSphericalLayerSource : public RadiationSourceInCylindrical {
protected:
	double my_rho;
	double my_rhoin;
	double my_minrho;
	double my_phi_sectoral;

	double my_z;

	double my_drho;
	double my_dz;
	double my_dphi;

	double*** my_area;
	double*** my_length;
	bool my_geometryCashed;
	double my_velocity;
	
	double*** my_v;
	double*** my_vtheta;
	double*** my_vphi;

	virtual double evaluateLength(int irho, int iz, int iphi);
	virtual double evaluateArea(int irho, int iz, int iphi);
	virtual void evaluateLengthAndArea();
public:
	SectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& velocity = 0, const double& redShift = 0);
	SectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	~SectoralSphericalLayerSource();

	virtual double getRho(int irho);
	virtual double getZ(int iz);
	virtual double getPhi(int iphi);
	virtual int getRhoIndex(const double& rho);
	double getMaxRho();
	double getMinRho();
	double getRhoin();
	double getMinZ();
	double getMaxZ();
	double getPhiSectoral();
	double getTotalVolume();

	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
	virtual double getArea(int irho, int iz, int iphi);
	virtual double getLength(int irho, int iz, int iphi);
	virtual double getCrossSectionArea(int irho, int iphi);
	virtual double getB(int irho, int iz, int iphi) = 0;
	virtual double getSinTheta(int irho, int iz, int iphi) = 0;
	virtual double getBTheta(int irho, int iz, int iphi) = 0;
	virtual double getBPhi(int irho, int iz, int iphi) = 0;
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi) = 0;

};

class TabulatedSectoralSphericalLayerSource : public SectoralSphericalLayerSource {
protected:
	double*** my_B;
	double*** my_theta;
	double*** my_phi;
	double*** my_concentration;
	double my_velocity;
	bool** my_isSource;
	MassiveParticleDistribution* my_distribution;

	bool rayTraceToNextCell(const double& rho0, const double& z0, int iphi, const double& theta, double& rho1, double& z1, double& lB2);
	double evaluateTotalLB2fromPoint(const double& rho0, const double& z0, int iphi, const double& theta);
public:
	TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& velocity = 0, const double& redShift = 0);
	TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& velocity = 0, const double& redShift = 0);
	TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	TabulatedSectoralSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	virtual ~TabulatedSectoralSphericalLayerSource();

	//virtual double getLength(int irho, int iz, int iphi);
	virtual void setMask(bool** mask);
	virtual bool isSource(int irho, int iphi);
	virtual double getB(int irho, int iz, int iphi);
	virtual double getMaxB();
	virtual double getMaxOuterB();
	virtual double getAverageSigma();
	virtual double getAverageBsquared();
	virtual double getAverageConcentration();
	double getConcentration(int irho, int iz, int iphi);
	virtual double getSinTheta(int irho, int iz, int iphi);
	virtual double getBTheta(int irho, int iz, int iphi);
	virtual double getBPhi(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class TabulatedSectoralSLSourceWithSynchCutoff : public TabulatedSectoralSphericalLayerSource {
protected:
	double my_downstreamVelocity;
	double my_meanB;
	double my_defaultCutoff;
	MassiveParticlePowerLawCutoffDistribution* my_cutoffDistribution;
	double*** my_LB2;
	void updateLB2();
public:
	TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, const double& velocity = 0, const double& redShift = 0);
	TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, const double& velocity = 0, const double& redShift = 0);
	TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, double*** B, double*** theta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	TabulatedSectoralSLSourceWithSynchCutoff(int Nrho, int Nz, int Nphi, MassiveParticleDistribution* electronDistribution, const double& B, const double& theta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& minrho, const double& phi_sectoral, const double& distance, const double& downstreamVelocity, double*** velocity, double*** vtheta, double*** vphi, const double& redShift = 0);
	virtual ~TabulatedSectoralSLSourceWithSynchCutoff();

	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual MassiveParticleDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class RadiationTimeDependentSource {
protected:
	double my_t0;
	RadiationSourceInCylindrical* my_radiationSource;
public:
	RadiationTimeDependentSource(RadiationSourceInCylindrical* source, const double& t0) {
		my_radiationSource = source;
		my_t0 = t0;
	}
    virtual ~RadiationTimeDependentSource(){

    }
	//note that number of parameters and they sence are on your responsibility
	virtual void resetParameters(const double* parameters, const double* normalizationUnits) = 0;
	virtual RadiationSourceInCylindrical* getRadiationSource(double& time, const double* normalizationUnits) = 0;
};

class ExpandingRemnantSource : public RadiationTimeDependentSource {
protected:
	double my_R0;
	double my_B0;
	double my_concentration0;
	double my_v;
	double my_widthFraction;

	double my_velocityPower;
	double my_concentrationPower;
	double my_Bpower;
	double my_widthPower;
public:
	ExpandingRemnantSource(const double& R0, const double& B0, const double& concentration0, const double& v, const double& widthFraction, RadiationSourceInCylindrical* source, const double& t0, const double& velocityPower = 1.0, const double& Bpower = 1.0, const double& concentrationPower = 2.0, const double& widthPower = 1.0);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual RadiationSourceInCylindrical* getRadiationSource(double& time, const double* normalizationUnits);
};


#endif
