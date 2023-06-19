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
	virtual double getMaxB() = 0;
	virtual double getMaxOuterB() = 0;
	virtual double getAverageSigma() = 0;
	virtual double getAverageConcentration() = 0;
	int getNrho();
	int getNz();
	int getNphi();
	double getDistance();
	virtual double getArea(int irho, int iz, int iphi)=0;
	virtual double getVolume(int irho, int iz, int iphi);
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi) = 0;

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
	double my_theta;
	double my_concentration;
	double my_velocity;
	MassiveParticleIsotropicDistribution* my_distribution;
public:
	SimpleFlatSource(MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& theta, const double& rho, const double& z, const double& distance, const double& velocity = 0);
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getMaxB();
	double getMaxOuterB();
	double getAverageSigma();
	double getAverageConcentration();
	double getConcentration(int irho, int iz, int iphi);
	void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
	double getSinTheta(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class TabulatedDiskSource : public DiskSource {
protected:
	double*** my_B;
	double*** my_theta;
	double*** my_concentration;
	double my_velocity;
	MassiveParticleIsotropicDistribution* my_distribution;
public:
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, double*** B, double*** theta, double*** concentration, const double& rho, const double& z, const double& distance, const double& velocity = 0);
	TabulatedDiskSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& theta, const double& concentration , const double& rho, const double& z, const double& distance, const double& velocity = 0);
    virtual ~TabulatedDiskSource();
	double getLength(int irho, int iz, int iphi);
	double getB(int irho, int iz, int iphi);
	double getConcentration(int irho, int iz, int iphi);
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
	double getSinTheta(int irho, int iz, int iphi);
	//void resetConcentration(const double& concentration);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	MassiveParticleIsotropicDistribution* getParticleDistribution(int irho, int iz, int iphi);
};

class SphericalLayerSource : public RadiationSource {
protected:
	double my_rho;
	double my_rhoin;

	double*** my_area;
	double*** my_length;
	bool my_geometryCashed;

	virtual double evaluateLength(int irho, int iz, int iphi);
	virtual double evaluateArea(int irho, int iz, int iphi);
	virtual void evaluateLengthAndArea();
public:
	SphericalLayerSource(int Nrho, int Nz, int Nphi, const double& rho, const double& rhoin, const double& distance);
	virtual ~SphericalLayerSource();
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
	double*** my_theta;
	double*** my_concentration;
	double my_velocity;
	MassiveParticleIsotropicDistribution* my_distribution;
public:
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, double*** B, double*** theta, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0);
	TabulatedSphericalLayerSource(int Nrho, int Nz, int Nphi, MassiveParticleIsotropicDistribution* electronDistribution, const double& B, const double& concentration, const double& theta, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0);
    virtual ~TabulatedSphericalLayerSource();

	//virtual double getLength(int irho, int iz, int iphi);
	virtual double getB(int irho, int iz, int iphi);
	virtual double getMaxB();
	virtual double getMaxOuterB();
	virtual double getAverageSigma();
	virtual double getAverageConcentration();
	double getConcentration(int irho, int iz, int iphi);
	virtual void getVelocity(int irho, int iz, int iphi, double& velocity, double& theta, double& phi);
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
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleIsotropicDistribution** electronDistributions, double*** B, double*** sinTheta, double*** phi, double*** concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0);
	AngleDependentElectronsSphericalSource(int Nrho, int Nz, int Nphi, int Ntheta, MassiveParticleIsotropicDistribution** electronDistributions, const double& B, const double& sinTheta, const double& phi, const double& concentration, const double& rho, const double& rhoin, const double& distance, const double& velocity = 0);
    virtual ~AngleDependentElectronsSphericalSource();

	//void resetConcentration(const double& concentration);
	double getShockWaveAngle(int irho, int iz, int iphi);
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

	double my_radiusPower;
	double my_concentrationPower;
	double my_Bpower;
public:
	ExpandingRemnantSource(const double& R0, const double& B0, const double& concentration0, const double& v, const double& widthFraction, RadiationSource* source, const double& t0);
	virtual void resetParameters(const double* parameters, const double* normalizationUnits);
	virtual RadiationSource* getRadiationSource(double& time, const double* normalizationUnits);
};

class RadiationSourceFactory {
private:
	const static int minModeNumber = 5;
	static double evaluateTurbulenceAmplitude(const double& k, const double& turbulenceKoef, const double& index, const double& L0);
	static double evaluateAnisotropicTurbulenceAmplitude(const double& kx, const double& ky, const double& kz, const double& turbulenceKoef, const double& index, const double& L0, const double& thetaN, const double& phiN, const double& anisotropy);
	static void normalizeTurbulenceKoef(double& turbulenceKoef, const double& index, const double& L0, const int Nmodes, const double& fraction, const double& B0);
	static void normalizeAnisotropicTurbulenceKoef(double& turbulenceKoef, const double& index, const double& L0, const int Nmodes, const double& fraction, const double& B0, const double& thetaN, const double& phiN, const double& anisotropy);
public:
	static void sumFields(double*** B, double*** theta, double*** phi, double*** B1, double*** theta1, double*** phi1, double*** B2, double*** theta2, double*** phi2, int Nrho, int Nz, int Nphi);

	static void initializeTurbulentField(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R);
	static void initializeAnisotropicTurbulentField(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& thetaN, const double& phiN, const double& anisotropy);

	static void initializeParkerField(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R);
	static void initializeParkerFieldWithRotation(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R, const double& thetaRot);

	static AngleDependentElectronsSphericalSource* createSourceWithTurbulentField(MassiveParticleIsotropicDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double n0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& rho, const double& rhoin, const double& distance);

	static AngleDependentElectronsSphericalSource* createSourceWithParkerField(MassiveParticleIsotropicDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& rho, const double& rhoin, const double& distance);
	static AngleDependentElectronsSphericalSource* createSourceWithParkerFieldWithRotation(MassiveParticleIsotropicDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& thetaRot, const double& rho, const double& rhoin, const double& distance);
};


#endif
