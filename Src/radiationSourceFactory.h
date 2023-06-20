#ifndef RADIATION_SOURCE_FACTORY
#define RADIATION_SOURCE_FACTORY

#include "radiation.h"
#include "radiationSource.h"

class RadiationSourceFactory {
private:
	const static int minModeNumber = 10;
	static double evaluateTurbulenceAmplitude(const double& k, const double& turbulenceKoef, const double& index, const double& L0);
	static double evaluateAnisotropicTurbulenceAmplitude(const double& kx, const double& ky, const double& kz, const double& turbulenceKoef, const double& index, const double& L0, const double& anisotropy);
	static void normalizeTurbulenceKoef(double& turbulenceKoef, const double& index, const double& L0, const int Nmodes, const double& fraction, const double& B0);
	static void normalizeAnisotropicTurbulenceKoef(double& turbulenceKoef, const double& index, const double& L0, const int Nmodes, const double& fraction, const double& B0, const double& anisotropy);
public:
	static void sumFields(double*** B, double*** theta, double*** phi, double*** B1, double*** theta1, double*** phi1, double*** B2, double*** theta2, double*** phi2, int Nrho, int Nz, int Nphi);

	static void initializeTurbulentField(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R);
	static void initializeAnisotropicLocalTurbulentFieldInSphericalSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& anisotropy);
	static void initializeAnisotropicLocalTurbulentFieldInDiskSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& anisotropy);

	static void initializeParkerField(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R);
	static void initializeParkerFieldWithRotation(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R, const double& thetaRot);

	static AngleDependentElectronsSphericalSource* createSourceWithTurbulentField(MassiveParticleIsotropicDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double n0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& rho, const double& rhoin, const double& distance);

	static AngleDependentElectronsSphericalSource* createSourceWithParkerField(MassiveParticleIsotropicDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& rho, const double& rhoin, const double& distance);
	static AngleDependentElectronsSphericalSource* createSourceWithParkerFieldWithRotation(MassiveParticleIsotropicDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& thetaRot, const double& rho, const double& rhoin, const double& distance);
};

#endif
