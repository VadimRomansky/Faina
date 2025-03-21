#ifndef RADIATION_SOURCE_FACTORY
#define RADIATION_SOURCE_FACTORY

#include "radiation.h"
#include "radiationSource.h"

enum SourceInputGeometry { CARTESIAN, CYLINDRICAL, SPHERICAL };

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
	static void initializeAnisotropicLocalTurbulentFieldInSectoralSphericalSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& Rmin, const double& phiR, const double& anisotropy);
	static void initializeAnisotropicLocalTurbulentFieldInDiskSource(double*** B, double*** theta, double*** phi, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& R, const double& anisotropy);

	static void initializeParkerField(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R);
	static void initializeParkerFieldWithRotation(double*** B, double*** theta, double*** phi, double*** concentration, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& R, const double& thetaRot);

	static void initializeAngularMask(bool** mask, int Nrho, int Nphi, const double& angle);
	static void initialize3dAngularMask(bool** mask, int Nrho, int Nphi, const double& angle);
	static void initializeRhoMask(bool** mask, int Nrho, int Nphi, const double& fraction);

	static AngleDependentElectronsSphericalSource* createSourceWithTurbulentField(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& theta0, const double& phi0, const double n0, const double& fraction, const double& index, const double& L0, int Nmodes, const double& rho, const double& rhoin, const double& distance);

	static AngleDependentElectronsSphericalSource* createSourceWithParkerField(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& rho, const double& rhoin, const double& distance);
	static AngleDependentElectronsSphericalSource* createSourceWithParkerFieldWithRotation(MassiveParticleDistribution** electronDistributions, int Ntheta, int Nrho, int Nz, int Nphi, const double& B0, const double& n0, const double& v, const double& d, const double& omega, const double& thetaRot, const double& rho, const double& rhoin, const double& distance);

	static TabulatedDiskSource* readDiskSourceFromFile(MassiveParticleDistribution* electronDistribution, const double& rho, const double& zmin, const double zmax, const int Nrho, const int Nz, const int Nphi, const double& distance, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir);
	static void readRectangularSourceArraysFromFile(double*** & B, double***& Btheta, double***& Bphi, double*** & concentration, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir);
	static void readRectangularSourceArraysWithVFromFile(double*** & B, double***& Btheta, double***& Bphi, double***& Vx, double***& Vy, double***& Vz, double*** & concentration, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, SourceInputGeometry geometry, const char* BFileName, const char* VFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir);
	static RectangularSource* readRectangularSourceFromFile(MassiveParticleDistribution* electronDistribution, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, const double& distance, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const double& thetar, const double& phir, const double& psir);
	static ThermalRectangularSource* readThermalRectangularSourceFromFile(const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const int Nx, const int Nz, const int Ny, const double& distance, SourceInputGeometry geometry, const char* BFileName, const char* concentrationFileName, const char* temperatureFileName, const double& thetar, const double& phir, const double& psir);

	static void transformScalarArrayToCartesian(double*** inputArray, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, int Nx, int Nz, int Ny, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const double& thetar, const double& phir, const double& psir);
	static void transformVectorArraysToCartesian(double*** inputArray1, double*** inputArray2, double*** inputArray3, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, double*** outputArrayTheta, double*** outputArrayPhi, int Nx, int Nz, int Ny, const double& minX, const double& maxX, const double& minZ, const double& maxZ, const double& minY, const double& maxY, const double& thetar, const double& phir, const double& psir);

	static void transformScalarArrayToCylindrical(double*** inputArray, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, int Nrho, int Nz, int Nphi, const double& rhomin, const double& rhomax, const double& zmin, const double& zmax, const double& phimin, const double& phimax, const double& thetar, const double& phir, const double& psir);
	static void transformVectorArraysToCylindrical(double*** inputArray1, double*** inputArray2, double*** inputArray3, int N1, int N2, int N3, const double& xmin1, const double& xmax1, const double& xmin2, const double& xmax2, const double& xmin3, const double& xmax3, SourceInputGeometry geometry, double*** outputArray, double*** outputArrayTheta, double*** outputArrayPhi, int Nrho, int Nz, int Nphi, const double& rhomin, const double& rhomax, const double& zmin, const double& zmax, const double& phimin, const double& phimax, const double& thetar, const double& phir, const double& psir);

	static bool*** breadthFirstSearchOfRegion(double*** inputArray, const int N1, const int N2, const int N3, const int start1, const int start2, const int start3, const double& threshold);
	static void correctSourceConcentration(double*** concentration, const int Nx, const int Ny, const int Nz, const int start1, const int start2, const int start3, const double& threshold);
	static void correctSourceTemperature(double*** temperature, const int Nx, const int Ny, const int Nz, const double& value);

	static void evaluateDistributionAfterDiffusion(double* energy, double* distribution, const int Ne, double* D, const double& x, const double& y, const double& z, const double& time, const int Nt, double* sourceDistribution, void (*sourceCoordinates)(const double& time, double& x1, double& x2, double& x3), double (*sourcePower)(const double& time));
	static RectangularSourceInhomogenousDistribution* createRectangularSourceFromDiffusion(const double& mass, double* energy, double* sourceDistribution, const int Ne, double* D, const int Nx, const int Ny, const int Nz, const double& minX, const double& maxX, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& time, const int Nt, void (*sourceCoordinates)(const double& time, double& x1, double& x2, double& x3), double (*sourcePower)(const double& time), const double& B, const double& theta, const double& phi, const double& distance, const double& velocity = 0, const double& redShift = 0);
	static void createRectangularSourceArraysFromDiffusion(const double& mass, double* energy, double* sourceDistribution, const int Ne, double* D, const int Nx, const int Ny, const int Nz, const double& minX, const double& maxX, const double& minY, const double& maxY, const double& minZ, const double& maxZ, const double& time, const int Nt, void (*sourceCoordinates)(const double& time, double& x1, double& x2, double& x3), double (*sourcePower)(const double& time), MassiveParticleDistribution**** & distributions, double*** & concentration);
};

#endif
