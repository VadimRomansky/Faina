#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "util.h"
#include "inverseCompton.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "optimization.h"
#include "pionDecay.h"
#include "bremsstrahlung.h"
#include "coordinateTransform.h"
#include "examples.h"

void evaluateFluxSNRtoWind() {
	double theta = pi/2;
	double rsource = 1.4E17;
	double index = 3.5;
	double electronConcentration = 1E5;
	double Tstar = 50 * 1000;
	double luminosity = 510000 * 4 * 1E33;
	double rsun = 7.5E10;
	double rstar = rsun*sqrt(510000.0/pow(Tstar/5500,4));
	double B = 0.6;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;

	double Emin = me_c2;
	//double Emax = 1E12 * me_c2;
	double Emax = 1E4 * me_c2;
	int Ne = 200;
	int Nmu = 50;
	int Nphi = 4;

	int Nrho = 20;
	int Nz = 20;

	//initializing mean galactic photon field
	double Ephmin = 0.01 * Tstar * kBoltzman;
	double Ephmax = 100 * Tstar * kBoltzman;
	//PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	//PhotonIsotropicDistribution* photonDistribution = PhotonPlankDistribution::getCMBradiation();
	PhotonIsotropicDistribution* photonDistribution = new PhotonPlankDistribution(Tstar, sqr(rstar/rsource));
	double*** Bturb = create3dArray(Nrho, Nz, Nphi);
	double*** thetaTurb = create3dArray(Nrho, Nz, Nphi);
	double*** phiTurb = create3dArray(Nrho, Nz, Nphi);
	double*** concentration = create3dArray(Nrho, Nz, Nphi, electronConcentration);
	RadiationSourceFactory::initializeTurbulentField(Bturb, thetaTurb, phiTurb, Nrho, Nz, Nphi, B, pi / 2, 0, 0.1, 11.0 / 6.0, rsource, 10, rsource);
	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);
	//MassiveParticleTabulatedIsotropicDistribution* electronsFromSmilei = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./input/Ee3.dat", "./input/Fs3.dat", 200, electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	//electronsFromSmilei->addPowerLaw(20 * massElectron * speed_of_light2, 3.5);
	//electronsFromSmilei->rescaleDistribution(sqrt(18));
	int Ndistributions = 10;
	//reading electron distributions from files
	MassiveParticleIsotropicDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(50 * me_c2, 3.5);
	}

	
	//creating radiation source
	//RadiationSource* source = new SimpleFlatSource(electronsFromSmilei, B, theta, rsource, 0.1*rsource, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(20, 20, 4, electronsFromSmilei, B, theta, electronConcentration, rsource, 0.9*rsource, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(20, 20, 4, electronsFromSmilei, Bturb, thetaTurb, concentration, rsource, 0.9*rsource, distance);
	RadiationSource* source = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, Bturb, thetaTurb, phiTurb, concentration, rsource, 0.9*rsource, distance);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_THOMSON);
	//SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true);
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true);
	//comptonEvaluator->outputDifferentialFlux("output1.dat");
	//return;

	//initializing photon energy grid for output
	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double EphFinalmin = 0.01 * kBoltzman * Tstar;
	double EphFinalmax = 2 * Emax + Emin;
	//photonDistribution->writeDistribution("output3.dat", 200, Ephmin, Ephmax);
	double factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
	E[0] = EphFinalmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}

	double minEev = 0.3 * 1000 * 1.6E-12;
	double maxEev = 10 * 1000 * 1.6E-12;
	double kevFlux = comptonEvaluator->evaluateTotalFluxInEnergyRange(minEev, maxEev, 10, source);
	double totalLuminosity = kevFlux * 4 * pi * distance * distance;
	double synchrotronKevFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(minEev, maxEev, 10, source);
	double synchrotronFluxGHz3 = synchrotronEvaluator->evaluateFluxFromSource(3E9 * hplank, source)*1E26*hplank;
	
	FILE* outFile = fopen("SNRtoWindData.dat", "w");
	printf("total luminosity = %g erg/s \n", totalLuminosity);
	fprintf(outFile, "total luminosity = %g erg/s \n", totalLuminosity);
	printf("total flux = %g erg/s cm^2 \n", kevFlux);
	fprintf(outFile, "total flux = %g erg/s cm^2 \n", kevFlux);
	printf("synchrotron total keV flux = %g erg/s cm^2 \n", synchrotronKevFlux);
	fprintf(outFile, "synchrotron total keV flux = %g erg/s cm^2 \n", synchrotronKevFlux);
	printf("99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	fprintf(outFile, "99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	printf("synchrotron radio flux at 3 GHz = %g mJy\n", synchrotronFluxGHz3);
	fprintf(outFile, "synchrotron radio flux at 3 GHz = %g mJy\n", synchrotronFluxGHz3);
	printf("99 days F(3GHz) = 4.3 mJy\n");
	fprintf(outFile, "99 days F(3GHz) = 4.3 mJy\n");
	fclose(outFile);
	//CSS161010
	//99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39
	//130 days F = 1.94+-0.74 10^-15 L = 5.0+-2.5 10^39
	//291 days F << 1.31 10^-15 L << 3.4 10^39
	//H density 4.7*10^20 cm^-2 what does it means?
	double theyRatio = 3.4E39 / 1.33E-15;
	double myRatio = totalLuminosity / kevFlux;

	//evaluating radiation flux
	printLog("evaluating\n");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		printLog("%d\n", i);
		F[i] = comptonEvaluator->evaluateFluxFromSource(E[i], source);
		//F[i] = comptonEvaluator->evaluateFluxFromSourceAnisotropic(E[i], 0, 0, CMBradiation, source);
	}

	//outputing
	FILE* output_ev_EFE = fopen("output.dat", "w");
	//FILE* output_GHz_Jansky = fopen("output.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i] / (1.6E-9), F[i]);
		//fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_ev_EFE);
	//fclose(output_GHz_Jansky);

	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch.dat", source, Ephmin, Ephmax, 200);
	FILE* output_GHz_Jansky = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = 0.0001*E[i] / hplank;
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(0.0001*E[i], source);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_GHz_Jansky);

	delete[] E;
	delete[] F;
	delete electrons;
	delete source;
	delete comptonEvaluator;
}


int main() {
	//evaluateSimpleSynchrotron();
	//evaluateComtonWithPowerLawDistribution();
	//fitCSS161010withPowerLawDistribition();
	//fitCSS161010withTabulatedDistributions();
	//fitTimeDependentCSS161010();
	//evaluatePionDecayWithPowerLawDistribution();
	//evaluateBremsstrahlung();
	//compareComptonSynchrotron();


	evaluateFluxSNRtoWind();

	return 0;
}
