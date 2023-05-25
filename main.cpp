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
	FILE* logFile = fopen("log.dat", "w");
	fclose(logFile);

	double theta = pi/2;
	double index = 3.5;
	double Tstar = 50 * 1000;
	double luminosity = 510000 * 4 * 1E33;
	double rsun = 7.5E10;
	double rstar = rsun*sqrt(510000.0/pow(Tstar/5500,4));
	

	double electronConcentration = 25;
	double B = 0.003;
	double rmax = 1.4E17;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	double fraction = 0.5;
	//sigma = 0.002;

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

	int Nrho = 10;
	int Nz = 20;
	int Nphi = 4;

	//initializing mean galactic photon field
	double Ephmin = 0.01 * Tstar * kBoltzman;
	double Ephmax = 100 * Tstar * kBoltzman;
	
	srand(100);
	double*** Bturb = create3dArray(Nrho, Nz, Nphi);
	double*** thetaTurb = create3dArray(Nrho, Nz, Nphi);
	double*** phiTurb = create3dArray(Nrho, Nz, Nphi);
	double*** concentration = create3dArray(Nrho, Nz, Nphi, electronConcentration);
	RadiationSourceFactory::initializeTurbulentField(Bturb, thetaTurb, phiTurb, Nrho, Nz, Nphi, B, pi / 2, 0, 0.9, 11.0 / 6.0, rmax, 10, rmax);
	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);
	//MassiveParticleTabulatedIsotropicDistribution* electronsFromSmilei = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./input/Ee3.dat", "./input/Fs3.dat", 200, electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	//electronsFromSmilei->addPowerLaw(20 * massElectron * speed_of_light2, 3.5);
	//electronsFromSmilei->rescaleDistribution(sqrt(18));
	int Ndistributions = 10;
	//reading electron distributions from files
	//MassiveParticleIsotropicDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./examples_data/gamma0.3_theta0-90/Ee", "./examples_data/gamma0.3_theta0-90/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	MassiveParticleIsotropicDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./examples_data/gamma0.5_theta0-90/Ee", "./examples_data/gamma0.5_theta0-90/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	double newEmax = 1E6 * me_c2;
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->prolongEnergyRange(newEmax, 100);
		if (i < 4) {
			//(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(100 * me_c2, 3.5);
			(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(300 * me_c2, 3.5);
			(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(500 * me_c2, 2);
		}
		
	}

	angleDependentDistributions[1]->writeDistribution("dist1.dat", 300, Emin, newEmax);
	angleDependentDistributions[4]->writeDistribution("dist4.dat", 300, Emin, newEmax);
	angleDependentDistributions[8]->writeDistribution("dist8.dat", 300, Emin, newEmax);

	//creating radiation source
	//RadiationSource* source = new SimpleFlatSource(electronsFromSmilei, B, theta, rsource, 0.1*rsource, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronsFromSmilei, B, theta, electronConcentration, rsource, 0.9*rsource, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronsFromSmilei, Bturb, thetaTurb, concentration, rsource, 0.9*rsource, distance);
	RadiationSource* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, Bturb, thetaTurb, phiTurb, concentration, rmax, 0.9*rmax, distance, 0.5*speed_of_light);
	int nangle = 0;
	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double angle = (dynamic_cast<AngleDependentElectronsSphericalSource*>(source))->getShockWaveAngle(irho, iz, iphi);
				if ((angle < 2*pi / 9)||(angle > 7*pi/9)) {
					nangle++;
				}
			}
		}
	}
	//SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true);
	//SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, newEmax, true, true);
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne+100, Emin, newEmax, true, true);
	//comptonEvaluator->outputDifferentialFlux("output1.dat");
	//return;

	//number of parameters of the source
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1.0E17, 0.000001, 20, 0.001, 0.3*speed_of_light };
	double maxParameters[Nparams] = { 1.5E17, 0.05, 2E6, 0.2, 0.5*speed_of_light };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, fraction, 0.5*speed_of_light };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}

	const int Nenergy1 = 4;
	double energy1[Nenergy1] = { 1.5E9, 3.0E9 , 6.1E9, 9.87E9 };
	double observedFlux[Nenergy1] = { 1.5, 4.3, 6.1, 4.2 };
	double observedError[Nenergy1] = { 0.1, 0.2 , 0.3, 0.2 };
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank;
		observedFlux[i] = observedFlux[i] / (hplank * 1E26);
		observedError[i] = observedError[i] / (hplank * 1E26);
	}

	printf("start optimization\n");
	printLog("start optimization\n");
	bool optPar[Nparams] = { true, true, true, true, false };
	int Niterations = 2;
	//creating gradient descent optimizer
	RadiationOptimizer* synchrotronOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations);
	//RadiationOptimizer* synchrotronOptimizer = new CoordinateRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 5,5,5,5,2 };
	//creating grid enumeration optimizer
	RadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints);
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints);
	//grid enumeration optimization, finding best starting point for gradien descent
	double error = synchrotronOptimizer->evaluateOptimizationFunction(vector, energy1, observedFlux, observedError, Nenergy1, source);
	printf("starting error = %g\n", error);
	printLog("starting error = %g\n", error);
	//enumOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source);
	//gradient descent optimization
	//synchrotronOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source);
	//reseting source parameters to found values
	//synchrotronOptimizer->outputProfileDiagrams(vector, energy1, observedFlux, observedError, Nenergy1, source, 10);
	//synchrotronOptimizer->outputOptimizedProfileDiagram(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source, 10, 1, 2);
	combinedOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source);
	combinedOptimizer->outputOptimizedProfileDiagram(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source, 10, 1, 2);
	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	error = synchrotronOptimizer->evaluateOptimizationFunction(vector, energy1, observedFlux, observedError, Nenergy1, source);
	printf("error = %g\n", error);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	printf("parameters:\n");
	fprintf(paramFile, "parameters:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("average sigma = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "average sigma = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	printf("velocity/c = %g\n", vector[4] * maxParameters[4]/speed_of_light);
	fprintf(paramFile, "celocity/c = %g\n", vector[4] * maxParameters[4]/speed_of_light);

	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("average B = %g\n", B);
	fprintf(paramFile, "average B = %g\n", B);
	fclose(paramFile);

	//initialization arrays for full synchrotron spectrum
	const int Nnu = 200;
	double* Nu = new double[Nnu];
	double* F = new double[Nnu];

	double Numin = 1E8;
	double Numax = 1E19;
	double factor = pow(Numax / Numin, 1.0 / (Nnu - 1));
	Nu[0] = Numin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		Nu[i] = Nu[i - 1] * factor;
		F[i] = 0;
	}

	//evaluationg full spectrum of the source
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nu[i], source);
	}

	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_GZ_Jansky, "%g %g\n", Nu[i] / 1E9, hplank * F[i] * 1E26);
	}
	fclose(output_GZ_Jansky);

	double synchrotronFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(hplank * 1E8, hplank * 1E11, 100, source);
	printf("total synchrotron flux = %g erg/cm^2 s\n", synchrotronFlux);

	//PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	//PhotonIsotropicDistribution* photonDistribution = PhotonPlankDistribution::getCMBradiation();
	PhotonIsotropicDistribution* photonDistribution = new PhotonPlankDistribution(Tstar, sqr(rstar / rmax));

	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_THOMSON);
	//initializing photon energy grid for output

	double* E = new double[Nnu];

	double EphFinalmin = 0.01 * kBoltzman * Tstar;
	double EphFinalmax = 2 * Emax + Emin;
	//photonDistribution->writeDistribution("output3.dat", 200, Ephmin, Ephmax);
	factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
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
	FILE* output_GHz_Jansky = fopen("outputCompt.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i] / (1.6E-9), F[i]);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch.dat", source, Ephmin, Ephmax, 200);
	/*FILE* output_GHz_Jansky = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = 0.0001*E[i] / hplank;
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(0.0001*E[i], source);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_GHz_Jansky);*/

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
