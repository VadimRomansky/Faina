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

#include "examples.h"

//example 0 evaluation simple synchrotron
void evaluateSimpleSynchrotron() {
	double B = 1.0;
	double electronConcentration = 1.0;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 3.0, me_c2, 1.0);
	RadiationSource* source = new SimpleFlatSource(distribution, B, 1.0, parsec, parsec, 1000 * parsec);
	RadiationEvaluator* evaluator = new SynchrotronEvaluator(1000, me_c2, 1000 * me_c2, true);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	evaluator->writeFluxFromSourceToFile("output.dat", source, 10 * hplank * cyclotronOmega, 100000 * hplank * cyclotronOmega, 1000);
}

// example 1. Evaluating inverse compton flux of powerlaw distributed electrons on CMB radiation
void evaluateComtonWithPowerLawDistribution() {
	double sinTheta = 1.0;
	//double rmax = 1.3E17;
	double rmax = 1.0 / sqrt(pi);
	double B = 0.0;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//const double distance = 150 * 1000000 * parsec;
	const double distance = 1.0;

	double Emin = 652.317 * me_c2 * 1;
	double Emax = 1E12 * me_c2;
	int Ne = 200;
	int Nmu = 50;
	int Nphi = 4;
	double index = 2.5;
	double KK = 24990.8;
	double electronConcentration = KK / (pow(652.317, index - 1) * (index - 1));

	//initializing mean galactic photon field
	//PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	PhotonIsotropicDistribution* photonDistribution = PhotonPlankDistribution::getCMBradiation();
	//PhotonIsotropicDistribution* photonDistribution = new PhotonPlankDistribution(1E5, 1.0);

	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);
	//creating radiation source
	RadiationSource* source = new SimpleFlatSource(electrons, B, sinTheta, rmax, rmax, distance);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, photonDistribution, ComptonSolverType::ISOTROPIC_THOMSON);

	//comptonEvaluator->outputDifferentialFlux("output1.dat");
	//return;

	//initializing photon energy grid for output
	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double Ephmin = 0.0001 * kBoltzman * 2.7;
	double Ephmax = 10 * Emax + Emin;
	//photonDistribution->writeDistribution("output3.dat", 200, Ephmin, Ephmax);
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}

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

	delete[] E;
	delete[] F;
	delete electrons;
	delete source;
	delete comptonEvaluator;
}

//example 2. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with simple flat disk source and powerlaw distribution
void fitCSS161010withPowerLawDistribition() {
	//initial parameters of the source
	double electronConcentration = 1.9;
	double B = 0.29;
	double R = 1.4E17;
	double fraction = 0.5 * 4 / 3;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	//sigma = 0.0002;
	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//distance to source
	const double distance = 150 * 1E6 * parsec;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = 5.05 * me_c2;
	double Emax = 10000 * me_c2;
	double index = 3.0;
	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//creating electrons powerlaw distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);
	//creating radiation source
	SimpleFlatSource* source = new SimpleFlatSource(electrons, B, 1.0, R, fraction * R, distance);
	//number of parameters of the source
	const int Nparams = 4;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E17, 0.0001, 0.5, 0.1 };
	double maxParameters[Nparams] = { 2E17, 1.0, 2000, 1.0 };
	//starting point of optimization and normalization
	double vector[Nparams] = { R, sigma, electronConcentration, fraction };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}

	double magneticDensity = B * B / (8 * pi);
	double particleDensity = (index - 1.0) * electronConcentration * Emin / (index - 2.0);

	double EminCop = ((B * B / (8 * pi)) * ((index - 2) / (index - 1)) / electronConcentration) / me_c2;
	double conc = (B * B / (8 * pi)) * ((index - 2) / (index - 1)) / Emin;

	//at2018 t = 15
	//nu1[0] = 35E9;
	//observedInu[0] = 8;
	//nu1[1] = 225E9;
	//observedInu[1] = 30.8;
	//nu1[2] = 233E9;
	//observedInu[2] = 28.6;
	//nu1[3] = 241E9;
	//observedInu[3] = 27.4;
	//Time = 16.5*24*3600;
	//at2018 t = 7.7
	//nu1[0] = 15.5E9;
	//observedInu[0] = 0.489;
	//nu1[1] = 214E9;
	//observedInu[1] = 36.425;
	//nu1[2] = 326E9;
	//observedInu[2] = 32.705;
	//Time = 7.7*24*3600;
	//css161010 t = 357
	/*nu1[0] = 0.33E9;
	observedInu[0] = 0.357;
	observedError[0] = 0.09;
	nu1[1] = 0.61E9;
	observedInu[1] = 0.79;
	observedError[1] = 0.09;
	nu1[2] = 1.5E9;
	observedInu[2] = 0.27;
	observedError[2] = 0.07;
	nu1[3] = 3.0E9;
	observedInu[3] = 0.17;
	observedError[3] = 0.03;
	nu1[4] = 6.05E9;
	observedInu[4] = 0.07;
	observedError[4] = 0.01;
	nu1[5] = 10.0E9;
	observedInu[5] = 0.032;
	observedError[5] = 0.008;
	Time = 357*24*3600;*/
	//css1610101 t = 98
	//observed parameters of the source in units of erg ana cm^-2 s^-1
	const int Nenergy1 = 4;
	double energy1[Nenergy1] = { 1.5E9 * hplank, 3.0E9 * hplank, 6.1E9 * hplank, 9.8E9 * hplank };
	double observedFlux[Nenergy1] = { 1.5 / (hplank * 1E26), 4.3 / (hplank * 1E26), 6.1 / (hplank * 1E26), 4.2 / (hplank * 1E26) };
	double observedError[Nenergy1] = { 0.1 / (hplank * 1E26), 0.2 / (hplank * 1E26), 0.3 / (hplank * 1E26), 0.2 / (hplank * 1E26) };


	//picking parameters to be optimized
	bool optPar[Nparams] = { true, false, true, false };
	int Niterations = 20;
	//creating gradient descent optimizer
	RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 10,10,10,10 };
	//creating grid enumeration optimizer
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints);

	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error0 = gradientOptimizer->evaluateOptimizationFunction(vector, energy1, observedFlux, observedError, Nenergy1, source);

	//grid enumeration optimization, finding best starting point for gradien descent
	enumOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source);
	//gradient descent optimization
	gradientOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, source);
	//reseting source parameters to found values
	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = gradientOptimizer->evaluateOptimizationFunction(vector, energy1, observedFlux, observedError, Nenergy1, source);

	//initialization arrays for full spectrum
	const int Nnu = 200;
	double* Nu = new double[Nnu];
	double* F = new double[Nnu];

	double Numin = 1E8;
	double Numax = 1E18;
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
	FILE* output_synchr = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_synchr, "%g %g\n", Nu[i] / 1E9, F[i] * hplank * 1E26);
	}
	fclose(output_synchr);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	printf("parameters:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("sigma = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "sigma = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	fclose(paramFile);

	//deleting arrays
	delete[] Nu;
	delete[] F;

	delete synchrotronEvaluator;
	delete electrons;
	delete gradientOptimizer;
	delete enumOptimizer;
}

//example 3. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with electron distributions read from files
void fitCSS161010withTabulatedDistributions() {
	//initial parameters of the source
	double electronConcentration = 150;
	double B = 0.6;
	double rmax = 1.3E17;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	sigma = 0.002;
	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//distance to source
	const double distance = 150 * 1E6 * parsec;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//number of different distributions depending on inclination angle, wich will be read from files
	int Ndistributions = 10;
	//reading electron distributions from files
	MassiveParticleIsotropicDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	}
	//angleDependentDistributions[4]->writeDistribution("output1.dat", 200, Emin, Emax);
	//creating radiation source
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, 1.0, 0, electronConcentration, rmax, 0.5 * rmax, distance);
	//number of parameters of the source
	const int Nparams = 4;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E16, 0.0001, 0.01, 0.1 };
	double maxParameters[Nparams] = { 2E17, 1, 1000, 1.0 };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, 0.5 };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//at2018 t = 15
	//nu1[0] = 35E9;
	//observedInu[0] = 8;
	//nu1[1] = 225E9;
	//observedInu[1] = 30.8;
	//nu1[2] = 233E9;
	//observedInu[2] = 28.6;
	//nu1[3] = 241E9;
	//observedInu[3] = 27.4;
	//Time = 16.5*24*3600;
	//at2018 t = 7.7
	//nu1[0] = 15.5E9;
	//observedInu[0] = 0.489;
	//nu1[1] = 214E9;
	//observedInu[1] = 36.425;
	//nu1[2] = 326E9;
	//observedInu[2] = 32.705;
	//Time = 7.7*24*3600;
	//css161010 t = 357
	/*nu1[0] = 0.33E9;
	observedInu[0] = 0.357;
	observedError[0] = 0.09;
	nu1[1] = 0.61E9;
	observedInu[1] = 0.79;
	observedError[1] = 0.09;
	nu1[2] = 1.5E9;
	observedInu[2] = 0.27;
	observedError[2] = 0.07;
	nu1[3] = 3.0E9;
	observedInu[3] = 0.17;
	observedError[3] = 0.03;
	nu1[4] = 6.05E9;
	observedInu[4] = 0.07;
	observedError[4] = 0.01;
	nu1[5] = 10.0E9;
	observedInu[5] = 0.032;
	observedError[5] = 0.008;
	Time = 357*24*3600;*/
	//css1610101 t = 98
	//observed parameters of the source in units of GHz and mJansky
	const int Nenergy1 = 4;
	double energy1[Nenergy1] = { 1.5E9 * hplank, 3.0E9 * hplank, 6.1E9 * hplank, 0.97E9 * hplank };
	double observedFlux[Nenergy1] = { 1.5 / (hplank * 1E26), 4.3 / (hplank * 1E26), 6.1 / (hplank * 1E26), 4.2 / (hplank * 1E26) };
	double observedError[Nenergy1] = { 0.1 / (hplank * 1E26), 0.2 / (hplank * 1E26), 0.3 / (hplank * 1E26), 0.2 / (hplank * 1E26) };
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, false, true, false };

	//creating gradient descent optimizer
	RadiationOptimizer* synchrotronOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 20);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 5,5,5,5 };
	//creating grid enumeration optimizer
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints);
	//grid enumeration optimization, finding best starting point for gradien descent
	enumOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, angleDependentSource);
	//gradient descent optimization
	synchrotronOptimizer->optimize(vector, optPar, energy1, observedFlux, observedError, Nenergy1, angleDependentSource);
	//reseting source parameters to found values
	angleDependentSource->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = synchrotronOptimizer->evaluateOptimizationFunction(vector, energy1, observedFlux, observedError, Nenergy1, angleDependentSource);

	//initialization arrays for full spectrum
	const int Nnu = 200;
	double* Nu = new double[Nnu];
	double* F = new double[Nnu];

	double Numin = 1E8;
	double Numax = 1E18;
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
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nu[i], angleDependentSource);
	}

	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_GZ_Jansky, "%g %g\n", Nu[i] / 1E9, hplank * F[i] * 1E26);
	}
	fclose(output_GZ_Jansky);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	printf("parameters:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("sigma = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "sigma = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);

	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	fclose(paramFile);

	//deleting arrays
	delete[] Nu;
	delete[] F;

	delete synchrotronEvaluator;
	for (int i = 0; i < Ndistributions; ++i) {
		delete angleDependentDistributions[i];
	}
	delete[] angleDependentDistributions;
	delete synchrotronOptimizer;
	delete enumOptimizer;
}

//example 4. Fitting observed synchrotron radio fluxes from CSS161010 at 3 time moments
void fitTimeDependentCSS161010() {
	//observed data at 99, 162 and 357 days after explosion in units erg and cm^-2 s^-2
	const double cssx1[4] = { 1.5 * hplank * 1E9, 3.0 * hplank * 1E9, 6.1 * hplank * 1E9, 9.87 * hplank * 1E9 };
	const double cssy1[4] = { 1.5 / (hplank * 1E26), 4.3 / (hplank * 1E26), 6.1 / (hplank * 1E26), 4.2 / (hplank * 1E26) };
	const double cssError1[4] = { 0.1 / (hplank * 1E26), 0.2 / (hplank * 1E26), 0.3 / (hplank * 1E26), 0.2 / (hplank * 1E26) };

	//remove first point because The measurements below 2 GHz at 162 days post explosion were strongly affected by radio frequency interference
	//const double cssx2[5] = { 1.5 * hplank * 1E9, 2.94 * hplank * 1E9, 6.1 * hplank * 1E9, 9.74 * hplank * 1E9, 22.0 * hplank * 1E9 };
	//const double cssy2[5] = { 4.7 / (hplank * 1E26), 2.9 / (hplank * 1E26), 2.3 / (hplank * 1E26), 1.74 / (hplank * 1E26), 0.56 / (hplank * 1E26) };
	//const double cssError2[5] = { 0.6 / (hplank * 1E26), 0.2 / (hplank * 1E26), 0.1 / (hplank * 1E26), 0.09 / (hplank * 1E26), 0.03 / (hplank * 1E26) };

	const double cssx2[4] = { 2.94 * hplank * 1E9, 6.1 * hplank * 1E9, 9.74 * hplank * 1E9, 22.0 * hplank * 1E9 };
	const double cssy2[4] = { 2.9 / (hplank * 1E26), 2.3 / (hplank * 1E26), 1.74 / (hplank * 1E26), 0.56 / (hplank * 1E26) };
	const double cssError2[4] = { 0.2 / (hplank * 1E26), 0.1 / (hplank * 1E26), 0.09 / (hplank * 1E26), 0.03 / (hplank * 1E26) };

	//todo first point
	const double cssx3[6] = { 0.33 * hplank * 1E9, 0.61 * hplank * 1E9, 1.5 * hplank * 1E9, 3.0 * hplank * 1E9, 6.05 * hplank * 1E9, 10.0 * hplank * 1E9 };
	const double cssy3[6] = { 0.0 / (hplank * 1E26), 0.79 / (hplank * 1E26), 0.27 / (hplank * 1E26), 0.17 / (hplank * 1E26), 0.07 / (hplank * 1E26), 0.032 / (hplank * 1E26) };
	const double cssError3[6] = { 0.375 / (hplank * 1E26), 0.09 / (hplank * 1E26), 0.07 / (hplank * 1E26), 0.03 / (hplank * 1E26), 0.01 / (hplank * 1E26), 0.008 / (hplank * 1E26) };

	//initializing time moments
	const int Ntimes = 3;
	double times[Ntimes] = { 99 * 24 * 3600, 162 * 24 * 3600, 357 * 24 * 3600 };


	//putting observed data into 2d arrays for optimizer
	int Nenergy[Ntimes];
	Nenergy[0] = 4;
	Nenergy[1] = 4;
	Nenergy[2] = 6;

	double** energy = new double* [Ntimes];
	double** F = new double* [Ntimes];
	double** Error = new double* [Ntimes];
	for (int m = 0; m < Ntimes; ++m) {
		energy[m] = new double[Nenergy[m]];
		F[m] = new double[Nenergy[m]];
		Error[m] = new double[Nenergy[m]];
	}

	for (int i = 0; i < Nenergy[0]; ++i) {
		energy[0][i] = cssx1[i];
		F[0][i] = cssy1[i];
		Error[0][i] = cssError1[i];
	}

	for (int i = 0; i < Nenergy[1]; ++i) {
		energy[1][i] = cssx2[i];
		F[1][i] = cssy2[i];
		Error[1][i] = cssError2[i];
	}

	for (int i = 0; i < Nenergy[2]; ++i) {
		energy[2][i] = cssx3[i];
		F[2][i] = cssy3[i];
		Error[2][i] = cssError3[i];
	}


	//distance to source
	const double distance = 150 * 1E6 * parsec;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;

	//number of different distributions depending on inclination angle, wich will be read from files
	const int Ndistributions = 10;

	//initial parameters of source
	double electronConcentration = 150;
	double rmax = 1.3E17;
	double B = 0.6;
	double widthFraction = 0.5;
	double v = 0.3 * speed_of_light;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	//number of optimized parameters
	const int Nparams = 8;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E16, 0.0001, 0.01, 0.1, 0.01 * speed_of_light, 1.1, 1.0, 1.0 };
	double maxParameters[Nparams] = { 2E17, 1, 1000, 1.0, 0.6 * speed_of_light, 2.0, 3.5, 3.5 };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, widthFraction, v, 2.0, 2.0, 3.0 };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, true, true, false, true, true };
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 3,3,3,3,3, 3, 3, 3 };
	//number of iterations in gradient descent optimizer
	int Niterations = 5;

	//reading electron distributions from files
	//ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	MassiveParticleIsotropicDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200, 20 * me_c2, 3.5);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	}
	angleDependentDistributions[9]->writeDistribution("output1.dat", 200, Emin, Emax);
	//creating radiation source, which does not depend on time
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, 1.0, 0, electronConcentration, rmax, 0.5 * rmax, distance);
	//creating time dependent radiation source
	RadiationTimeDependentSource* source = new ExpandingRemnantSource(rmax, B, electronConcentration, 0.3 * speed_of_light, 0.5, angleDependentSource, times[0]);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	RadiationTimeOptimizer* gridEnumOptimizer = new GridEnumRadiationTimeOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints);
	//gridEnumOptimizer->optimize(vector, optPar, energy, F, Error, Nenergy, Ntimes, times, source);
	/*vector[0] = 9.457E16 / maxParameters[0];
	vector[1] = 0.1 / maxParameters[1];
	vector[2] = 560.2341 / maxParameters[2];
	vector[3] = 0.1 / maxParameters[3];
	vector[4] = 1.798E10 / maxParameters[4];
	vector[5] = 2.0 / maxParameters[5];
	vector[6] = 2.0 / maxParameters[6];
	vector[7] = 3.0 / maxParameters[7];*/
	//creating gradient descent optimizer and optimizing
	RadiationTimeOptimizer* gradientOptimizer = new GradientDescentRadiationTimeOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations);
	gradientOptimizer->optimize(vector, optPar, energy, F, Error, Nenergy, Ntimes, times, source);
	//reset parameters of source to the found values
	source->resetParameters(vector, maxParameters);
	//evaluating final error
	double error = gradientOptimizer->evaluateOptimizationFunction(vector, energy, F, Error, Nenergy, Ntimes, times, source);

	//initializing arrays for evaluationg full spectrum of source with found values
	int Nout = 200;
	double* Nuout = new double[Nout];
	double** Fout = new double* [Ntimes];

	double Numin = 1E8;
	double Numax = 1E18;
	double factor = pow(Numax / Numin, 1.0 / (Nout - 1));
	Nuout[0] = Numin;
	for (int j = 0; j < Ntimes; ++j) {
		Fout[j] = new double[Nout];
		Fout[j][0] = 0;
	}

	for (int i = 1; i < Nout; ++i) {
		Nuout[i] = Nuout[i - 1] * factor;
		for (int j = 0; j < Ntimes; ++j) {
			Fout[j][i] = 0;
		}
	}

	//evaluating full spectrum at given time moments
	for (int j = 0; j < Ntimes; ++j) {
		RadiationSource* source1 = source->getRadiationSource(times[j], maxParameters);
		for (int i = 0; i < Nout; ++i) {
			Fout[j][i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nuout[i], source1);
		}
	}

	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("css161010.dat", "w");
	for (int i = 0; i < Nout; ++i) {
		//to GHz and mJansky
		fprintf(output_GZ_Jansky, "%g", Nuout[i] * 1E-9);
		for (int j = 0; j < Ntimes; ++j) {
			fprintf(output_GZ_Jansky, " %g", hplank * Fout[j][i] * 1E26);
		}
		fprintf(output_GZ_Jansky, "\n");
	}
	fclose(output_GZ_Jansky);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	printf("parameters at first time moment:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("B = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "sigma = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	printf("v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	fprintf(paramFile, "v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	printf("r power = %g\n", vector[5] * maxParameters[5] - 1.0);
	fprintf(paramFile, "r power = %g\n", vector[5] * maxParameters[5] - 1.0);
	printf("concentration power = %g\n", vector[6] * maxParameters[6] - 1.0);
	fprintf(paramFile, "concentration power = %g\n", vector[6] * maxParameters[6] - 1.0);
	printf("B power = %g\n", vector[7] * maxParameters[7] - 1.0);
	fprintf(paramFile, "B power = %g\n", vector[7] * maxParameters[7] - 1.0);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	fclose(paramFile);

	//deleting arrays
	for (int i = 0; i < Ntimes; ++i) {
		delete[] energy[i];
		delete[] F[i];
		delete[] Error[i];
		delete[] Fout[i];
	}
	delete[] energy;
	delete[] F;
	delete[] Error;
	delete[] Fout;
	delete[] Nuout;

	for (int i = 0; i < Ndistributions; ++i) {
		delete angleDependentDistributions[i];
	}
	delete[] angleDependentDistributions;
	delete angleDependentSource;
	delete synchrotronEvaluator;
	delete source;
	delete gridEnumOptimizer;
	delete gradientOptimizer;

}

// example 5. Evaluating pion decay gamma flux of powerlaw distributed protons in cygnus cocoon
void evaluatePionDecay() {
	double protonConcentration = 150;
	double rmax = 55 * parsec;
	double B = 0;
	double sinTheta = 1.0;

	//Cygnus
	const double distance = 1400 * parsec;
	double Emin = massProton * speed_of_light2 + 0.01E9 * 1.6E-12;
	double Emax = 1E13 * 1.6E-12;
	double Etrans = 2.2E12 * 1.6E-12;

	MassiveParticleBrokenPowerLawDistribution* protons = new MassiveParticleBrokenPowerLawDistribution(massProton, 2.1, 2.64, Emin, Etrans, protonConcentration);
	//MassiveParticlePowerLawDistribution* protons = new MassiveParticlePowerLawDistribution(massProton, 2.0, Emin, protonConcentration);
	//MassiveParticlePowerLawCutoffDistribution* protons = new MassiveParticlePowerLawCutoffDistribution(massProton, 2.0, Emin, 1.0, Emax, protonConcentration);
	//protons->writeDistribution("outputProtons.dat", 200, Emin, Emax);
	RadiationSource* source = new SimpleFlatSource(protons, B, sinTheta, rmax, rmax, distance);
	double protonAmbientConcentration = 20;
	PionDecayEvaluator* pionDecayEvaluator = new PionDecayEvaluator(200, Emin, Emax, protonAmbientConcentration);

	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double Ephmin = 0.01 * Emin;
	double Ephmax = 1E16 * 1.6E-12;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}

	/*double* sigmaGamma = new double[Nnu];
	double protonKineticEnergy = 100000E9 * 1.6E-12;
	for (int i = 0; i < Nnu; ++i) {
		sigmaGamma[i] = pionDecayEvaluator->sigmaGamma(E[i], protonKineticEnergy);
	}*/

	//double* sigmaInel = new double[Nnu];
	double* Ep = new double[Nnu];

	/*double Epmin = 0.01 * Emin;
	double Epmax = 1E17 * 1.6E-12;
	factor = pow(Epmax / Epmin, 1.0 / (Nnu - 1));
	Ep[0] = Epmin;
	for (int i = 1; i < Nnu; ++i) {
		Ep[i] = Ep[i - 1] * factor;
	}
	FILE* output_sigma_inel = fopen("outputSigmaInel.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_sigma_inel, "%g %g\n", Ep[i], pionDecayEvaluator->sigmaInelastic(Ep[i]));
	}
	fclose(output_sigma_inel);
	delete[] Ep;*/

	printLog("evaluating\n");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		printLog("%d\n", i);
		F[i] = pionDecayEvaluator->evaluateFluxFromSource(E[i], source);
		//F[i] = pionDecayEvaluator->evaluatePionDecayKelnerIsotropicFluxFromSource(E[i], source);
	}

	FILE* output_ev_dNdE = fopen("outputPionE.dat", "w");
	//FILE* output_GHz_Jansky = fopen("outputPionNu.dat", "w");
	//FILE* output_sigma = fopen("outputSigma.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_dNdE, "%g %g\n", E[i] / (1.6E-12), F[i] / E[i]);
		//fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
		//fprintf(output_sigma, "%g %g\n", E[i] / (1.6E-12), sigmaGamma[i]);
	}
	fclose(output_ev_dNdE);
	//fclose(output_GHz_Jansky);
	//fclose(output_sigma);

	delete[] E;
	delete[] F;

	delete protons;
	delete source;
	delete pionDecayEvaluator;
}

// example 6. Evaluating bremsstrahlung radiation from hot gas
void evaluateBremsstrahlung() {
	double electronConcentration = 150;
	double rmax = 1.3E17;
	double B = 0.6;
	double temperature = 1E6;
	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 3.08 * 1.0E24;
	double Emin = me_c2;
	double Emax = 10000 * me_c2;

	MassiveParticleMaxwellDistribution* electrons = new MassiveParticleMaxwellDistribution(massElectron, temperature, electronConcentration);
	RadiationSource* source = new SimpleFlatSource(electrons, 0, 0, rmax, rmax, distance);
	BremsstrahlungThermalEvaluator* bremsstrahlungEvaluator = new BremsstrahlungThermalEvaluator();

	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double Ephmin = 0.001 * kBoltzman * temperature;
	double Ephmax = 100 * kBoltzman * temperature;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}

	printLog("evaluating\n");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		printLog("%d\n", i);
		F[i] = bremsstrahlungEvaluator->evaluateFluxFromSource(E[i], source);
	}

	FILE* output_ev_EFE = fopen("outputBremE.dat", "w");
	FILE* output_GHz_Jansky = fopen("outputBremNu.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i] / (1.6E-12), E[i] * F[i]);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	delete[] E;
	delete[] F;
	delete electrons;
	delete source;
	delete bremsstrahlungEvaluator;
}

//example 7 compare compton and synchrotron
void compareComptonSynchrotron() {
	double B = 1.0E-5;
	double electronConcentration = 1.0;
	double Emin = 20 * massElectron * speed_of_light2;
	double Emax = 10000 * me_c2;
	double distance = 1000 * parsec;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 4.0, Emin, 1.0);
	RadiationSource* source = new SimpleFlatSource(distribution, B, 1.0, parsec, parsec, distance);
	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(400, Emin, Emax, false);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	synchrotronEvaluator->writeFluxFromSourceToFile("output1.dat", source, 0.001 * hplank * cyclotronOmega, 1000000 * hplank * cyclotronOmega, 100);
	PhotonIsotropicDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	RadiationEvaluator* inverseComptonEvaluator = new InverseComptonEvaluator(400, 50, 4, Emin, Emax, photons, ComptonSolverType::ISOTROPIC_JONES);
	inverseComptonEvaluator->writeFluxFromSourceToFile("output2.dat", source, 1E11 * hplank * cyclotronOmega, 1E16 * hplank * cyclotronOmega, 100);

	double meanFactor = 0;
	double currentE = Emin;
	double factor = pow(Emax / Emin, 1.0 / (1000 - 1));
	for (int i = 0; i < 1000; ++i) {
		double dE = currentE * (factor - 1);
		double gamma = currentE / (me_c2);
		double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
		meanFactor += gamma * gamma * beta * beta * distribution->distributionNormalized(currentE) * dE;
		currentE *= factor;
	}


	double magneticEnergy = B * B / (8 * pi);
	double photonEnergy = photons->getConcentration() * photons->getMeanEnergy();
	double ratioEnergy = magneticEnergy / photonEnergy;

	double sigmaT = 8.0 * pi * re2 / 3.0;
	double Ptheor = (4.0 / 3.0) * sigmaT * speed_of_light * photonEnergy * meanFactor * electronConcentration * source->getTotalVolume() / (4 * pi * distance * distance);

	double synchrotronFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(0.001 * hplank * cyclotronOmega, 1000000 * hplank * cyclotronOmega, 100, source);
	double inverseComptonFlux = inverseComptonEvaluator->evaluateTotalFluxInEnergyRange(1E11 * hplank * cyclotronOmega, 1E16 * hplank * cyclotronOmega, 100, source);
	double ratioFlux = synchrotronFlux / inverseComptonFlux;

	double ratioGhiselini = Ptheor / inverseComptonFlux;

	double ratioRatio = ratioFlux / ratioEnergy;

	printf("energy ratio %g\n", ratioEnergy);
	printf("flux ratio = %g\n", ratioFlux);
	printf("ratio ratio = %g\n", ratioRatio);
	printf("ghiselini ratio %g\n", ratioGhiselini * 4 * pi);
}