#include "stdio.h"
#include "math.h"
#include "time.h"

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
#include "radiationSourceFactory.h"
#include "coordinateTransform.h"

#include "examples.h"

//example 0 evaluation simple synchrotron
void evaluateSimpleSynchrotron() {
	double B = 1.0;
	double electronConcentration = 1.0;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 3.0, me_c2, 1.0);
	//MassiveParticleIsotropicDistribution* distribution = new MassiveParticleMonoenergeticDistribution(massElectron, 1000*me_c2, me_c2, 1.0);
	RadiationSource* source = new SimpleFlatSource(distribution, B, pi/2, 0, parsec, parsec, 1000 * parsec);
	RadiationEvaluator* evaluator = new SynchrotronEvaluator(10000, me_c2, 10000 * me_c2, true);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	evaluator->writeFluxFromSourceToFile("outputSynch.dat", source, 10 * hplank * cyclotronOmega, 100000 * hplank * cyclotronOmega, 1000);
}

// example 1. Evaluating inverse compton flux of powerlaw distributed electrons on CMB radiation
void evaluateComptonWithPowerLawDistribution() {
	double theta = pi/2;
	//double rmax = 2E14;
	double rmax = 1.0 / sqrt(pi);
	rmax = 1E16;
	double B = 0.6;
	double electronConcentration = 150;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;
	//const double distance = 1.0;

	//double Emin = 652.317 * me_c2 * 1;
	double Emin = me_c2;
	double Emax = 1E4 * me_c2;
	int Ne = 200;
	int Nmu = 20;
	int Nphi = 4;
	double index =3.5;

	//initializing mean galactic photon field
	
    PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);


	RadiationSource* source = new SimpleFlatSource(electrons, B, theta, 0, rmax, rmax, distance);

	double Ephmin = 0.1 * 2.7 * kBoltzman;
	double Ephmax = 10 * 10000 * kBoltzman;
	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);

	//initializing photon energy grid for output
	int Nnu = 100;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double EphFinalmin = 0.01 * kBoltzman * 2.725;
	double EphFinalmax = 2 * Emax;
	double factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
	E[0] = EphFinalmin;
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
	}

	//outputing
	FILE* output_ev_EFE1 = fopen("output1.dat", "w");


	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE1, "%g %g\n", E[i] / (1.6E-12), E[i]*F[i]);
	}
	fclose(output_ev_EFE1);

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
	double Emin = 1.0*me_c2;
	double Emax = 10000 * me_c2;
	double index = 3.0;
	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//creating electrons powerlaw distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);
	//creating radiation source
	SimpleFlatSource* source = new SimpleFlatSource(electrons, B, pi/2, 0, R, fraction * R, distance);
	//number of parameters of the source
	const int Nparams = 4;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E17, 0.0001, 0.5, 0.001 };
	double maxParameters[Nparams] = { 2E17, 1.0, 2E5, 0.1 };
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
	double energy1[Nenergy1] = { 1.5E9 * hplank, 3.0E9 * hplank, 6.1E9 * hplank, 9.87E9 * hplank };
	double observedFlux[Nenergy1] = { 1.5 / (hplank * 1E26), 4.3 / (hplank * 1E26), 6.1 / (hplank * 1E26), 4.2 / (hplank * 1E26) };
	double observedError[Nenergy1] = { 0.1 / (hplank * 1E26), 0.2 / (hplank * 1E26), 0.3 / (hplank * 1E26), 0.2 / (hplank * 1E26) };


	//picking parameters to be optimized
	bool optPar[Nparams] = { false, true, true, true };
	int Niterations = 20;
	//creating KPIevaluator
	LossEvaluator* KPIevaluator = new SpectrumLossEvaluator(energy1, observedFlux, observedError, Nenergy1, source);
	//creating gradient descent optimizer
	RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 10,10,10,10 };
	//creating grid enumeration optimizer
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);

	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error0 = gradientOptimizer->evaluateOptimizationFunction(vector);

	//grid enumeration optimization, finding best starting point for gradien descent
	enumOptimizer->optimize(vector, optPar);
	//gradient descent optimization
	gradientOptimizer->optimize(vector, optPar);
	//reseting source parameters to found values
	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = gradientOptimizer->evaluateOptimizationFunction(vector);

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

	double synchrotronFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(hplank * 0.5E9, hplank * 2E10, 100, source);
	printf("total synchrotron flux = %g erg/cm^2 s\n", synchrotronFlux);

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
	double rmax = 1.4E17;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	double fraction = 0.5;
	//sigma = 0.002;
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
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax, true, true);
	//number of different distributions depending on inclination angle, wich will be read from files
	int Ndistributions = 10;
	//reading electron distributions from files
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(50*me_c2, 3.5);
	}
	(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[4]))->writeDistribution("output4.dat", 200, Emin, Emax);
	//creating radiation source
	int Nrho = 20;
	int Nz = 20;
	int Nphi = 4;
	double*** Bturb = create3dArray(Nrho, Nz, Nphi);
	double*** thetaTurb = create3dArray(Nrho, Nz, Nphi);
	double*** phiTurb = create3dArray(Nrho, Nz, Nphi);
	double*** concentration = create3dArray(Nrho, Nz, Nphi, electronConcentration);
	RadiationSourceFactory::initializeTurbulentField(Bturb, thetaTurb, phiTurb, Nrho, Nz, Nphi, B, pi / 2, 0, 0.1, 11.0 / 6.0, rmax, 10, rmax);
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, Bturb, thetaTurb, phiTurb, concentration, rmax, 0.5 * rmax, distance, 0.3*speed_of_light);
	//AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, B, pi/2, 0, electronConcentration, rmax, 0.5 * rmax, distance);
	//number of parameters of the source
	const int Nparams = 4;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1.3E17, 0.0000001, 100, 0.01 };
	double maxParameters[Nparams] = { 1.4E17, 0.1, 2E6, 0.2 };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, fraction };
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
	double energy1[Nenergy1] = { 1.5E9, 3.0E9 , 6.1E9, 9.87E9 };
	double observedFlux[Nenergy1] = { 1.5, 4.3, 6.1, 4.2 };
	double observedError[Nenergy1] = { 0.1, 0.2 , 0.3, 0.2 };
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank;
		observedFlux[i] = observedFlux[i] / (hplank * 1E26);
		observedError[i] = observedError[i] / (hplank * 1E26);
	}
	//css1610101 t = 357
	//observed parameters of the source in units of GHz and mJansky
	/*const int Nenergy1 = 6;
	double energy1[Nenergy1] = { 0.33E9 * hplank, 0.61E9*hplank, 1.5E9 * hplank, 3.0E9 * hplank, 6.05E9 * hplank, 10E9 * hplank };
	double observedFlux[Nenergy1] = { 0.357 / (hplank * 1E26), 0.79 / (hplank * 1E26), 0.27 / (hplank * 1E26), 0.17 / (hplank * 1E26), 0.07 / (hplank * 1E26), 0.032 / (hplank * 1E26) };
	double observedError[Nenergy1] = { 0.09 / (hplank * 1E26), 0.09 / (hplank * 1E26), 0.07 / (hplank * 1E26), 0.03 / (hplank * 1E26), 0.01 / (hplank * 1E26), 0.008 / (hplank * 1E26) };*/
	//picking parameters to be optimized
	bool optPar[Nparams] = { false, true, true, false };
	int Niterations = 10;
	//creating KPIevaluator
	LossEvaluator* KPIevaluator = new SpectrumLossEvaluator(energy1, observedFlux, observedError, Nenergy1, angleDependentSource);
	//creating gradient descent optimizer
	RadiationOptimizer* synchrotronOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 5,5,5,5 };
	//creating grid enumeration optimizer
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	//grid enumeration optimization, finding best starting point for gradien descent
	enumOptimizer->optimize(vector, optPar);
	//gradient descent optimization
	synchrotronOptimizer->optimize(vector, optPar);
	//reseting source parameters to found values
	angleDependentSource->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = synchrotronOptimizer->evaluateOptimizationFunction(vector);
	printf("error = %g\n", error);
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

	double synchrotronFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(hplank * 1E8, hplank * 1E11, 100, angleDependentSource);
	printf("total synchrotron flux = %g erg/cm^2 s\n", synchrotronFlux);

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

	delete3dArray(Bturb, Nrho, Nz, Nphi);
	delete3dArray(thetaTurb, Nrho, Nz, Nphi);
	delete3dArray(phiTurb, Nrho, Nz, Nphi);
	delete3dArray(concentration, Nrho, Nz, Nphi);
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
	const double cssy3[6] = { 0.375 / (hplank * 1E26), 0.79 / (hplank * 1E26), 0.27 / (hplank * 1E26), 0.17 / (hplank * 1E26), 0.07 / (hplank * 1E26), 0.032 / (hplank * 1E26) };
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

	//initial parameters of source
	double electronConcentration = 150;
	double rmax = 1.3E17;
	//rmax = 1.0;
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
	bool optPar[Nparams] = { true, true, true, true, true, true, true, true };
	

	//number of different distributions depending on inclination angle, wich will be read from files
	const int Ndistributions = 10;

	//reading electron distributions from files
	//ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200, 20 * me_c2, 3.5);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	}
	//angleDependentDistributions[9]->writeDistribution("output1.dat", 200, Emin, Emax);
	//creating radiation source, which does not depend on time
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, pi/2, 0, electronConcentration, rmax, 0.5 * rmax, distance, 0.3*speed_of_light);
	//creating time dependent radiation source
	RadiationTimeDependentSource* source = new ExpandingRemnantSource(rmax, B, electronConcentration, 0.3 * speed_of_light, 0.5, angleDependentSource, times[0]);
	
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 3,3,3,3,3, 3, 3, 3 };
	//number of iterations in gradient descent optimizer
	int Niterations = 5;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating KPI evaluator
	LossEvaluator* KPIevaluator = new TimeDependentSpectrumLossEvaluator(energy, F, Error, Nenergy, times, Ntimes, source);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax, true, true);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	RadiationOptimizer* gridEnumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
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
	RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	gradientOptimizer->optimize(vector, optPar);
	//reset parameters of source to the found values
	source->resetParameters(vector, maxParameters);
	//evaluating final error
	double error = gradientOptimizer->evaluateOptimizationFunction(vector);

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
	printf("sigma = %g\n", vector[1] * maxParameters[1]);
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
	double protonConcentration = 1E-9;
	double rmax = 55 * parsec;
	double B = 0;
	double theta = pi/2;

	//Cygnus
	const double distance = 1400 * parsec;
	double Emin = massProton * speed_of_light2 + 0.01E9 * 1.6E-12;
	double Emax = 1E15 * 1.6E-12;
	double Etrans = 2.2E12 * 1.6E-12;

	MassiveParticleBrokenPowerLawDistribution* protons = new MassiveParticleBrokenPowerLawDistribution(massProton, 2.1, 2.64, Emin, Etrans, protonConcentration);
	//MassiveParticlePowerLawDistribution* protons = new MassiveParticlePowerLawDistribution(massProton, 2.0, Emin, protonConcentration);
	//MassiveParticlePowerLawCutoffDistribution* protons = new MassiveParticlePowerLawCutoffDistribution(massProton, 2.0, Emin, 1.0, Emax, protonConcentration);
	//protons->writeDistribution("outputProtons.dat", 200, Emin, Emax);
	RadiationSource* source = new SimpleFlatSource(protons, B, theta, 0, rmax, rmax, distance);
	double protonAmbientConcentration = 20;
	//PionDecayEvaluator* pionDecayEvaluator = new PionDecayEvaluator(200, Emin, Emax, protonAmbientConcentration);
	PionDecayEvaluatorBase* pionDecayEvaluator = new PionDecayEvaluatorKelner(200, Emin, Emax, protonAmbientConcentration);

	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double Ephmin = 1E9 * 1.6E-12;
	double Ephmax = 2E14 * 1.6E-12;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}

	for (int i = 0; i < Nnu; ++i) {
		F[i] = pionDecayEvaluator->evaluateFluxFromSource(E[i], source);
	}

	FILE* output_ev_dNdE = fopen("outputPionE.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_dNdE, "%g %g\n", E[i] / (1E9*1.6E-12), E[i]*F[i]);
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
	double electronConcentration = 1;
	double rmax = 1.3E17;
	double B = 0.6;
	double temperature = 1E10;
	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 3.08 * 1.0E24;
	double Emin = me_c2;
	double Emax = me_c2 + 100 * kBoltzman * temperature;
	int Ne = 500;

	MassiveParticleMaxwellJuttnerDistribution* electrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, temperature, electronConcentration);
	//MassiveParticleMaxwellDistribution* electrons = new MassiveParticleMaxwellDistribution(massElectron, temperature, electronConcentration);
	RadiationSource* source = new SimpleFlatSource(electrons, 0, 0, 0, rmax, rmax, distance);
	BremsstrahlungThermalEvaluator* bremsstrahlungEvaluator1 = new BremsstrahlungThermalEvaluator();
	BremsstrahlungEvaluator* bremsstrahlungEvaluator2 = new BremsstrahlungEvaluator(Ne, Emin, Emax, 1.0);

	int Nnu = 200;
	double* E = new double[Nnu];

	double Ephmin = 0.001 * kBoltzman * temperature;
	double Ephmax = 100 * kBoltzman * temperature;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
	}

	printLog("evaluating\n");

	FILE* output_ev_EFE = fopen("outputBremE.dat", "w");
	FILE* output_GHz_Jansky = fopen("outputBremNu.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		printf("%d\n", i);
		printLog("%d\n", i);
		double F1 = bremsstrahlungEvaluator1->evaluateFluxFromSource(E[i], source);
		double F2 = bremsstrahlungEvaluator2->evaluateFluxFromSource(E[i], source);
		fprintf(output_ev_EFE, "%g %g %g\n", E[i] / (1.6E-12), E[i] * F1, E[i] * F2);
		fprintf(output_GHz_Jansky, "%g %g %g\n", nu / 1E9, 1E26 * hplank * F1, 1E26 * hplank * F2);
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	delete[] E;
	delete electrons;
	delete source;
	delete bremsstrahlungEvaluator1;
	delete bremsstrahlungEvaluator2;
}

//example 7 compare compton and synchrotron
void compareComptonSynchrotron() {
	double B = 1.0E-5;
	double theta = pi / 2;
	double electronConcentration = 1.0;
	double Emin = 1 * massElectron * speed_of_light2;
	double Emax = 1E8 * me_c2;
	double distance = 1000 * parsec;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 4.0, Emin, 1.0);
	RadiationSource* source = new SimpleFlatSource(distribution, B, theta, 0, parsec, parsec, distance);
	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(400, Emin, Emax, false);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	synchrotronEvaluator->writeFluxFromSourceToFile("output1.dat", source, 0.001 * hplank * cyclotronOmega, 10000000 * hplank * cyclotronOmega, 200);
	double T = 1E3;
	double Ephmin = 0.01 * kBoltzman * T;
	double Ephmax = 100 * kBoltzman * T;
	//PhotonIsotropicDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonIsotropicDistribution* photons = new PhotonPlankDistribution(T, 1.0);
	RadiationEvaluator* inverseComptonEvaluator = new InverseComptonEvaluator(100, 50, 4, Emin, Emax, Ephmin, Ephmax, photons, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* inverseComptonEvaluator2 = new InverseComptonEvaluator(100, 50, 4, Emin, Emax, Ephmin, Ephmax, photons, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	double EminCompton = 0.1*kBoltzman*T;
	double EmaxCompton = 20*sqrt(kBoltzman * T * Emax);
	inverseComptonEvaluator->writeFluxFromSourceToFile("output2.dat", source, EminCompton, EmaxCompton, 100);
	inverseComptonEvaluator2->writeFluxFromSourceToFile("output3.dat", source, EminCompton, EmaxCompton, 100);


	double meanFactor = 0;
	double currentE = Emin;
	double factor = pow(Emax / Emin, 1.0 / (1000 - 1));

	double dzeta3 = 1.202056903;
	double intPlank2 = 2 * dzeta3;
	double dzeta5 = 1.036927755;
	double intPlank4 = 24 * dzeta5;
	double meanE2 = sqr(kBoltzman * T) * intPlank4 / intPlank2;

	for (int i = 0; i < 1000; ++i) {
		double dE = currentE * (factor - 1);
		double gamma = currentE / (me_c2);
		double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
		//double correction = (1.0 - 63.0 * gamma * meanE2 / (me_c2 * photons->getMeanEnergy()));
		double correction = 1.0/(1 + 4*pi*gamma*photons->getMeanEnergy()/me_c2);
		//correction = 1.0;
		meanFactor += 4*pi*gamma * gamma * beta * beta * distribution->distributionNormalized(currentE) * dE * correction;
		currentE *= factor;
	}


	double magneticEnergy = B * B / (8 * pi);
	double photonEnergy = photons->getConcentration() * photons->getMeanEnergy();
	double ratioEnergy = magneticEnergy / photonEnergy;


	double sigmaT = 8.0 * pi * re2 / 3.0;
	double Ptheor = (4.0 / 3.0) * sigmaT * speed_of_light * photonEnergy * meanFactor * electronConcentration * source->getTotalVolume() / (4 * pi * distance * distance);
	double PsynchTheor = (2*sin(theta)*sin(theta)) * sigmaT * speed_of_light * (B*B/(8*pi)) * meanFactor * electronConcentration * source->getTotalVolume() / (4 * pi * distance * distance);

	double synchrotronFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(0.001 * hplank * cyclotronOmega, 10000000 * hplank * cyclotronOmega, 200, source);
	double inverseComptonFlux = inverseComptonEvaluator->evaluateTotalFluxInEnergyRange(EminCompton, EmaxCompton, 200, source);
	double inverseComptonFlux2 = inverseComptonEvaluator2->evaluateTotalFluxInEnergyRange(EminCompton, EmaxCompton, 200, source);
	double ratioFlux = synchrotronFlux / inverseComptonFlux;

	double ratioGhiselini = Ptheor / inverseComptonFlux;
	double ratioKN = Ptheor / inverseComptonFlux2;

	double ratioRatio = ratioFlux / ratioEnergy;

	double ratioComptonEvaluators = inverseComptonFlux / inverseComptonFlux2;

	printf("magneticEnergy/photonEnergy = %g\n", ratioEnergy);
	printf("synchrotron / compton = %g\n", ratioFlux);
	printf("(synchrotron/magnetic) / (compton/photon) = %g\n", ratioRatio);
	printf("theoretical compton/ compton = %g\n", ratioGhiselini);
	printf("theoretical compton/ compton 2 = %g\n", ratioKN);
	printf("theoretical synchrotron/ synchrotron = %g\n", PsynchTheor/synchrotronFlux);
	printf("Jones flux/klein-nishina flux = %g\n", ratioComptonEvaluators);
}

//example 8 evaluation synchrotron Image
void evaluateSynchrotronImage() {
	double B = 1.0;
	double electronConcentration = 1.0;
	int Nrho = 20;
	int Nz = 20;
	int Nphi = 20;
	int Ne = 200;
	double R = 1.4E17;
	double Rin = 0.8 * R;
	const double distance = 150 * 3.08 * 1.0E24;
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 3.0, me_c2, 1.0);
	//MassiveParticleIsotropicDistribution* distribution = new MassiveParticleMonoenergeticDistribution(massElectron, 1000 * me_c2, me_c2, 1.0);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, distribution, B, electronConcentration, pi / 2, R, Rin, distance);
	MassiveParticleDistribution** distributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200, 20 * me_c2, 3.5);
	RadiationSource* source = RadiationSourceFactory::createSourceWithTurbulentField(distributions, 10, Nrho, Nz, Nphi, B, pi / 2, 0, electronConcentration, 0.01, 3.5, 0.5 * R, 10, R, Rin, distance);
	RadiationEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	//evaluator->writeFluxFromSourceToFile("outputSynch.dat", source, 10 * hplank * cyclotronOmega, 10000000 * hplank * cyclotronOmega, 1000);
	evaluator->writeImageFromSourceToFile("image.dat", source, 10 * hplank * cyclotronOmega, 1000000 * hplank * cyclotronOmega, 100);
}

//example 9 test rotation
void testRotation() {
	srand(time(NULL));
	double thetaR = pi * uniformDistribution();
	double phiR = 2 * pi * uniformDistribution();
	//thetaR = 3.1415926535897931;

	//thetaR = pi - 1E-14;
	//phiR = 0.15707963267948966;

	double theta0 = pi * uniformDistribution();
	double phi0 = 2 * pi * uniformDistribution();

	theta0 = 1E-7;
	//phi0 = 0.15707963267948966;

	double theta1;
	double phi1;

	double theta2;
	double phi2;

	resetLog();

	rotationSphericalCoordinates(thetaR, phiR, theta0, phi0, theta1, phi1);
	thetaR = 3.1415926535897931;
	phiR = 0;
	theta1 = 4.7121609281219889E-8;
	phi1 = 0.31415926535897931;
	inverseRotationSphericalCoordinates(thetaR, phiR, theta1, phi1, theta2, phi2);

	printf("theta rotation = %g\n", thetaR);
	printLog("theta rotation = %g\n", thetaR);
	printf("phi rotation = %g\n", phiR);
	printLog("phi rotation = %g\n", phiR);

	printf("theta0 = %g\n", theta0);
	printLog("theta0 = %g\n", theta0);
	printf("phi0 = %g\n", phi0);
	printLog("phi0 = %g\n", phi0);

	printf("theta1 = %g\n", theta1);
	printLog("theta1 = %g\n", theta1);
	printf("phi1 = %g\n", phi1);
	printLog("phi1 = %g\n", phi1);

	printf("theta2 = %g\n", theta2);
	printLog("theta2 = %g\n", theta2);
	printf("phi2 = %g\n", phi2);
	printLog("phi2 = %g\n", phi2);
}

//example 10 test anisotropic compton
void testAnisotropicCompton() {
	double rmax = 1E16;
	double B = 0.3;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;
	//const double distance = 1.0;

	//double Emin = 652.317 * me_c2 * 1;
	double Emin = me_c2;
	double Emax = 1E4 * me_c2;
	Emax = 1E9 * (1.6E-12);
	//Emax = 1000* me_c2;
	//Emin = Emax - me_c2;
	double dE = me_c2;
	int Ne = 50;
	int Nmu = 40;
	int Nrho = 2;
	int Nz = 4;
	int Nphi = 20;
	double index = 3.0;
	double electronConcentration = 5E5;


	double Tstar = 50 * 1000;
	//Tstar = 2.7;
	double Ephmin = 0.1 * Tstar * kBoltzman;
	double Ephmax = 100 * Tstar * kBoltzman;
	double luminosity = 510000 * 4 * 1E33;
	double rsun = 7.5E10;
	double rstar = rsun * sqrt(510000.0 / pow(Tstar / 5500, 4));

	int Nangles = 50;

	//MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, 3.0, Emin, electronConcentration);
	//MassiveParticleIsotropicDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, Emax, dE, electronConcentration);
	//MassiveParticleDistribution* electrons = new MassiveParticleMonoenergeticDirectedDistribution(massElectron, Emax, dE, electronConcentration, 0, 0, pi/100);
	MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee9.dat", "./examples_data/gamma0.5_theta0-90/Fs9.dat", 200, electronConcentration, GAMMA_KIN_FGAMMA);
	electrons->rescaleDistribution(sqrt(18));
	//electrons->addPowerLaw(300 * me_c2, 3.5);
	RadiationSource* source = new SimpleFlatSource(electrons, B, pi/2, 0, rmax, rmax, distance);
	PhotonIsotropicDistribution* photonDummyDistribution = new PhotonPlankDistribution(Tstar, sqr(rstar / rmax));
	InverseComptonEvaluator* comptonEvaluator1 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDummyDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator0 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDummyDistribution, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator1->outputDifferentialFlux("differentialFLux.dat");

	double minEev = 0.3 * 0.001 * 1.6E-12;
	double maxEev = 10 * 0.001 * 1.6E-12;
	int Nph = 10;
	//double kevFlux = comptonEvaluator2->evaluateTotalFluxInEnergyRange(minEev, maxEev, Nph, source);
	//double kevJonesFlux = comptonEvaluator1->evaluateTotalFluxInEnergyRange(minEev, maxEev, Nph, source);
	//printf("isotropic flux = %g\n", kevFlux);
	//printf("isotropic jones flux = %g\n", kevJonesFlux);

	FILE* outFile = fopen("anisotropicCompton.dat", "w");

	for (int i = 0; i < Nangles; ++i) {
		double theta = (i + 0.5) * pi / Nangles;
		PhotonPlankDirectedDistribution* photonDirectedDistribution = new PhotonPlankDirectedDistribution(Tstar, sqr(rstar / rmax), 0, 0, pi / 4);
		
		double kevAnisotropicFlux1 = 0;
		double kevAnisotropicFlux2 = 0;
		double kevAnisotropicFlux3 = 0;
		double factor = pow(maxEev / minEev, 1.0 / (Nph - 1));
		double currentE = minEev;
		double flux = 0;
		/*for (int j = 0; j < Nph; ++j) {
			printf("%d\n", j);
			double dE = currentE * (factor - 1.0);
			kevAnisotropicFlux += comptonEvaluator2->evaluateFluxFromSourceAnisotropic(currentE, 0, 0, photonDirectedDistribution, source) * dE;
			currentE = currentE * factor;
		}*/
		kevAnisotropicFlux1 = comptonEvaluator1->evaluateFluxFromSourceAnisotropic(currentE, theta, 0, photonDirectedDistribution, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
		kevAnisotropicFlux2 = comptonEvaluator1->evaluateFluxFromSourceAnisotropic(currentE, theta, 0, photonDirectedDistribution, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA2);
		kevAnisotropicFlux3 = comptonEvaluator1->evaluateFluxFromSourceAnisotropic(currentE, theta, 0, photonDirectedDistribution, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA3);
		printf("theta = %g flux1 = %g flux2 = %g flux3 = %g\n", theta, kevAnisotropicFlux1, kevAnisotropicFlux2, kevAnisotropicFlux3);
		fprintf(outFile, "%g %g %g %g\n", theta, kevAnisotropicFlux1, kevAnisotropicFlux2, kevAnisotropicFlux3);

		delete photonDirectedDistribution;
	}
	fclose(outFile);
}

//example 11 compare different compton evaluators
void compareComptonWithPowerLawDistribution() {
	double theta = pi / 2;
	//double rmax = 2E14;
	double rmax = 1.0 / sqrt(pi);
	rmax = 1E16;
	double B = 0.6;
	double electronConcentration = 150;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;
	//const double distance = 1.0;

	//double Emin = 652.317 * me_c2 * 1;
	double Emin = me_c2;
	double Emax = 1E10 * me_c2;
	int Ne = 200;
	int Nmu = 40;
	int Nphi = 4;
	double index = 3.5;

	//initializing mean galactic photon field
	double photonEnergy = 1.6E-12;
	double halfWidth = 0.1 * photonEnergy;
	PhotonIsotropicDistribution* photonDistribution = new PhotonMonoenergeticDistribution(photonEnergy, halfWidth, 1.0);
	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);


	RadiationSource* source = new SimpleFlatSource(electrons, B, theta, 0, rmax, rmax, distance);

	double Ephmin = photonEnergy - halfWidth;
	double Ephmax = photonEnergy + halfWidth;
	int evaluatorsNumber = 5;
	InverseComptonEvaluator** evaluators = new InverseComptonEvaluator * [evaluatorsNumber];
	evaluators[0] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES);
	evaluators[1] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	evaluators[2] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	evaluators[3] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA2);
	evaluators[4] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA3);

	//initializing photon energy grid for output
	int Nnu = 50;
	double* E = new double[Nnu];
	double** F = new double* [evaluatorsNumber];
	for (int i = 0; i < evaluatorsNumber; ++i) {
		F[i] = new double[Nnu];
	}

	double EphFinalmin = 0.1*Ephmin;
	double EphFinalmax = 2 * Emax;
	double factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
	E[0] = EphFinalmin;
	for (int i = 0; i < evaluatorsNumber; ++i) {
		F[i][0] = 0;
    }
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		for (int j = 0; j < evaluatorsNumber; ++j) {
			F[j][i] = 0;
		}
	}

	//evaluating radiation flux
	printLog("evaluating\n");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		printLog("%d\n", i);
		for (int j = 0; j < evaluatorsNumber; ++j) {
			F[j][i] = evaluators[j]->evaluateFluxFromSource(E[i], source);
		}
	}

	//outputing
	FILE* output_ev_EFE1 = fopen("outputCompt.dat", "w");


	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE1, "%g", E[i] / (1.6E-12));
		for (int j = 0; j < evaluatorsNumber; ++j) {
			fprintf(output_ev_EFE1, " %g", E[i] * F[j][i]);
		}
		fprintf(output_ev_EFE1, "\n");
	}
	fclose(output_ev_EFE1);

	for (int j = 0; j < evaluatorsNumber; ++j) {
		delete[] F[j];
		delete evaluators[j];
	}
	delete[] E;
	delete[] F;
	delete electrons;
	delete source;
	delete evaluators;
}

//example 12 test reading source from file
void testReadingSource() {
	const char* concentrationFileName = "../Pluto/Tools/pyPLUTO/density.dat";
	const char* BFileName = "../Pluto/Tools/pyPLUTO/B.dat";

	double rmax = 1.5E19;
	double f = 2;
	double B = 0.6;
	double electronConcentration = 150;

	const double distance = 150 * 1000000 * parsec;

	double Emin = me_c2;
	double Emax = 1E10 * me_c2;
	int Ne = 100;
	double index = 3.5;

	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);

	RadiationSource* source = RadiationSourceFactory::readSourceFromFile(electrons, rmax, -f*rmax/2, f * rmax/2, 300, 600, 300, distance, SourceInputGeometry::CYLINDRICAL, BFileName, concentrationFileName, pi/2, 0, pi/2);
	printf("finish creating source\n");
	SynchrotronEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax);

	evaluator->writeImageFromSourceAtEToFile(hplank * 1E9, "image.dat", source);
}

//example fit angle dependent flux

void fitAngleDependentFlux() {
	FILE* logFile = fopen("log.dat", "w");
	fclose(logFile);
	srand(time(NULL));


	//sigma = 0.0002;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;

	double electronConcentration = 25;
	double B = 0.29;
	double theta = pi / 2;
	double index = 3.5;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	double fraction = 0.5;

	double velocity = 0.75 * speed_of_light;
	double R = velocity * 99 * 24 * 3600;

	int Nrho = 10;
	int Nz = 10;
	int Nphi = 20;
	//double*** Bturb = create3dArray(Nrho, Nz, Nphi);
	//double*** thetaTurb = create3dArray(Nrho, Nz, Nphi);
	//double*** phiTurb = create3dArray(Nrho, Nz, Nphi);
	//double*** concentration = create3dArray(Nrho, Nz, Nphi, electronConcentration);
	//RadiationSourceFactory::initializeTurbulentField(Bturb, thetaTurb, phiTurb, Nrho, Nz, Nphi, B, pi / 2, 0, 0.9, 11.0 / 6.0, rmax, 10, rmax);
	//initializing electrons distribution
	//MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, me_c2, electronConcentration);
	//MassiveParticleBrokenPowerLawDistribution* electrons = new MassiveParticleBrokenPowerLawDistribution(massElectron, index, 2.001, 2*me_c2, 1000*me_c2, electronConcentration);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", 259, electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", 259, electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);


	
	//reading electron distributions from files
	int Ndistributions = 10;
	double Emin = me_c2;
	//double Emax = 1E12 * me_c2;
	double Emax = 1E4 * me_c2;
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./examples_data/gamma1.5_theta0-90/Ee", "./examples_data/gamma1.5_theta0-90/Fs", ".dat", Ndistributions, DistributionInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);

	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(1.2));
		(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->prolongEnergyRange(Emax, 100);
		if (i < 4) {
			//(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(100 * me_c2, 3.5);
			(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(300 * me_c2, 3.5);
			//(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->addPowerLaw(500 * me_c2, 2);
		}

	}

	//RadiationSource* source = RadiationSourceFactory::createSourceWithTurbulentField(angleDependentDistributions, 10, Nrho, Nz, Nphi, B, pi / 2, 0, electronConcentration, 0.00000000000001, 3.5, 0.5 * R, 10, R, fraction*R, distance);
	RadiationSource* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, B, theta, 0, electronConcentration, R, (1.0 - fraction) * R, distance);
	int Ne = 50;
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);



	//number of parameters of the source
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.5 * R, 1E-6, 1, 0.001, 0.5 * speed_of_light };
	double maxParameters[Nparams] = { 1.1 * R, 5E-1, 2E4, 0.5, 0.75 * speed_of_light };
	bool optPar[Nparams] = { false, true, true, false, false };
	//starting point of optimization and normalization
	//fraction = 2E13 / rmax;

	//fraction = (0.4-0.04 + 0.004/3);
	fraction = 0.2;
	int Ndays = 99;
	double denseFactor = 0.5 / fraction;
	electronConcentration = sqr(99.0 / Ndays) * 17 * denseFactor;
	//sigma = 0.5*(0.01 * 0.5 / (0.012))/denseFactor;
	//sigma = 0.34 * 0.34 / (4 * pi * massProton * speed_of_light2 * electronConcentration);
	sigma = 0.08333;
	//electronConcentration = 3167 * sqr(69.0 / Ndays);
	//sigma = 2.33E-5;
	double vector[Nparams] = { R, sigma, electronConcentration, fraction, velocity };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}

	//vector[1] = 0.00115637 / maxParameters[1];
	//vector[2] = 788.285 / maxParameters[2];
	//vector[3] = 0.1 / maxParameters[3];

	printf("new sigma = %g\n", sigma);
	printLog("new sigma = %g\n", sigma);
	printf("new concentration = %g\n", electronConcentration);
	printLog("new concentration = %g\n", electronConcentration);

	double* energy1;
	double* observedFlux1;
	double* observedError1;
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans99.txt");
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank * 1E9;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}


	printf("start optimization\n");
	printLog("start optimization\n");
	int Niterations = 5;
	//creating KPIevaluator
	LossEvaluator* lossEvaluator = new SpectrumLossEvaluator(energy1, observedFlux1, observedError1, Nenergy1, source);
	//creating gradient descent optimizer
	RadiationOptimizer* synchrotronOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, lossEvaluator);
	//RadiationOptimizer* synchrotronOptimizer = new CoordinateRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 2,50,50,2,2 };
	//creating grid enumeration optimizer
	RadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints, lossEvaluator);
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, lossEvaluator);
	//grid enumeration optimization, finding best starting point for gradien descent

	//double error = synchrotronOptimizer->evaluateOptimizationFunction(vector);
	//printf("starting error = %g\n", error);
	//printLog("starting error = %g\n", error);

	//enumOptimizer->optimize(vector, optPar, source);
	//gradient descent optimization
	//synchrotronOptimizer->optimize(vector, optPar);
	//reseting source parameters to found values
	//synchrotronOptimizer->outputProfileDiagrams(vector, source, 10);
	//synchrotronOptimizer->outputOptimizedProfileDiagram(vector, optPar, source, 10, 1, 2);

	combinedOptimizer->optimize(vector, optPar);
	

	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = combinedOptimizer->evaluateOptimizationFunction(vector);
	printf("resulting error = %g\n", error);
	printLog("resulting error = %g\n", error);

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
	printf("R/t*c = %g\n", vector[0] * maxParameters[0] / (Ndays * 24 * 3600 * speed_of_light));
	fprintf(paramFile, "R/t*c = %g\n", vector[0] * maxParameters[0] / (Ndays * 24 * 3600 * speed_of_light));
	printf("velocity/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	fprintf(paramFile, "celocity/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("average B = %g\n", B);
	fprintf(paramFile, "average B = %g\n", B);
	fclose(paramFile);



	R = vector[0] * maxParameters[0];
	electronConcentration = vector[2] * maxParameters[2];
	fraction = vector[3] * maxParameters[3];
	velocity = vector[4] * maxParameters[4] / speed_of_light;
	double gamma = 1.0 / sqrt(1.0 - velocity * velocity);

	//double totalEnergy = electronConcentration * massProton * speed_of_light2 * (gamma - 1.0)*(4*pi*rmax*rmax*rmax/3)*(1.0 - cube(1.0 - fraction));
	double totalEnergy = electronConcentration * massProton * speed_of_light2 * (gamma - 1.0) * source->getTotalVolume();
	//double energyInRadioElectrons = electronConcentration * (electrons->getMeanEnergy() - me_c2) * source->getTotalVolume();
	//double energyInRadioProtons = electronConcentration * (protons->getMeanEnergy() - massProton * speed_of_light2) * source->getTotalVolume();
	double magneticEnergy = (B * B / (8 * pi)) * source->getTotalVolume();
	printf("total kinetik energy = %g\n", totalEnergy);
	printLog("total kinetik energy = %g\n", totalEnergy);
	//printf("energy in radio electrons = %g\n", energyInRadioElectrons);
	//printLog("energy in radio electrons = %g\n", energyInRadioElectrons);
	//printf("energy in radio protons = %g\n", energyInRadioProtons);
	//printLog("energy in radio protons = %g\n", energyInRadioProtons);
	printf("energy in magnetic field = %g\n", magneticEnergy);
	printLog("energy in magnetic field = %g\n", magneticEnergy);
	double totalMass = electronConcentration * massProton * source->getTotalVolume();
	printf("total mass = %g\n", totalMass);
	printLog("total mass = %g\n", totalMass);

	
	//initialization arrays for full synchrotron spectrum
	const int Nnu = 50;
	double* Nu = new double[Nnu];
	double* F = new double[Nnu];
	double* F2 = new double[Nnu];

	double Numin = 1E8;
	double Numax = 1E19;
	double factor = pow(Numax / Numin, 1.0 / (Nnu - 1));
	Nu[0] = Numin;
	F[0] = 0;
	F2[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		Nu[i] = Nu[i - 1] * factor;
		F[i] = 0;
		F2[i] = 0;
	}

	//evaluationg full spectrum of the source
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nu[i], source);
	}

	synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch.dat", source, hplank * 1E8, hplank * 1E11, 100);

	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("outputSynch.dat", "w");
	//FILE* output_GZ_Jansky2 = fopen("outputSynch2.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_GZ_Jansky, "%g %g\n", Nu[i] / 1E9, hplank * F[i] * 1E26);
		//fprintf(output_GZ_Jansky2, "%g %g\n", Nu[i] / 1E9, hplank * F2[i] * 1E26);
	}
	fclose(output_GZ_Jansky);
	//fclose(output_GZ_Jansky2);
	
	synchrotronEvaluator->writeImageFromSourceAtEToFile(10*1E9 * hplank, "image.dat", source);
	synchrotronEvaluator->writeImageFromSourceAtEToFile(0.5*1E9 * hplank, "image1.dat", source);
	synchrotronEvaluator->writeImageFromSourceInRangeToFile(0.5*1E9 * hplank, 20*1E9*hplank, 20, "image_array.dat", source);

	//return;

}