#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "electronDistribution.h"
#include "photonDistribution.h"
#include "util.h"
#include "inverseCompton.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "optimization.h"



// example 1. Evaluating inverse compton flux of powerlaw distributed electrons on CMB radiation
void evaluateComtonWithPowerLawDistribution() {
	double electronConcentration = 150;
	double sinTheta = sin(pi / 6);
	double rmax = 1.3E17;
	double B = 0.6;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 3.08 * 1.0E24;
	double Emin = me_c2;
	double Emax = 10000 * me_c2;

	PhotonPlankDistribution* CMBradiation = PhotonPlankDistribution::getCMBradiation();
	ElectronPowerLawDistribution* electrons = new ElectronPowerLawDistribution(3.5, Emin, electronConcentration);
	RadiationSource* source = new SimpleFlatSource(electrons, B, sinTheta, electronConcentration, rmax, rmax, distance);
	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(200, 20, 20, Emin, Emax);

	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double Ephmin = 0.0001 * kBoltzman * 2.7;
	double Ephmax = 2 * Emax + Emin;
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
		F[i] = comptonEvaluator->evaluateComptonIsotropicFluxFromSource(E[i], CMBradiation, source);
	}

	FILE* output_ev_EFE = fopen("outputE.dat", "w");
	FILE* output_GHz_Jansky = fopen("outputNu.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i] / (1.6E-12), E[i] * E[i] * F[i]);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * E[i] * F[i]);
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	delete[] E;
	delete[] F;
	delete electrons;
	delete source;
	delete comptonEvaluator;
}

//example 2. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with simple flat disk source and powerlaw distribution
void fitCSS161010withPowerLawDistribition() {
	//initial parameters of the source
	double electronConcentration = 150;
	double B = 0.6;
	double rmax = 1.3E17;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//distance to source
	const double distance = 150 * 3.08 * 1.0E24;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//number of different distributions depending on inclination angle, wich will be read from files
	int Ndistributions = 10;
	//creating electrons powerlaw distribution
	ElectronPowerLawDistribution* electrons = new ElectronPowerLawDistribution(3.5, Emin, electronConcentration);
	//creating radiation source
	SimpleFlatSource* source = new SimpleFlatSource(electrons, B, 1.0, electronConcentration, rmax, rmax, distance);
	//number of parameters of the source
	const int Nparams = 4;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E16, 0.01, 0.01, 0.1 };
	double maxParameters[Nparams] = { 2E17, 10, 1000, 1.0 };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, B, electronConcentration, 0.5 };
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
	const int Nnu1 = 4;
	double nu1[Nnu1] = { 1.5E9, 3.0E9, 6.1E9, 0.97E9 };
	double observedFlux[Nnu1] = { 1.5, 4.3, 6.1, 4.2 };
	double observedError[Nnu1] = { 0.1, 0.2, 0.3, 0.2 };
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, true };

	//creating gradient descent optimizer
	SynchrotronOptimizer* synchrotronOptimizer = new GradientDescentSynchrotronOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 20, ErrorScale::LINEAR);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 5,5,5,5 };
	//creating grid enumeration optimizer
	SynchrotronOptimizer* enumOptimizer = new GridEnumSynchrotronOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, ErrorScale::LINEAR, Npoints);
	//grid enumeration optimization, finding best starting point for gradien descent
	enumOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, source);
	//gradient descent optimization
	synchrotronOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, source);
	//reseting source parameters to found values
	source->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = synchrotronOptimizer->evaluateOptimizationFunction(vector, nu1, observedFlux, observedError, Nnu1, source);

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
		F[i] = synchrotronEvaluator->evaluateSynchrotronFluxFromSource(source, Nu[i]);
	}

	//outputing spectrum
	FILE* output_synchr = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_synchr, "%g %g\n", Nu[i] / 1E9, F[i] * 1E26);
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
	printf("B = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "B = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	fclose(paramFile);

	//deleting arrays
	delete[] Nu;
	delete[] F;

	delete synchrotronEvaluator;
	delete electrons;
	delete synchrotronOptimizer;
	delete enumOptimizer;
}

//example 3. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with electron distributions read from files
void fitCSS161010withTabulatedDistributions() {
	//initial parameters of the source
	double electronConcentration = 150;
	double B = 0.6;
	double rmax = 1.3E17;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//distance to source
	const double distance = 150 * 3.08 * 1.0E24;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//number of different distributions depending on inclination angle, wich will be read from files
	int Ndistributions = 10;
	//reading electron distributions from files
	ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<ElectronTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	}
	//creating radiation source
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, 1.0, 0, electronConcentration, rmax, 0.5 * rmax, distance);
	//number of parameters of the source
	const int Nparams = 4;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E16, 0.01, 0.01, 0.1 };
	double maxParameters[Nparams] = { 2E17, 10, 1000, 1.0 };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, B, electronConcentration, 0.5 };
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
	const int Nnu1 = 4;
	double nu1[Nnu1] = { 1.5E9, 3.0E9, 6.1E9, 0.97E9 };
	double observedFlux[Nnu1] = { 1.5, 4.3, 6.1, 4.2 };
	double observedError[Nnu1] = { 0.1, 0.2, 0.3, 0.2 };
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, true };

	//creating gradient descent optimizer
	SynchrotronOptimizer* synchrotronOptimizer = new GradientDescentSynchrotronOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 20, ErrorScale::LINEAR);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 5,5,5,5 };
	//creating grid enumeration optimizer
	SynchrotronOptimizer* enumOptimizer = new GridEnumSynchrotronOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, ErrorScale::LINEAR, Npoints);
	//grid enumeration optimization, finding best starting point for gradien descent
	enumOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, angleDependentSource);
	//gradient descent optimization
	synchrotronOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, angleDependentSource);
	//reseting source parameters to found values
	angleDependentSource->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = synchrotronOptimizer->evaluateOptimizationFunction(vector, nu1, observedFlux, observedError, Nnu1, angleDependentSource);

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
		F[i] = synchrotronEvaluator->evaluateSynchrotronFluxFromSource(angleDependentSource, Nu[i]);
	}

	//outputing spectrum
	FILE* output_synchr = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_synchr, "%g %g\n", Nu[i] / 1E9, F[i] * 1E26);
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
	printf("B = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "B = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
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
	//observed data at 99, 162 and 357 days after explosion in units GHz and mJansky
	const double cssx1[4] = { 1.5, 3.0, 6.1, 9.87 };
	const double cssy1[4] = { 1.5, 4.3, 6.1, 4.2 };
	const double cssError1[4] = { 0.1, 0.2, 0.3, 0.2 };

	const double cssx2[5] = { 1.5, 2.94, 6.1, 9.74, 22.0 };
	const double cssy2[5] = { 4.7, 2.9, 2.3, 1.74, 0.56 };
	const double cssError2[5] = { 0.6, 0.2, 0.1, 0.09, 0.03 };

	const double cssx3[6] = { 0.33, 0.61, 1.5, 3.0, 6.05, 10.0 };
	const double cssy3[6] = { 0.1, 0.79, 0.27, 0.17, 0.07, 0.032 };
	const double cssError3[6] = { 0.375, 0.09, 0.07, 0.03, 0.01, 0.008 };

	//initializing time moments
	const int Ntimes = 3;
	double times[Ntimes] = { 99 * 24 * 3600, 162 * 24 * 3600, 357 * 24 * 3600 };


	//putting observed data into 2d arrays for optimizer
	int Nnu[Ntimes];
	Nnu[0] = 4;
	Nnu[1] = 5;
	Nnu[2] = 6;

	double** Nu = new double* [Ntimes];
	double** F = new double* [Ntimes];
	double** Error = new double* [Ntimes];
	for (int m = 0; m < Ntimes; ++m) {
		Nu[m] = new double[Nnu[m]];
		F[m] = new double[Nnu[m]];
		Error[m] = new double[Nnu[m]];
	}

	for (int i = 0; i < Nnu[0]; ++i) {
		Nu[0][i] = cssx1[i] * 1E9;
		F[0][i] = cssy1[i];
		Error[0][i] = cssError1[i];
	}

	for (int i = 0; i < Nnu[1]; ++i) {
		Nu[1][i] = cssx2[i] * 1E9;
		F[1][i] = cssy2[i];
		Error[1][i] = cssError2[i];
	}

	for (int i = 0; i < Nnu[2]; ++i) {
		Nu[2][i] = cssx3[i] * 1E9;
		F[2][i] = cssy3[i];
		Error[2][i] = cssError3[i];
	}


	//distance to source
	const double distance = 150 * 3.08 * 1.0E24;
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
	//number of optimized parameters
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 1E16, 0.01, 0.01, 0.1, 0.01*speed_of_light};
	double maxParameters[Nparams] = { 2E17, 10, 1000, 1.0, 0.9*speed_of_light};
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, B, electronConcentration, widthFraction, v};
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, true, true };
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 5,5,5,5,5};
	//number of iterations in gradient descent optimizer
	int Niterations = 20;

	//reading electron distributions from files
	//ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200, 50*me_c2, 3.5);
	for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
		(dynamic_cast<ElectronTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	}
	//creating radiation source, which does not depend on time
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, 1.0, 0, electronConcentration, rmax, 0.5 * rmax, distance);
	//creating time dependent radiation source
	RadiationTimeDependentSource* source = new ExpandingRemnantSource(rmax, B, electronConcentration, 0.3 * speed_of_light, 0.5, angleDependentSource, times[0]);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	SynchrotronTimeOptimizer* gridEnumOptimizer = new GridEnumSynchrotronTimeOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, ErrorScale::LINEAR, Npoints);
	gridEnumOptimizer->optimize(vector, optPar, Nu, F, Error, Nnu, Ntimes, times, source);
	/*vector[0] = 2E17 / maxParameters[0];
	vector[1] = 10 / maxParameters[1];
	vector[2] = 0.01 / maxParameters[2];
	vector[3] = 0.325 / maxParameters[3];
	vector[4] = 2.69813E10 / maxParameters[4];*/
	//creating gradient descent optimizer and optimizing
	SynchrotronTimeOptimizer* gradientOptimizer = new GradientDescentSynchrotronTimeOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, ErrorScale::LINEAR);
	gradientOptimizer->optimize(vector, optPar, Nu, F, Error, Nnu, Ntimes, times, source);
	//reset parameters of source to the found values
	source->resetParameters(vector, maxParameters);
	//evaluating final error
	double error = gradientOptimizer->evaluateOptimizationFunction(vector, Nu, F, Error, Nnu, Ntimes, times, source);

	//initializing arrays for evaluationg full spectrum of source with found values
	int Nout = 200;
	double* Nuout = new double[Nout];
	double** Fout = new double*[Ntimes];

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
			Fout[j][i] = synchrotronEvaluator->evaluateSynchrotronFluxFromSource(source1, Nuout[i]);
		}
	}

	//outputing spectrum
	FILE* outFile = fopen("css161010.dat", "w");
	for (int i = 0; i < Nout; ++i) {
		//to GHz and mJansky
		fprintf(outFile, "%g", Nuout[i]*1E-9);
		for (int j = 0; j < Ntimes; ++j) {
			fprintf(outFile, " %g", Fout[j][i]*1E26);
		}
		fprintf(outFile, "\n");
	}
	fclose(outFile);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile,"hi^2 = %g\n", error);
	printf("parameters at first time moment:\n");
	fprintf(paramFile,"parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("B = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "B = %g\n", vector[1] * maxParameters[1]);
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	printf("v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	fprintf(paramFile, "v/c = %g\n", vector[4] * maxParameters[4]/speed_of_light);
	fclose(paramFile);

	//deleting arrays
	for (int i = 0; i < Ntimes; ++i) {
		delete[] Nu[i];
		delete[] F[i];
		delete[] Error[i];
		delete[] Fout[i];
	}
	delete[] Nu;
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


int main() {
	//evaluateComtonWithPowerLawDistribution();
	//fitCSS161010withPowerLawDistribition();
	//fitCSS161010withTabulatedDistributions();
	fitTimeDependentCSS161010();
	return 0;
}