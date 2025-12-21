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
#include "./Math/specialmath.h"
#include "./Math/largeVectorBasis.h"
#include "./Math/matrixElement.h"

#include "examples.h"

//example 0 evaluation simple synchrotron
void evaluateSimpleSynchrotron() {
	double B = 1.0;
	double electronConcentration = 1.0;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 3.0, me_c2);
	//MassiveParticleIsotropicDistribution* distribution = new MassiveParticleMonoenergeticDistribution(massElectron, 1000*me_c2, me_c2);
	RadiationSourceInCylindrical* source = new SimpleFlatSource(distribution, B, pi/2, 0, electronConcentration, parsec, parsec, 1000 * parsec);
	RadiationEvaluator* evaluator = new SynchrotronEvaluator(10000, me_c2, 10000 * me_c2, true);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	evaluator->writeFluxFromSourceToFile("outputSynch.dat", source, 10 * hplank * cyclotronOmega, 100000 * hplank * cyclotronOmega, 1000);
}

// example 1. Evaluating inverse compton flux of powerlaw distributed electrons on CMB radiation
void evaluateComptonWithPowerLawDistribution() {
	double theta = pi/2;
	//double rmax = 2E14;
	double electronConcentration = 1.9;
	double B = 0.052;
	double rmax = 1E16;
	double fraction = 1.0;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;
	//const double distance = 1.0;

	//double Emin = 652.317 * me_c2 * 1;
	double Emin = me_c2;
	double Emax = 1E17 * 1.6E-12;
	int Ne = 400;
	int Nmu = 20;
	int Nphi = 10;
	double index =2.0;

	int Nph = 100;

	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);


	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, electronConcentration, rmax, fraction * rmax, distance);

	double T = 11000;
	//initializing mean galactic photon field
	double luminosity = 510000 * 4 * 1E33;
	double rsun = 7.5E10;
	double rstar = rsun * sqrt(510000.0 / pow(T / 5772, 4));

	T = 2.7;

	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();


	double Ephmin = 0.1 * T * kBoltzman;
	double Ephmax = 10 * T * kBoltzman;
	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photons, photons->getConcentration(), ComptonSolverType::ISOTROPIC_JONES);

	comptonEvaluator->writeEFEFromSourceToFile("./outputCompton.dat", source, 1.6E-12, 1E16*1.6E-12, 500);

	return;
	
	PhotonPlankDirectedDistribution* photonDistribution1 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4*rmax)), 0, 0, 0.24);
	double photonConcentration1 = photonDistribution1->getConcentration();
	PhotonPlankDirectedDistribution* photonDistribution2 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4 * rmax)), pi/6, 0, 0.24);
	double photonConcentration2 = photonDistribution2->getConcentration();
	PhotonPlankDirectedDistribution* photonDistribution3 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4 * rmax)), pi/3, 0, 0.24);
	double photonConcentration3 = photonDistribution3->getConcentration();
	PhotonPlankDirectedDistribution* photonDistribution4 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4 * rmax)), pi/2, 0, 0.24);
	double photonConcentration4 = photonDistribution4->getConcentration();
	PhotonPlankDirectedDistribution* photonDistribution5 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4 * rmax)), 2*pi/3, 0, 0.24);
	double photonConcentration5 = photonDistribution5->getConcentration();
	PhotonPlankDirectedDistribution* photonDistribution6 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4 * rmax)), 5*pi/6, 0, 0.24);
	double photonConcentration6 = photonDistribution6->getConcentration();
	PhotonPlankDirectedDistribution* photonDistribution7 = new PhotonPlankDirectedDistribution(T, sqr(rstar / (4 * rmax)), pi, 0, 0.24);
	double photonConcentration7 = photonDistribution7->getConcentration();
	PhotonPlankDistribution* photonDistribution8 = new PhotonPlankDistribution(T, 1.0);
	double photonConcentration8 = photonDistribution8->getConcentration();

	InverseComptonEvaluator* comptonEvaluator1 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution1, photonConcentration1, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator2 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution2, photonConcentration2, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator3 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution3, photonConcentration3, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator4 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution4, photonConcentration4, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator5 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution5, photonConcentration5, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator6 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution6, photonConcentration6, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator7 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution7, photonConcentration7, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator8 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution8, photonConcentration8, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator9 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution8, photonConcentration8, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator10 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution8, photonConcentration8, ComptonSolverType::ISOTROPIC_JONES);

	//initializing photon energy grid for output
	int Nnu = 50;
	double* E = new double[Nnu];

	double EphFinalmin = 0.01 * kBoltzman * T;
	double EphFinalmax = 2 * Emax;
	double factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
	E[0] = EphFinalmin;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
	}

	double** F = new double* [10];
	for (int i = 0; i < 10; ++i) {
		F[i] = new double[Nnu];
		for (int j = 0; j < Nnu; ++j) {
			F[i][j] = 0;
		}
	}

	//evaluating radiation flux
	printLog("evaluating\n");

	int i;
#pragma omp parallel for private(i) shared(source, E, F)
	for (i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		printLog("%d\n", i);
		F[0][i] = comptonEvaluator1->evaluateFluxFromSource(E[i], source);
		F[1][i] = comptonEvaluator2->evaluateFluxFromSource(E[i], source);
		F[2][i] = comptonEvaluator3->evaluateFluxFromSource(E[i], source);
		F[3][i] = comptonEvaluator4->evaluateFluxFromSource(E[i], source);
		F[4][i] = comptonEvaluator5->evaluateFluxFromSource(E[i], source);
		F[5][i] = comptonEvaluator6->evaluateFluxFromSource(E[i], source);
		F[6][i] = comptonEvaluator7->evaluateFluxFromSource(E[i], source);
		F[7][i] = comptonEvaluator8->evaluateFluxFromSource(E[i], source);
		F[8][i] = comptonEvaluator9->evaluateFluxFromSource(E[i], source);
		F[9][i] = comptonEvaluator10->evaluateFluxFromSource(E[i], source);
	}

	FILE* output_ev_EFE1 = fopen("outputCompt.dat", "w");
	for (i = 0; i < Nnu; ++i) {
		fprintf(output_ev_EFE1, "%g", E[i] / (1.6E-12));
		for(int j = 0; j < 10; ++j) {
			fprintf(output_ev_EFE1, " %g", E[i] * F[j][i]);
		}
		fprintf(output_ev_EFE1, "\n");
	}
	fclose(output_ev_EFE1);

	for (i = 0; i < Nnu; ++i) {
		delete[] F[i];
	}
	delete[] F;

	delete[] E;
	delete electrons;
	delete source;
	delete comptonEvaluator1;
	delete comptonEvaluator2;
	delete comptonEvaluator3;
	delete comptonEvaluator4;
	delete comptonEvaluator5;
	delete comptonEvaluator6;
	delete comptonEvaluator7;
	delete comptonEvaluator8;
	delete comptonEvaluator9;
	delete comptonEvaluator10;
}

// example 1. Evaluating inverse compton flux of powerlaw distributed electrons on CMB radiation
void evaluateComptonOneParticle() {
	double theta = pi / 2;
	//double rmax = 2E14;
	double electronConcentration = 1.9;
	double B = 0.052;
	double rmax = 1E16;
	double fraction = 1.0;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;
	//const double distance = 1.0;

	int Ne = 8;
	int Nmu = 500;
	int Nphi = 4;
	double index = 2.0;

	int Nph = 100;

	double Ee = 1E10 * me_c2;
	double Tph = 10000.0;
	double Eph = kBoltzman * Tph;
	double dEe = 0.1 * Ee;
	double dEph = 0.1 * Eph;

	//initializing electrons distribution
	MassiveParticleMonoenergeticDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, Ee, dEe);


	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, electronConcentration, rmax, fraction * rmax, distance);

	PhotonMonoenergeticDistribution* photons = new PhotonMonoenergeticDistribution(Eph, dEph);
	double photonConcentration = 1.0;

	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Ee - dEe, Ee + dEe, Nph, Eph - dEph, Eph + dEph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);

	comptonEvaluator->writeEFEFromSourceToFile("./outputCompton.dat", source, 1E6*1.6E-12, 1E16 * 1.6E-12, 100);
}

//example 2. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with simple flat disk source and powerlaw distribution
void fitCSS161010withPowerLawDistribition() {
	//initial parameters of the source
	/*double electronConcentration = 25;
	double B = 0.29;
	double R = 1.41E17;*/
	double electronConcentration = 1.9;
	double B = 0.052;
	double R = 3.5E17;
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
	double Emin = 8*me_c2;
	double Emax = 10000 * me_c2;
	double index = 3.5;

	double EminCop = ((B * B / (8 * pi)) * ((index - 2) / (index - 1)) / electronConcentration) / me_c2;
	double conc = (B * B / (8 * pi)) * ((index - 2) / (index - 1)) / Emin;

	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);
	//creating electrons powerlaw distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);
	//creating radiation source
	SimpleFlatSource* source = new SimpleFlatSource(electrons, B, pi/2, 0, electronConcentration, R, fraction * R, distance);
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
	//enumOptimizer->optimize(vector, optPar);
	//gradient descent optimization
	//gradientOptimizer->optimize(vector, optPar);
	//reseting source parameters to found values
	//source->resetParameters(vector, maxParameters);
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
	double electronConcentration = 50;
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
	double Emax = 1000 * me_c2;
	//creating synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax, true, true);
	//number of different distributions depending on inclination angle, wich will be read from files
	int Ndistributions = 10;
	//reading electron distributions from files
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, 200);
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
	double* energy1;
	double* observedFlux1;
	double* observedError1;
	//int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans69.txt");
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans99.txt");
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank * 1E9;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}

	double* energy2;
	double* observedFlux2;
	double* observedError2;
	//int Nenergy2 = readRadiationFromFile(energy2, observedFlux2, observedError2, "./examples_data/css_data/coppejans99.txt");
	int Nenergy2 = readRadiationFromFile(energy2, observedFlux2, observedError2, "./examples_data/css_data/coppejans162.txt");
	//int Nenergy2 = readRadiationFromFile(energy2, observedFlux2, observedError2, "./examples_data/css_data/coppejans357.txt");
	for (int i = 0; i < Nenergy2; ++i) {
		energy2[i] = energy2[i] * hplank * 1E9;
		observedFlux2[i] = observedFlux2[i] / (hplank * 1E26);
		observedError2[i] = observedError2[i] / (hplank * 1E26);
	}

	double* energy3;
	double* observedFlux3;
	double* observedError3;
	int Nenergy3 = readRadiationFromFile(energy3, observedFlux3, observedError3, "./examples_data/css_data/coppejans357.txt");
	for (int i = 0; i < Nenergy3; ++i) {
		energy3[i] = energy3[i] * hplank * 1E9;
		observedFlux3[i] = observedFlux3[i] / (hplank * 1E26);
		observedError3[i] = observedError3[i] / (hplank * 1E26);
	}

	printf("finish reading data\n");
	printLog("finish reading data\n");


	//initializing time moments
	const int Ntimes = 3;
	//const int Ntimes = 2;
	const int Ntimes3 = 3;
	double times[Ntimes] = { 99 * 24 * 3600, 162 * 24 * 3600, 357 * 24 * 3600 };
	//double times[Ntimes] = {99 * 24 * 3600, 357 * 24 * 3600 };
	double times3[Ntimes3] = { 99 * 24 * 3600, 162*24*3600, 357 * 24 * 3600 };


	//putting observed data into 2d arrays for optimizer
	int Nenergy[Ntimes];
	Nenergy[0] = Nenergy1;
	Nenergy[1] = Nenergy2;
	Nenergy[2] = Nenergy3;

	double** energy = new double* [Ntimes];
	double** F = new double* [Ntimes];
	double** Error = new double* [Ntimes];
	for (int m = 0; m < Ntimes; ++m) {
		energy[m] = new double[Nenergy[m]];
		F[m] = new double[Nenergy[m]];
		Error[m] = new double[Nenergy[m]];
	}

	for (int i = 0; i < Nenergy[0]; ++i) {
		energy[0][i] = energy1[i];
		F[0][i] = observedFlux1[i];
		Error[0][i] = observedError1[i];
	}

	for (int i = 0; i < Nenergy[1]; ++i) {
		energy[1][i] = energy2[i];
		F[1][i] = observedFlux2[i];
		Error[1][i] = observedError2[i];
	}

	for (int i = 0; i < Nenergy[2]; ++i) {
		energy[2][i] = energy3[i];
		F[2][i] = observedFlux3[i];
		Error[2][i] = observedError3[i];
	}


	//distance to source
	const double distance = 150 * 1E6 * parsec;

	//initial parameters of source
	double electronConcentration = 50;
	double rmax = times[0]*0.75*speed_of_light;
	//rmax = 1.0;
	double B = 0.6;
	double widthFraction = 0.25;
	double v = 0.75 * speed_of_light;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	//number of optimized parameters
	const int Nparams = 9;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.9*rmax, 0.000001, 1, 0.1, 0.5 * speed_of_light, 1.0, 1.0, 1.0, 1.0 };
	double maxParameters[Nparams] = { times[0] * 0.9 * speed_of_light, 2, 50, 0.2, 0.8 * speed_of_light, 3.0, 4.0, 4.0, 3.0 };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, widthFraction, v, 1.01, 2.0, 3.0, 1.0 };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, true, true, true, true, true, true };

	int numberOfOptpar = 0;
	for (int i = 0; i < Nparams; ++i) {
		if (optPar[i]) {
			numberOfOptpar++;
		}
	}

	int numberOfPoints = 0;
	for (int i = 0; i < Ntimes; ++i) {
		numberOfPoints += Nenergy[i];
	}
	
	printf("creating distributions\n");
	printLog("creating distributions\n");

	//number of different distributions depending on inclination angle, wich will be read from files
	const int Ndistributions = 10;

	//reading electron distributions from files
	//ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, 200, 20 * me_c2, 3.5);
	//for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
	//	(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	//}
	//angleDependentDistributions[9]->writeDistribution("output1.dat", 200, Emin, Emax);
	//creating radiation source, which does not depend on time
	//MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_theta0-90/Ee3.dat", "./examples_data/gamma1.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);

	printf("distribution 1\n");
	printLog("distribution 1\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution1 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_theta0-90/Ee3.dat", "./examples_data/gamma0.2_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 2\n");
	printLog("distribution 2\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution2 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.3_theta0-90/Ee3.dat", "./examples_data/gamma0.3_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 3\n");
	printLog("distribution 3\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution3 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 4\n");
	printLog("distribution 4\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution4 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	electronDistribution1->addPowerLaw(10 * me_c2, 3.5);
	electronDistribution1->rescaleDistribution(1.2);
	electronDistribution2->addPowerLaw(20 * me_c2, 3.5);
	electronDistribution2->rescaleDistribution(1.2);
	electronDistribution3->addPowerLaw(50 * me_c2, 3.5);
	electronDistribution3->rescaleDistribution(1.2);

	double velocities[4] = { 0.2 * speed_of_light, 0.3 * speed_of_light, 0.5 * speed_of_light, 0.75 * speed_of_light };

	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution*[4];
	distributions[0] = electronDistribution1;
	distributions[1] = electronDistribution2;
	distributions[2] = electronDistribution3;
	distributions[3] = electronDistribution4;

	//electronDistribution1->writeDistribution("distribution1.dat", 1000, me_c2, 1E9*me_c2);
	//electronDistribution2->writeDistribution("distribution2.dat", 1000, me_c2, 1E9 * me_c2);
	//electronDistribution3->writeDistribution("distribution3.dat", 1000, me_c2, 1E9 * me_c2);
	electronDistribution4->writeDistribution("distribution.dat", 1000, me_c2, 1E9 * me_c2);


	printf("creating sources\n");
	printLog("creating sources\n");

	//MassiveParticleIsotropicDistribution* electronDistribution = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 10*me_c2, electronConcentration);
	//electronDistribution->rescaleDistribution(1.2);
	//SimpleFlatSource* source1 = new SimpleFlatSource(electronDistribution, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);
	//SimpleFlatSource2* source1 = new SimpleFlatSource2(4, velocities, distributions, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);
	int Nrho = 100;
	int Nz = 200;
	int Nphi = 1;
	SphericalLayerSource* source1 = new TabulatedSphericalLayerSource2(4, velocities, distributions, Nrho, Nz, Nphi, B, pi / 2, 0, electronConcentration, rmax, (1.0 - widthFraction) * rmax, distance);
	//SimpleFlatSource2* source1 = new SimpleFlatSource2(4, velocities, distributions, B, pi / 2, 0, electronConcentration, rmax, widthFraction * rmax, distance);
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, pi/2, 0, electronConcentration, rmax, 0.5 * rmax, distance, 0.3*speed_of_light);
	//creating time dependent radiation source
	//RadiationTimeDependentSource* source = new ExpandingRemnantSource(rmax, B, electronConcentration, 0.3 * speed_of_light, 0.5, angleDependentSource, times[0]);
	RadiationTimeDependentSource* source = new ExpandingRemnantSource(rmax, B, electronConcentration, 0.75 * speed_of_light, 0.1, source1, times[0]);
	
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 4,4,4,4,4, 4, 4, 4,4 };
	//number of iterations in gradient descent optimizer
	int Niterations = 4;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating KPI evaluator
	LossEvaluator* KPIevaluator = new TimeDependentSpectrumLossEvaluator(energy, F, Error, Nenergy, times, Ntimes, source);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax, true, false);
	CombinedRadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints, KPIevaluator);
	SequentCoordinateEnumOptimizer* sequentOptimizer = new SequentCoordinateEnumOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 200, 2, KPIevaluator);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	//RadiationOptimizer* gridEnumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	//gridEnumOptimizer->optimize(vector, optPar, energy, F, Error, Nenergy, Ntimes, times, source);
	vector[0] = 1.76253e+17 / maxParameters[0];
	vector[1] = 0.1 / maxParameters[1];
	vector[2] = 50 / maxParameters[2];
	vector[3] = 0.100539 / maxParameters[3];
	vector[4] = 0.755575 *speed_of_light / maxParameters[4];
	vector[5] = 1.300079 / maxParameters[5];
	vector[6] = 1.155571 / maxParameters[6];
	vector[7] = 3.02757 / maxParameters[7];
	vector[8] = 1.00 / maxParameters[8];

	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 0, "error0.dat");
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 1, "error1.dat");*/
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 2, "error2.dat");
	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 3, "error3.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 4, "error4.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 5, "error5.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 6, "error6.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 7, "error7.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 8, "error8.dat");*/

	printf("optimization\n");
	printLog("optimization\n");
	//creating gradient descent optimizer and optimizing
	//RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//gradientOptimizer->optimize(vector, optPar);
	//combinedOptimizer->optimize(vector, optPar);
	//sequentOptimizer->optimize(vector, optPar);
	//reset parameters of source to the found values
	source->resetParameters(vector, maxParameters);
	//evaluating final error
	//double error = combinedOptimizer->evaluateOptimizationFunction(vector);
	double error = 100000;

	//initializing arrays for evaluationg full spectrum of source with found values
	printf("evaluationg radiation\n");
	printLog("evaluationg radiation\n");
	int Nout = 200;
	double* Nuout = new double[Nout];
	double** Fout = new double* [Ntimes3];

	double Numin = 1E8;
	double Numax = 1E12;
	double factor = pow(Numax / Numin, 1.0 / (Nout - 1));
	Nuout[0] = Numin;
	for (int j = 0; j < Ntimes3; ++j) {
		Fout[j] = new double[Nout];
		Fout[j][0] = 0;
	}

	for (int i = 1; i < Nout; ++i) {
		Nuout[i] = Nuout[i - 1] * factor;
		for (int j = 0; j < Ntimes3; ++j) {
			Fout[j][i] = 0;
		}
	}

	//evaluating full spectrum at given time moments
	for (int j = 0; j < Ntimes3; ++j) {
		RadiationSourceInCylindrical* source1 = source->getRadiationSource(times3[j], maxParameters);
		for (int i = 0; i < Nout; ++i) {
			Fout[j][i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nuout[i], source1);
			//Fout[j][i] = 0;
		}
	}

	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("css161010.dat", "w");
	for (int i = 0; i < Nout; ++i) {
		//to GHz and mJansky
		fprintf(output_GZ_Jansky, "%g", Nuout[i] * 1E-9);
		for (int j = 0; j < Ntimes3; ++j) {
			fprintf(output_GZ_Jansky, " %g", hplank * Fout[j][i] * 1E26);
		}
		fprintf(output_GZ_Jansky, "\n");
	}
	fclose(output_GZ_Jansky);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	int degreesOfFreedom = numberOfPoints - numberOfOptpar - 1;
	printf("relative hi^2 = %g\n", error/degreesOfFreedom);
	fprintf(paramFile, "relative hi^2 = %g\n", error / degreesOfFreedom);
	double normalDist = (error - degreesOfFreedom) / sqrt(2 * degreesOfFreedom);
	printf("degrees of freedom = %d\n", degreesOfFreedom);
	fprintf(paramFile, "degrees of freedom = %d\n", degreesOfFreedom);
	printf("normal distribution = %g\n", normalDist);
	fprintf(paramFile,"normal distribution = %g\n", normalDist);
	printf("parameters at first time moment:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("R/ct = %g\n", vector[0] * maxParameters[0]/(speed_of_light*times[0]));
	fprintf(paramFile, "R/ct = %g\n", vector[0] * maxParameters[0]/(speed_of_light*times[0]));
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
	printf("width power = %g\n", vector[8] * maxParameters[8] - 1.0);
	fprintf(paramFile, "width power = %g\n", vector[8] * maxParameters[8] - 1.0);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	fclose(paramFile);

	printf("evaluating wide-range radiation\n");
	printLog("evaluating wide-range radiation\n");


	Nrho = 2000;
	Nz = 4000;
	Nphi = 1;

	double R = vector[0] * maxParameters[0];
	//B
	electronConcentration = vector[2] * maxParameters[2];
	double f = vector[3] * maxParameters[3];
	double velocity = vector[4] * maxParameters[4];
	double downstreamV = 0.25 * velocity;

	TabulatedSLSourceWithSynchCutoff* source2 = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance, downstreamV, velocity);
	//TabulatedDiskSource* source2 = new TabulatedDiskSource(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance);
	TabulatedDiskSourceWithSynchAndComptCutoff* source3 = new TabulatedDiskSourceWithSynchAndComptCutoff(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, velocity);

	int Ne = 1000;
	Emin = me_c2;
	Emax = me_c2 * 2E9;

	electronDistribution4->writeDistribution("distribution.dat",Ne, Emin, Emax);

	MassiveParticleIsotropicDistribution* distribution1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source2->getParticleDistribution(Nrho-1, Nz/2, 0));
	distribution1->writeDistribution("distribution1.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution2 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source3->getParticleDistribution(0, Nz-1, 0));
	distribution2->writeDistribution("distribution2.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution3 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source2->getParticleDistribution(0, Nz-10, 0));
	distribution3->writeDistribution("distribution3.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution4 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source3->getParticleDistribution(0, Nz - 10, 0));
	distribution4->writeDistribution("distribution4.dat", Ne, Emin, Emax);

	SynchrotronEvaluator* evaluator2 = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);

	double kevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source2);

	double mevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source2);

	double kevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source3);

	double mevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source3);

	printf("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);
	printLog("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);

	printf("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);
	printLog("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);

	printf("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);
	printLog("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);

	printf("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);
	printLog("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);

	double totalKineticEnergy = source2->getTotalVolume() * massProton * electronConcentration * 0.5 * speed_of_light2;
	printf("totalKineticEnergy = %g\n", totalKineticEnergy);
	printLog("totalKineticEnergy = %g\n", totalKineticEnergy);

	evaluator2->writeFluxFromSourceToFile("wideRangeSynch.dat", source2, 1E8 * hplank, 200 * MeV, 500);
	evaluator2->writeFluxFromSourceToFile("wideRangeSynchDisk.dat", source3, 1E8 * hplank, 200 * MeV, 500);

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
	//delete gridEnumOptimizer;
	//delete gradientOptimizer;

}

// example 5. Evaluating pion decay gamma flux of powerlaw distributed protons in cygnus cocoon
void evaluatePionDecay() {
	double protonConcentration = 1E-9;
	double rmax = 55 * parsec;
	double B = 0;
	double theta = pi/2;

	//Cygnus
	const double distance = 1400 * parsec;
	double Emin = massProton * speed_of_light2;
	//Emin = 1E15 * 1.6E-12;
	double Emax = 1E16 * 1.6E-12;
	double Etrans = 2.2E12 * 1.6E-12;

	//MassiveParticleBrokenPowerLawDistribution* protons = new MassiveParticleBrokenPowerLawDistribution(massProton, 2.1, 2.64, Emin, Etrans);
	//MassiveParticlePowerLawDistribution* protons = new MassiveParticlePowerLawDistribution(massProton, 2.0, Emin);
	MassiveParticlePowerLawCutoffDistribution* protons = new MassiveParticlePowerLawCutoffDistribution(massProton, 2.0, Emin, 1.0, 1E15*1.6E-12);
	//protons->writeDistribution("outputProtons.dat", 200, Emin, Emax);
	RadiationSourceInCylindrical* source = new SimpleFlatSource(protons, B, theta, 0, protonConcentration, rmax, rmax, distance);
	double protonAmbientConcentration = 12;
	PionDecayEvaluator* pionDecayEvaluator = new PionDecayEvaluator(2000, Emin, Emax, protonAmbientConcentration);
	PionDecayEvaluatorBase* pionDecayEvaluator2 = new PionDecayEvaluatorKelner(2000, Emin, Emax, protonAmbientConcentration);

	int Nnu = 2000;
	double* E = new double[Nnu];
	double* F = new double[Nnu];
	double* F2 = new double[Nnu];

	double Ephmin = 1E5 * 1.6E-12;
	double Ephmax = 2E15 * 1.6E-12;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	F[0] = 0;
	F2[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
		F2[i] = 0;
	}

	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		F[i] = pionDecayEvaluator->evaluateFluxFromSource(E[i], source);
		F2[i] = pionDecayEvaluator2->evaluateFluxFromSource(E[i], source);
	}

	FILE* output_ev_dNdE = fopen("outputPionE.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_dNdE, "%g %g %g\n", E[i] / (1E9*1.6E-12), E[i]*F[i], E[i] * F2[i]);
	}
	fclose(output_ev_dNdE);
	//fclose(output_GHz_Jansky);
	//fclose(output_sigma);

	delete[] E;
	delete[] F;
	delete[] F2;

	delete protons;
	delete source;
	delete pionDecayEvaluator;
}

// example 6. Evaluating bremsstrahlung radiation from hot gas
void evaluateBremsstrahlung() {
	double electronConcentration = 0.2E8;
	double rmax = 0.75*1000000*1000*100;
	double B = 0.6;
	double temperature = 5772;
	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//const double distance = 150 * 3.08 * 1.0E24;
	const double distance = 150.0 * 1000000.0 * 1000.0 * 100.0;
	double Emin = me_c2;
	double Emax = me_c2 + 100 * kBoltzman * temperature;
	int Ne = 500;

	//MassiveParticleMaxwellJuttnerDistribution* electrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, temperature, electronConcentration);
	//MassiveParticleMaxwellDistribution* electrons = new MassiveParticleMaxwellDistribution(massElectron, temperature, electronConcentration);
	//MassiveParticleMaxwellDistribution* electrons = new MassiveParticleMaxwellDistribution(massElectron, temperature, electronConcentration);
	RadiationSource* source = new ThermalRectangularSource(1, 1, 1, massElectron, B, pi/2, 0, electronConcentration, temperature, -rmax, rmax, -rmax, rmax, -rmax, rmax, distance);
	BremsstrahlungThermalEvaluator* bremsstrahlungEvaluator1 = new BremsstrahlungThermalEvaluator(true, false);
	BremsstrahlungEvaluator* bremsstrahlungEvaluator2 = new BremsstrahlungEvaluator(Ne, Emin, Emax, 1.0, true, false);

	int Nnu = 200;
	double* E = new double[Nnu];

	double Ephmin = 0.00001 * kBoltzman * temperature;
	double Ephmax = 100 * kBoltzman * temperature;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
	}

	printLog("evaluating\n");

	int N = 5;

	double vector[6];
	double normalizationUnits[6] = { 1,1,1,1,1,1 };
	vector[0] = 2*rmax;
	vector[1] = 2 * rmax;
	vector[2] = 2 * rmax;
	vector[3] = B * B / (4 * pi * electronConcentration * massProton * speed_of_light2);
	vector[4] = electronConcentration;
	vector[5] = 0;

	FILE* output_ev_EFE = fopen("outputBremE.dat", "w");
	FILE* output_GHz_Jansky = fopen("outputBremNu.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		printf("%d\n", i);
		printLog("%d\n", i);
		fprintf(output_ev_EFE, "%g", E[i] / (1.6E-12));
		fprintf(output_GHz_Jansky, "%g", nu / 1E9);
		for(int j = 0; j < N; ++j){
			vector[4] = electronConcentration;
			source->resetParameters(vector, normalizationUnits);
			electronConcentration *= 20;
			double F1 = bremsstrahlungEvaluator1->evaluateFluxFromSource(E[i], source);
			double F2 = bremsstrahlungEvaluator2->evaluateFluxFromSource(E[i], source);
			fprintf(output_ev_EFE, " %g %g", E[i] * F1, E[i] * F2);
			fprintf(output_GHz_Jansky, " %g %g", 1E26 * hplank * F1, 1E26 * hplank * F2);
		}
		for (int j = 0; j < N; ++j) {
			electronConcentration /= 20;
		}
		double photonFinalFrequency = E[i] / hplank;
		double B = (1/hplank)*(2 * hplank * cube(photonFinalFrequency) / (speed_of_light2)) / (exp(E[i] / (kBoltzman * temperature)) - 1.0);
		double F3 = B * 4 * rmax * rmax / (distance * distance);
		fprintf(output_ev_EFE, " %g\n", E[i]*F3);
		fprintf(output_GHz_Jansky, " %g\n", 1E26*hplank*F3);
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	double totalFlux = bremsstrahlungEvaluator2->evaluateTotalFluxInEnergyRange(Ephmin, Ephmax, 1000, source);
	printf("total flux = %g erg /s cm^2\n", totalFlux);

	delete[] E;
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
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 4.0, Emin);
	RadiationSourceInCylindrical* source = new SimpleFlatSource(distribution, B, theta, 0, electronConcentration, parsec, parsec, distance);
	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(400, Emin, Emax, false);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	synchrotronEvaluator->writeFluxFromSourceToFile("output1.dat", source, 0.001 * hplank * cyclotronOmega, 10000000 * hplank * cyclotronOmega, 200);
	double T = 1E3;
	double Ephmin = 0.01 * kBoltzman * T;
	double Ephmax = 100 * kBoltzman * T;
	//PhotonIsotropicDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photons = new PhotonPlankDistribution(T, 1.0);
	double photonConcentration = photons->getConcentration();
	RadiationEvaluator* inverseComptonEvaluator = new InverseComptonEvaluator(100, 50, 4, Emin, Emax, 100, Ephmin, Ephmax, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* inverseComptonEvaluator2 = new InverseComptonEvaluator(100, 50, 4, Emin, Emax, 100, Ephmin, Ephmax, photons, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
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
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, 3.0, me_c2);
	//MassiveParticleIsotropicDistribution* distribution = new MassiveParticleMonoenergeticDistribution(massElectron, 1000 * me_c2, me_c2, 1.0);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, distribution, B, electronConcentration, pi / 2, R, Rin, distance);
	MassiveParticleDistribution** distributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, 200, 20 * me_c2, 3.5);
	RadiationSourceInCylindrical* source = RadiationSourceFactory::createSourceWithTurbulentField(distributions, 10, Nrho, Nz, Nphi, B, pi / 2, 0, electronConcentration, 0.01, 3.5, 0.5 * R, 10, R, Rin, distance);
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
	int Nph = 100;
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
	MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee9.dat", "./examples_data/gamma0.5_theta0-90/Fs9.dat", GAMMA_KIN_FGAMMA);
	electrons->rescaleDistribution(sqrt(18));
	//electrons->addPowerLaw(300 * me_c2, 3.5);
	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, pi/2, 0, electronConcentration, rmax, rmax, distance);
	PhotonPlankDistribution* photonDummyDistribution = new PhotonPlankDistribution(Tstar, sqr(rstar / rmax));
	double photonConcentration = photonDummyDistribution->getConcentration();
	InverseComptonEvaluator* comptonEvaluator1 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDummyDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator0 = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDummyDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator1->outputDifferentialFlux("differentialFLux.dat");

	double minEev = 0.3 * 0.001 * 1.6E-12;
	double maxEev = 10 * 0.001 * 1.6E-12;
	//int Nph = 10;
	//double kevFlux = comptonEvaluator2->evaluateTotalFluxInEnergyRange(minEev, maxEev, Nph, source);
	//double kevJonesFlux = comptonEvaluator1->evaluateTotalFluxInEnergyRange(minEev, maxEev, Nph, source);
	//printf("isotropic flux = %g\n", kevFlux);
	//printf("isotropic jones flux = %g\n", kevJonesFlux);

	FILE* outFile = fopen("anisotropicCompton.dat", "w");

	for (int i = 0; i < Nangles; ++i) {
		double theta = (i + 0.5) * pi / Nangles;
		PhotonPlankDirectedDistribution* photonDirectedDistribution = new PhotonPlankDirectedDistribution(Tstar, sqr(rstar / rmax), 0, 0, pi / 4);
		photonConcentration = photonDirectedDistribution->getConcentration();
		
		double kevAnisotropicFlux1 = 0;
		double kevAnisotropicFlux2 = 0;
		double kevAnisotropicFlux3 = 0;
		//double factor = pow(maxEev / minEev, 1.0 / (Nph - 1));
		double currentE = minEev;
		double flux = 0;
		/*for (int j = 0; j < Nph; ++j) {
			printf("%d\n", j);
			double dE = currentE * (factor - 1.0);
			kevAnisotropicFlux += comptonEvaluator2->evaluateFluxFromSourceAnisotropic(currentE, 0, 0, photonDirectedDistribution, source) * dE;
			currentE = currentE * factor;
		}*/
		kevAnisotropicFlux1 = comptonEvaluator1->evaluateFluxFromSourceAnisotropic(currentE, theta, 0, photonDirectedDistribution, photonConcentration, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
		kevAnisotropicFlux2 = comptonEvaluator1->evaluateFluxFromSourceAnisotropic(currentE, theta, 0, photonDirectedDistribution, photonConcentration, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
		kevAnisotropicFlux3 = comptonEvaluator1->evaluateFluxFromSourceAnisotropic(currentE, theta, 0, photonDirectedDistribution, photonConcentration, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
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
	int Nph = 100;
	double index = 3.5;

	//initializing mean galactic photon field
	double photonEnergy = 1.6E-12;
	double halfWidth = 0.1 * photonEnergy;
	PhotonIsotropicDistribution* photonDistribution = new PhotonMonoenergeticDistribution(photonEnergy, halfWidth);
	double photonConcentration = 1.0;
	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);


	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, electronConcentration, rmax, rmax, distance);

	double Ephmin = photonEnergy - halfWidth;
	double Ephmax = photonEnergy + halfWidth;
	int evaluatorsNumber = 5;
	InverseComptonEvaluator** evaluators = new InverseComptonEvaluator * [evaluatorsNumber];
	evaluators[0] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	evaluators[1] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	evaluators[2] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	evaluators[3] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	evaluators[4] = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Nph, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);

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

	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);

	RadiationSourceInCylindrical* source = RadiationSourceFactory::readDiskSourceFromFile(electrons, rmax, -f*rmax/2, f * rmax/2, 300, 600, 300, distance, SourceInputGeometry::CYLINDRICAL, BFileName, concentrationFileName, pi/2, 0, pi/2);
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

	double electronConcentration = 1.9;
	double B = 0.052;
	double theta = pi / 2;
	double index = 3.5;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	double fraction = 0.5;

	double velocity = 0.75 * speed_of_light;
	double R = velocity * 357 * 24 * 3600;

	int Nrho = 10;
	int Nz = 10;
	int Nphi = 20;
	//double*** Bturb = create3dArray(Nrho, Nz, Nphi);
	//double*** thetaTurb = create3dArray(Nrho, Nz, Nphi);
	//double*** phiTurb = create3dArray(Nrho, Nz, Nphi);
	//double*** concentration = create3dArray(Nrho, Nz, Nphi, electronConcentration);
	//RadiationSourceFactory::initializeTurbulentField(Bturb, thetaTurb, phiTurb, Nrho, Nz, Nphi, B, pi / 2, 0, 0.9, 11.0 / 6.0, rmax, 10, rmax);
	//initializing electrons distribution
	//MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, me_c2);
	//MassiveParticleBrokenPowerLawDistribution* electrons = new MassiveParticleBrokenPowerLawDistribution(massElectron, index, 2.001, 2*me_c2, 1000*me_c2);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", 259, DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", 259, DistributionInputType::GAMMA_KIN_FGAMMA);


	
	//reading electron distributions from files
	int Ndistributions = 10;
	double Emin = me_c2;
	//double Emax = 1E12 * me_c2;
	double Emax = 1E4 * me_c2;
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./examples_data/gamma1.5_theta0-90/Ee", "./examples_data/gamma1.5_theta0-90/Fs", ".dat", Ndistributions, DistributionInputType::GAMMA_KIN_FGAMMA, 200);

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
	RadiationSourceInCylindrical* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, B, theta, 0, electronConcentration, R, (1.0 - fraction) * R, distance, velocity);
	int Ne = 50;
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);



	//number of parameters of the source
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.5 * R, 1E-6, 0.01, 0.001, 0.5 * speed_of_light };
	double maxParameters[Nparams] = { 1.1 * R, 5E-1, 2E4, 0.5, 0.75 * speed_of_light };
	bool optPar[Nparams] = { false, true, true, false, false };
	//starting point of optimization and normalization
	//fraction = 2E13 / rmax;

	//fraction = (0.4-0.04 + 0.004/3);
	fraction = 0.2;
	int Ndays = 357;
	//double denseFactor = 0.5 / fraction;
	//electronConcentration = 20*sqr(99.0 / Ndays) * 17 * denseFactor;
	//sigma = 0.5*(0.01 * 0.5 / (0.012))/denseFactor;
	//sigma = 0.34 * 0.34 / (4 * pi * massProton * speed_of_light2 * electronConcentration);
	//sigma = 0.8333;
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
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans357.txt");
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
	int Npoints[Nparams] = { 2,20,20,2,2 };
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
	double Numax = 1E11;
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

	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch.dat", source, hplank * 1E8, hplank * 1E11, 100);

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

//test vary parameters
void varyParameters() {
	FILE* logFile = fopen("log.dat", "w");
	fclose(logFile);
	srand(time(NULL));

	const double distance = 150 * 1000000 * parsec;

	double electronConcentration = 25;
	double B = 0.29;
	double theta = pi / 2;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	double fraction = 0.5;

	double velocity = 0.75 * speed_of_light;
	double R = velocity * 99 * 24 * 3600;

	double Emin = me_c2;
	double index = 3.5;
	double Emax = 1E4 * me_c2;

	//velocity = 0.0001 * velocity;

	/*int Ndistributions = 10;
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributions(massElectron, "./examples_data/gamma1.5_theta0-90/Ee", "./examples_data/gamma1.5_theta0-90/Fs", ".dat", Ndistributions, DistributionInputType::GAMMA_KIN_FGAMMA, 200);

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

	int Nrho = 20;
	int Nphi = 10;
	int Nz = 10;
	RadiationSource* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, B, theta, 0, electronConcentration, R, (1.0 - fraction) * R, distance, 0.0001*velocity);*/

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);
	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, pi / 2, 0, R, fraction * R, distance, velocity);

	int Ne = 50;
	SynchrotronEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true, true);

	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.5 * R, 1E-6, 1, 0.001, 0.5 * speed_of_light };
	double maxParameters[Nparams] = { 1.1 * R, 5E-1, 2E5, 0.5, 0.75 * speed_of_light };

	double vector[Nparams] = { R, sigma, electronConcentration, fraction, velocity };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}

	char buffer[1];

	for (int i = 0; i < 5; ++i) {
		source->resetParameters(vector, maxParameters);
		vector[1] *= 2;
		//_itoa(i, buffer, 10);
		sprintf(buffer, "%d",i);
		evaluator->writeFluxFromSourceToFile((std::string("outputSynch") + std::string(buffer) + std::string(".dat")).c_str(), source, hplank * 1E8, hplank * 1E11, 100);
	}
	for (int i = 0; i < 5; ++i) {
		vector[1] /= 2;
	}
	/*for (int i = 5; i < 10; ++i) {
		source->resetParameters(vector, maxParameters);
		vector[2] *= 2;
		//_itoa(i, buffer, 10);
		sprintf(buffer, "%d", i);
		evaluator->writeFluxFromSourceToFile((std::string("outputSynch") + std::string(buffer) + std::string(".dat")).c_str(), source, hplank * 1E8, hplank * 1E11, 100);
	}
	for (int i = 0; i < 5; ++i) {
		vector[2] /= 2;
	}*/
	for (int i = 5; i < 10; ++i) {
		source->resetParameters(vector, maxParameters);
		vector[0] *= 2;
		//_itoa(i, buffer, 10);
		sprintf(buffer, "%d", i);
		evaluator->writeFluxFromSourceToFile((std::string("outputSynch") + std::string(buffer) + std::string(".dat")).c_str(), source, hplank * 1E8, hplank * 1E11, 100);
	}

}

//example test versin and delta
void testVersin() {
	double x = 0;
	double v = versin(x);

	double gamma = 1E8;
	double delta = relativisticDelta(gamma);
}

//example test Bessel functions from numerical recepies
void testBessel() {
	printf("Bessel0(0.00001) = %g\n", bessj0(0.00001));
	printf("Bessel0(0.001) = %g\n", bessj0(0.001));
	printf("Bessel0(0.05) = %g\n", bessj0(0.05));
	printf("Bessel0(1.0) = %g\n", bessj0(1.0));
	printf("Bessel0(5.0) = %g\n", bessj0(5.0));
	printf("Bessel0(10.0) = %g\n", bessj0(10.0));
	printf("Bessel0(20.0) = %g\n", bessj0(20.0));
	printf("Bessel0(100.0) = %g\n", bessj0(100.0));

	printf("Bessel1(0.00001) = %g\n", bessj1(0.00001));
	printf("Bessel1(0.001) = %g\n", bessj1(0.001));
	printf("Bessel1(0.05) = %g\n", bessj1(0.05));
	printf("Bessel1(1.0) = %g\n", bessj1(1.0));
	printf("Bessel1(5.0) = %g\n", bessj1(5.0));
	printf("Bessel1(10.0) = %g\n", bessj1(10.0));
	printf("Bessel1(20.0) = %g\n", bessj1(20.0));
	printf("Bessel1(100.0) = %g\n", bessj1(100.0));

	double I, K, Ip, Kp;
	bessik(0.00001, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3,0.00001) = %g\n", K);
	bessik(0.001, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 0.001) = %g\n", K);
	bessik(0.05, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 0.05) = %g\n", K);
	bessik(1.0, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 1.0) = %g\n", K);
	bessik(5.0, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 5.0) = %g\n", K);
	bessik(10.0, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 10.0) = %g\n", K);
	bessik(20.0, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 20.0) = %g\n", K);
	bessik(100.0, 2.0 / 3.0, &I, &K, &Ip, &Kp);
	printf("BesselK(2/3, 100.0) = %g\n", K);
}

//example test Chevalier model
void testChevalier() {
	double B = 10.0;
	double electronConcentration = 10.0;
	double E0 = me_c2;
	double index = 3.5;
	double c1 = 0.6265E19;
	double c5 = 0.5013E-23;
	double c6 = 0.49697E-40;
	double N0 = electronConcentration * (index - 1) * pow(E0, index - 1);
	double D = 1000 * parsec;
	double R = parsec;
	double f = 0.5;
	double s = (4.0 / 3.0)*f * R;
	MassiveParticleIsotropicDistribution* distribution = new MassiveParticlePowerLawDistribution(massElectron, index, E0);
	//MassiveParticleIsotropicDistribution* distribution = new MassiveParticleMonoenergeticDistribution(massElectron, 1000*me_c2, me_c2, 1.0);
	RadiationSourceInCylindrical* source = new SimpleFlatSource(distribution, B, pi / 2, 0, electronConcentration, R, s, D);
	RadiationEvaluator* evaluator = new SynchrotronEvaluator(10000, me_c2, 10000 * me_c2, true);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);

	const int Nnu = 1000;
	double* Nu = new double[Nnu];
	double* F = new double[Nnu];
	double* F1 = new double[Nnu];
	double* F2 = new double[Nnu];
	double Numin = cyclotronOmega;
	double Numax = 10000 * cyclotronOmega;
	double factor = pow(Numax / Numin, 1.0 / (Nnu - 1));
	Nu[0] = Numin;
	F[0] = 0;
	F1[0] = 0;
	F2[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		Nu[i] = Nu[i - 1] * factor;
		F[i] = 0;
		F1[i] = 0;
		F2[i] = 0;
	}

	FILE* file = fopen("outputSynch3.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		F[i] = hplank*evaluator->evaluateFluxFromSource(hplank * Nu[i], source);
		F1[i] = (pi * R * R / (D * D)) * (c5 / c6) * sqrt(1 / B) * pow(Nu[i] / (2 * c1), 5.0 / 2.0);
		F2[i] = (4 * pi * f * R * R * R / (3 * D * D)) * c5 * N0 * pow(B, (index + 1) / 2.0) * pow(Nu[i] / (2 * c1), -(index - 1) / 2.0);
		fprintf(file, "%g %g %g %g\n", Nu[i], F[i], F1[i], F2[i]);
	}
	fclose(file);

	delete[] Nu;
	delete[] F;
	delete[] F1;
	delete[] F2;
}

//example 4. Fitting observed synchrotron radio fluxes from CSS161010 at 3 time moments

void fitCSS161010() {
	//observed data at 99, 162 and 357 days after explosion in units erg and cm^-2 s^-2
	double* energy1;
	double* observedFlux1;
	double* observedError1;
	//int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans69.txt");
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans99.txt");
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank * 1E9;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}

	printf("finish reading data\n");
	printLog("finish reading data\n");



	double timeMoment = 99 * 24 * 3600;


	//distance to source
	const double distance = 150 * 1E6 * parsec;

	//initial parameters of source
	double electronConcentration = 17*5;
	double rmax = timeMoment * 0.75 * speed_of_light;
	//rmax = 1.0;
	double B = 0.34;
	double widthFraction = 0.1;
	double v = 0.75 * speed_of_light;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	//number of optimized parameters
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.6 * rmax, 0.000001, 1, 0.05, 0.5 * speed_of_light };
	double maxParameters[Nparams] = { timeMoment * 0.8 * speed_of_light, 0.05, 200, 0.2, 0.8 * speed_of_light };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, widthFraction, v };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//picking parameters to be optimized
	bool optPar[Nparams] = { false, true, true, false, false };

	int numberOfOptpar = 0;
	for (int i = 0; i < Nparams; ++i) {
		if (optPar[i]) {
			numberOfOptpar++;
		}
	}

	int numberOfPoints = Nenergy1;


	printf("creating distributions\n");
	printLog("creating distributions\n");

	int Nrho = 200;
	int Nz = 400;
	int Nphi = 1;
	const int Ndistributions = 10;

	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, 200, 20 * me_c2, 3.5);
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, pi / 2, 0, electronConcentration, rmax, 0.5 * rmax, distance, 0.3 * speed_of_light);





	MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//electronDistribution1->addPowerLaw(10 * me_c2, 3.5);
	//electronDistribution1->rescaleDistribution(1.2);


	electronDistribution->writeDistribution("distributionCSS161010.dat", 1000, me_c2, 1E9 * me_c2);


	printf("creating sources\n");
	printLog("creating sources\n");

	SimpleFlatSource* source1 = new SimpleFlatSource(electronDistribution, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);

	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 4,4,4,4,4 };
	//number of iterations in gradient descent optimizer
	int Niterations = 4;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating KPI evaluator
	LossEvaluator* KPIevaluator = new SpectrumLossEvaluator(energy1, observedFlux1, observedError1, Nenergy1, angleDependentSource);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(50, Emin, Emax, true, false);
	CombinedRadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints, KPIevaluator);
	SequentCoordinateEnumOptimizer* sequentOptimizer = new SequentCoordinateEnumOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 200, 2, KPIevaluator);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	//RadiationOptimizer* gridEnumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	//gridEnumOptimizer->optimize(vector, optPar, energy, F, Error, Nenergy, Ntimes, times, source);
	
	//vector[0] = 2.02253e+17 / maxParameters[0];
	//vector[1] = 0.0436792 / maxParameters[1];
	//vector[2] = 65.199 / maxParameters[2];
	//vector[3] = 0.104334 / maxParameters[3];
	//vector[4] = 0.755575 * speed_of_light / maxParameters[4];

	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 0, "error0.dat");
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 1, "error1.dat");*/
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 2, "error2.dat");
	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 3, "error3.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 4, "error4.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 5, "error5.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 6, "error6.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 7, "error7.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 8, "error8.dat");*/

	printf("optimization\n");
	printLog("optimization\n");
	//creating gradient descent optimizer and optimizing
	RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//combinedOptimizer->optimize(vector, optPar);
	//sequentOptimizer->optimize(vector, optPar);
	//gradientOptimizer->optimize(vector, optPar);
	//reset parameters of source to the found values
	source1->resetParameters(vector, maxParameters);
	//evaluating final error
	double error = sequentOptimizer->evaluateOptimizationFunction(vector);
	//double error = 100000;

	double* uncertainties = sequentOptimizer->evaluateUncertainties(vector, optPar);

	//initializing arrays for evaluationg full spectrum of source with found values
	printf("evaluationg radiation\n");
	printLog("evaluationg radiation\n");
	int Nout = 200;
	double* Nuout = new double[Nout];
	double* Fout = new double[Nout];

	double Numin = 1E8;
	double Numax = 1E12;
	double factor = pow(Numax / Numin, 1.0 / (Nout - 1));
	Nuout[0] = Numin;
	for (int j = 0; j < Nout; ++j) {
		Fout[j] = 0;
	}

	for (int i = 1; i < Nout; ++i) {
		Nuout[i] = Nuout[i - 1] * factor;
	}

	//evaluating full spectrum at given time moments
	for (int i = 0; i < Nout; ++i) {
		Fout[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nuout[i], angleDependentSource);
	}


	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("css161010.dat", "w");
	for (int i = 0; i < Nout; ++i) {
		//to GHz and mJansky
		fprintf(output_GZ_Jansky, "%g", Nuout[i] * 1E-9);
		fprintf(output_GZ_Jansky, " %g", hplank * Fout[i] * 1E26);
		fprintf(output_GZ_Jansky, "\n");
	}
	fclose(output_GZ_Jansky);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	int degreesOfFreedom = numberOfPoints - numberOfOptpar - 1;
	printf("relative hi^2 = %g\n", error / degreesOfFreedom);
	fprintf(paramFile, "relative hi^2 = %g\n", error / degreesOfFreedom);
	double normalDist = (error - degreesOfFreedom) / sqrt(2 * degreesOfFreedom);
	printf("degrees of freedom = %d\n", degreesOfFreedom);
	fprintf(paramFile, "degrees of freedom = %d\n", degreesOfFreedom);
	printf("normal distribution = %g\n", normalDist);
	fprintf(paramFile, "normal distribution = %g\n", normalDist);
	printf("parameters at first time moment:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("sigma R = %g\n", uncertainties[0]);
	fprintf(paramFile, "sigma R = %g\n", uncertainties[0]);
	printf("R/ct = %g\n", vector[0] * maxParameters[0] / (speed_of_light * timeMoment));
	fprintf(paramFile, "R/ct = %g\n", vector[0] * maxParameters[0] / (speed_of_light * timeMoment));
	printf("sigma R/ct = %g", uncertainties[0] / (speed_of_light * timeMoment));
	fprintf(paramFile, "sigma R/ct = %g", uncertainties[0] / (speed_of_light * timeMoment));
	sigma = vector[1] * maxParameters[1];
	printf("sigma = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "sigma = %g\n", vector[1] * maxParameters[1]);
	printf("sigma sigma = %g\n", uncertainties[1]);
	fprintf(paramFile, "sigma sigma = %g\n", uncertainties[1]);
	double n = vector[2] * maxParameters[2];
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("sigma n = %g\n", uncertainties[2]);
	fprintf(paramFile, "sigma n = %g\n", uncertainties[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	printf("sigma width fraction = %g\n", uncertainties[3]);
	fprintf(paramFile, "sigma width fraction = %g\n", uncertainties[3]);
	printf("v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	fprintf(paramFile, "v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	printf("sigma v/c = %g\n", uncertainties[5] / speed_of_light);
	fprintf(paramFile, "sigma v/c = %g\n", uncertainties[5] / speed_of_light);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	double sigmaB = 0.5 * B * (uncertainties[1] / sigma + uncertainties[2] / n);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	printf("sigma B = %g\n", sigmaB);
	fprintf(paramFile, "sigma B = %g\n", sigmaB);
	fclose(paramFile);

	printf("evaluating wide-range radiation\n");
	printLog("evaluating wide-range radiation\n");


	Nrho = 1;
	Nz = 4000;
	Nphi = 1;

	double R = vector[0] * maxParameters[0];
	//B
	electronConcentration = vector[2] * maxParameters[2];
	double f = vector[3] * maxParameters[3];
	double velocity = vector[4] * maxParameters[4];
	double downstreamV = 0.25 * velocity;

	double totalKineticEnergy = source1->getTotalVolume() * massProton * electronConcentration * 0.5 * speed_of_light2;
	printf("totalKineticEnergy = %g\n", totalKineticEnergy);
	printLog("totalKineticEnergy = %g\n", totalKineticEnergy);

	//return;

	//TabulatedSLSourceWithSynchCutoff* source2 = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance, downstreamV, velocity);
	//TabulatedDiskSource* source2 = new TabulatedDiskSource(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance);
	TabulatedDiskSourceWithSynchAndComptCutoff* source3 = new TabulatedDiskSourceWithSynchAndComptCutoff(1, Nz, Nphi, electronDistribution, B, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, 0, velocity);

	int Ne = 1000;
	Emin = me_c2;
	Emax = me_c2 * 2E9;

	SynchrotronEvaluator* evaluator2 = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);

	//double kevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source2);

	//double mevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source2);

	double kevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source3);

	double mevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source3);

	//printf("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);
	//printLog("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);

	//printf("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);
	//printLog("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);

	printf("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);
	printLog("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);

	printf("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);
	printLog("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);

	//evaluator2->writeFluxFromSourceToFile("wideRangeSynch.dat", source2, 1E8 * hplank, 200 * MeV, 500);
	evaluator2->writeFluxFromSourceToFile("wideRangeSynchDisk.dat", source3, 1E8 * hplank, 200 * MeV, 500);

	//deleting arrays

	delete[] Fout;
	delete[] Nuout;

	delete synchrotronEvaluator;
	//delete gridEnumOptimizer;
	//delete gradientOptimizer;

}

void fitCSS161010_2() {
	//observed data at 99, 162 and 357 days after explosion in units erg and cm^-2 s^-2
	double* energy1;
	double* observedFlux1;
	double* observedError1;
	//int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans69.txt");
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans99.txt");
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank * 1E9;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}

	printf("finish reading data\n");
	printLog("finish reading data\n");



	double timeMoment = 99 * 24 * 3600;


	//distance to source
	const double distance = 150 * 1E6 * parsec;

	//initial parameters of source
	double electronConcentration = 50;
	double rmax = timeMoment * 0.75 * speed_of_light;
	//rmax = 1.0;
	double B = 0.6;
	double widthFraction = 0.25;
	double v = 0.75 * speed_of_light;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	//number of optimized parameters
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.6 * rmax, 0.000001, 1, 0.05, 0.5 * speed_of_light};
	double maxParameters[Nparams] = { timeMoment * 0.8 * speed_of_light, 0.05, 200, 0.2, 0.8 * speed_of_light};
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, widthFraction, v};
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, true, false};

	int numberOfOptpar = 0;
	for (int i = 0; i < Nparams; ++i) {
		if (optPar[i]) {
			numberOfOptpar++;
		}
	}

	int numberOfPoints = Nenergy1;


	printf("creating distributions\n");
	printLog("creating distributions\n");

	//number of different distributions depending on inclination angle, wich will be read from files
	const int Ndistributions = 10;

	//reading electron distributions from files
	//ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, 200);
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, 200, 20 * me_c2, 3.5);
	//for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
	//	(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	//}
	//angleDependentDistributions[9]->writeDistribution("output1.dat", 200, Emin, Emax);
	//creating radiation source, which does not depend on time
	//MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_theta0-90/Ee3.dat", "./examples_data/gamma1.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);

	printf("distribution 1\n");
	printLog("distribution 1\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution1 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_theta0-90/Ee3.dat", "./examples_data/gamma0.2_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 2\n");
	printLog("distribution 2\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution2 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.3_theta0-90/Ee3.dat", "./examples_data/gamma0.3_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 3\n");
	printLog("distribution 3\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution3 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 4\n");
	printLog("distribution 4\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution4 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	electronDistribution1->addPowerLaw(10 * me_c2, 3.5);
	electronDistribution1->rescaleDistribution(1.2);
	electronDistribution2->addPowerLaw(20 * me_c2, 3.5);
	electronDistribution2->rescaleDistribution(1.2);
	electronDistribution3->addPowerLaw(50 * me_c2, 3.5);
	electronDistribution3->rescaleDistribution(1.2);

	double velocities[4] = { 0.2 * speed_of_light, 0.3 * speed_of_light, 0.5 * speed_of_light, 0.75 * speed_of_light };

	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution * [4];
	distributions[0] = electronDistribution1;
	distributions[1] = electronDistribution2;
	distributions[2] = electronDistribution3;
	distributions[3] = electronDistribution4;

	//electronDistribution1->writeDistribution("distribution1.dat", 1000, me_c2, 1E9*me_c2);
	//electronDistribution2->writeDistribution("distribution2.dat", 1000, me_c2, 1E9 * me_c2);
	//electronDistribution3->writeDistribution("distribution3.dat", 1000, me_c2, 1E9 * me_c2);
	electronDistribution4->writeDistribution("distribution.dat", 1000, me_c2, 1E9 * me_c2);


	printf("creating sources\n");
	printLog("creating sources\n");

	//MassiveParticleIsotropicDistribution* electronDistribution = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 10*me_c2, electronConcentration);
	//electronDistribution->rescaleDistribution(1.2);
	//SimpleFlatSource* source1 = new SimpleFlatSource(electronDistribution, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);
	//SimpleFlatSource2* source1 = new SimpleFlatSource2(4, velocities, distributions, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);
	int Nrho = 200;
	int Nz = 400;
	int Nphi = 1;
	SphericalLayerSource* source1 = new TabulatedSphericalLayerSource2(4, velocities, distributions, Nrho, Nz, Nphi, B, pi / 2, 0, electronConcentration, rmax, (1.0 - widthFraction) * rmax, distance);
	//SimpleFlatSource2* source1 = new SimpleFlatSource2(4, velocities, distributions, B, pi / 2, 0, electronConcentration, rmax, widthFraction * rmax, distance);
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, pi / 2, 0, electronConcentration, rmax, 0.5 * rmax, distance, 0.3 * speed_of_light);

	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 4,4,4,4,4};
	//number of iterations in gradient descent optimizer
	int Niterations = 4;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating KPI evaluator
	LossEvaluator* KPIevaluator = new SpectrumLossEvaluator(energy1, observedFlux1, observedError1, Nenergy1, source1);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(50, Emin, Emax, true, false);
	CombinedRadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints, KPIevaluator);
	SequentCoordinateEnumOptimizer* sequentOptimizer = new SequentCoordinateEnumOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 200, 2, KPIevaluator);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	//RadiationOptimizer* gridEnumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	//gridEnumOptimizer->optimize(vector, optPar, energy, F, Error, Nenergy, Ntimes, times, source);
	vector[0] = 2.02253e+17 / maxParameters[0];
	vector[1] = 0.0436792 / maxParameters[1];
	vector[2] = 65.199 / maxParameters[2];
	vector[3] = 0.104334 / maxParameters[3];
	vector[4] = 0.755575 * speed_of_light / maxParameters[4];

	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 0, "error0.dat");
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 1, "error1.dat");*/
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 2, "error2.dat");
	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 3, "error3.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 4, "error4.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 5, "error5.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 6, "error6.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 7, "error7.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 8, "error8.dat");*/

	printf("optimization\n");
	printLog("optimization\n");
	//creating gradient descent optimizer and optimizing
	RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//combinedOptimizer->optimize(vector, optPar);
	sequentOptimizer->optimize(vector, optPar);
	gradientOptimizer->optimize(vector, optPar);
	//reset parameters of source to the found values
	source1->resetParameters(vector, maxParameters);
	//evaluating final error
	double error = sequentOptimizer->evaluateOptimizationFunction(vector);
	//double error = 100000;

	double* uncertainties = sequentOptimizer->evaluateUncertainties(vector, optPar);

	//initializing arrays for evaluationg full spectrum of source with found values
	printf("evaluationg radiation\n");
	printLog("evaluationg radiation\n");
	int Nout = 200;
	double* Nuout = new double[Nout];
	double* Fout = new double [Nout];

	double Numin = 1E8;
	double Numax = 1E12;
	double factor = pow(Numax / Numin, 1.0 / (Nout - 1));
	Nuout[0] = Numin;
	for (int j = 0; j < Nout; ++j) {
		Fout[j] = 0;
	}

	for (int i = 1; i < Nout; ++i) {
		Nuout[i] = Nuout[i - 1] * factor;
	}

	//evaluating full spectrum at given time moments
		for (int i = 0; i < Nout; ++i) {
			Fout[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nuout[i], source1);
		}


	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("css161010.dat", "w");
	for (int i = 0; i < Nout; ++i) {
		//to GHz and mJansky
		fprintf(output_GZ_Jansky, "%g", Nuout[i] * 1E-9);
	    fprintf(output_GZ_Jansky, " %g", hplank * Fout[i] * 1E26);
		fprintf(output_GZ_Jansky, "\n");
	}
	fclose(output_GZ_Jansky);

	//outputing parameters
	FILE* paramFile = fopen("parametersCSS161010.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	int degreesOfFreedom = numberOfPoints - numberOfOptpar - 1;
	printf("relative hi^2 = %g\n", error / degreesOfFreedom);
	fprintf(paramFile, "relative hi^2 = %g\n", error / degreesOfFreedom);
	double normalDist = (error - degreesOfFreedom) / sqrt(2 * degreesOfFreedom);
	printf("degrees of freedom = %d\n", degreesOfFreedom);
	fprintf(paramFile, "degrees of freedom = %d\n", degreesOfFreedom);
	printf("normal distribution = %g\n", normalDist);
	fprintf(paramFile, "normal distribution = %g\n", normalDist);
	printf("parameters at first time moment:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("sigma R = %g\n", uncertainties[0]);
	fprintf(paramFile,"sigma R = %g\n", uncertainties[0]);
	printf("R/ct = %g\n", vector[0] * maxParameters[0] / (speed_of_light * timeMoment));
	fprintf(paramFile, "R/ct = %g\n", vector[0] * maxParameters[0] / (speed_of_light * timeMoment));
	printf("sigma R/ct = %g", uncertainties[0] / (speed_of_light * timeMoment));
	fprintf(paramFile,"sigma R/ct = %g", uncertainties[0] / (speed_of_light * timeMoment));
	sigma = vector[1] * maxParameters[1];
	printf("sigma = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "sigma = %g\n", vector[1] * maxParameters[1]);
	printf("sigma sigma = %g\n", uncertainties[1]);
	fprintf(paramFile, "sigma sigma = %g\n", uncertainties[1]);
	double n = vector[2] * maxParameters[2];
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("sigma n = %g\n", uncertainties[2]);
	fprintf(paramFile, "sigma n = %g\n", uncertainties[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	printf("sigma width fraction = %g\n", uncertainties[3]);
	fprintf(paramFile, "sigma width fraction = %g\n", uncertainties[3]);
	printf("v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	fprintf(paramFile, "v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	printf("sigma v/c = %g\n", uncertainties[5] / speed_of_light);
	fprintf(paramFile, "sigma v/c = %g\n", uncertainties[5] / speed_of_light);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	double sigmaB = 0.5 * B * (uncertainties[1] / sigma + uncertainties[2] / n);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	printf("sigma B = %g\n", sigmaB);
	fprintf(paramFile, "sigma B = %g\n", sigmaB);
	fclose(paramFile);

	printf("evaluating wide-range radiation\n");
	printLog("evaluating wide-range radiation\n");


	Nrho = 2000;
	Nz = 4000;
	Nphi = 1;

	double R = vector[0] * maxParameters[0];
	//B
	electronConcentration = vector[2] * maxParameters[2];
	double f = vector[3] * maxParameters[3];
	double velocity = vector[4] * maxParameters[4];
	double downstreamV = 0.25 * velocity;

    double totalKineticEnergy = source1->getTotalVolume() * massProton * electronConcentration * 0.5 * speed_of_light2;
	printf("totalKineticEnergy = %g\n", totalKineticEnergy);
	printLog("totalKineticEnergy = %g\n", totalKineticEnergy);

    //return;


	TabulatedSLSourceWithSynchCutoff* source2 = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance, downstreamV, velocity);
	//TabulatedDiskSource* source2 = new TabulatedDiskSource(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance);
	TabulatedDiskSourceWithSynchAndComptCutoff* source3 = new TabulatedDiskSourceWithSynchAndComptCutoff(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, velocity);

	int Ne = 1000;
	Emin = me_c2;
	Emax = me_c2 * 2E9;

	electronDistribution4->writeDistribution("distribution.dat", Ne, Emin, Emax);

	MassiveParticleIsotropicDistribution* distribution1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source2->getParticleDistribution(Nrho - 1, Nz / 2, 0));
	distribution1->writeDistribution("distribution1.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution2 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source3->getParticleDistribution(0, Nz - 1, 0));
	distribution2->writeDistribution("distribution2.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution3 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source2->getParticleDistribution(0, Nz - 10, 0));
	distribution3->writeDistribution("distribution3.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution4 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source3->getParticleDistribution(0, Nz - 10, 0));
	distribution4->writeDistribution("distribution4.dat", Ne, Emin, Emax);

	SynchrotronEvaluator* evaluator2 = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);

	double kevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source2);

	double mevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source2);

	double kevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source3);

	double mevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source3);

	printf("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);
	printLog("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);

	printf("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);
	printLog("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);

	printf("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);
	printLog("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);

	printf("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);
	printLog("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);

	evaluator2->writeFluxFromSourceToFile("wideRangeSynch.dat", source2, 1E8 * hplank, 200 * MeV, 500);
	evaluator2->writeFluxFromSourceToFile("wideRangeSynchDisk.dat", source3, 1E8 * hplank, 200 * MeV, 500);

	//deleting arrays

	delete[] Fout;
	delete[] Nuout;

	for (int i = 0; i < Ndistributions; ++i) {
		delete angleDependentDistributions[i];
	}
	delete[] angleDependentDistributions;
	delete angleDependentSource;
	delete synchrotronEvaluator;
	//delete gridEnumOptimizer;
	//delete gradientOptimizer;

}

void fitAT2025wpp_2() {
	//observed data at 19, 35 and 50 days after explosion in units erg and cm^-2 s^-2
	double* energy1;
	double* observedFlux1;
	double* observedError1;
	//int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/css_data/coppejans69.txt");
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/AT2025wpp_data/nayana78.dat");
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank * 1E9;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}

	printf("finish reading data\n");
	printLog("finish reading data\n");



	double timeMoment = 78 * 24 * 3600;


	//distance to source
	const double distance = 411 * 1E6 * parsec;

	//initial parameters of source
	double electronConcentration = 50;
	double rmax = timeMoment * 0.75 * speed_of_light;
	//rmax = 1.0;
	double B = 0.6;
	double widthFraction = 0.5;
	double v = 0.75 * speed_of_light;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	//number of optimized parameters
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.1 * rmax, 0.000001, 1, 0.5, 0.1 * speed_of_light };
	double maxParameters[Nparams] = { timeMoment * 0.8 * speed_of_light, 1.0, 5000, 0.5, 0.8 * speed_of_light };
	//starting point of optimization and normalization
	double vector[Nparams] = { rmax, sigma, electronConcentration, widthFraction, v };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	//picking parameters to be optimized
	bool optPar[Nparams] = { true, true, true, false, false };

	int numberOfOptpar = 0;
	for (int i = 0; i < Nparams; ++i) {
		if (optPar[i]) {
			numberOfOptpar++;
		}
	}

	int numberOfPoints = Nenergy1;


	printf("creating distributions\n");
	printLog("creating distributions\n");

	//number of different distributions depending on inclination angle, wich will be read from files
	const int Ndistributions = 10;

	//reading electron distributions from files
	//ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs", ".dat", 10, ElectronInputType::GAMMA_KIN_FGAMMA, 200);
	MassiveParticleDistribution** angleDependentDistributions = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionsAddPowerLawTail(massElectron, "./input/Ee", "./input/Fs", ".dat", 10, DistributionInputType::GAMMA_KIN_FGAMMA, 200, 20 * me_c2, 3.5);
	//for (int i = 0; i < Ndistributions; ++i) {
		//rescale distributions to real mp/me relation
	//	(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	//}
	//angleDependentDistributions[9]->writeDistribution("output1.dat", 200, Emin, Emax);
	//creating radiation source, which does not depend on time
	//MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_theta0-90/Ee3.dat", "./examples_data/gamma1.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);

	printf("distribution 1\n");
	printLog("distribution 1\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution1 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_theta0-90/Ee3.dat", "./examples_data/gamma0.2_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 2\n");
	printLog("distribution 2\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution2 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.3_theta0-90/Ee3.dat", "./examples_data/gamma0.3_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 3\n");
	printLog("distribution 3\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution3 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	printf("distribution 4\n");
	printLog("distribution 4\n");
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution4 = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	electronDistribution1->addPowerLaw(10 * me_c2, 3.5);
	electronDistribution1->rescaleDistribution(1.2);
	electronDistribution2->addPowerLaw(20 * me_c2, 3.5);
	electronDistribution2->rescaleDistribution(1.2);
	electronDistribution3->addPowerLaw(50 * me_c2, 3.5);
	electronDistribution3->rescaleDistribution(1.2);

	double velocities[4] = { 0.2 * speed_of_light, 0.3 * speed_of_light, 0.5 * speed_of_light, 0.75 * speed_of_light };

	MassiveParticleIsotropicDistribution** distributions = new MassiveParticleIsotropicDistribution * [4];
	distributions[0] = electronDistribution1;
	distributions[1] = electronDistribution2;
	distributions[2] = electronDistribution3;
	distributions[3] = electronDistribution4;

	//electronDistribution1->writeDistribution("distribution1.dat", 1000, me_c2, 1E9*me_c2);
	//electronDistribution2->writeDistribution("distribution2.dat", 1000, me_c2, 1E9 * me_c2);
	//electronDistribution3->writeDistribution("distribution3.dat", 1000, me_c2, 1E9 * me_c2);
	electronDistribution4->writeDistribution("distribution.dat", 1000, me_c2, 1E9 * me_c2);


	printf("creating sources\n");
	printLog("creating sources\n");

	//MassiveParticleIsotropicDistribution* electronDistribution = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 10*me_c2, electronConcentration);
	//electronDistribution->rescaleDistribution(1.2);
	//SimpleFlatSource* source1 = new SimpleFlatSource(electronDistribution, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);
	//SimpleFlatSource2* source1 = new SimpleFlatSource2(4, velocities, distributions, B, pi/2, 0, electronConcentration, rmax, widthFraction*rmax, distance);
	int Nrho = 200;
	int Nz = 400;
	int Nphi = 1;
	SphericalLayerSource* source1 = new TabulatedSphericalLayerSource2(4, velocities, distributions, Nrho, Nz, Nphi, B, pi / 2, 0, electronConcentration, rmax, (1.0 - widthFraction) * rmax, distance);
	//SimpleFlatSource2* source1 = new SimpleFlatSource2(4, velocities, distributions, B, pi / 2, 0, electronConcentration, rmax, widthFraction * rmax, distance);
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, pi / 2, 0, electronConcentration, rmax, 0.5 * rmax, distance, 0.3 * speed_of_light);

	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 4,4,4,4,4 };
	//number of iterations in gradient descent optimizer
	int Niterations = 4;
	//energies of electrons wich will be used for evaluatig radiation
	double Emin = me_c2;
	double Emax = 10000 * me_c2;
	//creating KPI evaluator
	LossEvaluator* KPIevaluator = new SpectrumLossEvaluator(energy1, observedFlux1, observedError1, Nenergy1, source1);
	//creating time dependent synchrotron evaluator
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(50, Emin, Emax, true, false);
	CombinedRadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints, KPIevaluator);
	SequentCoordinateEnumOptimizer* sequentOptimizer = new SequentCoordinateEnumOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 200, 2, KPIevaluator);
	//creating time depedent grid enumeration optimizer, which will chose the best starting poin for gradien descent
	//RadiationOptimizer* gridEnumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	//gridEnumOptimizer->optimize(vector, optPar, energy, F, Error, Nenergy, Ntimes, times, source);
	vector[0] = 3e+16 / maxParameters[0];
	vector[1] = 0.0436792 / maxParameters[1];
	vector[2] = 2500 / maxParameters[2];
	vector[3] = 0.1 / maxParameters[3];
	vector[4] = 0.2 * speed_of_light / maxParameters[4];

	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 0, "error0.dat");
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 1, "error1.dat");*/
	//combinedOptimizer->outputOneVariableProfile(vector, 100, 2, "error2.dat");
	/*combinedOptimizer->outputOneVariableProfile(vector, 100, 3, "error3.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 4, "error4.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 5, "error5.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 6, "error6.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 7, "error7.dat");
	combinedOptimizer->outputOneVariableProfile(vector, 100, 8, "error8.dat");*/

	printf("optimization\n");
	printLog("optimization\n");
	//creating gradient descent optimizer and optimizing
	RadiationOptimizer* gradientOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//combinedOptimizer->optimize(vector, optPar);
	sequentOptimizer->optimize(vector, optPar);
	gradientOptimizer->optimize(vector, optPar);
	//reset parameters of source to the found values
	source1->resetParameters(vector, maxParameters);
	//evaluating final error
	double error = sequentOptimizer->evaluateOptimizationFunction(vector);
	//double error = 100000;

	double* uncertainties = sequentOptimizer->evaluateUncertainties(vector, optPar);

	//initializing arrays for evaluationg full spectrum of source with found values
	printf("evaluationg radiation\n");
	printLog("evaluationg radiation\n");
	int Nout = 200;
	double* Nuout = new double[Nout];
	double* Fout = new double[Nout];

	double Numin = 1E8;
	double Numax = 1E12;
	double factor = pow(Numax / Numin, 1.0 / (Nout - 1));
	Nuout[0] = Numin;
	for (int j = 0; j < Nout; ++j) {
		Fout[j] = 0;
	}

	for (int i = 1; i < Nout; ++i) {
		Nuout[i] = Nuout[i - 1] * factor;
	}

	//evaluating full spectrum at given time moments
	for (int i = 0; i < Nout; ++i) {
		Fout[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nuout[i], source1);
	}


	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("AT2025wpp.dat", "w");
	for (int i = 0; i < Nout; ++i) {
		//to GHz and mJansky
		fprintf(output_GZ_Jansky, "%g", Nuout[i] * 1E-9);
		fprintf(output_GZ_Jansky, " %g", hplank * Fout[i] * 1E26);
		fprintf(output_GZ_Jansky, "\n");
	}
	fclose(output_GZ_Jansky);

	//outputing parameters
	FILE* paramFile = fopen("parametersAT2025wpp.dat", "w");
	printf("hi^2 = %g\n", error);
	fprintf(paramFile, "hi^2 = %g\n", error);
	int degreesOfFreedom = numberOfPoints - numberOfOptpar - 1;
	printf("relative hi^2 = %g\n", error / degreesOfFreedom);
	fprintf(paramFile, "relative hi^2 = %g\n", error / degreesOfFreedom);
	double normalDist = (error - degreesOfFreedom) / sqrt(2 * degreesOfFreedom);
	printf("number of optimized parameters = %d\n", numberOfOptpar);
	fprintf(paramFile, "number of optimized parameters = %d\n", numberOfOptpar);
	printf("degrees of freedom = %d\n", degreesOfFreedom);
	fprintf(paramFile, "degrees of freedom = %d\n", degreesOfFreedom);
	printf("normal distribution = %g\n", normalDist);
	fprintf(paramFile, "normal distribution = %g\n", normalDist);
	printf("parameters at first time moment:\n");
	fprintf(paramFile, "parameters at first time moment:\n");
	printf("R = %g\n", vector[0] * maxParameters[0]);
	fprintf(paramFile, "R = %g\n", vector[0] * maxParameters[0]);
	printf("sigma R = %g\n", uncertainties[0]);
	fprintf(paramFile, "sigma R = %g\n", uncertainties[0]);
	printf("R/ct = %g\n", vector[0] * maxParameters[0] / (speed_of_light * timeMoment));
	fprintf(paramFile, "R/ct = %g\n", vector[0] * maxParameters[0] / (speed_of_light * timeMoment));
	printf("sigma R/ct = %g", uncertainties[0] / (speed_of_light * timeMoment));
	fprintf(paramFile, "sigma R/ct = %g", uncertainties[0] / (speed_of_light * timeMoment));
	sigma = vector[1] * maxParameters[1];
	printf("sigma = %g\n", vector[1] * maxParameters[1]);
	fprintf(paramFile, "sigma = %g\n", vector[1] * maxParameters[1]);
	printf("sigma sigma = %g\n", uncertainties[1]);
	fprintf(paramFile, "sigma sigma = %g\n", uncertainties[1]);
	double n = vector[2] * maxParameters[2];
	printf("n = %g\n", vector[2] * maxParameters[2]);
	fprintf(paramFile, "n = %g\n", vector[2] * maxParameters[2]);
	printf("sigma n = %g\n", uncertainties[2]);
	fprintf(paramFile, "sigma n = %g\n", uncertainties[2]);
	printf("width fraction = %g\n", vector[3] * maxParameters[3]);
	fprintf(paramFile, "width fraction = %g\n", vector[3] * maxParameters[3]);
	printf("sigma width fraction = %g\n", uncertainties[3]);
	fprintf(paramFile, "sigma width fraction = %g\n", uncertainties[3]);
	printf("v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	fprintf(paramFile, "v/c = %g\n", vector[4] * maxParameters[4] / speed_of_light);
	printf("sigma v/c = %g\n", uncertainties[5] / speed_of_light);
	fprintf(paramFile, "sigma v/c = %g\n", uncertainties[5] / speed_of_light);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	double sigmaB = 0.5 * B * (uncertainties[1] / sigma + uncertainties[2] / n);
	printf("B = %g\n", B);
	fprintf(paramFile, "B = %g\n", B);
	printf("sigma B = %g\n", sigmaB);
	fprintf(paramFile, "sigma B = %g\n", sigmaB);
	fclose(paramFile);

	printf("evaluating wide-range radiation\n");
	printLog("evaluating wide-range radiation\n");


	Nrho = 2000;
	Nz = 4000;
	Nphi = 1;

	double R = vector[0] * maxParameters[0];
	//B
	electronConcentration = vector[2] * maxParameters[2];
	double f = vector[3] * maxParameters[3];
	double velocity = vector[4] * maxParameters[4];
	double downstreamV = 0.25 * velocity;

	double totalKineticEnergy = source1->getTotalVolume() * massProton * electronConcentration * 0.5 * speed_of_light2;
	printf("totalKineticEnergy = %g\n", totalKineticEnergy);
	printLog("totalKineticEnergy = %g\n", totalKineticEnergy);

	return;


	TabulatedSLSourceWithSynchCutoff* source2 = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance, downstreamV, velocity);
	//TabulatedDiskSource* source2 = new TabulatedDiskSource(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance);
	TabulatedDiskSourceWithSynchAndComptCutoff* source3 = new TabulatedDiskSourceWithSynchAndComptCutoff(1, Nz, Nphi, electronDistribution4, B, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, velocity);

	int Ne = 1000;
	Emin = me_c2;
	Emax = me_c2 * 2E9;

	electronDistribution4->writeDistribution("distribution.dat", Ne, Emin, Emax);

	MassiveParticleIsotropicDistribution* distribution1 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source2->getParticleDistribution(Nrho - 1, Nz / 2, 0));
	distribution1->writeDistribution("distribution1.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution2 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source3->getParticleDistribution(0, Nz - 1, 0));
	distribution2->writeDistribution("distribution2.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution3 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source2->getParticleDistribution(0, Nz - 10, 0));
	distribution3->writeDistribution("distribution3.dat", Ne, Emin, Emax);
	MassiveParticleIsotropicDistribution* distribution4 = dynamic_cast<MassiveParticleIsotropicDistribution*>(source3->getParticleDistribution(0, Nz - 10, 0));
	distribution4->writeDistribution("distribution4.dat", Ne, Emin, Emax);

	SynchrotronEvaluator* evaluator2 = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);

	double kevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source2);

	double mevFlux = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source2);

	double kevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source3);

	double mevFluxDisk = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source3);

	printf("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);
	printLog("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);

	printf("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);
	printLog("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);

	printf("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);
	printLog("keV flux disk= %g, luminosity = %g\n", kevFluxDisk, kevFluxDisk * 4 * pi * distance * distance);

	printf("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);
	printLog("MeV flux disk= %g, luminocity = %g\n", mevFluxDisk, mevFluxDisk * 4 * pi * distance * distance);

	evaluator2->writeFluxFromSourceToFile("wideRangeSynch.dat", source2, 1E8 * hplank, 200 * MeV, 500);
	evaluator2->writeFluxFromSourceToFile("wideRangeSynchDisk.dat", source3, 1E8 * hplank, 200 * MeV, 500);

	//deleting arrays

	delete[] Fout;
	delete[] Nuout;

	for (int i = 0; i < Ndistributions; ++i) {
		delete angleDependentDistributions[i];
	}
	delete[] angleDependentDistributions;
	delete angleDependentSource;
	delete synchrotronEvaluator;
	//delete gridEnumOptimizer;
	//delete gradientOptimizer;

}

void testMatrixInverse()
{
	int N = 2;

	printf("2x2 matrix\n");

	double** matrix1 = new double* [N];
	for (int i = 0; i < N; ++i) {
		matrix1[i] = new double[N];
		for (int j = 0; j < N; ++j) {
			matrix1[i][j] = 0;
		}
	}

	matrix1[0][0] = 1;
	matrix1[0][1] = 2;
	matrix1[1][0] = 3;
	matrix1[1][1] = 4;

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			printf("%g ", matrix1[i][j]);
		}
		printf("\n");
	}

	printf("\n");

	double** matrix1_inverse = inverseMatrix(matrix1, N);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			printf("%g ", matrix1_inverse[i][j]);
		}
		printf("\n");
	}

	printf("4x4 matrix\n");

	N = 4;

	double** matrix4 = new double* [N];
	for (int i = 0; i < N; ++i) {
		matrix4[i] = new double[N];
		for (int j = 0; j < N; ++j) {
			matrix4[i][j] = 0;
		}
	}

	matrix4[0][0] = 1;
	matrix4[0][1] = 2;
	matrix4[0][2] = 3;
	matrix4[0][3] = 4;
	matrix4[1][0] = 2;
	matrix4[1][1] = 1;
	matrix4[1][2] = 0;
	matrix4[1][3] = 1;
	matrix4[2][0] = 1;
	matrix4[2][1] = 0;
	matrix4[2][2] = 3;
	matrix4[2][3] = 2;
	matrix4[3][0] = 5;
	matrix4[3][1] = 1;
	matrix4[3][2] = 1;
	matrix4[3][3] = 1;

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			printf("%g ", matrix4[i][j]);
		}
		printf("\n");
	}

	printf("\n");

	double** matrix4_inverse = inverseMatrix(matrix4, N);

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			printf("%g ", matrix4_inverse[i][j]);
		}
		printf("\n");
	}
}

void testGMRES() {
	//with this example gmres does not converge until N iterations
	int N = 9;
	/*std::vector<MatrixElement>**** matrix = new std::vector<MatrixElement>***[N];
	for (int i = 0; i < N; ++i) {
		matrix[i] = new std::vector<MatrixElement>**[1];
		matrix[i][0] = new std::vector<MatrixElement>*[1];
		matrix[i][0][0] = new std::vector<MatrixElement>[1];
	}*/

	std::vector<MatrixElement>**** matrix = new std::vector<MatrixElement>***[1];
	matrix[0] = new std::vector<MatrixElement>**[1];
	matrix[0][0] = new std::vector<MatrixElement>*[1];
	matrix[0][0][0] = new std::vector<MatrixElement>[9];
	for (int i = 0; i < N; ++i) {
		//matrix[0][0][i] = new std::vector<MatrixElement>[1];
	}

	matrix[0][0][0][0].push_back(MatrixElement(1.0, 0, 0, 0, 0));
	matrix[0][0][0][0].push_back(MatrixElement(2.0, 0, 0, 0, 3));
	matrix[0][0][0][0].push_back(MatrixElement(5.0, 0, 0, 0, 6));

	matrix[0][0][0][1].push_back(MatrixElement(1.0, 0, 0, 0, 1));
	matrix[0][0][0][1].push_back(MatrixElement(4.0, 0, 0, 0, 4));

	matrix[0][0][0][2].push_back(MatrixElement(1.0, 0, 0, 0, 2));
	matrix[0][0][0][2].push_back(MatrixElement(-2.0, 0, 0, 0, 5));

	matrix[0][0][0][3].push_back(MatrixElement(1.0, 0, 0, 0, 3));

	matrix[0][0][0][4].push_back(MatrixElement(1.0, 0, 0, 0, 4));
	matrix[0][0][0][4].push_back(MatrixElement(7.0, 0, 0, 0, 2));

	matrix[0][0][0][5].push_back(MatrixElement(1.0, 0, 0, 0, 5));
	matrix[0][0][0][5].push_back(MatrixElement(1.0, 0, 0, 0, 1));

	matrix[0][0][0][6].push_back(MatrixElement(1.0, 0, 0, 0, 6));
	matrix[0][0][0][6].push_back(MatrixElement(1.0, 0, 0, 0, 0));

	matrix[0][0][0][7].push_back(MatrixElement(1.0, 0, 0, 0, 7));
	matrix[0][0][0][7].push_back(MatrixElement(1.0, 0, 0, 0, 5));

	matrix[0][0][0][8].push_back(MatrixElement(1.0, 0, 0, 0, 8));

	double**** R = new double***[1];
	R[0] = new double** [1];
	R[0][0] = new double* [1];
	R[0][0][0] = new double[N];
	for (int i = 0; i < N; ++i) {
		//R[i] = new double** [1];
		//R[0][i] = new double* [1];
		//R[0][0][i] = new double[1];
		R[0][0][0][i] = 0;
	}
	R[0][0][0][0] = 1.0;
	R[0][0][0][1] = 2.0;
	R[0][0][0][2] = 3.0;

	double**** X = new double***[1];
	X[0] = new double** [1];
	X[0][0] = new double* [1];
	X[0][0][0] = new double[N];
	for (int i = 0; i < N; ++i) {
		//X[i] = new double** [1];
		//X[0][i] = new double* [1];
		//X[0][0][i] = new double[1];
	}

	double**** initialVector = new double*** [1];
	initialVector[0] = new double** [1];
	initialVector[0][0] = new double* [1];
	initialVector[0][0][0] = new double[1];
	for (int i = 0; i < N; ++i) {
		initialVector[0][0][0][i] = 1.0;
	}

	LargeVectorBasis* gmresBasis = new LargeVectorBasis(5, 1, 1, 1, N);

	//for (int i = 0; i < 1000; ++i) {
		generalizedMinimalResidualMethod(matrix, R, X, initialVector, 1, 1, 1, N, 1E-5, N, 2, gmresBasis);
	//}

	for (int i = 0; i < N; ++i) {
		printf("%g\n", X[0][0][0][i]);
	}

	/*conjugateGradientMethod(matrix, R, X, N, 1E-10, N, 2);

	for (int i = 0; i < N; ++i) {
		printf("%g\n", X[i]);
	}

	gaussSeidelMethod(matrix, R, X, N, 1E-10, N, 2);

	for (int i = 0; i < N; ++i) {
		printf("%g\n", X[i]);
	}*/
}

void testNishinaLosses() {
	double distance = 1E18;
	double size = 1E16;
	double concentration = 1;
	double B = 1.0;
	double Ee = me_c2 * 1E12;
	double halfEe = me_c2 * 10;
	double gamma = Ee / me_c2;
	MassiveParticleMonoenergeticDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, Ee, halfEe);
	SimpleFlatSource* source = new SimpleFlatSource(electrons, B, pi / 2, 0, concentration, size, size, distance);
	double volume = source->getTotalVolume();

	int Neph = 50;
	double minEph = 0.01*kBoltzman*2.725;
	double maxEph = 1E10*minEph;
	double factor = pow(maxEph / minEph, 1.0 / (Neph - 1.0));
	double Eph = minEph;

	double* L = new double[Neph];

	double* L1 = new double[Neph];
	double* L2 = new double[Neph];
	double* L3 = new double[Neph];

	double r2 = sqr(electron_charge * electron_charge / me_c2);
	double sigmaT = 8 * pi * r2 / 3.0;

	for (int i = 0; i < Neph; ++i) {
		printf("i = %d\n", i);
		PhotonMonoenergeticDistribution* photons = new PhotonMonoenergeticDistribution(Eph, 0.01 * Eph);
		double photonConcentration = 1.0/Eph;
		double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();
		double dEph = 0.01 * Eph;
		InverseComptonEvaluator* evaluator = new InverseComptonEvaluator(100, 50, 4, Ee - halfEe, Ee + halfEe, 100, Eph - dEph, Eph+dEph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
		//InverseComptonEvaluator* evaluator1 = new InverseComptonEvaluator(100, 50, 4, Ee - halfEe, Ee + halfEe, 50, Eph - dEph, Eph + dEph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
		//InverseComptonEvaluator* evaluator2 = new InverseComptonEvaluator(100, 50, 4, Ee - halfEe, Ee + halfEe, 200, Eph - dEph, Eph + dEph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
		//InverseComptonEvaluator* evaluator3 = new InverseComptonEvaluator(10, 50, 4, Ee - halfEe, Ee + halfEe, 500, Eph - dEph, Eph + dEph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

		InverseComptonEvaluator* evaluator1 = new InverseComptonEvaluator(10, 10, 4, Ee - halfEe, Ee + halfEe, 10, Eph - dEph, Eph+dEph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);

		L[i] = evaluator->evaluateTotalFluxInEnergyRange(0.01 * Eph, 2 * Ee, 200, source) * 4 * pi * distance * distance;
		L1[i] = volume * concentration * (4.0 / 3.0) * speed_of_light * sigmaT * gamma * gamma * photonEnergyDensity / pow(1.0 + 4 * gamma * Eph / me_c2, 3.0 / 2.0);
		L3[i] = volume * concentration * (4.0 / 3.0) * speed_of_light * sigmaT * gamma * gamma * photonEnergyDensity / pow(1.0 + 4 * pi * gamma * Eph / me_c2, 3.0 / 2.0);

		L2[i] = evaluator1->evaluateTotalFluxInEnergyRange(0.01 * Eph, 2 * Ee, 200, source) * 4 * pi * distance * distance;
		//L2[i] = evaluator2->evaluateTotalFluxInEnergyRange(0.01 * Eph, 2 * Ee, 200, source) * 4 * pi * distance * distance;
		//L3[i] = evaluator3->evaluateTotalFluxInEnergyRange(0.01 * Eph, 2 * Ee, 200, source) * 4 * pi * distance * distance;

		Eph = Eph * factor;

		delete evaluator;
		//delete evaluator1;
		delete photons;
	}

	Eph = minEph;
	FILE* outFile = fopen("nishina_losses.dat", "w");
	for (int i = 0; i < Neph; ++i) {
		fprintf(outFile, "%g %g %g %g %g\n", Eph/1.6E-12 , L[i], L1[i], L2[i], L3[i]);
		Eph = Eph * factor;
	}
	fclose(outFile);
}

void testNishinaLosses2() {
	double distance = 1E18;
	double size = 1E16;
	double concentration = 1;
	double B = 1.0;

	int Nee = 70;
	double minEe = 1E5*me_c2;
	double maxEe = 1E12*me_c2;
	double factor = pow(maxEe / minEe, 1.0 / (Nee - 1.0));
	double E = minEe;

	double* L = new double[Nee];

	double* L1 = new double[Nee];
	double* L2 = new double[Nee];
	double* L3 = new double[Nee];

	double r2 = sqr(electron_charge * electron_charge / me_c2);
	double sigmaT = 8 * pi * r2 / 3.0;

	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();

	double Ee = minEe;
	for (int i = 0; i < Nee; ++i) {
		printf("i = %d\n", i);
		double halfEe = 0.01*Ee;
		double gamma = Ee / me_c2;
		MassiveParticleMonoenergeticDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, Ee, halfEe);
		SimpleFlatSource* source = new SimpleFlatSource(electrons, B, pi / 2, 0, concentration, size, size, distance);
		double volume = source->getTotalVolume();

		double Eph = 2.8 * kBoltzman * 2.7;
		
		InverseComptonEvaluator* evaluator = new InverseComptonEvaluator(10, 100, 4, Ee - halfEe, Ee + halfEe, 50, 0.1 * Eph, 10 * Eph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
		//InverseComptonEvaluator* evaluator2 = new InverseComptonEvaluator(10, 100, 4, Ee - halfEe, Ee + halfEe, 50, 0.1 * Eph, 10 * Eph, photons, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);

		L[i] = evaluator->evaluateTotalFluxInEnergyRange(0.1 * Eph, 2 * Ee, 2000, source) * 4 * pi * distance * distance;
		//L[i] = evaluator->evaluateTotalFluxInEnergyRangeOneParticle(0.1 * Eph, 2 * Ee, 20000, Ee, Eph) * 4 * pi * distance * distance;
		L1[i] = volume * concentration * (4.0 / 3.0) * speed_of_light * sigmaT * gamma * gamma * photonEnergyDensity / pow(1.0 + 4 * gamma * Eph / me_c2, 3.0 / 2.0);
		L2[i] = volume * concentration * (4.0 / 3.0) * speed_of_light * sigmaT * gamma * gamma * photonEnergyDensity / pow(1.0 + pow(4 * pi * gamma * Eph / me_c2, 0.6), 1.9/0.6);

		//L1[i] = evaluator2->evaluateTotalFluxInEnergyRange(0.1 * Eph, 2 * Ee, 200, source) * 4 * pi * distance * distance;

		Ee = Ee * factor;

		delete source;
		delete electrons;

		delete evaluator;
		//delete evaluator2;
	}

	Ee = minEe;
	FILE* outFile = fopen("nishina_losses2.dat", "w");
	for (int i = 0; i < Nee; ++i) {
		fprintf(outFile, "%g %g %g %g %g\n", Ee / 1.6E-12, L[i], L1[i], L2[i], L3[i]);
		Ee = Ee * factor;
	}
	fclose(outFile);
}

void testNishinaSpectrum() {
	double distance = 1E18;
	double size = 1E16;
	double concentration = 1;
	double B = 1.0;

	int Nee = 6;
	double minEe = 1E5 * me_c2;
	double maxEe = 511*1E12*1.6E-12;
	double factor = pow(maxEe / minEe, 1.0 / (Nee - 1.0));
	double E = minEe;

	double* Ee = new double[Nee];
	for (int i = 0; i < Nee; ++i) {
		Ee[i] = E;
		E = E * factor;
	}
	Ee[0] = 1E8 * me_c2;
	Ee[1] = 1E9 * me_c2;
	Ee[2] = 1E10 * me_c2;
	Ee[3] = 1E11 * me_c2;
	Ee[4] = 1E12 * me_c2;
	Ee[5] = 1E13 * me_c2;

	int Nph = 1000;
	double minEph = 10000 * kBoltzman * 2.7;
	double maxEph = 1.5*Ee[Nee-1];
	factor = pow(maxEph / minEph, 1.0 / (Nph - 1.0));

	double* Eph = new double[Nph];
	E = minEph;
	for (int i = 0; i < Nph; ++i) {
		Eph[i] = E;
		E = E * factor;
	}

	double** L = new double* [Nee];
	double** L2 = new double* [Nee];

	for (int i = 0; i < Nee; ++i) {
		L[i] = new double[Nph];
		L2[i] = new double[Nph];
		for (int j = 0; j < Nph; ++j) {
			L[i][j] = 0;
			L2[i][j] = 0;
		}
	}

	double r2 = sqr(electron_charge * electron_charge / me_c2);
	double sigmaT = 8 * pi * r2 / 3.0;

	PhotonMultiPlankDistribution* photons = PhotonMultiPlankDistribution::getGalacticField();
	//PhotonPlankDistribution* photons = new PhotonPlankDistribution(10000, 1.0);
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();

	RadiationSource** sources = new RadiationSource * [Nee];
	RadiationEvaluator** evaluators1 = new RadiationEvaluator * [Nee];
	RadiationEvaluator** evaluators2 = new RadiationEvaluator * [Nee];

	for (int i = 0; i < Nee; ++i) {
		double halfEe = 0.01 * Ee[i];
		double gamma = Ee[i] / me_c2;
		MassiveParticleMonoenergeticDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, Ee[i], halfEe);
		sources[i] = new SimpleFlatSource(electrons, B, pi / 2, 0, concentration, size, size, distance);
		evaluators1[i] = new InverseComptonEvaluator(10, 100, 4, Ee[i] - halfEe, Ee[i] + halfEe, 50, 0.1 * kBoltzman*2.7, 2*kBoltzman*5000, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
		evaluators2[i] = new InverseComptonEvaluator(10, 100, 4, Ee[i] - halfEe, Ee[i] + halfEe, 50, 0.1 * kBoltzman*2.7, 2*kBoltzman*5000, photons, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	}

	for (int i = 0; i < Nph; ++i) {
		printf("i = %d\n", i);
		
		for (int j = 0; j < Nee; ++j) {
			L[j][i] = Eph[i]*evaluators1[j]->evaluateFluxFromSource(Eph[i], sources[j]);			L2[j][i] = Eph[i]*evaluators2[j]->evaluateFluxFromSource(Eph[i], sources[j]);
		}
	}

	FILE* outFile = fopen("nishina_spectrum.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		fprintf(outFile, "%g", Eph[i] / 1.6E-12);
		for (int j = 0; j < Nee; ++j) {
			fprintf(outFile, " %g", L[j][i]);
		}
		for (int j = 0; j < Nee; ++j) {
			fprintf(outFile, " %g", L2[j][i]);
		}
		fprintf(outFile, "\n");
	}
	fclose(outFile);
}

void testBigSource() {
	double distance = (18000 / 3.26) * parsec;

	double* energy;
	double* xgrid;
	double* concentration;
	double** distributions;

	int Nenergy;
	int Nx;

	const char* xfileName = "./examples_data/W50/newdistribution/x_grid.dat";

	Nx = 0;
	FILE* xfile = fopen(xfileName, "r");
	while (!feof(xfile)) {
		double a;
		fscanf(xfile, "%lf", &a);
		Nx = Nx + 1;
	}
	fclose(xfile);
	Nx = Nx - 1;
	double* xgrid1 = new double[Nx];
	xfile = fopen(xfileName, "r");
	for (int i = 0; i < Nx; ++i) {
		fscanf(xfile, "%lf", &xgrid1[i]);
	}
	fclose(xfile);
	int zeroIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}

	Nx = Nx - zeroIndex;

	xgrid = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		xgrid[Nx - i - 1] = 1000*xgrid1[i + zeroIndex];
	}

	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	double size = 0.5 * fabs(headMaxX);
	double B0 = 3E-7;
	double magneticEnergyDensity = B0 * B0 / (8 * pi);

	//RadiationSourceInCylindrical* downstreamSource = new SimpleFlatSource(upstreamElectrons, downstreamB, pi / 2, 0, concentration, size, size, distance);
	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photonsIR = new PhotonPlankDistribution(140, 0.8 / 1800000);
	double photonIRconcentration = photonsIR->getConcentration();
	double photonIRenergyDensity = photonIRconcentration * photonsIR->getMeanEnergy();
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();
	PhotonMultiPlankDistribution* photonsTotal = new PhotonMultiPlankDistribution(2.725, 1.0, 140, 0.8 / 1800000);
	//PhotonMultiPlankDistribution* photonsTotal = PhotonMultiPlankDistribution::getGalacticField();

	double photonTotalConcentration = photonsTotal->getConcentration();
	double photonTotalEnergyDensity = photonTotalConcentration * photonsTotal->getMeanEnergy();

	int Nz = 1;
	int Ny = 1;


	for (int i = 0; i < Nx; ++i) {
		xgrid[i] = -xgrid[i];
	}
	const char* fileName = "./examples_data/W50/newdistribution/electrons.dat";

	MassiveParticleTabulatedIsotropicDistribution* electrons1;
	double concentration1;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons1, concentration1);

	double*** B = new double** [Nx];
	double*** Btheta = new double** [Nx];
	double*** Bphi = new double** [Nx];
	double*** concentrationArray = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		B[i] = new double* [Ny];
		Btheta[i] = new double* [Ny];
		Bphi[i] = new double* [Ny];
		concentrationArray[i] = new double* [Ny];
		for (int j = 0; j < Ny; ++j) {
			B[i][j] = new double[Nz];
			Btheta[i][j] = new double[Nz];
			Bphi[i][j] = new double[Nz];
			concentrationArray[i][j] = new double[Nz];
			for (int k = 0; k < Nz; ++k) {
				B[i][j][k] = B0;
				Btheta[i][j][k] = pi / 2;
				Bphi[i][j][k] = 0;
				concentrationArray[i][j][k] = concentration1;
			}
		}
	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	double u = 0.25 * 0.2 * speed_of_light;
	RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, xgrid, Ny, Nz, electrons1, B, Btheta, Bphi, concentrationArray, 0, size, 0, pi * size, distance, u, photonTotalEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(Nx, downstreamXgrid, Ny, Nz, electrons2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* downstreamSource = new RectangularSource(1, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, downstreamXgrid[0], downstreamXgrid[Nx - 1], 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	int Ne = 1000;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 1000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	comptonEvaluator->writeEFEFromSourceToFile("./testBigSource.dat", source, 1.6E-10, 1.6E4, 500);



	//frontElectrons->writeDistribution("./output/thinDistribution.dat", 2000, me_c2, 1E10 * me_c2);

	double norm = 0;
	electrons1->transformToThickRegime(photonTotalEnergyDensity + magneticEnergyDensity, norm);

	double E0 = 1.6E-1;
	concentration1 *= u * pi * size * size * norm;

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, downstreamXgrid, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(Nx, downstreamXgrid, Ny, Nz, electrons2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	RectangularSource* source1 = new RectangularSource(1, 1, 1, electrons1, B0, pi / 2, 0, concentration1, 0, 1, 0, 1, 0, 1, distance);


	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false, false);
	//RadiationEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false, false);

	comptonEvaluator->writeEFEFromSourceToFile("./thickCompton.dat", source1, 1.6E-12, 1.6E4, 500);

	return;

}