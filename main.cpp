#include "stdio.h"
#include "math.h"
#include <omp.h>
#include <time.h>

#include "./Src/constants.h"
#include "./Src/massiveParticleDistribution.h"
#include "./Src/photonDistribution.h"
#include "./Src/util.h"
#include "./Src/inverseCompton.h"
#include "./Src/radiationSource.h"
#include "./Src/synchrotron.h"
#include "./Src/optimization.h"
#include "./Src/pionDecay.h"
#include "./Src/bremsstrahlung.h"
#include "./Src/coordinateTransform.h"
#include "./Src/radiationSourceFactory.h"
#include "./Src/diffusion.h"
#include "./Src/examples.h"

void evaluateFluxSNRtoWind() {
	FILE* logFile = fopen("log.dat", "w");
	fclose(logFile);

	double theta = pi/2;
	double index = 3.5;
	double Tstar = 3 * 1000;
	double luminosity = 10000 * 4 * 1E33;
	double rsun = 6.9E10;
	double rstar = rsun*sqrt((luminosity/4E33)/ pow(Tstar / 5800, 4));
	

	double electronConcentration = 2500;
	double B = 0.29;
	double rmax = 1.4E17;
	double sigma = B * B / (4 * pi * massProton * electronConcentration * speed_of_light2);
	double fraction = 0.5;
	//sigma = 0.0002;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//const double distance = 150 * 1000000 * parsec;
	//AT2020xnd
	const double distance = 1232 * 1000000 * parsec;

	double Emin = me_c2;
	//double Emax = 1E12 * me_c2;
	double Emax = 1E4 * me_c2;
	int Ne = 40;
	int Nmu = 200;

	int Nrho = 50;
	int Nz = 50;
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
	//initializing upstreamElectrons distribution
	//MassiveParticlePowerLawDistribution* upstreamElectrons = new MassiveParticlePowerLawDistribution(massElectron, index, 10*me_c2, electronConcentration);
	//MassiveParticleBrokenPowerLawDistribution* upstreamElectrons = new MassiveParticleBrokenPowerLawDistribution(massElectron, index, 2.001, 2*me_c2, 1000*me_c2, electronConcentration);
    //MassiveParticleTabulatedIsotropicDistribution* upstreamElectrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_combined_cutoff/Ee3.dat", "./examples_data/gamma0.2_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
    MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* upstreamElectrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_theta0-90/Ee3.dat", "./examples_data/gamma0.2_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleMaxwellJuttnerDistribution* upstreamElectrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, 1.08 * me_c2 / kBoltzman);
	MassiveParticleTabulatedIsotropicDistribution* protons = new MassiveParticleTabulatedIsotropicDistribution(massProton, "./examples_data/gamma1.5_theta0-90_protons/Ee3.dat", "./examples_data/gamma1.5_theta0-90_protons/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	MassiveParticleMaxwellJuttnerDistribution* thermalElectrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, 4 * me_c2 / kBoltzman);
	double velocity = 0.2 * speed_of_light;
	//double velocity = 0.55 * speed_of_light;
	double gamma = 1.0 / sqrt(1.0 - velocity * velocity / speed_of_light2);
	//upstreamElectrons->rescaleDistribution(1.2);
	//upstreamElectrons->addPowerLaw(10 * massElectron * speed_of_light2, 3.0);
	//(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(upstreamElectrons))->prolongEnergyRange(1E6, 100);
	//upstreamElectrons->addPowerLaw(200 * massElectron * speed_of_light2, 3.5);
	electrons->writeDistribution("distribution.dat", 200, Emin, Emax);
	thermalElectrons->writeDistribution("distribution1.dat", 200, Emin, Emax);
	

	//creating radiation downstreamSource
	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, rmax, fraction*rmax, distance, 0, 0.2433);
	//RadiationSource* downstreamSource = new TabulatedSphericalLayerSource(20, 20, 1,upstreamElectrons, downstreamB, theta, electronConcentration, rmax, (1-fraction)*rmax, distance);
	RadiationSourceInCylindrical* thermalSource = new SimpleFlatSource(thermalElectrons, B, theta, 0, electronConcentration, rmax, fraction*rmax, distance);
	//RadiationSource* downstreamSource = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, electronConcentration, rmax, 0.5*rmax, distance);
	//RadiationSource* downstreamSource = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronsFromSmilei, Bturb, thetaTurb, concentration, rsource, 0.9*rsource, distance);
	//AngleDependentElectronsSphericalSource* downstreamSource = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, Bturb, thetaTurb, phiTurb, concentration, rmax, 0.9*rmax, distance, 0.5*speed_of_light);
	//counting of quai parallel
	/*int nangle = 0;
	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double angle = (dynamic_cast<AngleDependentElectronsSphericalSource*>(downstreamSource))->getShockWaveAngle(irho, iz, iphi);
				if ((angle < 2*pi / 9)||(angle > 7*pi/9)) {
					nangle++;
				}
			}
		}
	}*/



	//SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true);
	//SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, newEmax, true, true);
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true, false);
	SynchrotronEvaluator* synchrotronEvaluator2 = new SynchrotronEvaluator(Ne, Emin, Emax, false, false);
	//comptonEvaluator->outputDifferentialFlux("output1.dat");
	//return;
	int Ndays = 26;
	rmax = velocity * Ndays * 24 * 3600;
	//rmax = velocity * 99 * 24 * 3600 + 0.5 * speed_of_light * (Ndays - 99) * 24 * 3600;
	//number of parameters of the downstreamSource
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { 0.9*rmax, 1E-7, 2E4, 0.01, 0.1*speed_of_light };
	double maxParameters[Nparams] = { 1.5*rmax, 5E0, 5E6, 0.7, 0.3*speed_of_light };
	//starting point of optimization and normalization
	//fraction = 2E13 / rmax;

	//fraction = (0.4-0.04 + 0.004/3);
	fraction = 0.5*4.0/3.0;

	double denseFactor = 0.5 / fraction;
	electronConcentration = 1.0;
	//sigma = 0.5*(0.01 * 0.5 / (0.012))/denseFactor;
	//sigma = 0.34 * 0.34 / (4 * pi * massProton * speed_of_light2 * electronConcentration);
	sigma = 0.01;
	//electronConcentration = 3167 * sqr(69.0 / Ndays);
	//sigma = 2.33E-5;
	double vector[Nparams] = { rmax, sigma, electronConcentration, fraction, velocity };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
	printf("new sigma = %g\n", sigma);
	printLog("new sigma = %g\n", sigma);
	printf("new concentration = %g\n", electronConcentration);
	printLog("new concentration = %g\n", electronConcentration);

	/*const int Nenergy1 = 4;
	double energy1[Nenergy1] = { 1.5E9, 3.0E9 , 6.1E9, 9.87E9 };
	double observedFlux1[Nenergy1] = { 1.5, 4.3, 6.1, 4.2 };
	double observedError1[Nenergy1] = { 0.1, 0.2 , 0.3, 0.2 };*/

	/*const int Nenergy1 = 4;
	double energy1[Nenergy1] = { 2.94E9 , 6.1E9, 9.74E9, 22.0E9};
	double observedFlux1[Nenergy1] = {2.9, 2.3, 1.74, 0.56};
	double observedError1[Nenergy1] = {0.2 , 0.1, 0.09, 0.03 };*/

	/*const int Nenergy1 = 7;
	double energy1[Nenergy1] = { 0.33E9, 0.61E9 , 1.39E9, 1.5E9, 3.0E9, 6.05E9, 10E9 };
	double observedFlux1[Nenergy1] = { 0.375, 0.79, 0.38, 0.27, 0.17, 0.07, 0.032 };
	double observedError1[Nenergy1] = { 0.375, 0.09 , 0.05, 0.07, 0.03, 0.01, 0.008 };*/

	/*for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}*/

	double* energy1;
	double* observedFlux1;
	double* observedError1;
	int Nenergy1 = readRadiationFromFile(energy1, observedFlux1, observedError1, "./examples_data/AT2020xnd_data/bright26.dat");
	for (int i = 0; i < Nenergy1; ++i) {
		energy1[i] = energy1[i] * hplank * 1E9;
		observedFlux1[i] = observedFlux1[i] / (hplank * 1E26);
		observedError1[i] = observedError1[i] / (hplank * 1E26);
	}


	printf("start optimization\n");
	printLog("start optimization\n");
	bool optPar[Nparams] = { true, true, true, false, false };
	int Niterations = 5;
	//creating KPIevaluator
	LossEvaluator* KPIevaluator = new SpectrumLossEvaluator(energy1, observedFlux1, observedError1, Nenergy1, source);
	//creating gradient descent optimizer
	RadiationOptimizer* synchrotronOptimizer = new GradientDescentRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//RadiationOptimizer* synchrotronOptimizer = new CoordinateRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	//number of points per axis in gridEnumOptimizer
	int Npoints[Nparams] = { 50,50,50,10,2 };
	//creating grid enumeration optimizer
	RadiationOptimizer* combinedOptimizer = new CombinedRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Niterations, Npoints, KPIevaluator);
	RadiationOptimizer* enumOptimizer = new GridEnumRadiationOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);
	//grid enumeration optimization, finding best starting point for gradien descent
	
	//double error = synchrotronOptimizer->evaluateOptimizationFunction(vector);
	//printf("starting error = %g\n", error);
	//printLog("starting error = %g\n", error);
	
	//enumOptimizer->optimize(vector, optPar, downstreamSource);
	//gradient descent optimization
	//synchrotronOptimizer->optimize(vector, optPar);
	//reseting downstreamSource parameters to found values
	//synchrotronOptimizer->outputProfileDiagrams(vector, downstreamSource, 10);
	//synchrotronOptimizer->outputOptimizedProfileDiagram(vector, optPar, downstreamSource, 10, 1, 2);
	//combinedOptimizer->optimize(vector, optPar);
	B = 4.986;
	electronConcentration = 1.15E6;
	vector[0] = 1.782E16 / maxParameters[0];
	vector[1] = (B*B/(4*pi*massProton*electronConcentration*speed_of_light2)) / maxParameters[1];
	vector[2] = electronConcentration / maxParameters[2];
	//vector[3] = 1.33*0.5 / maxParameters[3];
	//vector[4] = 0.55*speed_of_light / maxParameters[4];
	source->resetParameters(vector, maxParameters);
	thermalSource->resetParameters(vector, maxParameters);
	//evaluating resulting error
	double error = combinedOptimizer->evaluateOptimizationFunction(vector);
	printf("resulting error = %g\n", error);
	printLog("resulting error = %g\n", error);

	/*double totalVolume = downstreamSource->getTotalVolume();
	double totalVolume2 = 0;
	for (int i = 0; i < 20; ++i) {
		for (int j = 0; j < 20; ++j) {
			for (int k = 0; k < 1; ++k) {
				totalVolume2 += downstreamSource->getVolume(i, j, k);
			}
		}
	}

	double area1 = pi * rmax * rmax * (1 - (1 - fraction) * (1 - fraction));
	double area = 0;
	for (int i = 0; i < 20; ++i) {
		for (int k = 0; k < 1; ++k) {
			area += downstreamSource->getArea(i, 10, k);
		}
	}


	printf("total volume = %g\n", totalVolume);
	printf("total volume2 = %g\n", totalVolume2);*/

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
	fprintf(paramFile,"R/t*c = %g\n", vector[0] * maxParameters[0] / (Ndays * 24 * 3600 * speed_of_light));
	printf("velocity/c = %g\n", vector[4] * maxParameters[4]/speed_of_light);
	fprintf(paramFile, "celocity/c = %g\n", vector[4] * maxParameters[4]/speed_of_light);
	B = sqrt(vector[1] * maxParameters[1] * 4 * pi * massProton * vector[2] * maxParameters[2] * speed_of_light2);
	printf("average B = %g\n", B);
	fprintf(paramFile, "average B = %g\n", B);
	fclose(paramFile);

	//TabulatedDiskSourceWithSynchCutoff* source3 = new TabulatedDiskSourceWithSynchCutoff(1, 10000, 1, upstreamElectrons, downstreamB, pi / 2, 0, electronConcentration, rmax, fraction * rmax, distance, 0.2 * speed_of_light);


	rmax = vector[0] * maxParameters[0];
	electronConcentration = vector[2] * maxParameters[2];
	fraction = vector[3] * maxParameters[3];
	velocity = vector[4] * maxParameters[4] / speed_of_light;
	gamma = 1.0 / sqrt(1.0 - velocity * velocity);

	double actualVelocity = rmax / (Ndays * 24 * 3600 * speed_of_light);
	double actualGamma = 1.0 / sqrt(1.0 - actualVelocity * actualVelocity);

	double energy2 = 0;
	double Estart = 10 * me_c2;
	double factor = pow(Emax / Estart, 1.0 / (200.0 - 1));
	double e1 = Estart;
	for (int i = 0; i < 200; ++i) {
		double de = e1 * (factor - 1);
		energy2 = energy2 + electronConcentration * electrons->distributionNormalized(e1) *(e1 - me_c2)* (pi * rmax * rmax * rmax * fraction)* de;
	}

	double totalEnergy = electronConcentration * massProton * speed_of_light2 * (actualGamma - 1.0)*(pi*rmax*rmax*rmax*fraction);
	//double totalEnergy = electronConcentration * massProton * speed_of_light2 * (gamma - 1.0)*downstreamSource->getTotalVolume();
	double energyInRadioElectrons = electronConcentration * (electrons->getMeanEnergy() - me_c2) * (pi * rmax * rmax * rmax * fraction);
	//double energyInRadioProtons = electronConcentration * (protons->getMeanEnergy() - massProton*speed_of_light2) * downstreamSource->getTotalVolume();
	double magneticEnergy = (B * B / (8 * pi)) * (pi * rmax * rmax * rmax * fraction);
	printf("total kinetik energy = %g\n", totalEnergy);
	printLog("total kinetik energy = %g\n", totalEnergy);
	printf("energy in radio electrons = %g\n", energyInRadioElectrons);
	printLog("energy in radio electrons = %g\n", energyInRadioElectrons);
	//printf("energy in radio protons = %g\n", energyInRadioProtons);
	//printLog("energy in radio protons = %g\n", energyInRadioProtons);
	printf("energy in magnetic field = %g\n", magneticEnergy);
	printLog("energy in magnetic field = %g\n", magneticEnergy);
	//double totalMass = electronConcentration * massProton * downstreamSource->getTotalVolume();
	//printf("total mass = %g\n", totalMass);
	//printLog("total mass = %g\n", totalMass);

	double epsilone = energyInRadioElectrons / (totalEnergy);
	double epsilonB = magneticEnergy / totalEnergy;
	printf("epsilon e = %g\n", epsilone);
	printLog("epsilon e = %g\n", epsilone);
	printf("epsilon B = %g\n", epsilonB);
	printLog("epsilon B = %g\n", epsilonB);

	BremsstrahlungThermalEvaluator* bremsstrahlungEvaluator = new BremsstrahlungThermalEvaluator();
	//bremsstrahlungEvaluator->writeFluxFromSourceToFile("bremsstrahlung.dat", thermalSource, 1.6E-9, 1.6E-5, 200);
	//double bremMevFlux = bremsstrahlungEvaluator->evaluateTotalFluxInEnergyRange(0.3 * 1000000 * 1.6E-12, 10 * 1000000 * 1.6E-12, 200, thermalSource);
	//printf("bremMevFlux = %g\n", bremMevFlux);
	//printLog("bremMevFlux = %g\n", bremMevFlux);
	//initialization arrays for full synchrotron spectrum
	const int Nnu = 50;
	double* Nu = new double[Nnu];
	double* F = new double[Nnu];
	double* F2 = new double[Nnu];

	double Numin = 1E8;
	double Numax = 1E14;
	factor = pow(Numax / Numin, 1.0 / (Nnu - 1));
	Nu[0] = Numin;
	F[0] = 0;
	F2[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		Nu[i] = Nu[i - 1] * factor;
		F[i] = 0;
		F2[i] = 0;
	}

	//evaluationg full spectrum of the downstreamSource
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(hplank * Nu[i], source);
		F2[i] = synchrotronEvaluator2->evaluateFluxFromSource(hplank * Nu[i], source);
	}

	//outputing spectrum
	FILE* output_GZ_Jansky = fopen("outputSynch1.dat", "w");
	//FILE* output_GZ_Jansky2 = fopen("outputSynch2.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_GZ_Jansky, "%g %g\n", Nu[i] / 1E9, hplank * F[i] * 1E26);
		//fprintf(output_GZ_Jansky2, "%g %g\n", Nu[i] / 1E9, hplank * F2[i] * 1E26);
	}
	fclose(output_GZ_Jansky);

	return;
	//fclose(output_GZ_Jansky2);
	//return;
	double synchrotronFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(hplank * 1E8, hplank * 1E11, 200, source);
	double synchrotronKevFlux1 = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(0.3*1000*1.6E-12, 10 * 1000 * 1.6E-12, 200, source);
	//double synchrotronKevFlux2 = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(0.3*1000*1.6E-12, 10 * 1000 * 1.6E-12, 200, source3);
	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch3.dat", downstreamSource, hplank * 1E8, 10 * 1000 * 1.6E-12, 1000);
	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch5.dat", source3, hplank * 1E8, 50*1000 * 1000 * 1.6E-12, 500);
	double synchrotronFlux2 = synchrotronEvaluator2->evaluateTotalFluxInEnergyRange(hplank * 1E8, hplank * 1E11, 20, source);
	printf("total synchrotron flux = %g erg/cm^2 s\n", synchrotronFlux);
	printLog("total synchrotron flux = %g erg/cm^2 s\n", synchrotronFlux);
	printf("total synchrotron flux without absorption = %g erg/cm^2 s\n", synchrotronFlux2);
	printLog("total synchrotron flux without absorption = %g erg/cm^2 s\n", synchrotronFlux2);
	printf("synchrotron keV flux = %g erg/cm^2 s\n", synchrotronKevFlux1);
	printLog("synchrotron keV flux = %g erg/cm^2 s\n", synchrotronKevFlux1);
	//printf("cutoff model synchrotron keV flux = %g erg/cm^2 s\n", synchrotronKevFlux2);
	//printLog("cutoff model synchrotron keV flux = %g erg/cm^2 s\n", synchrotronKevFlux2);

	//PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	//PhotonIsotropicDistribution* photonDistribution = PhotonPlankDistribution::getCMBradiation();
	double rcompton = 2E14;
	rmax = 2E15;
	fraction = 0.5;
	rcompton = rmax + 0.5E16;
	rcompton = 2E14;
	PhotonPlankDistribution* photonDistribution = new PhotonPlankDistribution(Tstar, 0.25*sqr(rstar / rcompton));
	double photonConcentration = photonDistribution->getConcentration();
	//PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	PhotonPlankDirectedDistribution* photonDirectedDistribution = new PhotonPlankDirectedDistribution(Tstar, 0.25*sqr(rstar / rcompton), pi*17/18, 0, atan2(1, 1));
	double photonDirectedConcentration = photonDirectedDistribution->getConcentration();
	PhotonMultiPlankDistribution* galacticField = PhotonMultiPlankDistribution::getGalacticField();
	double photonGalacticConcentration = galacticField->getConcentration();
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
    InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ne, Ephmin, Ephmax, photonDirectedDistribution, photonDirectedConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluatorWithSource(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES, rmax + 0.5E16, 0, 0);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluatorWithSource(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA, 0, -rmax - 0.5E16, 0);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_THOMSON);
	//initializing photon energy grid for output
	InverseComptonEvaluator* galacticEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ne, Ephmin, Ephmax, galacticField, photonGalacticConcentration, ComptonSolverType::ISOTROPIC_JONES);
	galacticEvaluator->writeFluxFromSourceToFile("outputCompton.dat", source, 0.1 * 1.6E-12, 1.6E-5, 1000);
	double galacticComptonFlux = galacticEvaluator->evaluateTotalFluxInEnergyRange(0.1 * 1.6E-12, 1.6E-5, 1000, source);
	printf("galactic Compton flux = %g\n", galacticComptonFlux);
	printLog("galactic Compton flux = %g\n", galacticComptonFlux);

	double* E = new double[Nnu];

	double EphFinalmin = 0.1*1000*1.6E-12;
	double EphFinalmax = 10000 * 1000 * 1.6E-12;
	//photonDistribution->writeDistribution("output3.dat", 200, Ephmin, Ephmax);
	factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
	E[0] = EphFinalmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}
	double jetelectronConcentration = 1E7;
	//MassiveParticleDistribution* jetelectrons = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 30 * me_c2);
	MassiveParticleDistribution* jetelectrons1 = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 1.1 * me_c2);
	//MassiveParticleDistribution* jetelectrons = new MassiveParticleMonoenergeticDistribution(massElectron, 20 * me_c2, me_c2);
	//MassiveParticleDistribution* jetelectrons1 = new MassiveParticleMonoenergeticDistribution(massElectron, 1.1 * me_c2, 0.1*me_c2);
	//MassiveParticleDistribution* jetelectrons = new MassiveParticleMonoenergeticDirectedDistribution(massElectron, 30 * me_c2, me_c2, 0, 0, 0.1);
	MassiveParticleDistribution* jetelectrons = new MassiveParticleMovingDistribution(jetelectrons1, sqrt(1.0 - 1.0/(30*30)) * speed_of_light);
	FILE* movingDistributionFile = fopen("distribution.dat", "w");
	factor = pow(1E6, 1.0 / (Ne - 1));
	double currentE = me_c2;
	for (int i = 0; i < Ne; ++i) {
		double distribution = jetelectrons->distributionNormalized(currentE, 1.0, 0.0);
		fprintf(movingDistributionFile, "%g %g\n", currentE, distribution);
		currentE = currentE * factor;
	}
	fclose(movingDistributionFile);

	FILE* angleDistribution = fopen("angledistribution.dat", "w");
	for (int i = 0; i < 100; ++i) {
		double mu = -1 + (i + 0.5) * 2.0 / 100;
		double distribution1 = jetelectrons->distributionNormalized(1.1 * me_c2, mu, 0.0);
		double distribution2 = jetelectrons->distributionNormalized(10 * me_c2, mu, 0.0);
		double distribution3 = jetelectrons->distributionNormalized(20 * me_c2, mu, 0.0);
		double distribution4 = jetelectrons->distributionNormalized(30 * me_c2, mu, 0.0);
		double distribution5 = jetelectrons->distributionNormalized(50 * me_c2, mu, 0.0);
		double distribution6 = jetelectrons->distributionNormalized(1000 * me_c2, mu, 0.0);
		fprintf(angleDistribution, "%g %g %g %g %g %g %g\n", mu, distribution1, distribution2, distribution3, distribution4, distribution5, distribution6);
	}
	fclose(angleDistribution);

	RadiationSource* source2 = new SimpleFlatSource(jetelectrons, B, pi / 2, 0, jetelectronConcentration, rmax, rmax*fraction, distance);

	double V = source2->getTotalVolume();
	double totalEnergy2 = jetelectronConcentration * massProton * speed_of_light2 * 10 * source2->getTotalVolume();
	double energyInComptonElectrons = jetelectronConcentration * (jetelectrons->getMeanEnergy()) * source2->getTotalVolume();
	double photonEnergyDensity = photonDistribution->getConcentration() * photonDistribution->getMeanEnergy();
	printf("source volume = %g\n", V);
	printLog("source volume = %g\n", V);
	printf("jet electron concentration = %g\n", jetelectronConcentration);
	printLog("jet electron concentration = %g\n", jetelectronConcentration);
	printf("photon energy density = %g\n", photonEnergyDensity);
	printLog("photon energy density = %g\n", photonEnergyDensity);
	printf("energy in compton electrons = %g\n", energyInComptonElectrons);
	printLog("energy in compton electrons = %g\n", energyInComptonElectrons);
	printf("energy in protons = %g\n", totalEnergy2);
	printLog("energy in protons = %g\n", totalEnergy2);

	FILE* output_GZ_Jansky3 = fopen("outputSynch6.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		//fprintf(output_GZ_Jansky3, "%g %g\n", Nu[i]/1E9, hplank*1E26*synchrotronEvaluator->evaluateFluxFromSource(hplank * Nu[i], source2));
	}
	fclose(output_GZ_Jansky3);

	double minEev = 0.3 * 1000 * 1.6E-12;
	double maxEev = 10 * 1000 * 1.6E-12;
	//double kevFlux = comptonEvaluator->evaluateTotalFluxInEnergyRange(minEev, maxEev, 10, downstreamSource);
	double kevFlux = 0;
	int Nph = 100;
	factor = pow(maxEev / minEev, 1.0 / (Nph - 1));
	currentE = minEev;
	FILE* output7 = fopen("output7.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		printf("%d\n", i);
		double dE = currentE * (factor - 1.0);
		//double dE = currentE * (1.0 - 1.0/factor);
		//kevFlux += comptonEvaluator->evaluateFluxFromSourceAnisotropic(currentE, 0, 0, photonDirectedDistribution, downstreamSource, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA) * dE;
		double flux = comptonEvaluator->evaluateFluxFromSource(currentE, source2);
		//double flux = galacticEvaluator->evaluateFluxFromSource(currentE, source2);
		fprintf(output7,"%g %g\n", currentE, flux);
		kevFlux += flux * dE;
		currentE = currentE * factor;
	}
	fclose(output7);
	double totalLuminosity = kevFlux * 4 * pi * distance * distance;
	//double synchrotronKevFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(minEev, maxEev, 10, downstreamSource);
	//double synchrotronFluxGHz3 = synchrotronEvaluator->evaluateFluxFromSource(3E9 * hplank, downstreamSource)*1E26*hplank;
	
	FILE* outFile = fopen("SNRtoWindData.dat", "w");
	printf("total luminosity = %g erg/s \n", totalLuminosity);
	printLog("total luminosity = %g erg/s \n", totalLuminosity);
	fprintf(outFile, "total luminosity = %g erg/s \n", totalLuminosity);
	printf("total flux = %g erg/s cm^2 \n", kevFlux);
	printLog("total flux = %g erg/s cm^2 \n", kevFlux);
	fprintf(outFile, "total flux = %g erg/s cm^2 \n", kevFlux);
	//printf("synchrotron total keV flux = %g erg/s cm^2 \n", synchrotronKevFlux);
	//printLog("synchrotron total keV flux = %g erg/s cm^2 \n", synchrotronKevFlux);
	//fprintf(outFile, "synchrotron total keV flux = %g erg/s cm^2 \n", synchrotronKevFlux);
	printf("99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	printLog("99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	fprintf(outFile, "99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	//printf("synchrotron radio flux at 3 GHz = %g mJy\n", synchrotronFluxGHz3);
	//printLog("synchrotron radio flux at 3 GHz = %g mJy\n", synchrotronFluxGHz3);
	//fprintf(outFile, "synchrotron radio flux at 3 GHz = %g mJy\n", synchrotronFluxGHz3);
	printf("99 days F(3GHz) = 4.3 mJy\n");
	printLog("99 days F(3GHz) = 4.3 mJy\n");
	fprintf(outFile, "99 days F(3GHz) = 4.3 mJy\n");
	fclose(outFile);

	//return;
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
		F[i] = comptonEvaluator->evaluateFluxFromSource(E[i], source2);
		//F[i] = comptonEvaluator->evaluateFluxFromSourceAnisotropic(E[i], 0, 0, CMBradiation, downstreamSource);
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

	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch.dat", downstreamSource, Ephmin, Ephmax, 200);
	/*FILE* output_GHz_Jansky = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = 0.0001*E[i] / hplank;
		F[i] = synchrotronEvaluator->evaluateFluxFromSource(0.0001*E[i], downstreamSource);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_GHz_Jansky);*/

	delete[] E;
	delete[] F;
	delete[] F2;
	delete electrons;
	delete source;
	delete comptonEvaluator;
}

void evaluateComtonFromWind() {
	double theta = pi / 2;
	double rmax = 2E15;
	double rtotal = 1.4E17;
	//double rmax = 1.0 / sqrt(pi);
	double B = 0.01;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 1000000 * parsec;
	//const double distance = 1.0;

	//double Emin = 652.317 * me_c2 * 1;
	double Emin = 10 * me_c2;
	double Emax = 1E6 * me_c2;
	int Ne = 200;
	int Nmu = 20;
	int Nrho = 20;
	int Nz = 20;
	int Nphi = 4;
	double index = 3.5;
	double KK = 24990.8;
	//double electronConcentration = KK / (pow(652.317, index - 1) * (index - 1));
	double electronConcentration = 1.0E7;
	double protonBulkConcentration = 1E4;

	double Tstar = 25 * 1000;
	double Ephmin = 0.01 * Tstar * kBoltzman;
	double Ephmax = 100 * Tstar * kBoltzman;
	double luminosity = 510000 * 4 * 1E33;
	double rsun = 7.5E10;
	double rstar = rsun * sqrt(510000.0 / pow(Tstar / 5772, 4));
	//PhotonIsotropicDistribution* photonDistribution = new PhotonPlankDistribution(Tstar, sqr(rstar / rmax));
	//PhotonIsotropicDistribution* photonDistribution = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDirectedDistribution* photonDistribution = new PhotonPlankDirectedDistribution(Tstar, sqr(rstar / rmax), 11*pi/12, 0, pi / 4);
	double photonConcentration = photonDistribution->getConcentration();

	//initializing upstreamElectrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);
	//MassiveParticleTabulatedIsotropicDistribution* upstreamElectrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.3_theta0-90/Ee3.dat", "./examples_data/gamma0.3_theta0-90/Fs3.dat", 200, GAMMA_KIN_FGAMMA);
	//MassiveParticleIsotropicDistribution* upstreamElectrons = new MassiveParticleMonoenergeticDistribution(massElectron, 2*Emin, me_c2 / 2);
	//upstreamElectrons->addPowerLaw(1.02 * me_c2, 3);
	//upstreamElectrons->rescaleDistribution(sqrt(18));
	//upstreamElectrons->addPowerLaw(100 * me_c2, 3.5);
	electrons->writeDistribution("dist1.dat", 2000, Emin, Emax);
	double electronsEnergy = electronConcentration * (electrons->getMeanEnergy() - me_c2);
	double v = 0.5 * speed_of_light;
	double gam = 1.0 / sqrt(1 - v * v / speed_of_light2);
	double ramPressure = protonBulkConcentration * massProton * speed_of_light2 * (gam - 1);
	double Etotal = ramPressure * 4 * pi * rtotal * rtotal * 0.01 * rtotal;
	double EtotalElectrons = electronsEnergy * pi * rmax * rmax * 0.2 * rmax;
	double availableEnergyFraction = rmax * rmax / (4 * rtotal * rtotal);
	printf("electronsEnergyDensity/ramPressure = %g\n", electronsEnergy / ramPressure);
	printLog("electronsEnergyDensity/ramPressure = %g\n", electronsEnergy / ramPressure);

	printf("total Energy = %g\n", Etotal);
	printLog("total Energy = %g\n", Etotal);

	printf("total electrons Energy = %g\n", EtotalElectrons);
	printLog("total electrons Energy = %g\n", EtotalElectrons);

	printf("total electrons/bulk = %g\n", EtotalElectrons/Etotal);
	printLog("total electrons/bulk = %g\n", EtotalElectrons/Etotal);

	printf("available fraction = %g\n", availableEnergyFraction);
	printLog("tavailable fraction = %g\n", availableEnergyFraction);
	//creating radiation downstreamSource
	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, electronConcentration, rmax, 0.2*rmax, distance);
	//RadiationSource* downstreamSource = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, electronConcentration, rmax, 0.8 * rmax, distance);
	//InverseComptonEvaluator* comptonEvaluator1 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator2 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator3 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, 100, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_THOMSON);

	//comptonEvaluator2->outputDifferentialFlux("output3.dat");
	//comptonEvaluator->outputDifferentialFluxJones("output2.dat", photonDistribution, upstreamElectrons);
	//return;

	double minEev = 0.3 * 1000 * 1.6E-12;
	double maxEev = 10 * 1000 * 1.6E-12;
	double kevFlux = comptonEvaluator3->evaluateTotalFluxInEnergyRange(minEev, maxEev, 20, source);
	double MevFlux = comptonEvaluator3->evaluateTotalFluxInEnergyRange(0.1*1.6E-6, 1.6E-6, 20, source);
	double totalLuminosity = kevFlux * 4 * pi * distance * distance;
	FILE* outFile = fopen("SNRtoWindData.dat", "w");
	printf("total luminosity = %g erg/s \n", totalLuminosity);
	fprintf(outFile, "total luminosity = %g erg/s \n", totalLuminosity);
	printf("total kev flux = %g erg/s cm^2 \n", kevFlux);
	printf("total Mev flux = %g erg/s cm^2 \n", MevFlux);
	printf("99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	fprintf(outFile, "99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	fclose(outFile);

	//initializing photon energy grid for output
	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double EphFinalmin = 0.1 * kBoltzman * Tstar;
	double EphFinalmax = 1.6E-6 + Emin;
	//photonDistribution->writeDistribution("output3.dat", 200, Ephmin, Ephmax);
	double factor = pow(EphFinalmax / EphFinalmin, 1.0 / (Nnu - 1));
	E[0] = EphFinalmin;
	F[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		F[i] = 0;
	}

	//evaluating radiation flux
	printLog("evaluating\n");
	/*for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		printLog("%d\n", i);
		F[i] = comptonEvaluator->evaluateFluxFromSource(E[i], downstreamSource);
		//F[i] = comptonEvaluator->evaluateFluxFromSourceAnisotropic(E[i], 0, 0, CMBradiation, downstreamSource);
	}*/

	//outputing
	FILE* output_ev_EFE1 = fopen("output1.dat", "w");
	FILE* output_ev_EFE2 = fopen("output2.dat", "w");
	FILE* output_ev_EFE3 = fopen("output3.dat", "w");
	//FILE* output_GHz_Jansky = fopen("output.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		double nu = E[i] / hplank;
		//fprintf(output_ev_EFE1, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator1->evaluateFluxFromSource(E[i], downstreamSource));
		//fprintf(output_ev_EFE2, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator2->evaluateFluxFromSource(E[i], downstreamSource));
		fprintf(output_ev_EFE3, "%g %g\n", E[i], comptonEvaluator3->evaluateFluxFromSource(E[i], source));
		//fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * F[i]);
	}
	fclose(output_ev_EFE1);
	fclose(output_ev_EFE2);
	fclose(output_ev_EFE3);
	//fclose(output_GHz_Jansky);

	delete[] E;
	delete[] F;
	delete electrons;
	delete source;
	//delete comptonEvaluator1;
	//delete comptonEvaluator2;
	delete comptonEvaluator3;
}

void evaluateTychoProfile() {
	int Nrho = 20;
	int Nphi = 8;
	int Nz = 40;

	double B0 = 30*3E-6;
	double magneticEnergy = B0 * B0 / (8 * pi);
	double theta0 = pi / 2;
	double phi0 = 0;
	double*** B = create3dArray(Nrho, Nz, Nphi, B0);
	double*** theta = create3dArray(Nrho, Nz, Nphi, theta0);
	double*** phi = create3dArray(Nrho, Nz, Nphi, phi0);
	double turbulentEnergy = 150*((280 + 2 * 218)*1E-12)/(8*pi);
	double fraction = turbulentEnergy/magneticEnergy;
	double index = 11.0 / 6.0;
	int Nmodes = 20;
	double anisotropy = 0.75; //0.75 to Bx^2/By^2 = 281/218 as uvarov fit
	double concentration = 1000;
	double*** concentrations = create3dArray(Nrho, Nz, Nphi, concentration);
	double Emin = me_c2;
	double Emax = me_c2 * 1E9;
	int Ne = 200;
	double distance = 2500 * parsec;
	double R = distance*258*pi/(180*60*60);
	double widthFraction = 0.25;
	double Rin = R * (1 - widthFraction);
	double Rmin = R * (1 - 2 * widthFraction);
	double lturb = 1E17;

	double Uupstream = 4E8;
	double Udownstream = 0.25 * Uupstream;
	double meanB = sqrt(8 * pi * (magneticEnergy + turbulentEnergy));

	double nuph = (5000 * 1.6E-12) / hplank;
	double Energy = sqrt(4 * pi * nuph * cube(massElectron) * speed_of_light * speed_of_light4 / (0.29 * 3 * electron_charge * meanB));
	double gamma = Energy / me_c2;
	double Rlosses = (2.0 / 3.0) * sqr(electron_charge * electron_charge / me_c2) * speed_of_light * (4.0 / 9.0) * meanB * meanB * sqr(Energy / me_c2);
	double tau = Energy / Rlosses;
	double L = tau * Udownstream;

	resetLog();

	printf("Tycho profile\n");
	printLog("Tycho profile\n");

	printf("Turbulent energy fraction = %g\n", fraction);
	printLog("Turbulent energy fraction = %g\n", fraction);

	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInDiskSource(downstreamB, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSphericalSource(downstreamB, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSectoralSphericalSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, Rmin, 2*pi, anisotropy);
	write3dArrayToFile(B, Nrho, Nz, Nphi, "B.dat");

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 2.0, 100*Energy);

	//RadiationSource* downstreamSource = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, concentrations, R, (1.0 - widthFraction) * R, distance);
	//RadiationSource* downstreamSource = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, concentrations, R, (1.0 - widthFraction) * R, distance, Udownstream);
	RadiationSourceInCylindrical* source = new TabulatedSectoralSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, 0, concentrations, R, Rin,Rmin, 2*pi, distance, Udownstream);
	//RadiationSource* downstreamSource = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, concentrations, R, R, distance, Udownstream);

	SynchrotronEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax, false);

	evaluator->writeImageFromSourceToFile("image.dat", source, 1000*1.6E-12, 10000*1.6E-12, 20);

	delete3dArray(B, Nrho, Nz, Nphi);
	delete3dArray(theta, Nrho, Nz, Nphi);
	delete3dArray(phi, Nrho, Nz, Nphi);
}

void fitTychoProfile() {
	int Nrho = 200;
	int Nphi = 8;
	int Nz = 400;

	double B0 = 100*3E-6;
	double magneticEnergy = B0 * B0 / (8 * pi);
	double theta0 = pi / 2;
	double phi0 = 0;
	double*** B = create3dArray(Nrho, Nz, Nphi, B0);
	double*** theta = create3dArray(Nrho, Nz, Nphi, theta0);
	double*** phi = create3dArray(Nrho, Nz, Nphi, phi0);
	double turbulentEnergy = 1E-20*((280 + 2 * 218) * 1E-12) / (8 * pi);
	double turbulentFraction = turbulentEnergy / magneticEnergy;
	double index = 11.0 / 6.0;
	int Nmodes = 20;
	double anisotropy = 0.75; //0.75 to Bx^2/By^2 = 281/218 as uvarov fit
	double concentration = 1000;
	double*** concentrations = create3dArray(Nrho, Nz, Nphi, concentration);
	double Emin = me_c2;
	double Emax = me_c2 * 1E10;
	int Ne = 200;
	double distance = 2500 * parsec;
	double R = distance * 258 * pi / (180 * 60 * 60);
	double widthFraction = 0.1;
	double Rin = R * (1 - widthFraction);
	double Rmin = R * (1 - 2 * widthFraction);
	double lturb = 1E17;

	double Uupstream = 4E8;
	double Udownstream = 0.25 * Uupstream;
	double meanB = sqrt(8 * pi * (magneticEnergy + turbulentEnergy));

	double nuph = (5000 * 1.6E-12) / hplank;
	double Energy = sqrt(4 * pi * nuph * cube(massElectron) * speed_of_light * speed_of_light4 / (0.29 * 3 * electron_charge * meanB));
	double gamma = Energy / me_c2;
	double Rlosses = (2.0 / 3.0) * sqr(electron_charge * electron_charge / me_c2) * speed_of_light * (4.0 / 9.0) * meanB * meanB * sqr(Energy / me_c2);
	double tau = Energy / Rlosses;
	double L = tau * Udownstream;
	
	resetLog();

	const int Ndata = 20;
	double rhoPoints[Ndata];
	double observedFlux[Ndata] = { 0.35, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.6, 0.75, 0.8, 0.9, 1.0, 0.95, 0.8, 0.55, 0.45, 0.25, 0.2, 0.15, 0.15 };
	double observedError[Ndata];
	double energyPoints[Ndata];

	for (int irho = 0; irho < Ndata; ++irho) {
		rhoPoints[irho] = (250 + 0.5 * irho) * distance * pi / (180 * 60 * 60);
		observedError[irho] = 0.1;
		energyPoints[irho] = 1000 * 1.6 * 1E-12;

		observedFlux[irho] = observedFlux[irho] * 1E-30;
		observedError[irho] = observedError[irho] * 1E-30;
	}

	printf("Tycho profile\n");
	printLog("Tycho profile\n");

	printf("Turbulent energy fraction = %g\n", turbulentFraction);
	printLog("Turbulent energy fraction = %g\n", turbulentFraction);

	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInDiskSource(downstreamB, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSphericalSource(downstreamB, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSectoralSphericalSource(downstreamB, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, turbulentFraction, index, lturb, Nmodes, R, Rmin, 2 * pi, anisotropy);
	write3dArrayToFile(B, Nrho, Nz, Nphi, "B.dat");

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 2.0, 100 * Energy);

	//RadiationSource* downstreamSource = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, concentrations, R, (1.0 - widthFraction) * R, distance);
	//RadiationSource* downstreamSource = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, concentrations, R, (1.0 - widthFraction) * R, distance, Udownstream);
	RadiationSourceInCylindrical* source = new TabulatedSectoralSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, 0, concentrations, R, Rin, Rmin, 2 * pi, distance, Udownstream);
	//RadiationSource* downstreamSource = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, upstreamElectrons, downstreamB, theta, concentrations, R, R, distance, Udownstream);

	//number of parameters of the downstreamSource
	const int Nparams = 5;
	//min and max parameters, which defind the region to find minimum. also max parameters are used for normalization of units
	double minParameters[Nparams] = { distance * 250.5 * pi / (180 * 60 * 60), 1E-22, 2E-12, 0.001, 0.3 * speed_of_light };
	double maxParameters[Nparams] = { 1.2E19, 1E-3, 2E10, 0.5, 0.5 * speed_of_light };
	//starting point of optimization and normalization
	double sigma = source->getAverageSigma();

	double vector[Nparams] = { R, sigma, concentration, widthFraction, 0.5 * speed_of_light };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}

	SynchrotronEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax, false);

	bool optPar[Nparams] = { true, true, true, false, false };
	int Niterations = 5;
	LossEvaluator* KPIevaluator = new RadialProfileLossEvaluator(energyPoints[0], observedFlux, observedError, rhoPoints, Ndata, source);
	GradientDescentRadiationOptimizer* optimizer = new GradientDescentRadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, Niterations, KPIevaluator);
	int Npoints[Nparams] = { 5,5,5,10,2 };
	GridEnumRadiationOptimizer* gridOptimizer = new GridEnumRadiationOptimizer(evaluator, minParameters, maxParameters, Nparams, Npoints, KPIevaluator);

	gridOptimizer->optimize(vector, optPar);
	optimizer->optimize(vector, optPar);

	source->resetParameters(vector, maxParameters);
	double error = optimizer->evaluateOptimizationFunction(vector);

	R = vector[0] * maxParameters[0];
	sigma = vector[1] * maxParameters[1];
	concentration = vector[2] * maxParameters[2];
	widthFraction = vector[3] * maxParameters[3];
	Rin = R * (1 - widthFraction);
	Rmin = R * (1 - 2 * widthFraction);
	double arcR = R * 180 * 60 * 60 / (pi * distance);
	double arcMinR = Rmin * 180 * 60 * 60 / (pi * distance);

	printf("optimization function = %g\n", error);
	printLog("optimization function = %g\n", error);
	printf("R = %g\n", R);
	printLog("R = %g\n", R);
	printf("sigma = %g\n", sigma);
	printLog("sigma = %g\n", sigma);
	printf("concentration = %g\n", concentration);
	printLog("concentration = %g\n", concentration);
	printf("width fraction = %g\n", widthFraction);
	printLog("width fraction = %g\n", widthFraction);
	printf("Rmin = %g\n", Rmin);
	printLog("Rmin = %g\n", Rmin);
	printf("arcR = %g\n", arcR);
	printLog("arcR = %g\n", arcR);
	printf("arcMinR = %g\n", arcMinR);
	printLog("arcMinR = %g\n", arcMinR);

	B0 = sqrt(sigma * 4 * pi * massProton * speed_of_light2 * concentration);
	printf("B = %g", B0);
	printLog("B = %g", B0);


	
	evaluator->writeImageFromSourceAtEToFile(energyPoints[0], "image.dat", source);

	FILE* outFile = fopen("outputRadial.dat", "w");
	for (int i = 0; i < Ndata; ++i) {
		double I = 0;
		double rho = rhoPoints[i];
		if (rho < source->getMaxRho() && rho > Rmin) {
			int irho = source->getRhoIndex(rho);
			int iphi = 0;
			double s = source->getCrossSectionArea(irho, iphi);
			I = evaluator->evaluateFluxFromSourceAtPoint(energyPoints[0], source, irho, iphi) / s;
		}
		fprintf(outFile, "%g, %g\n", rho*180*60*60/(pi*distance), I);
	}
	fclose(outFile);

	delete3dArray(B, Nrho, Nz, Nphi);
	delete3dArray(theta, Nrho, Nz, Nphi);
	delete3dArray(phi, Nrho, Nz, Nphi);
}

void evaluateSynchrotronInWideRange() {
	double R = 2.28784E17;
	double B = 0.601848;
	double f = 0.1;
	double electronConcentration = 48.0628;
	double distance = 150 * 1000000 * parsec;

	double velocity = 0.677527 * speed_of_light;
	double downstreamV = 0.25 * velocity;
	
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	
	int Nrho = 1;
	int Nz = 5000;
	int Nphi = 1;


	//TabulatedDiskSourceWithSynchCutoff* downstreamSource = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution, downstreamB, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, velocity);
	TabulatedSLSourceWithSynchCutoff* source = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution, B, pi / 2, 0, electronConcentration, R, (1.0 - f) * R, distance, downstreamV, velocity);

	int Ne = 5000;
	double Emin = me_c2;
	double Emax = me_c2 * 1E8;
	SynchrotronEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true, true);

	double kevFlux = evaluator->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source);

	double mevFlux = evaluator->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source);

	printf("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);
	printLog("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);

	printf("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);
	printLog("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);

	evaluator->writeFluxFromSourceToFile("wideRangeSynch.dat", source, 1E8 * hplank, 20 * MeV, 2000);
}

void evaluateW50bremsstrahlung() {
	double distance = (18000/3.26) * parsec;

	const char* concentrationFileName = "../PLUTO/Tools/pyPLUTO/density.dat";
	const char* BFileName = "../PLUTO/Tools/pyPLUTO/B.dat";
	const char* temperatureFileName = "../PLUTO/Tools/pyPLUTO/T.dat";

	FILE* concentrationFile = fopen(concentrationFileName, "r");
	int Nx, Ny, Nz;
	fscanf(concentrationFile, "%d %d %d", &Nz, &Nx, &Ny);
	Ny = Nz;
	double minX, maxX, minY, maxY, minZ, maxZ;
	fscanf(concentrationFile, "%lf %lf %lf", &minZ, &minX, &minY);
	fscanf(concentrationFile, "%lf %lf %lf", &maxZ, &maxX, &maxY);

	minZ = -maxZ;
	minY = 0;
	minY = -maxZ;
	maxY = maxZ;

	Nx = 200;
	Nz = 200;
	Ny = 200;


	ThermalRectangularSource* source = RadiationSourceFactory::readThermalRectangularSourceFromFile(minX, maxX, minZ, maxZ, minY, maxY, Nx, Nz, Ny, distance, SourceInputGeometry::CYLINDRICAL, BFileName, concentrationFileName, temperatureFileName, 0.8*pi/2, 0, 0);

	BremsstrahlungThermalEvaluator* evaluator = new BremsstrahlungThermalEvaluator(true, false);

	printf("start evaluating spectrum\n");
	evaluator->writeEFEFromSourceToFile("W50bremsstrahlung.dat", source, 1.6E-12, 1.6E-6, 50);
	printf("start writing image\n");
	evaluator->writeImageFromSourceToFile("W50bremsstrahlungImageeV.dat", source, 1.6E-11, 1.6E-10, 20);
	evaluator->writeImageFromSourceToFile("W50bremsstrahlungImageKeV.dat", source, 1.6E-9, 1.6E-8, 20);
	evaluator->writeImageFromSourceToFile("W50bremsstrahlungImageMeV.dat", source, 1.6E-7, 1.6E-6, 20);
}

void evaluateW50synchrotron() {
	double distance = (18000 / 3.26) * parsec;
	const char* fileName = "./examples_data/W50/electrons.dat";

	MassiveParticleTabulatedIsotropicDistribution* electrons;
	double concentration;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons, concentration);

	double size = 1E19;
	double B = 6E-5;

	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, pi / 2, 0, concentration, size, size, distance);
	RadiationEvaluator* evaluator = new SynchrotronEvaluator(100000, me_c2, 1E10 * me_c2, false);
	double cyclotronOmega = electron_charge * B / (massElectron * speed_of_light);
	//evaluator->writeFluxFromSourceToFile("outputSynch.dat", downstreamSource, 10 * hplank * cyclotronOmega, 100000 * hplank * cyclotronOmega, 1000);
	evaluator->writeEFEFromSourceToFile("W50synchrotron.dat", source, 1.6E-18, 1.6E-5, 2000);


}

double* getUvarovBpar(int Nx, double minX, double maxX, double L0) {
	double A0 = 41.0755*1E-6;
	double A1 = 9.27296;
	double A2 = 1.41262;
	double A3 = 3.26888;
	double A4 = 1.88243;
	double A5 = 0.0767564;
	double A6 = 0.00540849;
	double UNIT_LENGTH = L0;
	double* B = new double [Nx];
	double dx = (maxX - minX) / Nx;
	for (int i = 0; i < Nx; ++i) {
		double x = minX + (i + 0.5) * dx;
		double l = (maxX - x) / UNIT_LENGTH + 2;
		//downstreamB[i] = (A0 / (1 + pow(sqr((l - A3) / A1), A2))) * (1.0 / (1 + A6 * exp(-sqr((l - A4) / A5))));
		B[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) / (1.0 + A6 * exp(-(l - A4) / A5));
		if (B[i] < 1E-6 ) {
			B[i] = 1E-6;
		}
	}

	return B;
}

double* getUvarovBper(int Nx, double minX, double maxX, double L0) {
	double A0 = 2*98.3917*1E-6;
	double A1 = 30.533;
	double A2 = 2.33138;
	double A3 = -23.2141;
	double A4 = 18.847;
	double A5 = 1.05756;
	double A6 = 0.695847;
	double UNIT_LENGTH = L0;
	double* B = new double[Nx];
	double dx = (maxX - minX) / Nx;
	for (int i = 0; i < Nx; ++i) {
		double x = minX + (i + 0.5) * dx;
		double l = (maxX - x) / UNIT_LENGTH + 2;
		//downstreamB[i] = A0 / (1 + pow(sqr((l - A3) / A1), A2) + A4 * exp(-sqr((l - A5) / A6)));
		B[i] = (A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) + A4 * exp(-sqr((l - A5) / A6)));
		if (B[i] < 1E-5) {
			B[i] = 1E-5;
		}
	}

	return B;
}

double* getUvarovBpar2(int Nx, double* xgrid, double L0, double factor) {
	factor = factor * 1E-6;
	double A0 = 41.0755 * factor;
	double A1 = 9.27296;
	double A2 = 1.41262;
	double A3 = 3.26888;
	double A4 = 1.88243;
	double A5 = 0.0767564;
	double A6 = 0.00540849;
	double UNIT_LENGTH = L0;
	double* B = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		double x = xgrid[i];
		double shockx = 1.5;
		double l = x / UNIT_LENGTH + shockx;
		if (l >= shockx) {
			//downstreamB[i] = (A0 / (1 + pow(sqr((l - A3) / A1), A2))) * (1.0 / (1 + A6 * exp(-sqr((l - A4) / A5))));
			//downstreamB[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) / (1.0 + A6 * exp(-(l - A4) / A5));
			B[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) / (1.0 + A6 * exp(-(l - A4) / A5));
			if (B[i] < 1E-6) {
				B[i] = 1E-6;
			}
		}
		else {
			B[i] = 1E-6;
		}
	}

	return B;
}

double* getUvarovBper2(int Nx, double* xgrid, double L0, double factor) {
	factor = factor*1E-6;
	double A0 = 98.3917 * factor;
	double A1 = 30.533;
	double A2 = 2.33138;
	double A3 = -23.2141;
	double A4 = 18.847*factor;
	double A5 = 1.05756;
	double A6 = 0.695847;
	double UNIT_LENGTH = L0;
	double* B = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		double x = xgrid[i];
		double shockx = 1.5;
		double l = x / UNIT_LENGTH + shockx;
		if (l >= shockx) {
			//downstreamB[i] = A0 / (1 + pow(sqr((l - A3) / A1), A2) + A4 * exp(-sqr((l - A5) / A6)));
			//downstreamB[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2) + A4 * exp(-sqr((l - A5) / A6)));
			B[i] = (A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) + A4 * exp(-sqr((l - A5) / A6)));
			//downstreamB[i] = (A0 / ((pow(1.0 + sqr((l - A3) / A1), A2)) + A4 * exp(-sqr((l - A5) / A6))));
			if (B[i] < 1E-6) {
				B[i] = 1E-6;
			}
		}
		else {
			B[i] = 1E-6;
		}
	}

	return B;
}

double* getUvarovBpar2new(int Nx, double* xgrid, double L0, double factor) {
	factor = factor * 1E-6;
	double A0 = 7.66093 * factor;
	double A1 = 9.7218;
	double A2 = 2.14722;
	double A3 = 2.01985;
	double A4 = 1.9744;
	double A5 = 0.0736655;
	double A6 = 0.232626;
	double UNIT_LENGTH = L0;
	double* B = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		double x = xgrid[i];
		double shockx = 1.94;
		double l = x / UNIT_LENGTH + shockx;
		if (l >= shockx) {
			//downstreamB[i] = (A0 / (1 + pow(sqr((l - A3) / A1), A2))) * (1.0 / (1 + A6 * exp(-sqr((l - A4) / A5))));
			//downstreamB[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) / (1.0 + A6 * exp(-(l - A4) / A5));
			B[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) / (1.0 + A6 * exp(-(l - A4) / A5));
			if (B[i] < 1E-6) {
				B[i] = 1E-6;
			}
		}
		else {
			B[i] = 1E-6;
		}
	}

	return B;
}

double* getUvarovBper2new(int Nx, double* xgrid, double L0, double factor) {
	factor = factor * 1E-6;
	double A0 = 5.26842 * factor;
	double A1 = 1.05181;
	double A2 = 0.319873;
	double A3 = 2.87845;
	double A4 = 21.0461 * factor;
	double A5 = 0.042243;
	double A6 = 1.21898;
	double UNIT_LENGTH = L0;
	double* B = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		double x = xgrid[i];
		double shockx = 1.94;
		double l = x / UNIT_LENGTH + shockx;
		if (l >= shockx) {
			//downstreamB[i] = A0 / (1 + pow(sqr((l - A3) / A1), A2) + A4 * exp(-sqr((l - A5) / A6)));
			//downstreamB[i] = A0 / (pow(1.0 + sqr((l - A3) / A1), A2) + A4 * exp(-sqr((l - A5) / A6)));
			B[i] = (A0 / (pow(1.0 + sqr((l - A3) / A1), A2)) + A4 * exp(-sqr((l - A5) / A6)));
			//downstreamB[i] = (A0 / ((pow(1.0 + sqr((l - A3) / A1), A2)) + A4 * exp(-sqr((l - A5) / A6))));
			if (B[i] < 1E-6) {
				B[i] = 1E-6;
			}
		}
		else {
			B[i] = 1E-6;
		}
	}

	return B;
}

void evaluateW50comptonAndSynchrotron() {
	double distance = (18000 / 3.26) * parsec;
	const char* fileName = "./examples_data/W50/electrons.dat";

	MassiveParticleTabulatedIsotropicDistribution* electrons;
	double concentration;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons, concentration);

	double size = 5E20;
	double B0 = 6E-5;

	//RadiationSourceInCylindrical* downstreamSource = new SimpleFlatSource(upstreamElectrons, downstreamB, pi / 2, 0, concentration, size, size, distance);
	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photonsIR = new PhotonPlankDistribution(140, 0.8/1800000);
	double photonIRconcentration = photonsIR->getConcentration();
	double photonIRenergyDensity = photonIRconcentration * photonsIR->getMeanEnergy();
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();
	PhotonMultiPlankDistribution* photonsTotal = new PhotonMultiPlankDistribution(2.725, 1.0, 140, 0.8 / 1800000);
	double photonTotalConcentration = photonsTotal->getConcentration();
	double photonTotalEnergyDensity = photonTotalConcentration * photonsTotal->getMeanEnergy();

	int Nrho = 100000;
	int Nz = 1;
	int Ny = 1;

	double* Bpar = getUvarovBpar(Nrho, 0, size, 2.5E18);
	double* Bper = getUvarovBper(Nrho, 0, size, 2.5E18);
	double*** B = new double** [Nrho];
	double*** Btheta = new double** [Nrho];
	double*** Bphi = new double** [Nrho];
	double*** concentrationArray = new double** [Nrho];
	for (int i = 0; i < Nrho; ++i) {
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		concentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
			concentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				B[i][j][k] = sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]);
				//par - x, per - y and z
				Btheta[i][j][k] = atan2(Bper[i]/sqrt(2), sqrt(Bpar[i]*Bpar[i] + 0.5*Bper[i]*Bper[i]));
				Bphi[i][j][k] = pi / 4;
				concentrationArray[i][j][k] = concentration;
			}
		}
	}

	FILE* Bfile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < Nrho; ++i) {
		fprintf(Bfile, "%g %g %g\n", (i + 0.5) * size / Nrho, Bpar[i], Bper[i]);
	}

	fclose(Bfile);

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance, 0.25 * 0.1 * speed_of_light, photonTotalEnergyDensity);
	RectangularSource* source2 = new RectangularSource(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance);

	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(0, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	int Ne = 200;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 250, 0.1*kBoltzman*2.75, 10*kBoltzman*140, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", source, 1.6E-10, 1.6E3, 100);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", source, 1.6E-1, 1.6E3, 100);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", source, 1.6E-9, 50*1.6E-9, 100);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageeV.dat", source, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageKeV.dat", source, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageMeV.dat", source, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageGeV.dat", source, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageTeV.dat", source, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImagePeV.dat", source, 1.6E3, 1.6E4, 20);

	printf("start writing x-E diagram\n");
	printLog("start writing x-E diagram\n");

	double Emin = 1.6E-12;
	double Emax = 1.6E4;

	int Nnu = 200;
	
	double factor = pow(Emax / Emin, 1.0 / (Nnu - 1));
	double currentE = Emin;

	double** F = new double*[Nnu];
	for (int i = 0; i < Nnu; ++i) {
		F[i] = new double[Nrho];
	}

	FILE* Efile = fopen("./output/E_grid.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(Efile, "%g\n", currentE);
		currentE = currentE * factor;
	}
	fclose(Efile);

	currentE = Emin;

	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, source, currentE, factor, sumEvaluator, i)
		for (j = 0; j < Nrho; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0)/source->getCrossSectionArea(j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outFile = fopen("./output/xE.dat", "w");
	for (int i = 0; i < Nrho; ++i) {
		for (int j = 0; j < Nnu; ++j) {
			if (j == 0) {
				fprintf(outFile, "%g", F[j][i]);
			}
			else {
				fprintf(outFile, " %g", F[j][i]);
			}
		}
		fprintf(outFile, "\n");
	}
	fclose(outFile);

}

void evaluateW50comptonAndSynchrotron2() {
	double distance = (18000 / 3.26) * parsec;
	const char* fileName = "./examples_data/W50/electrons.dat";

	MassiveParticleTabulatedIsotropicDistribution* electrons;
	double concentration;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons, concentration);
	double* energy0 = electrons->getEnergyArray();
	double* distribution0 = electrons->getDistributionArray();
	int Nee = electrons->getN();
	double size = 5E20;
	double B0 = 6E-5;

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

	int Nrho = 100000;
	int Nz = 1;
	int Ny = 1;

	double* Bpar = getUvarovBpar(Nrho, 0, size, 2.5E18);
	double* Bper = getUvarovBper(Nrho, 0, size, 2.5E18);
	double*** B = new double** [Nrho];
	double*** Btheta = new double** [Nrho];
	double*** Bphi = new double** [Nrho];
	double*** concentrationArray = new double** [Nrho];
	for (int i = 0; i < Nrho; ++i) {
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		concentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
			concentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				B[i][j][k] = sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]);
				//downstreamB[i][j][k] = downstreamB[i][j][k] / 3;
				//par - x, per - y and z
				Btheta[i][j][k] = atan2(Bper[i] / sqrt(2), sqrt(Bpar[i] * Bpar[i] + 0.5 * Bper[i] * Bper[i]));
				Bphi[i][j][k] = pi / 4;
				concentrationArray[i][j][k] = concentration;
			}
		}
	}

	FILE* Bfile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < Nrho; ++i) {
		fprintf(Bfile, "%g %g %g\n", (i + 0.5) * size / Nrho, Bpar[i], Bper[i]);
	}

	fclose(Bfile);

	double* x = new double[Nrho];
	double dx = size / Nrho;

	double* inversedB = new double[Nrho];
	for (int i = 0; i < Nrho; ++i) {
		x[i] = dx * (i + 0.5);
		inversedB[i] = B[Nrho - i - 1][0][0];
	}

	double** energy = new double* [Nrho];
	double** distributions = new double* [Nrho];
	for (int i = 0; i < Nrho; ++i) {
		energy[i] = new double[Nee];
		distributions[i] = new double[Nee];
	}

	MassiveParticleDistributionFactory::evaluateDistributionAdvectionWithLosses(massElectron, energy0, distribution0, energy, distributions, Nee, Nrho, x, 0.25 * 0.1 * speed_of_light, inversedB, photonEnergyDensity, photons->getMeanEnergy(), photonIRenergyDensity, photonsIR->getMeanEnergy());

	MassiveParticleDistribution**** electrons2 = new MassiveParticleDistribution * **[Nrho];
	for (int i = 0; i < Nrho; ++i) {
		electrons2[i] = new MassiveParticleDistribution * *[Nz];
		for (int j = 0; j < Nz; ++j) {
			electrons2[i][j] = new MassiveParticleDistribution * [Ny];
			for (int k = 0; k < Ny; ++k) {
				electrons2[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energy[Nrho - i - 1], distributions[Nrho - i - 1], Nee, DistributionInputType::ENERGY_FE);
			}
		}
	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(Nrho, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, size, 0, pi * size, distance, 0.25 * 0.1 * speed_of_light, photonTotalEnergyDensity);
	RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nrho, Ny, Nz, electrons2, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance);
	RectangularSource* source2 = new RectangularSource(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(0, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	int Ne = 100;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 20 * kBoltzman * 140, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("./output/W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("./output/W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", source, 1.6E-18, 1.6E3, 200);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", source, 1.6E-1, 1.6E3, 200);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", source, 1.6E-9, 50*1.6E-9, 200);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageeV.dat", source, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageKeV.dat", source, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageMeV.dat", source, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageGeV.dat", source, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageTeV.dat", source, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImagePeV.dat", source, 1.6E3, 1.6E4, 20);

	printf("start writing x-E diagram\n");
	printLog("start writing x-E diagram\n");

	double Emin = 1.6E-12;
	double Emax = 1.6E4;

	int Nnu = 200;

	double factor = pow(Emax / Emin, 1.0 / (Nnu - 1));
	double currentE = Emin;

	double** F = new double* [Nnu];
	for (int i = 0; i < Nnu; ++i) {
		F[i] = new double[Nrho];
	}

	FILE* Efile = fopen("./output/E_grid.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(Efile, "%g\n", currentE);
		currentE = currentE * factor;
	}
	fclose(Efile);

	currentE = Emin;

	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, source, currentE, factor, sumEvaluator, i)
		for (j = 0; j < Nrho; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0) / source->getCrossSectionArea(j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outFile = fopen("./output/xE.dat", "w");
	for (int i = 0; i < Nrho; ++i) {
		for (int j = 0; j < Nnu; ++j) {
			if (j == 0) {
				fprintf(outFile, "%g", F[j][i]);
			}
			else {
				fprintf(outFile, " %g", F[j][i]);
			}
		}
		fprintf(outFile, "\n");
	}
	fclose(outFile);

}

void evaluateW50comptonAndSynchrotronMCfunctionUpstream() {
	double distance = (18000 / 3.26) * parsec;
	const char* distributionFileName = "./examples_data/W50/B15FEB6E18/pdf_sf.dat";
	const char* xfileName = "./examples_data/W50/B15FEB6E18/x_grid.dat";
	const char* pfileName = "./examples_data/W50/B15FEB6E18/p_grid.dat";
	const char* BfileName = "./examples_data/W50/B15FEB6E18/Beff.dat";

	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	double* energy;
	double* xgrid1;
	double* concentration;
	double** distributions;

	int Nenergy;
	int Nx1;
	int Nx;

	double electronToProtonCorrection = 3E-7;

	MassiveParticleDistributionFactory::readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(massElectron, xfileName, pfileName, distributionFileName, xgrid1, energy, distributions, concentration, Nenergy, Nx1);

	int zeroIndex = 0;
	for (int i = 0; i < Nx1; ++i) {
		if (xgrid1[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}
	double downstreamSize = 1E20;
	double upstreamSize = 1.6E20;
	int maxIndex = Nx1 - 1;
	for (int i = 0; i < Nx1; ++i) {
		if (xgrid1[i] >= downstreamSize) {
			maxIndex = i;
			break;
		}
	}
	//maxIndex = Nx1 - 1;
	int minIndex = 0;
	for (int i = 0; i < Nx1; ++i) {
		if (xgrid1[i] >= -upstreamSize) {
			minIndex = i;
			break;
		}
	}
	//minIndex = 0;

	int downstreamNx = maxIndex + 1 - zeroIndex;
	int upstreamNx = zeroIndex - minIndex;

	double* downstreamXgrid = new double[downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[downstreamNx - i - 1] = xgrid1[i + zeroIndex];
	}

	Nx = maxIndex - minIndex + 1;

	double* xgrid = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		xgrid[Nx - i - 1] = -xgrid1[i + minIndex];
	}

	double size = 0.5 * fabs(headMaxX);
	double B0 = 6E-5;

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

	double L0 = 0.5E18;
	double* Bpar1 = getUvarovBpar2(downstreamNx, downstreamXgrid, L0, 0.25);
	double* Bper1 = getUvarovBper2(downstreamNx, downstreamXgrid, L0, 0.25);

	for (int i = 0; i < downstreamNx; ++i) {
		if (sqrt(Bpar1[i] * Bpar1[i] + 2 * Bper1[i] * Bper1[i]) < 1E-5) {
			Bpar1[i] = 1E-5 / sqrt(3.0);
			Bper1[i] = 1E-5 / sqrt(3.0);
		}
	}

	//double* Bpar = getUvarovBpar2(Nx, xgrid, L0);
	//double* Bper = getUvarovBper2(Nx, xgrid, L0);
	double*** B = new double** [Nx];
	double*** Btheta = new double** [Nx];
	double*** Bphi = new double** [Nx];
	double*** concentrationArray = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		concentrationArray[Nx - i - 1] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
			concentrationArray[Nx - i - 1][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				//B[i][j][k] = sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]);
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				//Btheta[i][j][k] = atan2(Bper[i] / sqrt(2), sqrt(Bpar[i] * Bpar[i] + 0.5 * Bper[i] * Bper[i]));
				//downstreamBtheta[i][j][k] = pi / 2;
				//Bphi[i][j][k] = pi / 4;
				concentrationArray[Nx - i - 1][j][k] = concentration[i + minIndex]*electronToProtonCorrection;
			}
		}
	}

	double* Beff = new double[Nx1];
	FILE* Bfile = fopen(BfileName, "r");
	for (int i = 0; i < Nx1; ++i) {
		fscanf(Bfile, "%lf", &Beff[i]);
	}
	fclose(Bfile);

	double* Bpar = new double[Nx];
	double* Bper = new double[Nx];

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				if (i >= Nx - downstreamNx) {
					Bpar[Nx - i - 1] = Bpar1[Nx - i - 1];
					Bper[Nx - i - 1] = Bper1[Nx - i - 1];
					B[Nx - i - 1][j][k] = sqrt(Bpar1[Nx - i - 1] * Bpar1[Nx - i - 1] + 2 * Bper1[Nx - i - 1] * Bper1[Nx - i - 1]);
					Btheta[i][j][k] = atan2(sqrt(Bpar1[Nx - i - 1] * Bpar1[Nx - i - 1] + Bper1[Nx - i - 1] * Bper1[Nx - i - 1]), Bper1[Nx - i - 1]);
					Bphi[i][j][k] = atan2(Bpar1[Nx - i - 1], Bper1[Nx - i - 1]);
				}
				else {
					Bpar[Nx - i - 1] = 0;
					Bper[Nx - i - 1] = Beff[i + minIndex]/sqrt(2.0);
					B[Nx - i - 1][j][k] = Beff[i + minIndex];
					Btheta[i][j][k] = pi / 4;
					Bphi[i][j][k] = 0;
				}
				if ((B[Nx - i - 1][j][k] != B[Nx - i - 1][j][k]) || (0 * B[Nx - i - 1][j][k] != 0 * B[Nx - i - 1][j][k])) {
					printf("B = NaN\n");
					printLog("B = NaN\n");
					exit(0);
				}
			}
		}
	}

	FILE* BoutFile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < Nx; ++i) {
		fprintf(BoutFile, "%g %g %g\n", xgrid[i], Bpar[i], Bper[i]);
	}

	fclose(BoutFile);

	MassiveParticleDistribution**** electrons2 = new MassiveParticleDistribution * **[Nx];
	for (int i = 0; i < Nx; ++i) {
		electrons2[Nx - 1 - i] = new MassiveParticleDistribution * *[Nz];
		for (int j = 0; j < Nz; ++j) {
			electrons2[Nx - 1 - i][j] = new MassiveParticleDistribution * [Ny];
			for (int k = 0; k < Ny; ++k) {
				electrons2[Nx - 1 - i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energy, distributions[i + minIndex], Nenergy, DistributionInputType::ENERGY_FE);
				//electrons2[i][j][k] = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
			}
		}

	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(Nrho, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, size, 0, pi * size, distance, 0.25 * 0.1 * speed_of_light, photonTotalEnergyDensity);
	RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nx, xgrid, Ny, Nz, electrons2, B, Btheta, Bphi, concentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* source2 = new RectangularSource(Nrho, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx/2 + 200, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx/2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx/2 - 200, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(outXfile, "%g\n", xgrid[i]);
	}
	fclose(outXfile);

	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));
	FILE* outPfile = fopen("./output/p_grid.dat", "w");
	double p = pmin;
	for (int i = 0; i < Np; ++i) {
		fprintf(outPfile, "%g\n", p * massElectron / massProton);
		p = p * factorp;
	}
	fclose(outPfile);

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < Nx; ++i) {
	p = pmin;
	MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(i, 0, 0));
	//MassiveParticleIsotropicDistribution* distribution = frontElectrons;
	for (int j = 0; j < Np; ++j) {
		double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
		double F = distribution->distributionNormalized(E)*concentrationArray[i][0][0];
		F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
		fprintf(outDistributionFile, "%g %g\n", p, F);
		p = p * factorp;
	}
	}
	fclose(outDistributionFile);

	FILE* concentrationFile = fopen("./output/concentration.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(concentrationFile, "%g %g\n", xgrid[i], concentrationArray[i][0][0]);
	}
	fclose(concentrationFile);


	int Ne = 10;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", downstreamSource, 1.6E-12, 1.6E4, 1000);


	double Ephmin = 1.6E-12;
	double Ephmax = 1.6E4;
	int Nph = 1000;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	FILE* outFile = fopen("./output/W50synchandcompt.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		//omp_set_lock(&my_lock);
		printf("writeEFEFromSourceToFile iph = %d\n", i);
		printLog("writeEFEFromSourceToFile iph = %d\n", i);
		//omp_unset_lock(&my_lock);
		//double flux1 = sumEvaluator->evaluateFluxFromSource(currentE, downstreamSource);
		double flux1 = 0;
		double fluxHead = 0;
		double fluxCone = 0;
		int j;
#pragma omp parallel for private(j) shared(currentE, source, xgrid, Nx) reduction(+:flux1, fluxHead, fluxCone)
		for (j = 0; j < Nx; ++j) {
			double flux = sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0);
			flux1 += flux;
			if ((xgrid[j] >= headMaxX) && (xgrid[j] <= headMinX)) {
				fluxHead += flux;
			}
			if ((xgrid[j] >= coneMaxX) && (xgrid[j] <= coneMinX)) {
				fluxCone += flux;
			}
		}
		fprintf(outFile, "%g %g %g %g\n", currentE / 1.6E-12, currentE * (flux1), currentE * fluxHead, currentE * fluxCone);
		currentE = currentE * factor;
	}
	fclose(outFile);

	double* profileXMM = new double[Nx];
	double* profileNuSTAR = new double[Nx];
	int irho;

	Ephmin = 0.3 * 1000 * 1.6E-12;
	Ephmax = 10 * 1000 * 1.6E-12;
	Nph = 20;
	omp_lock_t lock;
	omp_init_lock(&lock);
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, source, Nx, sumEvaluator, Nph, profileXMM, lock)
	for (irho = 0; irho < Nx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = source->getCrossSectionArea(irho, 0);
		double d = source->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileXMM[irho] = localFlux;
	}

	FILE* xmmFile = fopen("./output/xmmprofile.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(xmmFile, "%g %g\n", xgrid[i], profileXMM[i]);
	}
	fclose(xmmFile);

	Ephmin = 10 * 1000 * 1.6E-12;
	Ephmax = 20 * 1000 * 1.6E-12;
	Nph = 20;
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, source,Nx, sumEvaluator, Nph, profileNuSTAR, lock)
	for (irho = 0; irho < Nx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating image irho = %d\n", irho);
		printLog("evaluating image irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = source->getCrossSectionArea(irho, 0);
		double d = source->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileNuSTAR[irho] = localFlux;
	}

	omp_destroy_lock(&lock);

	FILE* nustarFile = fopen("./output/nustarprofile.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(nustarFile, "%g %g\n", xgrid[i], profileNuSTAR[i]);
	}
	fclose(nustarFile);

}

void evaluateW50comptonAndSynchrotronAdvectionfunction() {
	double distance = (18000 / 3.26) * parsec;

	double* energy;
	double* xgrid;
	double* concentration;
	double** distributions;

	int Nenergy;
	int Nx;

	const char* xfileName = "./examples_data/W50/lowfield/x_grid.dat";

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
	int maxIndex = Nx - 1;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 1E20) {
			maxIndex = i;
			break;
		}
	}

	Nx = maxIndex + 1 - zeroIndex;

	xgrid = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		xgrid[Nx - i - 1] = xgrid1[i + zeroIndex];
	}

	double size = 5E20;
	double B0 = 6E-5;

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

	double L0 = 1E18;
	double* Bpar = getUvarovBpar2(Nx, xgrid, L0, 0.25);
	double* Bper = getUvarovBper2(Nx, xgrid, L0, 0.25);
	double*** B = new double** [Nx];
	double*** Btheta = new double** [Nx];
	double*** Bphi = new double** [Nx];
	double*** concentrationArray = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		concentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
			concentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				B[i][j][k] = sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]);
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				Btheta[i][j][k] = atan2(Bper[i] , sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]));
				//downstreamBtheta[i][j][k] = pi / 2;
				Bphi[i][j][k] = pi / 4;
				concentrationArray[i][j][k] =1.0;
			}
		}
	}

	for (int i = 0; i < Nx; ++i) {
		xgrid[i] = -xgrid[i];
	}

	FILE* Bfile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < Nx; ++i) {
		//fprintf(Bfile, "%g %g %g\n", -downstreamXgrid[i]/L0 + 1.5, Bpar[i], Bper[i]);
		fprintf(Bfile, "%g %g %g\n", xgrid[i], Bpar[i], Bper[i]);
	}

	fclose(Bfile);

	const char* fileName = "./examples_data/W50/lowfield/GLE_pdf_sf1.dat";

	MassiveParticleTabulatedIsotropicDistribution* electrons1;
	double concentration1;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons1, concentration1);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				concentrationArray[i][j][k] = concentration1;
			}
		}
	}

	MassiveParticleIsotropicDistribution* distribution0 = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
	MassiveParticleTabulatedIsotropicDistribution* electrons2 = new MassiveParticleTabulatedIsotropicDistribution(distribution0, me_c2, 1E16 * 1.6E-12, 1000);
	//upstreamElectrons = frontElectrons;
	double E0 = 1.6E-1;
	MassiveParticleIsotropicDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, E0, 0.01 * E0);

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				concentrationArray[i][j][k] = concentration1*electrons2->distributionNormalized(E0)*0.02*E0;
			}
		}
	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, xgrid, Ny, Nz, electrons1, B, Btheta, Bphi, concentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(Nx, downstreamXgrid, Ny, Nz, electrons2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* downstreamSource = new RectangularSource(1, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, downstreamXgrid[0], downstreamXgrid[Nx - 1], 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(outXfile, "%g\n", xgrid[i]);
	}
	fclose(outXfile);

	double pmin = 0.1*massProton/massElectron;
	double pmax = 5E6*massProton/massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));
	FILE* outPfile = fopen("./output/p_grid.dat", "w");
	double p = pmin;
	for (int i = 0; i < Np; ++i) {
		fprintf(outPfile, "%g\n", p*massElectron/massProton);
		p = p * factorp;
	}
	fclose(outPfile);

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - i - 1, 0, 0));
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E);
			F = (F * p * p * p * me_c2 * me_c2 / E)*massElectron/massProton;
			fprintf(outDistributionFile, "%g\n", F);
			p = p * factorp;
		}
	}
	fclose(outDistributionFile);


	int Ne = 10;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2*500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2*500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2*500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", source, 1.6E-12, 1.6E4, 1000);

	double* profileXMM = new double[Nx];
	double* profileNuSTAR = new double[Nx];
	int irho;

	double Ephmin = 0.3 * 1000 * 1.6E-12;
	double Ephmax = 10 * 1000 * 1.6E-12;
	int Nph = 20;
	omp_lock_t lock;
	omp_init_lock(&lock);
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, source, Nx, sumEvaluator, Nph, profileXMM, lock)
	for (irho = 0; irho < Nx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
			double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
			double currentE = Ephmin;
			double localFlux = 0;
			double s = source->getCrossSectionArea(irho, 0);
			double d = source->getDistance();
			for (int ie = 0; ie < Nph; ++ie) {
				double dE = currentE * (factor - 1.0);
				localFlux += sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, irho, 0) * dE * d * d / s;
				currentE = currentE * factor;
			}
			profileXMM[irho] = localFlux;
	}
	FILE* xmmFile = fopen("./output/xmmprofile.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(xmmFile, "%g %g\n", xgrid[i], profileXMM[i]);
	}
	fclose(xmmFile);

	Ephmin = 10 * 1000 * 1.6E-12;
	Ephmax = 20 * 1000 * 1.6E-12;
	Nph = 20;
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, source, Nx, sumEvaluator, Nph, profileNuSTAR, lock)
	for (irho = 0; irho < Nx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating image irho = %d\n", irho);
		printLog("evaluating image irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = source->getCrossSectionArea(irho, 0);
		double d = source->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileNuSTAR[irho] = localFlux;
	}
	omp_destroy_lock(&lock);


	FILE* nustarFile = fopen("./output/nustarprofile.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(nustarFile, "%g %g\n", xgrid[i], profileNuSTAR[i]);
	}
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", downstreamSource, 1.6E-1, 1.6E3, 300);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", downstreamSource, 1.6E-9, 50*1.6E-9, 300);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	return;

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageeV.dat", source, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageKeV.dat", downstreamSource, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageMeV.dat", downstreamSource, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageGeV.dat", downstreamSource, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageTeV.dat", downstreamSource, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("../output/W50scImagePeV.dat", downstreamSource, 1.6E3, 1.6E4, 20);

	printf("start writing x-E diagram\n");
	printLog("start writing x-E diagram\n");

	double Emin = 1.6E-12;
	double Emax = 1.6E4;

	int Nnu = 200;

	double factor = pow(Emax / Emin, 1.0 / (Nnu - 1));
	double currentE = Emin;

	double** F = new double* [Nnu];
	for (int i = 0; i < Nnu; ++i) {
		F[i] = new double[Nx];
	}

	FILE* Efile = fopen("./output/E_grid.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(Efile, "%g\n", currentE);
		currentE = currentE * factor;
	}
	fclose(Efile);

	currentE = Emin;
	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, source, currentE, factor, sumEvaluator, i)
		for (j = 0; j < Nx; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0) / source->getCrossSectionArea(j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outFile = fopen("./output/xE.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nnu; ++j) {
			if (j == 0) {
				fprintf(outFile, "%g", F[j][i]);
			}
			else {
				fprintf(outFile, " %g", F[j][i]);
			}
		}
		fprintf(outFile, "\n");
	}
	fclose(outFile);

}

void evaluateW50comptonThickRegime() {
	double distance = (18000 / 3.26) * parsec;

	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	double* energy;
	double* xgrid;
	double* concentration;
	double** distributions;

	int Nenergy;

	double size = 0.5 * fabs(headMaxX);
	double B0 = 3E-6;
	double magneticEnergyDensity = B0 * B0 / (8 * pi);

	//RadiationSourceInCylindrical* downstreamSource = new SimpleFlatSource(upstreamElectrons, downstreamB, pi / 2, 0, concentration, size, size, distance);
	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photonsIR = new PhotonPlankDistribution(30, 0.8 / 1000000000);
	double photonIRconcentration = photonsIR->getConcentration();
	double photonIRenergyDensity = photonIRconcentration * photonsIR->getMeanEnergy();
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();
	PhotonMultiPlankDistribution* photonsTotal = new PhotonMultiPlankDistribution(2.725, 1.0, 30, 0.8 / 100000000);
	//PhotonMultiPlankDistribution* photonsTotal = PhotonMultiPlankDistribution::getGalacticField();

	double photonTotalConcentration = photonsTotal->getConcentration();
	double photonTotalEnergyDensity = photonTotalConcentration * photonsTotal->getMeanEnergy();

	const char* fileName = "./examples_data/W50/newdistribution/electrons.dat";
	const char* farFileName = "./examples_data/w50/newdistribution/fardownstreamelectrons.dat";
	const char* farUpFileName = "./examples_data/w50/newdistribution/farupstreamelectrons.dat";
	const char* protonsFileName = "./examples_data/W50/newdistribution/protons.dat";




	MassiveParticleTabulatedIsotropicDistribution* electrons1;
	double concentration1;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons1, concentration1);

	MassiveParticleTabulatedIsotropicDistribution* frontProtons;
	double concentration3;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileName, frontProtons, concentration3);

	double electronToProtonCorrection = concentration3 * frontProtons->getDistributionArray()[70] / (concentration1 * electrons1->getDistributionArray()[70]);

	MassiveParticleTabulatedIsotropicDistribution* fardownstreamDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, farFileName, DistributionInputType::ENERGY_FE);
	double farUpstreamConcentration;
	MassiveParticleTabulatedIsotropicDistribution* farupstreamDistribution = MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionAndConcentration(massElectron, farUpFileName, DistributionInputType::ENERGY_FE, farUpstreamConcentration);

	//frontElectrons->writeDistribution("./output/thinDistribution.dat", 2000, me_c2, 1E10 * me_c2);
	FILE* outDistributionFile = fopen("./output/thinDistribution.dat", "w");
	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 2000;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));
	double p = pmin;
	for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = electrons1->distributionNormalized(E);
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g %g\n", p, F);
			p = p * factorp;
	}
	fclose(outDistributionFile);
	double norm = 0;
	electrons1->transformToThickRegime(photonTotalEnergyDensity+magneticEnergyDensity, norm);
	double norm2 = 0;
	double concentration2 = concentration1;
	fardownstreamDistribution->transformToThickRegime(photonTotalEnergyDensity+magneticEnergyDensity, norm2);
	double norm3 = 0;
	farupstreamDistribution->transformToThickRegime(photonTotalEnergyDensity+magneticEnergyDensity, norm3);
	//frontElectrons->writeDistribution("./output/thickDistribution.dat", 2000, me_c2, 1E10 * me_c2);
	outDistributionFile = fopen("./output/thickDistribution.dat", "w");
	p = pmin;
	for (int j = 0; j < Np; ++j) {
		double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
		double F = electrons1->distributionNormalized(E);
		F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
		fprintf(outDistributionFile, "%g %g\n", p, F);
		p = p * factorp;
	}
	fclose(outDistributionFile);
	double E0 = 1.6E-1;

	double u = 0.26 * 0.15 * speed_of_light;
	double outOfJetFactor = 0.1;
	concentration1 *= outOfJetFactor*u * pi * size * size*electronToProtonCorrection*norm;
	concentration2 *= u * pi * size * size * electronToProtonCorrection * norm2;

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, downstreamXgrid, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(Nx, downstreamXgrid, Ny, Nz, electrons2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	RectangularSource* source = new RectangularSource(1, 1, 1, electrons1, B0, pi / 2, 0, concentration1, 0, 1, 0, 1, 0, 1, distance);
	RectangularSource* source2 = new RectangularSource(1, 1, 1, fardownstreamDistribution, B0, pi / 2, 0, concentration2, 0, 1, 0, 1, 0, 1, distance);


	int Ne = 10000;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 20000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false, false);
	RadiationEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50thickCompton.dat", source, 1.6E-12, 1.6E4,1000);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50thickCompton2.dat", source2, 1.6E-12, 1.6E4, 1000);

	return;
}

void evaluateW50comptonAdvectionBigSource() {
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
		xgrid[Nx - i - 1] = xgrid1[i + zeroIndex];
	}

	double size = 5E20;
	double B0 = 0.0;

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

	MassiveParticleIsotropicDistribution* distribution0 = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
	MassiveParticleTabulatedIsotropicDistribution* electrons2 = new MassiveParticleTabulatedIsotropicDistribution(distribution0, me_c2, 1E16 * 1.6E-12, 1000);
	//upstreamElectrons = frontElectrons;
	//double E0 = 1.6E-1;
	//MassiveParticleIsotropicDistribution* upstreamElectrons = new MassiveParticleMonoenergeticDistribution(massElectron, E0, 0.01 * E0);

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
	RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, xgrid, Ny, Nz, electrons2, B, Btheta, Bphi, concentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(Nx, downstreamXgrid, Ny, Nz, electrons2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* downstreamSource = new RectangularSource(1, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, downstreamXgrid[0], downstreamXgrid[Nx - 1], 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(outXfile, "%g\n", xgrid[i]);
	}
	fclose(outXfile);

	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));
	FILE* outPfile = fopen("./output/p_grid.dat", "w");
	double p = pmin;
	for (int i = 0; i < Np; ++i) {
		fprintf(outPfile, "%g\n", p * massElectron / massProton);
		p = p * factorp;
	}
	fclose(outPfile);

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - i - 1, 0, 0));
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E);
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g\n", F);
			p = p * factorp;
		}
	}
	fclose(outDistributionFile);


	int Ne = 1000;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 1000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	comptonEvaluator->writeEFEFromSourceToFile("./output/W50comptonBigSource.dat", source, 1.6E-10, 1.6E4, 1000);

}

void evaluateW50comptonAndSynchrotronMCwithoutupstream() {
	double distance = (18000 / 3.26) * parsec;
	const char* distributionFileName = "./examples_data/W50/lowfield/pdf_sf.dat";
	const char* xfileName = "./examples_data/W50/lowfield/x_grid.dat";
	const char* pfileName = "./examples_data/W50/lowfield/p_grid.dat";

	double* energy;
	double* xgrid1;
	double* concentration1;
	double** distributions1;

	int Nenergy;
	int Nx;

	MassiveParticleDistributionFactory::readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(massElectron, xfileName, pfileName, distributionFileName, xgrid1, energy, distributions1, concentration1, Nenergy, Nx);

	int zeroIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}

	Nx = Nx - zeroIndex;

	double* xgrid = new double[Nx];
	double* concentration = new double[Nx];
	double** distributions = new double* [Nx];
	for (int i = 0; i < Nx; ++i) {
		xgrid[Nx - i - 1] = xgrid1[i + zeroIndex];
		concentration[Nx - i - 1] = concentration1[i + zeroIndex];
		distributions[Nx - i - 1] = distributions1[i + zeroIndex];
	}
	delete[] xgrid1;
	delete[] concentration1;
	delete[] distributions1;

	double size = 5E20;
	double B0 = 6E-5;

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

	double* Bpar = getUvarovBpar2(Nx, xgrid, 2.5E18, 0.25);
	double* Bper = getUvarovBper2(Nx, xgrid, 2.5E18, 0.25);
	double*** B = new double** [Nx];
	double*** Btheta = new double** [Nx];
	double*** Bphi = new double** [Nx];
	double*** concentrationArray = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		B[i] = new double* [Nz];
		Btheta[i] = new double* [Nz];
		Bphi[i] = new double* [Nz];
		concentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			B[i][j] = new double[Ny];
			Btheta[i][j] = new double[Ny];
			Bphi[i][j] = new double[Ny];
			concentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				B[i][j][k] = sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]);
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				Btheta[i][j][k] = atan2(Bper[i] / sqrt(2), sqrt(Bpar[i] * Bpar[i] + 0.5 * Bper[i] * Bper[i]));
				//downstreamBtheta[i][j][k] = pi / 2;
				Bphi[i][j][k] = pi / 4;
				concentrationArray[i][j][k] = concentration[i];
			}
		}
	}

	for (int i = 0; i < Nx; ++i) {
		xgrid[i] = -xgrid[i];
	}

	MassiveParticleDistribution**** electrons = new MassiveParticleDistribution * **[Nx];
	for (int i = 0; i < Nx; ++i) {
		electrons[i] = new MassiveParticleDistribution * *[Nz];
		for (int j = 0; j < Nz; ++j) {
			electrons[i][j] = new MassiveParticleDistribution * [Ny];
			for (int k = 0; k < Ny; ++k) {
				electrons[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energy, distributions[i], Nenergy, DistributionInputType::ENERGY_FE);
				//electrons2[i][j][k] = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
			}
		}
	}

	FILE* Bfile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < Nx; ++i) {
		fprintf(Bfile, "%g %g %g\n", xgrid[i], Bpar[i], Bper[i]);
	}

	fclose(Bfile);

	const char* fileName = "./examples_data/W50/electrons.dat";


	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, downstreamXgrid, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nx, xgrid, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* source2 = new RectangularSource(Nrho, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		fprintf(outXfile, "%g\n", xgrid[i]);
	}
	fclose(outXfile);

	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));
	FILE* outPfile = fopen("./output/p_grid.dat", "w");
	double p = pmin;
	for (int i = 0; i < Np; ++i) {
		fprintf(outPfile, "%g\n", p * massElectron / massProton);
		p = p * factorp;
	}
	fclose(outPfile);

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx - i - 1, 0, 0));
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E);
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g\n", F);
			p = p * factorp;
		}
	}
	fclose(outDistributionFile);


	int Ne = 200;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", source, 1.6E-12, 1.6E4, 300);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", source, 1.6E-1, 1.6E3, 300);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", source, 1.6E-9, 50 * 1.6E-9, 300);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageeV.dat", downstreamSource, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageKeV.dat", downstreamSource, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageMeV.dat", downstreamSource, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageGeV.dat", downstreamSource, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageTeV.dat", downstreamSource, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("../output/W50scImagePeV.dat", downstreamSource, 1.6E3, 1.6E4, 20);

	printf("start writing x-E diagram\n");
	printLog("start writing x-E diagram\n");

	double Emin = 1.6E-12;
	double Emax = 1.6E4;

	int Nnu = 200;

	double factor = pow(Emax / Emin, 1.0 / (Nnu - 1));
	double currentE = Emin;

	double** F = new double* [Nnu];
	for (int i = 0; i < Nnu; ++i) {
		F[i] = new double[Nx];
	}

	FILE* Efile = fopen("./output/E_grid.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(Efile, "%g\n", currentE);
		currentE = currentE * factor;
	}
	fclose(Efile);

	currentE = Emin;
	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, source, currentE, factor, sumEvaluator, i)
		for (j = 0; j < Nx; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0) / source->getCrossSectionArea(j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outFile = fopen("./output/xE.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nnu; ++j) {
			if (j == 0) {
				fprintf(outFile, "%g", F[j][i]);
			}
			else {
				fprintf(outFile, " %g", F[j][i]);
			}
		}
		fprintf(outFile, "\n");
	}
	fclose(outFile);

}

void evaluateW50comptonAndSynchrotronAdvectionfunctionWithUpstream() {
	double distance = (18000 / 3.26) * parsec;

	double* energy;
	double* downstreamXgrid;
	double* upstreamXgrid;
	double* concentration;
	double** distributions;

	int Nenergy;
	int Nx;

	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	const char* xfileName = "./examples_data/W50/newdistribution/x_grid.dat";
	const char* BfileName = "./examples_data/W50/B15FEB6E18/Beff.dat";

	const char* distributionFileName = "./examples_data/W50/newdistribution/electrons_full.dat";
	const char* pfileName = "./examples_data/W50/newdistribution/p_grid.dat";

	const char* fileName = "./examples_data/W50/newdistribution/electrons.dat";
	const char* protonsFileName = "./examples_data/W50/newdistribution/protons.dat";

	/*Nx = 0;
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
	fclose(xfile);*/

	double* xgrid1;

	double* concentration1;
	double** distributions1;

	//double electronToProtonCorrection = 3E-7;

	MassiveParticleTabulatedIsotropicDistribution* frontElectrons;
	double concentration2;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, frontElectrons, concentration2);
	//frontElectrons = new MassiveParticleTabulatedIsotropicDistribution(new MassiveParticlePowerLawDistribution(massElectron, 2.0, me_c2), me_c2, 1600, 1000);
	MassiveParticleTabulatedIsotropicDistribution* frontProtons;
	double concentration3;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileName, frontProtons, concentration3);

	double electronToProtonCorrection = concentration3 * frontProtons->getDistributionArray()[70] / (concentration2 * frontElectrons->getDistributionArray()[70]);

	MassiveParticleDistributionFactory::readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(massElectron, xfileName, pfileName, distributionFileName, xgrid1, energy, distributions1, concentration1, Nenergy, Nx);

	double* Beff = new double[Nx];
	FILE* Bfile = fopen(BfileName, "r");
	for (int i = 0; i < Nx; ++i) {
		fscanf(Bfile, "%lf", &Beff[i]);
	}
	fclose(Bfile);

	int zeroIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}
	double downstreamSize = 1E20;
	double upstreamSize = 1.6E20;
	int maxIndex = Nx - 1;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 1E20) {
			maxIndex = i;
			break;
		}
	}
	//maxIndex = Nx - 1;
	int minIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= -upstreamSize) {
			minIndex = i;
			break;
		}
	}
	//minIndex = 0;

	int downstreamNx = maxIndex + 1 - zeroIndex;
	int upstreamNx = zeroIndex - minIndex;

	downstreamXgrid = new double[downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[downstreamNx - i - 1] = xgrid1[i + zeroIndex];
	}

	upstreamXgrid = new double[upstreamNx];
	double* upstreamB1 = new double[upstreamNx];
	double* upstreamConcentration1 = new double[upstreamNx];
	double** upstreamDistributions1 = new double* [upstreamNx];
	for (int i = 0; i < upstreamNx; ++i) {
		upstreamXgrid[upstreamNx - 1 - i] = -xgrid1[i + minIndex];
		upstreamB1[upstreamNx - 1 - i] = Beff[i + minIndex];
		upstreamConcentration1[upstreamNx - 1 - i] = concentration1[i + minIndex];
		upstreamDistributions1[upstreamNx - 1 - i] = distributions1[i + minIndex];
	}

	double size = 0.5*fabs(headMaxX);
	double B0 = 6E-5;

	//RadiationSourceInCylindrical* downstreamSource = new SimpleFlatSource(upstreamElectrons, downstreamB, pi / 2, 0, concentration, size, size, distance);
	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photonsIR = new PhotonPlankDistribution(30, 0.8 / 10000);
	double photonIRconcentration = photonsIR->getConcentration();
	double photonIRenergyDensity = photonIRconcentration * photonsIR->getMeanEnergy();
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();
	PhotonMultiPlankDistribution* photonsTotal = new PhotonMultiPlankDistribution(2.725, 1, 30, 0.8 / 10000);
	//PhotonMultiPlankDistribution* photonsTotal = PhotonMultiPlankDistribution::getGalacticField();

	double photonTotalConcentration = photonsTotal->getConcentration();
	double photonTotalEnergyDensity = photonTotalConcentration * photonsTotal->getMeanEnergy();

	int Nz = 1;
	int Ny = 1;

	double L0 = 0.3E18;
	double* Bpar = getUvarovBpar2(downstreamNx, downstreamXgrid, L0, 0.125);
	double* Bper = getUvarovBper2(downstreamNx, downstreamXgrid, L0, 0.125);
	//double* Bpar = getUvarovBpar2new(downstreamNx, downstreamXgrid, L0, 0.6);
	//double* Bper = getUvarovBper2new(downstreamNx, downstreamXgrid, L0, 0.4);
	double L1 = 3E19;
	//double* Bpar1 = getUvarovBpar2(downstreamNx, downstreamXgrid, L1, 0.125);
	//double* Bper1 = getUvarovBper2(downstreamNx, downstreamXgrid, L1, 0.125);
	double* Bpar1 = getUvarovBpar2new(downstreamNx, downstreamXgrid, L1, 0.5);
	double* Bper1 = getUvarovBper2new(downstreamNx, downstreamXgrid, L1, 0.5);

	for (int i = 0; i < downstreamNx; ++i) {
		Bpar[i] = Bpar[i] + Bpar1[i];
		Bper[i] = Bper[i] + Bper1[i];
	}

	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[i] = -downstreamXgrid[i];
	}

	double minField = 0.0;
	double sinField = 0.0 * minField;

	for (int i = 0; i < downstreamNx; ++i) {
		/*if (downstreamXgrid[i] < -2E19) {
			Bpar[i] = 3E-5;
		}*/
                if(sqrt(Bpar[i]*Bpar[i] + 2*Bper[i]*Bper[i]) < minField){
                    Bpar[i] = minField / sqrt(3.0);
					Bper[i] = minField / sqrt(3.0);
                }
	}

	double* downstreamB1 = new double [downstreamNx];
	double*** downstreamB = new double** [downstreamNx];
	double*** downstreamBtheta = new double** [downstreamNx];
	double*** downstreamBphi = new double** [downstreamNx];
	double*** downstreamConcentrationArray = new double** [downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamB[i] = new double* [Nz];
		downstreamBtheta[i] = new double* [Nz];
		downstreamBphi[i] = new double* [Nz];
		downstreamConcentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			downstreamB[i][j] = new double[Ny];
			downstreamBtheta[i][j] = new double[Ny];
			downstreamBphi[i][j] = new double[Ny];
			downstreamConcentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				downstreamB[i][j][k] = sqrt(Bpar[i] * Bpar[i] + 2*Bper[i] * Bper[i]);
				downstreamB1[i] = downstreamB[i][j][k];
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				downstreamBtheta[i][j][k] = atan2(sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]), Bper[i]);
				//downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphi[i][j][k] = atan2(Bpar[i], Bper[i]);
				downstreamConcentrationArray[i][j][k] = 1.0;
			}
		}
	}

	double Bfactor = Bper[downstreamNx - 1] * sqrt(2.0) / upstreamB1[0];

	double*** upstreamB = new double** [upstreamNx];
	double*** upstreamBtheta = new double** [upstreamNx];
	double*** upstreamBphi = new double** [upstreamNx];
	double*** upstreamConcentrationArray = new double** [upstreamNx];
	for (int i = 0; i < upstreamNx; ++i) {
		upstreamB[i] = new double* [Nz];
		upstreamBtheta[i] = new double* [Nz];
		upstreamBphi[i] = new double* [Nz];
		upstreamConcentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			upstreamB[i][j] = new double[Ny];
			upstreamBtheta[i][j] = new double[Ny];
			upstreamBphi[i][j] = new double[Ny];
			upstreamConcentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				upstreamB[i][j][k] = upstreamB1[i]*Bfactor;
				upstreamBtheta[i][j][k] = pi / 4;
				upstreamBphi[i][j][k] = 0;
				upstreamConcentrationArray[i][j][k] = upstreamConcentration1[i]*electronToProtonCorrection;
			}
		}
	}

	/*for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				Bpar[i] = Beff[maxIndex - i - 1];
				Bper[i] = 0;
				downstreamB[i][j][k] = Beff[maxIndex - i - 1];
				downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphi[i][j][k] = 0;
			}
		}
	}*/

	FILE* BoutputFile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < downstreamNx; ++i) {
		//fprintf(Bfile, "%g %g %g\n", -downstreamXgrid[i]/L0 + 1.5, Bpar[i], Bper[i]);
		fprintf(BoutputFile, "%g %g %g\n", downstreamXgrid[i], Bpar[i], Bper[i]);
		//fprintf(BoutputFile, "%g %g %g\n", downstreamXgrid[i], 0.0, downstreamB[i][0][0]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(BoutputFile, "%g %g %g\n", upstreamXgrid[i], 0.0, upstreamB[i][0][0]/sqrt(2.0));
	}

	fclose(BoutputFile);
	
	int Nediff = frontElectrons->getN();
	double* energyGrid = frontElectrons->getEnergyArray();
	double* frontDistribution = frontElectrons->getDistributionArray();
	for (int i = 0; i < Nediff; ++i) {
		frontDistribution[i] = frontDistribution[i] * concentration2*electronToProtonCorrection;
	}
	double norm1 = 4*pi*MassiveParticleDistributionFactory::evaluateNorm(energyGrid, frontDistribution, Nediff);
	double** diffDistributions = NULL;
	double Uph[1];
	double Eph[1];
	Uph[0] = photonEnergyDensity;
	//Uph[1] = photonIRenergyDensity;
	Eph[0] = 2.8 * kBoltzman * 2.725;
	//Eph[1] = 2.8 * kBoltzman * 140;
	/*MassiveParticleDistributionFactory::evaluateDistributionDiffusionAdvectionWithLosses(massElectron, energyGrid, frontDistribution, diffDistributions, Nediff, downstreamNx, downstreamXgrid, 0.15 * 0.26 * speed_of_light, downstreamB1, 1, Uph, Eph);

	FILE* outXfile1 = fopen("./output/x_grid1.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(outXfile1, "%g\n", downstreamXgrid[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(outXfile1, "%g\n", upstreamXgrid[i]);
	}
	fclose(outXfile1);
	FILE* outPfile1 = fopen("./output/p_grid1.dat", "w");
	for (int i = 0; i < Nediff; ++i) {
		double p = sqrt(energyGrid[i] * energyGrid[i] - me_c2 * me_c2) / (speed_of_light * speed_of_light);
		fprintf(outPfile1, "%g\n", p / massProton);
	}
	fclose(outPfile1);

	FILE* outDistributionFile1 = fopen("./output/pdf1.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		double normD = 4 * pi * MassiveParticleDistributionFactory::evaluateNorm(energyGrid, diffDistributions[i], Nediff);
		for (int j = 0; j < Nediff; ++j) {
			double p = sqrt(energyGrid[j] * energyGrid[j] - me_c2 * me_c2) / (speed_of_light * speed_of_light);
			double F = diffDistributions[i][j];
			F = ((F * p * p * p * me_c2 * me_c2 / energyGrid[j]) * massElectron / massProton)/cube(massElectron);
			fprintf(outDistributionFile1, "%g %g\n", p/massProton, F);
		}
	}
	for (int i = 0; i < upstreamNx; ++i) {
		double normU = 4 * pi * MassiveParticleDistributionFactory::evaluateNorm(energy, upstreamDistributions1[i], Nenergy);
		for (int j = 0; j < Nediff; ++j) {
			double p = sqrt(energyGrid[j] * energyGrid[j] - me_c2 * me_c2) / (speed_of_light * speed_of_light);
			double F = MassiveParticleDistributionFactory::getDistribution(energyGrid[j], energy, upstreamDistributions1[i], Nenergy)*upstreamConcentrationArray[i][0][0];
			F = ((F * p * p * p * me_c2 * me_c2 / energyGrid[j]) * massElectron / massProton)/cube(massElectron);
			fprintf(outDistributionFile1, "%g %g\n", p/massProton, F);
		}
	}
	fclose(outDistributionFile1);

	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				double norm = 4*pi*MassiveParticleDistributionFactory::evaluateNorm(energyGrid, diffDistributions[i], Nediff);
				//downstreamConcentrationArray[i][j][k] = concentration2*electronToProtonCorrection;
				downstreamConcentrationArray[i][j][k] = norm;
			}
		}
	}*/

	//frontElectrons->rescaleDistribution(0.4);

	/*MassiveParticleDistribution**** downstreamElectrons = new MassiveParticleDistribution * **[downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamElectrons[i] = new MassiveParticleDistribution * *[Ny];
		for (int j = 0; j < Ny; ++j) {
			downstreamElectrons[i][j] = new MassiveParticleDistribution * [Nz];
			for (int k = 0; k < Nz; ++k) {
				downstreamElectrons[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energyGrid, diffDistributions[i], Nediff, DistributionInputType::ENERGY_FE);
			}
		}
	}*/

	MassiveParticleDistribution**** upstreamElectrons = new MassiveParticleDistribution * **[upstreamNx];
	for (int i = 0; i < upstreamNx; ++i) {
		upstreamElectrons[i] = new MassiveParticleDistribution * *[Nz];
		for (int j = 0; j < Nz; ++j) {
			upstreamElectrons[i][j] = new MassiveParticleDistribution * [Ny];
			for (int k = 0; k < Ny; ++k) {
				MassiveParticleTabulatedIsotropicDistribution* localDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energy, upstreamDistributions1[i], Nenergy, DistributionInputType::ENERGY_FE);
				//electrons2[i][j][k] = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
				//localDistribution->rescaleDistribution(0.4);
				upstreamElectrons[i][j][k] = localDistribution;
			}
		}
	}

	MassiveParticleIsotropicDistribution* distribution0 = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
	MassiveParticleTabulatedIsotropicDistribution* electrons2 = new MassiveParticleTabulatedIsotropicDistribution(distribution0, me_c2, 1E16 * 1.6E-12, 1000);
	//upstreamElectrons = frontElectrons;
	//double E0 = 1.6E-1;
	//MassiveParticleIsotropicDistribution* upstreamElectrons = new MassiveParticleMonoenergeticDistribution(massElectron, E0, 0.01 * E0);

	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				downstreamConcentrationArray[i][j][k] = concentration2 * electronToProtonCorrection;
			}
		}
	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.15 * 0.26 * speed_of_light, photonTotalEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.15 * 0.26 * speed_of_light, 0.15 * 0.26 * speed_of_light, photonEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(downstreamNx, downstreamXgrid, Ny, Nz, downstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	RadiationSource* upstreamSource = new RectangularSourceInhomogenousDistribution(upstreamNx, upstreamXgrid, Ny, Nz, upstreamElectrons, upstreamB, upstreamBtheta, upstreamBphi, upstreamConcentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* downstreamSource = new RectangularSource(1, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, downstreamXgrid[0], downstreamXgrid[Nx - 1], 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(downstreamNx - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(downstreamNx - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(0, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleTabulatedIsotropicDistribution* fardownstreamDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(0, 0, 0));
	int Nfardownstream = fardownstreamDistribution->getN();
	double* Efardownstream = fardownstreamDistribution->getEnergyArray();
	double* Ffardownstream = fardownstreamDistribution->getDistributionArray();
	FILE* outFarFile = fopen("fardownstreamelectrons.dat", "w");
	for (int i = 0; i < Nfardownstream; ++i) {
		fprintf(outFarFile, "%g %g\n", Efardownstream[i], Ffardownstream[i]);
	}
	fclose(outFarFile);

	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));

	FILE* outLeftFile = fopen("./output/electronsDownstream.dat", "w");
		double p = pmin;
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = fardownstreamDistribution->distributionNormalized(E) * downstreamConcentrationArray[0][0][0];
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outLeftFile, "%g %g\n", p, F);
			p = p * factorp;
		}
	fclose(outLeftFile);

	MassiveParticleTabulatedIsotropicDistribution* farupstreamDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(upstreamSource->getParticleDistribution(upstreamNx - 2, 0, 0));
	int Nfarupstream = farupstreamDistribution->getN();
	double* Efarupstream = farupstreamDistribution->getEnergyArray();
	double* Ffarupstream = farupstreamDistribution->getDistributionArray();
	FILE* outFarUpFile = fopen("farupstreamelectrons.dat", "w");
	for (int i = 0; i < Nfarupstream; ++i) {
		fprintf(outFarUpFile, "%g %g\n", Efarupstream[i], Ffarupstream[i] *upstreamConcentrationArray[upstreamNx - 2][0][0]);
	}
	fclose(outFarUpFile);


	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(outXfile, "%g\n", downstreamXgrid[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(outXfile, "%g\n", upstreamXgrid[i]);
	}
	fclose(outXfile);

	FILE* outPfile = fopen("./output/p_grid.dat", "w");
	p = pmin;
	for (int i = 0; i < Np; ++i) {
		fprintf(outPfile, "%g\n", p * massElectron / massProton);
		p = p * factorp;
	}
	fclose(outPfile);

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(i, 0, 0));
		//MassiveParticleIsotropicDistribution* distribution = frontElectrons;
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E) * downstreamConcentrationArray[i][0][0];
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g %g\n", p, F);
			p = p * factorp;
		}
	}
	for (int i = 0; i < upstreamNx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(upstreamSource->getParticleDistribution(i, 0, 0));
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E)*upstreamConcentrationArray[i][0][0];
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g %g\n", p, F);
			p = p * factorp;
		}
	}
	fclose(outDistributionFile);

	FILE* concentrationFile = fopen("./output/concentration.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(concentrationFile, "%g %g\n", downstreamXgrid[i], downstreamConcentrationArray[i][0][0]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(concentrationFile, "%g %g\n", upstreamXgrid[i], upstreamConcentrationArray[i][0][0]);
	}
	fclose(concentrationFile);

	/*FILE* outDiffusionConvectionFile = fopen("./output/diffusionConvection.dat", "w");
	double tempE[3] = { 48, 160, 480 };
	for (int i = 1; i < downstreamNx - 1; ++i) {
		MassiveParticleIsotropicDistribution* leftDistribution = new MassiveParticleTabulatedIsotropicDistribution(*dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i - 1, 0, 0)));
		MassiveParticleIsotropicDistribution* middleDistribution = new MassiveParticleTabulatedIsotropicDistribution(*dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i, 0, 0)));
		MassiveParticleIsotropicDistribution* rightDistribution = new MassiveParticleTabulatedIsotropicDistribution(*dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i + 1, 0, 0)));
		fprintf(outDiffusionConvectionFile, "%g", downstreamXgrid[i]);
		for (int j = 0; j < 3; ++j) {
			double D = tempE[j] * speed_of_light / (3 * electron_charge * downstreamB[i][0][0]);
			double Dd2f = fabs(2*D * ((rightDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i + 1][0][0] - middleDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i][0][0]) / (downstreamXgrid[i + 1] - downstreamXgrid[i]) -
				(middleDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i][0][0] - leftDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i - 1][0][0]) / (downstreamXgrid[i] - downstreamXgrid[i - 1])) / (downstreamXgrid[i + 1] - downstreamXgrid[i - 1]));
			double udf = fabs(0.15 * 0.2 * speed_of_light * (rightDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i + 1][0][0] - middleDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i][0][0]) / (downstreamXgrid[i + 1] - downstreamXgrid[i]));
			fprintf(outDiffusionConvectionFile, " %g %g", Dd2f, udf);
		}
		delete[] leftDistribution;
		delete[] middleDistribution;
		delete[] rightDistribution;
		fprintf(outDiffusionConvectionFile, "\n");
	}
	fclose(outDiffusionConvectionFile);*/


	int Ne = 10;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", downstreamSource, 1.6E-12, 1.6E4, 1000);


	double Ephmin = 1.6E-12;
	double Ephmax = 1.6E4;
	int Nph = 1000;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	FILE* outFile = fopen("./output/W50synchandcompt.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		//omp_set_lock(&my_lock);
		printf("writeEFEFromSourceToFile iph = %d\n", i);
		printLog("writeEFEFromSourceToFile iph = %d\n", i);
		//omp_unset_lock(&my_lock);
		//double flux1 = sumEvaluator->evaluateFluxFromSource(currentE, downstreamSource);
		double flux1 = 0;
		double flux2 = sumEvaluator->evaluateFluxFromSource(currentE, upstreamSource);
		double fluxHead = 0;
		double fluxCone = 0;
		int j;
#pragma omp parallel for private(j) shared(currentE, downstreamSource, downstreamXgrid, downstreamNx) reduction(+:flux1, fluxHead, fluxCone)
		for (j = 0; j < downstreamNx; ++j) {
			double flux = sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, j, 0);
			flux1 += flux;
			if ((downstreamXgrid[j] >= headMaxX) && (downstreamXgrid[j] <= headMinX)) {
				fluxHead += flux;
			}
			if ((downstreamXgrid[j] >= coneMaxX) && (downstreamXgrid[j] <= coneMinX)){
				fluxCone += flux;
			}
		}
		fprintf(outFile, "%g %g %g %g\n", currentE / 1.6E-12, currentE * (flux1 + flux2), currentE*fluxHead, currentE*fluxCone);
		currentE = currentE * factor;
	}
	fclose(outFile);

	double* profileXMM = new double[downstreamNx + upstreamNx];
	double* profileNuSTAR = new double[downstreamNx + upstreamNx];
	int irho;

	Ephmin = 0.3 * 1000 * 1.6E-12;
	Ephmax = 10 * 1000 * 1.6E-12;
	Nph = 20;
	omp_lock_t lock;
	omp_init_lock(&lock);
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, downstreamSource, downstreamNx, sumEvaluator, Nph, profileXMM, lock)
	for (irho = 0; irho < downstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = downstreamSource->getCrossSectionArea(irho, 0);
		double d = downstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0/currentE)*sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileXMM[irho] = localFlux;
	}
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, upstreamSource, downstreamNx, upstreamNx, sumEvaluator, Nph, profileXMM, lock)
	for (irho = 0; irho < upstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = upstreamSource->getCrossSectionArea(irho, 0);
		double d = upstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, upstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileXMM[downstreamNx + irho] = localFlux;
	}
	FILE* xmmFile = fopen("./output/xmmprofile.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(xmmFile, "%g %g\n", downstreamXgrid[i], profileXMM[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(xmmFile, "%g %g\n", upstreamXgrid[i], profileXMM[downstreamNx + i]);
	}
	fclose(xmmFile);

	Ephmin = 10 * 1000 * 1.6E-12;
	Ephmax = 20 * 1000 * 1.6E-12;
	Nph = 20;
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, downstreamSource, downstreamNx, sumEvaluator, Nph, profileNuSTAR, lock)
	for (irho = 0; irho < downstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating image irho = %d\n", irho);
		printLog("evaluating image irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = downstreamSource->getCrossSectionArea(irho, 0);
		double d = downstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileNuSTAR[irho] = localFlux;
	}
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, upstreamSource, downstreamNx, upstreamNx, sumEvaluator, Nph, profileNuSTAR, lock)
	for (irho = 0; irho < upstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating image irho = %d\n", irho);
		printLog("evaluating image irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = upstreamSource->getCrossSectionArea(irho, 0);
		double d = upstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, upstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileNuSTAR[downstreamNx + irho] = localFlux;
	}

	omp_destroy_lock(&lock);

	FILE* nustarFile = fopen("./output/nustarprofile.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(nustarFile, "%g %g\n", downstreamXgrid[i], profileNuSTAR[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(nustarFile, "%g %g\n", upstreamXgrid[i], profileNuSTAR[downstreamNx + i]);
	}
	fclose(nustarFile);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", downstreamSource, 1.6E-1, 1.6E3, 300);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", downstreamSource, 1.6E-9, 50*1.6E-9, 300);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	return;

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	sumEvaluator->writeImageFromSourceToFile("./output/W50scImageeV.dat", downstreamSource, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageKeV.dat", downstreamSource, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageMeV.dat", downstreamSource, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageGeV.dat", downstreamSource, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("./output/W50scImageTeV.dat", downstreamSource, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	//sumEvaluator->writeImageFromSourceToFile("../output/W50scImagePeV.dat", downstreamSource, 1.6E3, 1.6E4, 20);

	printf("start writing x-E diagram\n");
	printLog("start writing x-E diagram\n");

	double Emin = 1.6E-12;
	double Emax = 1.6E4;

	int Nnu = 200;

	factor = pow(Emax / Emin, 1.0 / (Nnu - 1));
	currentE = Emin;

	double** F = new double* [Nnu];
	for (int i = 0; i < Nnu; ++i) {
		F[i] = new double[downstreamNx];
	}

	FILE* Efile = fopen("./output/E_grid.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(Efile, "%g\n", currentE);
		currentE = currentE * factor;
	}
	fclose(Efile);

	currentE = Emin;
	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, downstreamSource, currentE, factor, sumEvaluator, i)
		for (j = 0; j < downstreamNx; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, j, 0) / downstreamSource->getCrossSectionArea(j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outxEfile = fopen("./output/xE.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nnu; ++j) {
			if (j == 0) {
				fprintf(outxEfile, "%g", F[j][i]);
			}
			else {
				fprintf(outxEfile, " %g", F[j][i]);
			}
		}
		fprintf(outxEfile, "\n");
	}
	fclose(outxEfile);

}

void evaluateW50comptonAndSynchrotronAdvectionfunctionWithBrinkmann() {
	double distance = (18000 / 3.26) * parsec;

	double* energy;
	double* downstreamXgrid;
	double* upstreamXgrid;
	double* concentration;
	double** distributions;

	int Nenergy;
	int Nx;

	double* energyBrinkmann;
	double* downstreamXgridBrinkmann;
	double* upstreamXgridBrinkmann;
	double* concentrationBrinkmann;
	double** distributionsBrinkmann;

	int NenergyBrinkmann;
	int NxBrinkmann;

	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	const char* xfileName = "./examples_data/W50/newdistribution/x_grid.dat";
	const char* BfileName = "./examples_data/W50/B15FEB6E18/Beff.dat";

	const char* distributionFileName = "./examples_data/W50/newdistribution/electrons_full.dat";
	const char* pfileName = "./examples_data/W50/newdistribution/p_grid.dat";

	const char* fileName = "./examples_data/W50/newdistribution/electrons.dat";
	const char* protonsFileName = "./examples_data/W50/newdistribution/protons.dat";

	const char* xfileNameBrinkmann = "./examples_data/W50/Brinkmann2/x_grid.dat";
	//const char* BfileNameBrinkmann = "./examples_data/W50/B15FEB6E18/Beff.dat";

	const char* distributionFileNameBrinkmann = "./examples_data/W50/Brinkmann/electrons_full.dat";
	const char* pfileNameBrinkmann = "./examples_data/W50/Brinkmann2/p_grid.dat";

	const char* fileNameBrinkmann = "./examples_data/W50/Brinkmann2/electrons.dat";
	const char* protonsFileNameBrinkmann = "./examples_data/W50/Brinkmann2/protons.dat";

	/*Nx = 0;
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
	fclose(xfile);*/

	double* xgrid1;

	double* concentration1;
	double** distributions1;

	//double electronToProtonCorrection = 3E-7;

	MassiveParticleTabulatedIsotropicDistribution* frontElectrons;
	double concentration2;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, frontElectrons, concentration2);
	//frontElectrons = new MassiveParticleTabulatedIsotropicDistribution(new MassiveParticlePowerLawDistribution(massElectron, 2.0, me_c2), me_c2, 1600, 1000);
	MassiveParticleTabulatedIsotropicDistribution* frontProtons;
	double concentration3;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileName, frontProtons, concentration3);

	double electronToProtonCorrection = concentration3 * frontProtons->getDistributionArray()[70] / (concentration2 * frontElectrons->getDistributionArray()[70]);

	MassiveParticleTabulatedIsotropicDistribution* frontElectronsBrinkmann;
	double concentration2Brinkmann;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileNameBrinkmann, frontElectronsBrinkmann, concentration2Brinkmann);
	//frontElectrons = new MassiveParticleTabulatedIsotropicDistribution(new MassiveParticlePowerLawDistribution(massElectron, 2.0, me_c2), me_c2, 1600, 1000);
	MassiveParticleTabulatedIsotropicDistribution* frontProtonsBrinkmann;
	double concentration3Brinkmann;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileNameBrinkmann, frontProtonsBrinkmann, concentration3Brinkmann);

	double electronToProtonCorrectionBrinkmann = concentration3Brinkmann * frontProtonsBrinkmann->getDistributionArray()[70] / (concentration2Brinkmann * frontElectronsBrinkmann->getDistributionArray()[70]);

	MassiveParticleDistributionFactory::readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(massElectron, xfileName, pfileName, distributionFileName, xgrid1, energy, distributions1, concentration1, Nenergy, Nx);

	double* xgrid1Brinkmann;

	double* concentration1Brinkmann;
	double** distributions1Brinkmann;

	MassiveParticleDistributionFactory::readInhomogenousTabulatedIsotropicDistributionFromMonteCarlo(massElectron, xfileNameBrinkmann, pfileNameBrinkmann, distributionFileNameBrinkmann, xgrid1Brinkmann, energyBrinkmann, distributions1Brinkmann, concentration1Brinkmann, NenergyBrinkmann, NxBrinkmann);

	double* Beff = new double[Nx];
	FILE* Bfile = fopen(BfileName, "r");
	for (int i = 0; i < Nx; ++i) {
		fscanf(Bfile, "%lf", &Beff[i]);
	}
	fclose(Bfile);

	int zeroIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}
	double downstreamSize = 1E20;
	double upstreamSize = 1.6E20;
	int maxIndex = Nx - 1;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= downstreamSize) {
			maxIndex = i;
			break;
		}
	}
	//maxIndex = Nx - 1;
	int minIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= -upstreamSize) {
			minIndex = i;
			break;
		}
	}
	//minIndex = 0;

	int downstreamNx = maxIndex + 1 - zeroIndex;
	int upstreamNx = zeroIndex - minIndex;

	downstreamXgrid = new double[downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[downstreamNx - i - 1] = xgrid1[i + zeroIndex];
	}

	upstreamXgrid = new double[upstreamNx];
	double* upstreamB1 = new double[upstreamNx];
	double* upstreamConcentration1 = new double[upstreamNx];
	double** upstreamDistributions1 = new double* [upstreamNx];
	for (int i = 0; i < upstreamNx; ++i) {
		upstreamXgrid[upstreamNx - 1 - i] = -xgrid1[i + minIndex];
		upstreamB1[upstreamNx - 1 - i] = Beff[i + minIndex];
		upstreamConcentration1[upstreamNx - 1 - i] = concentration1[i + minIndex];
		upstreamDistributions1[upstreamNx - 1 - i] = distributions1[i + minIndex];
	}

	int zeroIndexBrinkmann = 0;
	for (int i = 0; i < NxBrinkmann; ++i) {
		if (xgrid1Brinkmann[i] >= 0) {
			zeroIndexBrinkmann = i;
			break;
		}
	}
	double downstreamSizeBrinkmann = 2E20;
	int maxIndexBrinkmann = NxBrinkmann - 1;
	for (int i = 0; i < NxBrinkmann; ++i) {
		if (xgrid1Brinkmann[i] >= downstreamSizeBrinkmann) {
			maxIndexBrinkmann = i;
			break;
		}
	}
	//maxIndex = Nx - 1;

	int downstreamNxBrinkmann = maxIndexBrinkmann + 1 - zeroIndexBrinkmann;

	downstreamXgridBrinkmann = new double[downstreamNxBrinkmann];
	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		downstreamXgridBrinkmann[downstreamNxBrinkmann - i - 1] = xgrid1Brinkmann[i + zeroIndexBrinkmann];
	}

	double size = 0.5 * fabs(headMaxX);
	double B0 = 6E-5;
	double sizeBrinkmann = 1.5*size;

	//RadiationSourceInCylindrical* downstreamSource = new SimpleFlatSource(upstreamElectrons, downstreamB, pi / 2, 0, concentration, size, size, distance);
	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photonsIR = new PhotonPlankDistribution(30, 0.8 / 10000);
	double photonIRconcentration = photonsIR->getConcentration();
	double photonIRenergyDensity = photonIRconcentration * photonsIR->getMeanEnergy();
	double photonConcentration = photons->getConcentration();
	double photonEnergyDensity = photonConcentration * photons->getMeanEnergy();
	PhotonMultiPlankDistribution* photonsTotal = new PhotonMultiPlankDistribution(2.725, 1, 30, 0.8 / 10000);
	//PhotonMultiPlankDistribution* photonsTotal = PhotonMultiPlankDistribution::getGalacticField();

	double photonTotalConcentration = photonsTotal->getConcentration();
	double photonTotalEnergyDensity = photonTotalConcentration * photonsTotal->getMeanEnergy();

	int Nz = 1;
	int Ny = 1;

	double L0 = 0.3E18;
	double* Bpar = getUvarovBpar2(downstreamNx, downstreamXgrid, L0, 0.125);
	double* Bper = getUvarovBper2(downstreamNx, downstreamXgrid, L0, 0.125);
	//double* Bpar = getUvarovBpar2new(downstreamNx, downstreamXgrid, L0, 0.6);
	//double* Bper = getUvarovBper2new(downstreamNx, downstreamXgrid, L0, 0.4);
	double L1 = 3E19;
	//double* Bpar1 = getUvarovBpar2(downstreamNx, downstreamXgrid, L1, 0.125);
	//double* Bper1 = getUvarovBper2(downstreamNx, downstreamXgrid, L1, 0.125);
	double* Bpar1 = getUvarovBpar2new(downstreamNx, downstreamXgrid, L1, 0.5);
	double* Bper1 = getUvarovBper2new(downstreamNx, downstreamXgrid, L1, 0.5);

	for (int i = 0; i < downstreamNx; ++i) {
		Bpar[i] = Bpar[i] + Bpar1[i];
		Bper[i] = Bper[i] + Bper1[i];
	}

	double L0Brinkmann = 1.0E18;
	double* BparBrinkmann = getUvarovBpar2(downstreamNxBrinkmann, downstreamXgridBrinkmann, L0Brinkmann, 9.3/60.0);
	double* BperBrinkmann = getUvarovBper2(downstreamNxBrinkmann, downstreamXgridBrinkmann, L0Brinkmann, 9.3/60.0);

	/*for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		BparBrinkmann[i] = 30E-6;
		BperBrinkmann[i] = 0.0;
	}*/

	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[i] = -downstreamXgrid[i];
	}

	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		downstreamXgridBrinkmann[i] = -downstreamXgridBrinkmann[i];
	}

	double minField = 0.0;
	double sinField = 0.0 * minField;

	for (int i = 0; i < downstreamNx; ++i) {
		/*if (downstreamXgrid[i] < -2E19) {
			Bpar[i] = 3E-5;
		}*/
		if (sqrt(Bpar[i] * Bpar[i] + 2 * Bper[i] * Bper[i]) < minField) {
			Bpar[i] = minField / sqrt(3.0);
			Bper[i] = minField / sqrt(3.0);
		}
	}

	double minFieldBrinkmann = 0.0;
	double sinFieldBrinkmann = 0.0 * minFieldBrinkmann;

	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		/*if (downstreamXgrid[i] < -2E19) {
			Bpar[i] = 3E-5;
		}*/
		if (sqrt(BparBrinkmann[i] * BparBrinkmann[i] + 2 * BperBrinkmann[i] * BperBrinkmann[i]) < minFieldBrinkmann) {
			BparBrinkmann[i] = minFieldBrinkmann / sqrt(3.0);
			BperBrinkmann[i] = minFieldBrinkmann / sqrt(3.0);
		}
	}

	double* downstreamB1 = new double[downstreamNx];
	double*** downstreamB = new double** [downstreamNx];
	double*** downstreamBtheta = new double** [downstreamNx];
	double*** downstreamBphi = new double** [downstreamNx];
	double*** downstreamConcentrationArray = new double** [downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamB[i] = new double* [Nz];
		downstreamBtheta[i] = new double* [Nz];
		downstreamBphi[i] = new double* [Nz];
		downstreamConcentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			downstreamB[i][j] = new double[Ny];
			downstreamBtheta[i][j] = new double[Ny];
			downstreamBphi[i][j] = new double[Ny];
			downstreamConcentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				downstreamB[i][j][k] = sqrt(Bpar[i] * Bpar[i] + 2 * Bper[i] * Bper[i]);
				downstreamB1[i] = downstreamB[i][j][k];
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				downstreamBtheta[i][j][k] = atan2(sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]), Bper[i]);
				//downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphi[i][j][k] = atan2(Bpar[i], Bper[i]);
				downstreamConcentrationArray[i][j][k] = 1.0;
			}
		}
	}

	double* downstreamB1Brinkmann = new double[downstreamNxBrinkmann];
	double*** downstreamBBrinkmann = new double** [downstreamNxBrinkmann];
	double*** downstreamBthetaBrinkmann = new double** [downstreamNxBrinkmann];
	double*** downstreamBphiBrinkmann = new double** [downstreamNxBrinkmann];
	double*** downstreamConcentrationArrayBrinkmann = new double** [downstreamNxBrinkmann];
	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		downstreamBBrinkmann[i] = new double* [Nz];
		downstreamBthetaBrinkmann[i] = new double* [Nz];
		downstreamBphiBrinkmann[i] = new double* [Nz];
		downstreamConcentrationArrayBrinkmann[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			downstreamBBrinkmann[i][j] = new double[Ny];
			downstreamBthetaBrinkmann[i][j] = new double[Ny];
			downstreamBphiBrinkmann[i][j] = new double[Ny];
			downstreamConcentrationArrayBrinkmann[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				downstreamBBrinkmann[i][j][k] = sqrt(BparBrinkmann[i] * BparBrinkmann[i] + 2 * BperBrinkmann[i] * BperBrinkmann[i]);
				downstreamB1Brinkmann[i] = downstreamBBrinkmann[i][j][k];
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				downstreamBthetaBrinkmann[i][j][k] = atan2(sqrt(BparBrinkmann[i] * BparBrinkmann[i] + BperBrinkmann[i] * BperBrinkmann[i]), BperBrinkmann[i]);
				//downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphiBrinkmann[i][j][k] = atan2(BparBrinkmann[i], BperBrinkmann[i]);
				downstreamConcentrationArrayBrinkmann[i][j][k] = 1.0;
			}
		}
	}

	double Bfactor = Bper[downstreamNx - 1] * sqrt(2.0) / upstreamB1[0];

	double*** upstreamB = new double** [upstreamNx];
	double*** upstreamBtheta = new double** [upstreamNx];
	double*** upstreamBphi = new double** [upstreamNx];
	double*** upstreamConcentrationArray = new double** [upstreamNx];
	for (int i = 0; i < upstreamNx; ++i) {
		upstreamB[i] = new double* [Nz];
		upstreamBtheta[i] = new double* [Nz];
		upstreamBphi[i] = new double* [Nz];
		upstreamConcentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			upstreamB[i][j] = new double[Ny];
			upstreamBtheta[i][j] = new double[Ny];
			upstreamBphi[i][j] = new double[Ny];
			upstreamConcentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				upstreamB[i][j][k] = upstreamB1[i] * Bfactor;
				upstreamBtheta[i][j][k] = pi / 4;
				upstreamBphi[i][j][k] = 0;
				upstreamConcentrationArray[i][j][k] = upstreamConcentration1[i] * electronToProtonCorrection;
			}
		}
	}

	/*for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				Bpar[i] = Beff[maxIndex - i - 1];
				Bper[i] = 0;
				downstreamB[i][j][k] = Beff[maxIndex - i - 1];
				downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphi[i][j][k] = 0;
			}
		}
	}*/

	FILE* BoutputFile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < downstreamNx; ++i) {
		//fprintf(Bfile, "%g %g %g\n", -downstreamXgrid[i]/L0 + 1.5, Bpar[i], Bper[i]);
		fprintf(BoutputFile, "%g %g %g\n", downstreamXgrid[i], Bpar[i], Bper[i]);
		//fprintf(BoutputFile, "%g %g %g\n", downstreamXgrid[i], 0.0, downstreamB[i][0][0]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(BoutputFile, "%g %g %g\n", upstreamXgrid[i], 0.0, upstreamB[i][0][0] / sqrt(2.0));
	}

	fclose(BoutputFile);

	int Nediff = frontElectrons->getN();
	double* energyGrid = frontElectrons->getEnergyArray();
	double* frontDistribution = frontElectrons->getDistributionArray();
	for (int i = 0; i < Nediff; ++i) {
		frontDistribution[i] = frontDistribution[i] * concentration2 * electronToProtonCorrection;
	}
	double norm1 = 4 * pi * MassiveParticleDistributionFactory::evaluateNorm(energyGrid, frontDistribution, Nediff);
	double** diffDistributions = NULL;
	double Uph[1];
	double Eph[1];
	Uph[0] = photonEnergyDensity;
	//Uph[1] = photonIRenergyDensity;
	Eph[0] = 2.8 * kBoltzman * 2.725;
	//Eph[1] = 2.8 * kBoltzman * 140;
	/*MassiveParticleDistributionFactory::evaluateDistributionDiffusionAdvectionWithLosses(massElectron, energyGrid, frontDistribution, diffDistributions, Nediff, downstreamNx, downstreamXgrid, 0.15 * 0.26 * speed_of_light, downstreamB1, 1, Uph, Eph);

	FILE* outXfile1 = fopen("./output/x_grid1.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(outXfile1, "%g\n", downstreamXgrid[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(outXfile1, "%g\n", upstreamXgrid[i]);
	}
	fclose(outXfile1);
	FILE* outPfile1 = fopen("./output/p_grid1.dat", "w");
	for (int i = 0; i < Nediff; ++i) {
		double p = sqrt(energyGrid[i] * energyGrid[i] - me_c2 * me_c2) / (speed_of_light * speed_of_light);
		fprintf(outPfile1, "%g\n", p / massProton);
	}
	fclose(outPfile1);

	FILE* outDistributionFile1 = fopen("./output/pdf1.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		double normD = 4 * pi * MassiveParticleDistributionFactory::evaluateNorm(energyGrid, diffDistributions[i], Nediff);
		for (int j = 0; j < Nediff; ++j) {
			double p = sqrt(energyGrid[j] * energyGrid[j] - me_c2 * me_c2) / (speed_of_light * speed_of_light);
			double F = diffDistributions[i][j];
			F = ((F * p * p * p * me_c2 * me_c2 / energyGrid[j]) * massElectron / massProton)/cube(massElectron);
			fprintf(outDistributionFile1, "%g %g\n", p/massProton, F);
		}
	}
	for (int i = 0; i < upstreamNx; ++i) {
		double normU = 4 * pi * MassiveParticleDistributionFactory::evaluateNorm(energy, upstreamDistributions1[i], Nenergy);
		for (int j = 0; j < Nediff; ++j) {
			double p = sqrt(energyGrid[j] * energyGrid[j] - me_c2 * me_c2) / (speed_of_light * speed_of_light);
			double F = MassiveParticleDistributionFactory::getDistribution(energyGrid[j], energy, upstreamDistributions1[i], Nenergy)*upstreamConcentrationArray[i][0][0];
			F = ((F * p * p * p * me_c2 * me_c2 / energyGrid[j]) * massElectron / massProton)/cube(massElectron);
			fprintf(outDistributionFile1, "%g %g\n", p/massProton, F);
		}
	}
	fclose(outDistributionFile1);

	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				double norm = 4*pi*MassiveParticleDistributionFactory::evaluateNorm(energyGrid, diffDistributions[i], Nediff);
				//downstreamConcentrationArray[i][j][k] = concentration2*electronToProtonCorrection;
				downstreamConcentrationArray[i][j][k] = norm;
			}
		}
	}*/

	//frontElectrons->rescaleDistribution(0.4);

	/*MassiveParticleDistribution**** downstreamElectrons = new MassiveParticleDistribution * **[downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamElectrons[i] = new MassiveParticleDistribution * *[Ny];
		for (int j = 0; j < Ny; ++j) {
			downstreamElectrons[i][j] = new MassiveParticleDistribution * [Nz];
			for (int k = 0; k < Nz; ++k) {
				downstreamElectrons[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energyGrid, diffDistributions[i], Nediff, DistributionInputType::ENERGY_FE);
			}
		}
	}*/

	MassiveParticleDistribution**** upstreamElectrons = new MassiveParticleDistribution * **[upstreamNx];
	for (int i = 0; i < upstreamNx; ++i) {
		upstreamElectrons[i] = new MassiveParticleDistribution * *[Nz];
		for (int j = 0; j < Nz; ++j) {
			upstreamElectrons[i][j] = new MassiveParticleDistribution * [Ny];
			for (int k = 0; k < Ny; ++k) {
				MassiveParticleTabulatedIsotropicDistribution* localDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energy, upstreamDistributions1[i], Nenergy, DistributionInputType::ENERGY_FE);
				//electrons2[i][j][k] = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
				//localDistribution->rescaleDistribution(0.4);
				upstreamElectrons[i][j][k] = localDistribution;
			}
		}
	}

	MassiveParticleIsotropicDistribution* distribution0 = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
	MassiveParticleTabulatedIsotropicDistribution* electrons2 = new MassiveParticleTabulatedIsotropicDistribution(distribution0, me_c2, 1E16 * 1.6E-12, 1000);
	//upstreamElectrons = frontElectrons;
	//double E0 = 1.6E-1;
	//MassiveParticleIsotropicDistribution* upstreamElectrons = new MassiveParticleMonoenergeticDistribution(massElectron, E0, 0.01 * E0);

	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				downstreamConcentrationArray[i][j][k] = concentration2 * electronToProtonCorrection;
			}
		}
	}

	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				downstreamConcentrationArrayBrinkmann[i][j][k] = concentration2Brinkmann * electronToProtonCorrectionBrinkmann;
			}
		}
	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.15 * 0.26 * speed_of_light, photonTotalEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.15 * 0.26 * speed_of_light, 0.15 * 0.26 * speed_of_light, photonEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSourceBrinkmann = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNxBrinkmann, downstreamXgridBrinkmann, Ny, Nz, frontElectronsBrinkmann, downstreamBBrinkmann, downstreamBthetaBrinkmann, downstreamBphiBrinkmann, downstreamConcentrationArrayBrinkmann, 0, sizeBrinkmann, 0, pi * sizeBrinkmann, distance, 0.093E10, 0.093E10, photonEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(downstreamNx, downstreamXgrid, Ny, Nz, downstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	RadiationSource* upstreamSource = new RectangularSourceInhomogenousDistribution(upstreamNx, upstreamXgrid, Ny, Nz, upstreamElectrons, upstreamB, upstreamBtheta, upstreamBphi, upstreamConcentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* downstreamSource = new RectangularSource(1, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, downstreamXgrid[0], downstreamXgrid[Nx - 1], 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(downstreamNx - 1, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(downstreamNx - 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(0, 0, 0));
	distributionLeft->writeDistribution("./output/distributionLeft.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleTabulatedIsotropicDistribution* fardownstreamDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(0, 0, 0));
	int Nfardownstream = fardownstreamDistribution->getN();
	double* Efardownstream = fardownstreamDistribution->getEnergyArray();
	double* Ffardownstream = fardownstreamDistribution->getDistributionArray();
	FILE* outFarFile = fopen("fardownstreamelectrons.dat", "w");
	for (int i = 0; i < Nfardownstream; ++i) {
		fprintf(outFarFile, "%g %g\n", Efardownstream[i], Ffardownstream[i]);
	}
	fclose(outFarFile);

	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));

	FILE* outLeftFile = fopen("./output/electronsDownstream.dat", "w");
	double p = pmin;
	for (int j = 0; j < Np; ++j) {
		double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
		double F = fardownstreamDistribution->distributionNormalized(E) * downstreamConcentrationArray[0][0][0];
		F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
		fprintf(outLeftFile, "%g %g\n", p, F);
		p = p * factorp;
	}
	fclose(outLeftFile);

	MassiveParticleTabulatedIsotropicDistribution* farupstreamDistribution = dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(upstreamSource->getParticleDistribution(upstreamNx - 2, 0, 0));
	int Nfarupstream = farupstreamDistribution->getN();
	double* Efarupstream = farupstreamDistribution->getEnergyArray();
	double* Ffarupstream = farupstreamDistribution->getDistributionArray();
	FILE* outFarUpFile = fopen("farupstreamelectrons.dat", "w");
	for (int i = 0; i < Nfarupstream; ++i) {
		fprintf(outFarUpFile, "%g %g\n", Efarupstream[i], Ffarupstream[i] * upstreamConcentrationArray[upstreamNx - 2][0][0]);
	}
	fclose(outFarUpFile);


	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(outXfile, "%g\n", downstreamXgrid[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(outXfile, "%g\n", upstreamXgrid[i]);
	}
	fclose(outXfile);

	FILE* outXfileBrinkmann = fopen("./output/x_grid_Brinkmann.dat", "w");
	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		fprintf(outXfileBrinkmann, "%g\n", downstreamXgridBrinkmann[i]);
	}
	fclose(outXfileBrinkmann);

	FILE* outPfile = fopen("./output/p_grid.dat", "w");
	p = pmin;
	for (int i = 0; i < Np; ++i) {
		fprintf(outPfile, "%g\n", p * massElectron / massProton);
		p = p * factorp;
	}
	fclose(outPfile);

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(downstreamSource->getParticleDistribution(i, 0, 0));
		//MassiveParticleIsotropicDistribution* distribution = frontElectrons;
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E) * downstreamConcentrationArray[i][0][0];
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g %g\n", p, F);
			p = p * factorp;
		}
	}
	for (int i = 0; i < upstreamNx; ++i) {
		p = pmin;
		MassiveParticleIsotropicDistribution* distribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(upstreamSource->getParticleDistribution(i, 0, 0));
		for (int j = 0; j < Np; ++j) {
			double E = sqrt(p * p * me_c2 * me_c2 + me_c2 * me_c2);
			double F = distribution->distributionNormalized(E) * upstreamConcentrationArray[i][0][0];
			F = (F * p * p * p * me_c2 * me_c2 / E) * massElectron / massProton;
			fprintf(outDistributionFile, "%g %g\n", p, F);
			p = p * factorp;
		}
	}
	fclose(outDistributionFile);

	FILE* concentrationFile = fopen("./output/concentration.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(concentrationFile, "%g %g\n", downstreamXgrid[i], downstreamConcentrationArray[i][0][0]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(concentrationFile, "%g %g\n", upstreamXgrid[i], upstreamConcentrationArray[i][0][0]);
	}
	fclose(concentrationFile);

	FILE* concentrationFileBrinkmann = fopen("./output/concentration_Brinkmann.dat", "w");
	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		fprintf(concentrationFileBrinkmann, "%g %g\n", downstreamXgridBrinkmann[i], downstreamConcentrationArrayBrinkmann[i][0][0]);
	}
	fclose(concentrationFileBrinkmann);

	/*FILE* outDiffusionConvectionFile = fopen("./output/diffusionConvection.dat", "w");
	double tempE[3] = { 48, 160, 480 };
	for (int i = 1; i < downstreamNx - 1; ++i) {
		MassiveParticleIsotropicDistribution* leftDistribution = new MassiveParticleTabulatedIsotropicDistribution(*dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i - 1, 0, 0)));
		MassiveParticleIsotropicDistribution* middleDistribution = new MassiveParticleTabulatedIsotropicDistribution(*dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i, 0, 0)));
		MassiveParticleIsotropicDistribution* rightDistribution = new MassiveParticleTabulatedIsotropicDistribution(*dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i + 1, 0, 0)));
		fprintf(outDiffusionConvectionFile, "%g", downstreamXgrid[i]);
		for (int j = 0; j < 3; ++j) {
			double D = tempE[j] * speed_of_light / (3 * electron_charge * downstreamB[i][0][0]);
			double Dd2f = fabs(2*D * ((rightDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i + 1][0][0] - middleDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i][0][0]) / (downstreamXgrid[i + 1] - downstreamXgrid[i]) -
				(middleDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i][0][0] - leftDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i - 1][0][0]) / (downstreamXgrid[i] - downstreamXgrid[i - 1])) / (downstreamXgrid[i + 1] - downstreamXgrid[i - 1]));
			double udf = fabs(0.15 * 0.2 * speed_of_light * (rightDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i + 1][0][0] - middleDistribution->distributionNormalized(tempE[j]) * downstreamConcentrationArray[i][0][0]) / (downstreamXgrid[i + 1] - downstreamXgrid[i]));
			fprintf(outDiffusionConvectionFile, " %g %g", Dd2f, udf);
		}
		delete[] leftDistribution;
		delete[] middleDistribution;
		delete[] rightDistribution;
		fprintf(outDiffusionConvectionFile, "\n");
	}
	fclose(outDiffusionConvectionFile);*/


	int Ne = 10;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", downstreamSource, 1.6E-12, 1.6E4, 1000);


	double Ephmin = 1.6E-12;
	double Ephmax = 1.6E4;
	int Nph = 1000;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	FILE* outFile = fopen("./output/W50synchandcompt.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		//omp_set_lock(&my_lock);
		printf("writeEFEFromSourceToFile iph = %d\n", i);
		printLog("writeEFEFromSourceToFile iph = %d\n", i);
		//omp_unset_lock(&my_lock);
		//double flux1 = sumEvaluator->evaluateFluxFromSource(currentE, downstreamSource);
		double flux1 = 0;
		double flux2 = sumEvaluator->evaluateFluxFromSource(currentE, upstreamSource);
		double fluxBrinkmann = 0;
		double fluxHead = 0;
		double fluxCone = 0;
		int j;
#pragma omp parallel for private(j) shared(currentE, downstreamSource, downstreamXgrid, downstreamNx) reduction(+:flux1, fluxHead, fluxCone)
		for (j = 0; j < downstreamNx; ++j) {
			double flux = sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, j, 0);
			flux1 += flux;
			if ((downstreamXgrid[j] >= headMaxX) && (downstreamXgrid[j] <= headMinX)) {
				fluxHead += flux;
			}
			if ((downstreamXgrid[j] >= coneMaxX) && (downstreamXgrid[j] <= coneMinX)) {
				fluxCone += flux;
			}
		}
#pragma omp parallel for private(j) shared(currentE, downstreamSourceBrinkmann, downstreamXgridBrinkmann, downstreamNxBrinkmann) reduction(+:fluxBrinkmann)
		for (j = 0; j < downstreamNxBrinkmann; ++j) {
			double flux = sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSourceBrinkmann, j, 0);
			fluxBrinkmann += flux;
		}

		fprintf(outFile, "%g %g %g %g %g\n", currentE / 1.6E-12, currentE * (flux1 + flux2), currentE * fluxHead, currentE * fluxCone, currentE*fluxBrinkmann);
		currentE = currentE * factor;
	}
	fclose(outFile);

	double* profileXMM = new double[downstreamNx + upstreamNx];
	double* profileNuSTAR = new double[downstreamNx + upstreamNx];
	double* profileXMMBrinkmann = new double[downstreamNxBrinkmann];
	int irho;

	Ephmin = 0.3 * 1000 * 1.6E-12;
	Ephmax = 10 * 1000 * 1.6E-12;
	Nph = 20;
	omp_lock_t lock;
	omp_init_lock(&lock);
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, downstreamSource, downstreamNx, sumEvaluator, Nph, profileXMM, lock)
	for (irho = 0; irho < downstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = downstreamSource->getCrossSectionArea(irho, 0);
		double d = downstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileXMM[irho] = localFlux;
	}
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, upstreamSource, downstreamNx, upstreamNx, sumEvaluator, Nph, profileXMM, lock)
	for (irho = 0; irho < upstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = upstreamSource->getCrossSectionArea(irho, 0);
		double d = upstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, upstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileXMM[downstreamNx + irho] = localFlux;
	}

#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, downstreamSourceBrinkmann, downstreamNxBrinkmann, sumEvaluator, Nph, profileXMMBrinkmann, lock)
	for (irho = 0; irho < downstreamNxBrinkmann; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating profile irho = %d\n", irho);
		printLog("evaluating profile irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = downstreamSourceBrinkmann->getCrossSectionArea(irho, 0);
		double d = downstreamSourceBrinkmann->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSourceBrinkmann, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileXMMBrinkmann[irho] = localFlux;
	}
	FILE* xmmFile = fopen("./output/xmmprofile.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(xmmFile, "%g %g\n", downstreamXgrid[i], profileXMM[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(xmmFile, "%g %g\n", upstreamXgrid[i], profileXMM[downstreamNx + i]);
	}
	fclose(xmmFile);

	FILE* xmmFileBrinkmann = fopen("./output/xmmprofile_Brinkmann.dat", "w");
	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		fprintf(xmmFileBrinkmann, "%g %g\n", downstreamXgridBrinkmann[i], profileXMMBrinkmann[i]);
	}
	fclose(xmmFileBrinkmann);

	Ephmin = 10 * 1000 * 1.6E-12;
	Ephmax = 20 * 1000 * 1.6E-12;
	Nph = 20;
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, downstreamSource, downstreamNx, sumEvaluator, Nph, profileNuSTAR, lock)
	for (irho = 0; irho < downstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating image irho = %d\n", irho);
		printLog("evaluating image irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = downstreamSource->getCrossSectionArea(irho, 0);
		double d = downstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileNuSTAR[irho] = localFlux;
	}
#pragma omp parallel for private(irho) shared(Ephmin, Ephmax, upstreamSource, downstreamNx, upstreamNx, sumEvaluator, Nph, profileNuSTAR, lock)
	for (irho = 0; irho < upstreamNx; ++irho) {
		omp_set_lock(&lock);
		printf("evaluating image irho = %d\n", irho);
		printLog("evaluating image irho = %d\n", irho);
		omp_unset_lock(&lock);
		double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
		double currentE = Ephmin;
		double localFlux = 0;
		double s = upstreamSource->getCrossSectionArea(irho, 0);
		double d = upstreamSource->getDistance();
		for (int ie = 0; ie < Nph; ++ie) {
			double dE = currentE * (factor - 1.0);
			localFlux += (1.0 / currentE) * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, upstreamSource, irho, 0) * dE * d * d / s;
			currentE = currentE * factor;
		}
		profileNuSTAR[downstreamNx + irho] = localFlux;
	}

	omp_destroy_lock(&lock);

	FILE* nustarFile = fopen("./output/nustarprofile.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(nustarFile, "%g %g\n", downstreamXgrid[i], profileNuSTAR[i]);
	}
	for (int i = 0; i < upstreamNx; ++i) {
		fprintf(nustarFile, "%g %g\n", upstreamXgrid[i], profileNuSTAR[downstreamNx + i]);
	}
	fclose(nustarFile);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", downstreamSource, 1.6E-1, 1.6E3, 300);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", downstreamSource, 1.6E-9, 50*1.6E-9, 300);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	return;
}

void evaluateW50comptonDiffusion() {
	double distance = (18000 / 3.26) * parsec;

	const char* concentrationFileName = "../PLUTO/Tools/pyPLUTO/density.dat";
	const char* BFileName = "../PLUTO/Tools/pyPLUTO/B.dat";
	const char* VFileName = "../PLUTO/Tools/pyPLUTO/velocity.dat";

	FILE* concentrationFile = fopen(concentrationFileName, "r");
	int Nx, Ny, Nz;
	fscanf(concentrationFile, "%d %d %d", &Nz, &Nx, &Ny);
	Ny = Nz;
	double minX, maxX, minY, maxY, minZ, maxZ;
	fscanf(concentrationFile, "%lf %lf %lf", &minZ, &minX, &minY);
	fscanf(concentrationFile, "%lf %lf %lf", &maxZ, &maxX, &maxY);

	minZ = -maxZ;
	//minY = 0;
	minY = -maxZ;
	maxY = maxZ;

	minX = -3E20;
	maxX = 3E20;
	minZ = -2E20;
	maxZ = 2E20;
	minY = -2E20;
	maxY = 2E20;

	Nx = 90;
	Nz = 60;
	Ny = 60;

	double dx = (maxX - minX) / Nx;
	double dy = (maxY - minY) / Ny;
	double dz = (maxZ - minZ) / Nz;

	double* xgrid = new double[Nx];
	for (int i = 0; i < Nx; ++i) {
		xgrid[i] = minX + (i + 0.5) * dx;
	}

	double* ygrid = new double[Ny];
	for (int i = 0; i < Ny; ++i) {
		ygrid[i] = minY + (i + 0.5) * dy;
	}

	double* zgrid = new double[Nz];
	for (int i = 0; i < Nz; ++i) {
		zgrid[i] = minZ + (i + 0.5) * dz;
	}

	double time = 30000 * pi * 1E7;

	double*** B;
	double*** Btheta;
	double*** Bphi;
	double*** V;
	double*** Vtheta;
	double*** Vphi;
	double*** ambientConcentration;

	printf("reading source arrays\n");
	printLog("reading source arrays\n");

	RadiationSourceFactory::readRectangularSourceArraysWithVFromFile(B, Btheta, Bphi, V, Vtheta, Vphi, ambientConcentration, minX, maxX, minZ, maxZ, minY, maxY, Nx, Nz, Ny, SourceInputGeometry::CYLINDRICAL, BFileName, VFileName, concentrationFileName, pi / 2, 0.0, 0.0);

	double*** Vx = new double** [Nx];
	double*** Vy = new double** [Nx];
	double*** Vz = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		Vx[i] = new double* [Nz];
		Vy[i] = new double* [Nz];
		Vz[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			Vx[i][j] = new double[Ny];
			Vy[i][j] = new double[Ny];
			Vz[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				double v = V[i][j][k];
				double theta = Vtheta[i][j][k];
				double phi = Vphi[i][j][k];
				Vx[i][j][k] = v * sin(theta) * cos(phi);
				Vy[i][j][k] = v * sin(theta) * sin(phi);
				Vz[i][j][k] = v * cos(theta);
			}
		}
	}

	printf("reading source distribution\n");

	const char* fileName = "./examples_data/W50/lowfield/GLE_pdf_sf1.dat";

	//to remove zeros
	int Nmomentum = 100;
	/*FILE* file = fopen(fileName, "r");
	while (!feof(file)) {
		double a;
		double b;
		fscanf(file, "%lf %lf", &a, &b);
		Nmomentum = Nmomentum + 1;
	}
	fclose(file);
	Nmomentum = Nmomentum - 1;*/


	double* p = new double[Nmomentum];
	double* distribution = new double[Nmomentum];

	FILE* file = fopen(fileName, "r");
	for (int i = 0; i < Nmomentum; ++i) {
		double x;
		double y;
		fscanf(file, "%lf", &x);
		fscanf(file, "%lf", &y);
		if (y < 0) {
			printf("input distribution < 0\n");
			printLog("input distribution < 0\n");
			exit(0);
		}
		double m_c2 = massElectron * speed_of_light2;


		//energy[i] = sqrt(x * x * speed_of_light2 + m_c2 * m_c2);
		p[i] = x*massProton/massElectron;

		distribution[i] = y * cube(massElectron * speed_of_light) /(pow(massProton*speed_of_light,3)*x*x*x*x);

	}
	fclose(file);

	double norm = distribution[0] * p[0] * p[0] * (p[1] - p[0]);
	for (int i = 1; i < Nmomentum; ++i) {
		norm += distribution[i] *p[i]*p[i]*(p[i] - p[i - 1]);
	}
	norm *= 4 * pi;
	for (int i = 0; i < Nmomentum; ++i) {
		//substitution f = p^3 f
		distribution[i] *= p[i]*p[i]*p[i];
	}
	double concentration1 = norm;

	double size = 5 * parsec;
	double u = 0.25 * 0.2 * speed_of_light;
	double sourcePower = 0.25 * pi * size * size * u;

	double**** F = create4dArray(Nx, Nz, Ny, Nmomentum, 0.0);

	double** sourceCoords = new double* [2];
	sourceCoords[0] = new double[3];
	sourceCoords[1] = new double[3];
	sourceCoords[0][0] = -1E20;
	sourceCoords[0][1] = 0.0;
	sourceCoords[0][2] = 0.0;
	sourceCoords[1][0] = 1E20;
	sourceCoords[1][1] = 0.0;
	sourceCoords[1][2] = 0.0;

	printf("solving diffusion\n");
	printLog("solving diffusion\n");

	solveDiffusion(F, B, Vx, Vz, Vy, xgrid, zgrid, ygrid, p, Nx, Nz, Ny, Nmomentum, distribution, sourcePower, time, sourceCoords, 2);

	printf("outputing distribution\n");
	printLog("outputing distribution\n");

	FILE* distributionFile = fopen("./output/diffusionDistribution.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					fprintf(distributionFile, "%g\n", F[i][j][k][l]);
				}
			}
		}
	}
	fclose(distributionFile);

	/*FILE* distributionFile = fopen("./output/diffusionDistribution.dat", "r");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				for (int l = 0; l < Nmomentum; ++l) {
					fscanf(distributionFile, "%lf", &F[i][j][k][l]);
				}
			}
		}
	}
	fclose(distributionFile);*/

	printf("creating distributions\n");
	printLog("creating distributions\n");

	for (int i = 0; i < Nmomentum; ++i) {
		p[i] = p[i] * massElectron*speed_of_light;
	}

	MassiveParticleDistribution**** electrons = new MassiveParticleDistribution***[Nx];
	double*** concentrationArray = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		electrons[i] = new MassiveParticleDistribution**[Nz];
		concentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			electrons[i][j] = new MassiveParticleDistribution * [Ny];
			concentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				concentrationArray[i][j][k] = 4 * pi * F[i][j][k][0] * (p[1] - p[0]) / p[0];
				F[i][j][k][0] /= (p[0] * p[0] * p[0]);
				for (int l = 1; l < Nmomentum; ++l) {
					concentrationArray[i][j][k] += 4*pi*F[i][j][k][l] * (p[l] - p[l - 1]) / p[l];
					F[i][j][k][l] /= (p[l] * p[l] * p[l]);
				}
				if (concentrationArray[i][j][k] <= 0) {
					F[i][j][k][0] = 1.0;
				}
				electrons[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, p, F[i][j][k], Nmomentum, DistributionInputType::MOMENTUM_FP);
			}
		}
	}

	FILE* concentrationFile2 = fopen("concentration.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				fprintf(concentrationFile2, "%g\n", concentrationArray[i][j][k]);
			}
		}
	}
	fclose(concentrationFile2);

	FILE* distributionFile2 = fopen("F.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				fprintf(distributionFile2, "%g\n", concentrationArray[i][j][k] * F[i][j][k][Nmomentum-1]);
			}
		}
	}
	fclose(distributionFile2);

	printf("creating source\n");
	printLog("creating source\n");

	RadiationSource* source = new RectangularSourceInhomogenousDistribution(Nx, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, minX, maxX, minY, maxY, minZ, maxZ, distance);

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

	int Ne = 1000;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 1000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	printf("evaluating radiation\n");
	printLog("evaluating radiation\n");

	comptonEvaluator->writeEFEFromSourceToFile("./output/W50comptonDiffusion.dat", source, 1.6E-10, 1.6E4, 500);
}

void sourceCoordinates(const double& time, double& x, double& y, double& z) {
	x = 1.5E19;
	y = 0;
	z = 0;
}

void sourceCoordinates2(const double& time, double& x, double& y, double& z) {
	x = -1.5E19;
	y = 0;
	z = 0;
}

double sourcePower(const double& time) {
	double u = 0.25 * 0.1 * speed_of_light;
	double r = 1E19;
	double concentration = 5E-4;
	return pi*r*r*concentration*u;
}

double* getDiffusionCoefficient(double* energy, const int Ne) {
	//todo
	double* D = new double[Ne];
	double E0 = 1.6E-12 * 1E10;

	for (int i = 0; i < Ne; ++i) {
		D[i] = 5E28 * pow(energy[i] / E0, 1.0 / 3.0);
	}

	return D;
}

void evaluateW50pion() {
	double distance = (18000 / 3.26) * parsec;

	const char* concentrationFileName = "../PLUTO/Tools/pyPLUTO/density.dat";
	const char* BFileName = "../PLUTO/Tools/pyPLUTO/B.dat";

	FILE* concentrationFile = fopen(concentrationFileName, "r");
	int Nx, Ny, Nz;
	fscanf(concentrationFile, "%d %d %d", &Nz, &Nx, &Ny);
	Ny = Nz;
	double minX, maxX, minY, maxY, minZ, maxZ;
	fscanf(concentrationFile, "%lf %lf %lf", &minZ, &minX, &minY);
	fscanf(concentrationFile, "%lf %lf %lf", &maxZ, &maxX, &maxY);

	minZ = -maxZ;
	minY = 0;
	//minY = -maxZ;
	maxY = maxZ;

	Nx = 200;
	Nz = 100;
	Ny = 50;

	double dx = (maxX - minX) / Nx;
	double dy = (maxY - minY) / Ny;
	double dz = (maxZ - minZ) / Nz;

	double time = 30000 * pi * 1E7;
	int Nt = 100;

	double*** B;
	double*** Btheta;
	double*** Bphi;
	double*** ambientConcentration;

	RadiationSourceFactory::readRectangularSourceArraysFromFile(B, Btheta, Bphi, ambientConcentration, minX, maxX, minZ, maxZ, minY, maxY, Nx, Nz, Ny, SourceInputGeometry::CYLINDRICAL, BFileName, concentrationFileName, pi / 2, 0.0, 0.0);
	//fix unphysical concentration in the center
	for (int i = 0; i < Nx; ++i) {
		double x = minX + (i + 0.5) * dx;
		for (int j = 0; j < Nz; ++j) {
			double z = minZ + (j + 0.5) * dz;
			for (int k = 0; k < Ny; ++k) {
				double y = minY + (k + 0.5) * dy;
				double r = sqrt(x * x + y * y + z * z);
				if (r < 0.9E19) {
					ambientConcentration[i][j][k] = 0.0;
				}
			}
		}
	}

	const char* protonsFileName = "./examples_data/W50/protons.dat";

	MassiveParticleTabulatedIsotropicDistribution* protons;
	double sourceConcentration;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileName, protons, sourceConcentration);


	int Ne = protons->getN();
	double* energy = protons->getEnergyArray();
	double* sourceDistribution = protons->getDistributionArray();

	MassiveParticleDistribution**** distributions;
	double*** concentration;

	MassiveParticleDistribution**** distributions2;
	double*** concentration2;

	double* D = getDiffusionCoefficient(energy, Ne);

	RadiationSourceFactory::createRectangularSourceArraysFromDiffusion(massProton, energy, sourceDistribution, Ne, D, Nx, Ny, Nz, minX, maxX, minY, maxY, minZ, maxZ, time, Nt, sourceCoordinates, sourcePower, distributions, concentration);
	RadiationSourceFactory::createRectangularSourceArraysFromDiffusion(massProton, energy, sourceDistribution, Ne, D, Nx, Ny, Nz, minX, maxX, minY, maxY, minZ, maxZ, time, Nt, sourceCoordinates2, sourcePower, distributions2, concentration2);

	FILE* concentrationFile2 = fopen("concentration.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				fprintf(concentrationFile2, "%g\n", concentration[i][j][k]);
			}
		}
	}
	fclose(concentrationFile2);

	FILE* distributionFile2 = fopen("distribution.dat", "w");
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				double* d1 = (dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(distributions[i][j][k]))->getDistributionArray();
				fprintf(distributionFile2, "%g\n", concentration[i][j][k]*d1[100]);
				delete[] d1;
			}
		}
	}
	fclose(distributionFile2);

	//hack - ambient and cr concentrations participate as multiplication so we use namb = 1 ncr = ncr*namb
	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				double* d1 = (dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(distributions[i][j][k]))->getDistributionArray();
				double* d2 = (dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(distributions2[i][j][k]))->getDistributionArray();
				for (int l = 0; l < Ne; ++l) {
					d1[l] = concentration[i][j][k]*d1[l] + concentration2[i][j][k]*d2[l];
				}

				delete distributions[i][j][k];
				delete distributions2[i][j][k];

				distributions[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massProton, energy, d1, Ne, DistributionInputType::ENERGY_FE);

				delete[] d1;
				delete[] d2;

				concentration[i][j][k] = concentration[i][j][k] + concentration2[i][j][k];
				concentration[i][j][k] = concentration[i][j][k] * ambientConcentration[i][j][k];
			}
		}
	}

	RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nx, Ny, Nz, distributions, B, Btheta, Bphi, concentration, minX, maxX, minY, maxY, minZ, maxZ, distance);

	PionDecayEvaluatorKelner* evaluator = new PionDecayEvaluatorKelner(1000, massProton * speed_of_light2, 1E8 * massProton * speed_of_light2, 1.0);

	//evaluator->writeEFEFromSourceToFile("W50pion.dat", downstreamSource, 1.6E-7, 1.6E4, 500);

	printf("start writing image\n");
	evaluator->writeImageFromSourceToFile("W50pionImageGeV.dat", source, 1.6E-3, 1.6E-2, 20);
	evaluator->writeImageFromSourceToFile("W50pionImageTeV.dat", source, 1.6E0, 1.6E1, 20);
	evaluator->writeImageFromSourceToFile("W50pionImagePeV.dat", source, 1.6E3, 1.6E4, 20);
}

int main() {
	resetLog();
	srand(time(NULL));
	//evaluateSimpleSynchrotron();
	//evaluateComptonWithPowerLawDistribution();
	//fitCSS161010withPowerLawDistribition();
	//fitCSS161010withTabulatedDistributions();
	//fitTimeDependentCSS161010();
    //evaluatePionDecay();
	//evaluateBremsstrahlung();
	//compareComptonSynchrotron();
	//evaluateSynchrotronImage();
	//testRotation();
	//testAnisotropicCompton();
	//compareComptonWithPowerLawDistribution();
	//testReadingSource();
	//fitAngleDependentFlux();
	//varyParameters();
	//testVersin();
	//testBessel();
	//testChevalier();
	//fitCSS161010_2();
	//testMatrixInverse();
	//testGMRES();
	//testNishinaLosses();

	//evaluateFluxSNRtoWind();
	//evaluateComtonFromWind();
	//evaluateTychoProfile();
	//fitTychoProfile();
	//evaluateSynchrotronInWideRange();
	//evaluateW50bemsstrahlung();
	//evaluateW50synchrotron();
	//evaluateW50comptonAndSynchrotron();
	//evaluateW50comptonAndSynchrotron2();
	//evaluateW50comptonAndSynchrotronMCfunctionUpstream();
	//evaluateW50comptonAndSynchrotronAdvectionfunction();
	//evaluateW50comptonThickRegime();
	//evaluateW50comptonAdvectionBigSource();
	//evaluateW50comptonAndSynchrotronMCwithoutupstream();
	//evaluateW50comptonAndSynchrotronAdvectionfunctionWithUpstream();
	evaluateW50comptonAndSynchrotronAdvectionfunctionWithBrinkmann();
	//evaluateW50comptonDiffusion();
	//evaluateW50pion();

	return 0;
}
