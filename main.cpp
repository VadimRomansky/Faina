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
	//initializing electrons distribution
	//MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, 10*me_c2, electronConcentration);
	//MassiveParticleBrokenPowerLawDistribution* electrons = new MassiveParticleBrokenPowerLawDistribution(massElectron, index, 2.001, 2*me_c2, 1000*me_c2, electronConcentration);
    //MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_combined_cutoff/Ee3.dat", "./examples_data/gamma0.2_combined_cutoff/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
    MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_theta0-90/Ee3.dat", "./examples_data/gamma0.2_theta0-90/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleMaxwellJuttnerDistribution* electrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, 1.08 * me_c2 / kBoltzman);
	MassiveParticleTabulatedIsotropicDistribution* protons = new MassiveParticleTabulatedIsotropicDistribution(massProton, "./examples_data/gamma1.5_theta0-90_protons/Ee3.dat", "./examples_data/gamma1.5_theta0-90_protons/Fs3.dat", DistributionInputType::GAMMA_KIN_FGAMMA);
	MassiveParticleMaxwellJuttnerDistribution* thermalElectrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, 4 * me_c2 / kBoltzman);
	double velocity = 0.2 * speed_of_light;
	//double velocity = 0.55 * speed_of_light;
	double gamma = 1.0 / sqrt(1.0 - velocity * velocity / speed_of_light2);
	//electrons->rescaleDistribution(1.2);
	//electrons->addPowerLaw(10 * massElectron * speed_of_light2, 3.0);
	//(dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(electrons))->prolongEnergyRange(1E6, 100);
	//electrons->addPowerLaw(200 * massElectron * speed_of_light2, 3.5);
	electrons->writeDistribution("distribution.dat", 200, Emin, Emax);
	thermalElectrons->writeDistribution("distribution1.dat", 200, Emin, Emax);
	

	//creating radiation source
	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, rmax, fraction*rmax, distance, 0, 0.2433);
	//RadiationSource* source = new TabulatedSphericalLayerSource(20, 20, 1,electrons, B, theta, electronConcentration, rmax, (1-fraction)*rmax, distance);
	RadiationSourceInCylindrical* thermalSource = new SimpleFlatSource(thermalElectrons, B, theta, 0, electronConcentration, rmax, fraction*rmax, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, electronConcentration, rmax, 0.5*rmax, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electronsFromSmilei, Bturb, thetaTurb, concentration, rsource, 0.9*rsource, distance);
	//AngleDependentElectronsSphericalSource* source = new AngleDependentElectronsSphericalSource(Nrho, Nz, Nphi, Ndistributions, angleDependentDistributions, Bturb, thetaTurb, phiTurb, concentration, rmax, 0.9*rmax, distance, 0.5*speed_of_light);
	//counting of quai parallel
	/*int nangle = 0;
	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double angle = (dynamic_cast<AngleDependentElectronsSphericalSource*>(source))->getShockWaveAngle(irho, iz, iphi);
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
	//number of parameters of the source
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
	
	//enumOptimizer->optimize(vector, optPar, source);
	//gradient descent optimization
	//synchrotronOptimizer->optimize(vector, optPar);
	//reseting source parameters to found values
	//synchrotronOptimizer->outputProfileDiagrams(vector, source, 10);
	//synchrotronOptimizer->outputOptimizedProfileDiagram(vector, optPar, source, 10, 1, 2);
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

	/*double totalVolume = source->getTotalVolume();
	double totalVolume2 = 0;
	for (int i = 0; i < 20; ++i) {
		for (int j = 0; j < 20; ++j) {
			for (int k = 0; k < 1; ++k) {
				totalVolume2 += source->getVolume(i, j, k);
			}
		}
	}

	double area1 = pi * rmax * rmax * (1 - (1 - fraction) * (1 - fraction));
	double area = 0;
	for (int i = 0; i < 20; ++i) {
		for (int k = 0; k < 1; ++k) {
			area += source->getArea(i, 10, k);
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

	//TabulatedDiskSourceWithSynchCutoff* source3 = new TabulatedDiskSourceWithSynchCutoff(1, 10000, 1, electrons, B, pi / 2, 0, electronConcentration, rmax, fraction * rmax, distance, 0.2 * speed_of_light);


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
	//double totalEnergy = electronConcentration * massProton * speed_of_light2 * (gamma - 1.0)*source->getTotalVolume();
	double energyInRadioElectrons = electronConcentration * (electrons->getMeanEnergy() - me_c2) * (pi * rmax * rmax * rmax * fraction);
	//double energyInRadioProtons = electronConcentration * (protons->getMeanEnergy() - massProton*speed_of_light2) * source->getTotalVolume();
	double magneticEnergy = (B * B / (8 * pi)) * (pi * rmax * rmax * rmax * fraction);
	printf("total kinetik energy = %g\n", totalEnergy);
	printLog("total kinetik energy = %g\n", totalEnergy);
	printf("energy in radio electrons = %g\n", energyInRadioElectrons);
	printLog("energy in radio electrons = %g\n", energyInRadioElectrons);
	//printf("energy in radio protons = %g\n", energyInRadioProtons);
	//printLog("energy in radio protons = %g\n", energyInRadioProtons);
	printf("energy in magnetic field = %g\n", magneticEnergy);
	printLog("energy in magnetic field = %g\n", magneticEnergy);
	//double totalMass = electronConcentration * massProton * source->getTotalVolume();
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

	//evaluationg full spectrum of the source
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
	//synchrotronEvaluator->writeFluxFromSourceToFile("outputSynch3.dat", source, hplank * 1E8, 10 * 1000 * 1.6E-12, 1000);
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
    InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDirectedDistribution, photonDirectedConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluatorWithSource(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES, rmax + 0.5E16, 0, 0);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluatorWithSource(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA, 0, -rmax - 0.5E16, 0);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_THOMSON);
	//initializing photon energy grid for output
	InverseComptonEvaluator* galacticEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, galacticField, photonGalacticConcentration, ComptonSolverType::ISOTROPIC_JONES);
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
	//double kevFlux = comptonEvaluator->evaluateTotalFluxInEnergyRange(minEev, maxEev, 10, source);
	double kevFlux = 0;
	int Nph = 100;
	factor = pow(maxEev / minEev, 1.0 / (Nph - 1));
	currentE = minEev;
	FILE* output7 = fopen("output7.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		printf("%d\n", i);
		double dE = currentE * (factor - 1.0);
		//double dE = currentE * (1.0 - 1.0/factor);
		//kevFlux += comptonEvaluator->evaluateFluxFromSourceAnisotropic(currentE, 0, 0, photonDirectedDistribution, source, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA) * dE;
		double flux = comptonEvaluator->evaluateFluxFromSource(currentE, source2);
		//double flux = galacticEvaluator->evaluateFluxFromSource(currentE, source2);
		fprintf(output7,"%g %g\n", currentE, flux);
		kevFlux += flux * dE;
		currentE = currentE * factor;
	}
	fclose(output7);
	double totalLuminosity = kevFlux * 4 * pi * distance * distance;
	//double synchrotronKevFlux = synchrotronEvaluator->evaluateTotalFluxInEnergyRange(minEev, maxEev, 10, source);
	//double synchrotronFluxGHz3 = synchrotronEvaluator->evaluateFluxFromSource(3E9 * hplank, source)*1E26*hplank;
	
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

	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.3_theta0-90/Ee3.dat", "./examples_data/gamma0.3_theta0-90/Fs3.dat", 200, GAMMA_KIN_FGAMMA);
	//MassiveParticleIsotropicDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, 2*Emin, me_c2 / 2);
	//electrons->addPowerLaw(1.02 * me_c2, 3);
	//electrons->rescaleDistribution(sqrt(18));
	//electrons->addPowerLaw(100 * me_c2, 3.5);
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
	//creating radiation source
	RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, theta, 0, electronConcentration, rmax, 0.2*rmax, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, electronConcentration, rmax, 0.8 * rmax, distance);
	//InverseComptonEvaluator* comptonEvaluator1 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator2 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator3 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, photonConcentration, ComptonSolverType::ISOTROPIC_THOMSON);

	//comptonEvaluator2->outputDifferentialFlux("output3.dat");
	//comptonEvaluator->outputDifferentialFluxJones("output2.dat", photonDistribution, electrons);
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
		F[i] = comptonEvaluator->evaluateFluxFromSource(E[i], source);
		//F[i] = comptonEvaluator->evaluateFluxFromSourceAnisotropic(E[i], 0, 0, CMBradiation, source);
	}*/

	//outputing
	FILE* output_ev_EFE1 = fopen("output1.dat", "w");
	FILE* output_ev_EFE2 = fopen("output2.dat", "w");
	FILE* output_ev_EFE3 = fopen("output3.dat", "w");
	//FILE* output_GHz_Jansky = fopen("output.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		double nu = E[i] / hplank;
		//fprintf(output_ev_EFE1, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator1->evaluateFluxFromSource(E[i], source));
		//fprintf(output_ev_EFE2, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator2->evaluateFluxFromSource(E[i], source));
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

	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInDiskSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSphericalSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSectoralSphericalSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, Rmin, 2*pi, anisotropy);
	write3dArrayToFile(B, Nrho, Nz, Nphi, "B.dat");

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 2.0, 100*Energy);

	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance);
	//RadiationSource* source = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance, Udownstream);
	RadiationSourceInCylindrical* source = new TabulatedSectoralSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, 0, concentrations, R, Rin,Rmin, 2*pi, distance, Udownstream);
	//RadiationSource* source = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, R, distance, Udownstream);

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

	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInDiskSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSphericalSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, fraction, index, lturb, Nmodes, R, anisotropy);
	//RadiationSourceFactory::initializeAnisotropicLocalTurbulentFieldInSectoralSphericalSource(B, theta, phi, Nrho, Nz, Nphi, B0, theta0, phi0, turbulentFraction, index, lturb, Nmodes, R, Rmin, 2 * pi, anisotropy);
	write3dArrayToFile(B, Nrho, Nz, Nphi, "B.dat");

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 2.0, 100 * Energy);

	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance);
	//RadiationSource* source = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance, Udownstream);
	RadiationSourceInCylindrical* source = new TabulatedSectoralSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, 0, concentrations, R, Rin, Rmin, 2 * pi, distance, Udownstream);
	//RadiationSource* source = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, R, distance, Udownstream);

	//number of parameters of the source
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


	//TabulatedDiskSourceWithSynchCutoff* source = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution, B, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, velocity);
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
	//evaluator->writeFluxFromSourceToFile("outputSynch.dat", source, 10 * hplank * cyclotronOmega, 100000 * hplank * cyclotronOmega, 1000);
	evaluator->writeEFEFromSourceToFile("W50synchrotron.dat", source, 1.6E-18, 1.6E-5, 2000);


}

double* getUvarovBpar(int Nx, double minX, double maxX) {
	double A0 = 41.0755*1E-6;
	double A1 = 9.27296;
	double A2 = 1.41262;
	double A3 = 3.26888;
	double A4 = 1.88243;
	double A5 = 0.0767564;
	double A6 = 0.00540849;
	double UNIT_LENGTH = 2.5E18;
	double* B = new double [Nx];
	double dx = (maxX - minX) / Nx;
	for (int i = 0; i < Nx; ++i) {
		double x = minX + (i + 0.5) * dx;
		double l = (maxX - x) / UNIT_LENGTH + 2;
		B[i] = (A0 / (1 + pow(sqr((l - A3) / A1), A2))) * (1.0 / (1 + A6 * exp(-sqr((l - A4) / A5))));
	}

	return B;
}

double* getUvarovBper(int Nx, double minX, double maxX) {
	double A0 = 98.3917*1E-6;
	double A1 = 30.533;
	double A2 = 2.33138;
	double A3 = -23.2141;
	double A4 = 18.847;
	double A5 = 1.05756;
	double A6 = 0.695847;
	double UNIT_LENGTH = 2.5E18;
	double* B = new double[Nx];
	double dx = (maxX - minX) / Nx;
	for (int i = 0; i < Nx; ++i) {
		double x = minX + (i + 0.5) * dx;
		double l = (maxX - x) / UNIT_LENGTH + 2;
		B[i] = A0 / (1 + pow(sqr((l - A3) / A1), A2) + A4 * exp(-sqr((l - A5) / A6)));
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

	//RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, pi / 2, 0, concentration, size, size, distance);
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

	double* Bpar = getUvarovBpar(Nrho, 0, size);
	double* Bper = getUvarovBper(Nrho, 0, size);
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

	FILE* Bfile = fopen("Bturb.dat", "w");

	for (int i = 0; i < Nrho; ++i) {
		fprintf(Bfile, "%g %g %g\n", (i + 0.5) * size / Nrho, Bpar[i], Bper[i]);
	}

	fclose(Bfile);

	//TabulatedDiskSourceWithSynchAndComptCutoff* source = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, electrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance, 0.25 * 0.1 * speed_of_light, photonTotalEnergyDensity);
	RectangularSource* source2 = new RectangularSource(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance);

	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 1, 0, 0));
	distributionRight->writeDistribution("distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 2, 0, 0));
	distributionMiddle->writeDistribution("distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(0, 0, 0));
	distributionLeft->writeDistribution("distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	int Ne = 1000;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 0.1*kBoltzman*2.75, 50*kBoltzman*2.75, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", source, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", source, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("W50synchandcompt.dat", source, 1.6E-18, 1.6E3, 200);
	sumEvaluator->writeEFEFromSourceToFile("W50highenergy.dat", source, 1.6E-1, 1.6E3, 200);

	sumEvaluator->writeEFEFromSourceToFile("W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	sumEvaluator->writeEFEFromSourceToFile("W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageeV.dat", source, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageKeV.dat", source, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageMeV.dat", source, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageGeV.dat", source, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageTeV.dat", source, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImagePeV.dat", source, 1.6E3, 1.6E4, 20);

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


	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, source, currentE, factor, sumEvaluator, i)
		for (j = 0; j < Nrho; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outFile = fopen("xE.dat", "w");
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

	//RadiationSourceInCylindrical* source = new SimpleFlatSource(electrons, B, pi / 2, 0, concentration, size, size, distance);
	PhotonPlankDistribution* photons = PhotonPlankDistribution::getCMBradiation();
	PhotonPlankDistribution* photonsIR = new PhotonPlankDistribution(140, 0.8 / 1800000);
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

	double* Bpar = getUvarovBpar(Nrho, 0, size);
	double* Bper = getUvarovBper(Nrho, 0, size);
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
				Btheta[i][j][k] = atan2(Bper[i] / sqrt(2), sqrt(Bpar[i] * Bpar[i] + 0.5 * Bper[i] * Bper[i]));
				Bphi[i][j][k] = pi / 4;
				concentrationArray[i][j][k] = concentration;
			}
		}
	}

	FILE* Bfile = fopen("Bturb.dat", "w");

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

	MassiveParticleDistribution**** electrons2 = new MassiveParticleDistribution ***[Nrho];
	for (int i = 0; i < Nrho; ++i) {
		electrons2[i] = new MassiveParticleDistribution ** [Nz];
		for (int j = 0; j < Nz; ++j) {
			electrons2[i][j] = new MassiveParticleDistribution * [Ny];
			for (int k = 0; k < Ny; ++k) {
				electrons2[i][j][k] = new MassiveParticleTabulatedIsotropicDistribution(massElectron, energy[Nrho - i - 1], distributions[Nrho - i - 1], Nee, DistributionInputType::ENERGY_FE);
			}
		}
	}

	//TabulatedDiskSourceWithSynchAndComptCutoff* source = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, electrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance, 0.25 * 0.1 * speed_of_light, photonTotalEnergyDensity);
	RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nrho, Ny, Nz, electrons2, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance);
	RectangularSource* source2 = new RectangularSource(Nrho, Ny, Nz, electrons, B, Btheta, Bphi, concentrationArray, 0, size, 0, size, 0, pi * size, distance);
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 1, 0, 0));
	distributionRight->writeDistribution("distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nrho - 2, 0, 0));
	distributionMiddle->writeDistribution("distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(0, 0, 0));
	distributionLeft->writeDistribution("distributionLeft.dat", 200, me_c2, 1E10 * me_c2);

	int Ne = 1000;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 0.1 * kBoltzman * 2.75, 50 * kBoltzman * 2.75, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", source, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", source, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("W50synchandcompt.dat", source, 1.6E-18, 1.6E3, 200);
	sumEvaluator->writeEFEFromSourceToFile("W50highenergy.dat", source, 1.6E-1, 1.6E3, 200);

	//sumEvaluator->writeEFEFromSourceToFile("W50synchandcompt2.dat", source2, 1.6E-18, 1.6E3, 200);
	//sumEvaluator->writeEFEFromSourceToFile("W50highenergy2.dat", source2, 1.6E-1, 1.6E3, 200);

	printf("start writing images\n");
	printLog("start writing images\n");

	printf("start writing ev image\n");
	printLog("start writing ev image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageeV.dat", source, 1.6E-12, 1.6E-11, 20);

	printf("start writing keV image\n");
	printLog("start writing keV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageKeV.dat", source, 1.6E-9, 1.6E-8, 20);

	printf("start writing MeV image\n");
	printLog("start writing MeV images\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageMeV.dat", source, 1.6E-6, 1.6E-5, 20);

	printf("start writing GeV image\n");
	printLog("start writing GeV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageGeV.dat", source, 1.6E-3, 1.6E-2, 20);

	printf("start writing TeV image\n");
	printLog("start writing TeV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImageTeV.dat", source, 1.6E0, 1.6E1, 20);

	printf("start writin PeV image\n");
	printLog("start writing PeV image\n");
	sumEvaluator->writeImageFromSourceToFile("W50scImagePeV.dat", source, 1.6E3, 1.6E4, 20);

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


	for (int i = 0; i < Nnu; ++i) {
		printf("inu = %d\n", i);
		printLog("inu = %d\n", i);
		int j;
#pragma omp parallel for private(j) shared(F, source, currentE, factor, sumEvaluator, i)
		for (j = 0; j < Nrho; ++j) {
			F[i][j] = currentE * sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, source, j, 0);
		}
		currentE = factor * currentE;
	}

	FILE* outFile = fopen("xE.dat", "w");
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

	//evaluator->writeEFEFromSourceToFile("W50pion.dat", source, 1.6E-7, 1.6E4, 500);

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
	//testNishinaLosses();

	//evaluateFluxSNRtoWind();
	//evaluateComtonFromWind();
	//evaluateTychoProfile();
	//fitTychoProfile();
	//evaluateSynchrotronInWideRange();
	//evaluateW50bemsstrahlung();
	//evaluateW50synchrotron();
	//evaluateW50comptonAndSynchrotron();
	evaluateW50comptonAndSynchrotron2();
	//evaluateW50pion();

	return 0;
}
