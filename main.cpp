#include "stdio.h"
#include "math.h"
#include <omp.h>

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
    //MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_combined_cutoff/Ee3.dat", "./examples_data/gamma0.2_combined_cutoff/Fs3.dat", electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
    MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.5_theta0-90/Ee3.dat", "./examples_data/gamma0.5_theta0-90/Fs3.dat", electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.2_theta0-90/Ee3.dat", "./examples_data/gamma0.2_theta0-90/Fs3.dat", electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	//MassiveParticleMaxwellJuttnerDistribution* electrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, 1.08 * me_c2 / kBoltzman, electronConcentration);
	MassiveParticleTabulatedIsotropicDistribution* protons = new MassiveParticleTabulatedIsotropicDistribution(massProton, "./examples_data/gamma1.5_theta0-90_protons/Ee3.dat", "./examples_data/gamma1.5_theta0-90_protons/Fs3.dat", electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	MassiveParticleMaxwellJuttnerDistribution* thermalElectrons = new MassiveParticleMaxwellJuttnerDistribution(massElectron, 4 * me_c2 / kBoltzman, electronConcentration);
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
	RadiationSource* source = new SimpleFlatSource(electrons, B, theta, 0, rmax, fraction*rmax, distance, 0, 0.2433);
	//RadiationSource* source = new TabulatedSphericalLayerSource(20, 20, 1,electrons, B, theta, electronConcentration, rmax, (1-fraction)*rmax, distance);
	RadiationSource* thermalSource = new SimpleFlatSource(thermalElectrons, B, theta, 0, electronConcentration, rmax, fraction*rmax, distance);
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
	PhotonIsotropicDistribution* photonDistribution = new PhotonPlankDistribution(Tstar, 0.25*sqr(rstar / rcompton));
	//PhotonIsotropicDistribution* photonDistribution = PhotonMultiPlankDistribution::getGalacticField();
	PhotonPlankDirectedDistribution* photonDirectedDistribution = new PhotonPlankDirectedDistribution(Tstar, 0.25*sqr(rstar / rcompton), pi*17/18, 0, atan2(1, 1));
	PhotonIsotropicDistribution* galacticField = PhotonMultiPlankDistribution::getGalacticField();
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
    InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDirectedDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluatorWithSource(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES, rmax + 0.5E16, 0, 0);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluatorWithSource(Ne, Nmu, Nphi, Emin, newEmax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA, 0, -rmax - 0.5E16, 0);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_THOMSON);
	//initializing photon energy grid for output
	InverseComptonEvaluator* galacticEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, galacticField, ComptonSolverType::ISOTROPIC_JONES);
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
	electrons->resetConcentration(jetelectronConcentration);
	//MassiveParticleDistribution* jetelectrons = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 30 * me_c2, jetelectronConcentration);
	MassiveParticleDistribution* jetelectrons1 = new MassiveParticlePowerLawDistribution(massElectron, 3.5, 1.1 * me_c2, jetelectronConcentration);
	//MassiveParticleDistribution* jetelectrons = new MassiveParticleMonoenergeticDistribution(massElectron, 20 * me_c2, me_c2, jetelectronConcentration);
	//MassiveParticleDistribution* jetelectrons1 = new MassiveParticleMonoenergeticDistribution(massElectron, 1.1 * me_c2, 0.1*me_c2, jetelectronConcentration);
	//MassiveParticleDistribution* jetelectrons = new MassiveParticleMonoenergeticDirectedDistribution(massElectron, 30 * me_c2, me_c2, jetelectronConcentration, 0, 0, 0.1);
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

	jetelectrons->resetConcentration(jetelectronConcentration);
	RadiationSource* source2 = new SimpleFlatSource(jetelectrons, B, pi / 2, 0, electronConcentration, rmax, rmax*fraction, distance);

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
	//double rmax = 3E16;
	double rtotal = 1.4E17;
	double rmax = 1.0 / sqrt(pi);
	double B = 0.01;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	//const double distance = 150 * 1000000 * parsec;
	const double distance = 1.0;

	double Emin = 652.317 * me_c2 * 1;
	//double Emin = 1 * me_c2;
	double Emax = 1E8 * me_c2;
	int Ne = 200;
	int Nmu = 20;
	int Nrho = 20;
	int Nz = 20;
	int Nphi = 4;
	double index = 2.5;
	double KK = 24990.8;
	double electronConcentration = KK / (pow(652.317, index - 1) * (index - 1));
	//double electronConcentration = 5E6;
	double protonBulkConcentration = 1E4;

	double Tstar = 50 * 1000;
	double Ephmin = 0.01 * Tstar * kBoltzman;
	double Ephmax = 100 * Tstar * kBoltzman;
	double luminosity = 510000 * 4 * 1E33;
	double rsun = 7.5E10;
	double rstar = rsun * sqrt(510000.0 / pow(Tstar / 5772, 4));
	//PhotonIsotropicDistribution* photonDistribution = new PhotonPlankDistribution(Tstar, sqr(rstar / rmax));
	PhotonIsotropicDistribution* photonDistribution = PhotonPlankDistribution::getCMBradiation();

	//initializing electrons distribution
	MassiveParticlePowerLawDistribution* electrons = new MassiveParticlePowerLawDistribution(massElectron, index, Emin, electronConcentration);
	//MassiveParticleTabulatedIsotropicDistribution* electrons = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma0.3_theta0-90/Ee3.dat", "./examples_data/gamma0.3_theta0-90/Fs3.dat", 200, electronConcentration, GAMMA_KIN_FGAMMA);
	//MassiveParticleIsotropicDistribution* electrons = new MassiveParticleMonoenergeticDistribution(massElectron, 2*Emin, me_c2 / 2, electronConcentration);
	//electrons->addPowerLaw(1.02 * me_c2, 3);
	//electrons->rescaleDistribution(sqrt(18));
	//electrons->addPowerLaw(100 * me_c2, 3.5);
	electrons->writeDistribution("dist1.dat", 2000, Emin, Emax);
	double electronsEnergy = electrons->getConcentration() * (electrons->getMeanEnergy() - me_c2);
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
	RadiationSource* source = new SimpleFlatSource(electrons, B, theta, 0, electronConcentration, rmax, 0.2*rmax, distance);
	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, electronConcentration, rmax, 0.8 * rmax, distance);
	InverseComptonEvaluator* comptonEvaluator1 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_JONES);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator2 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_KLEIN_NISHINA);
	InverseComptonEvaluator* comptonEvaluator3 = new InverseComptonEvaluator(100, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ANISOTROPIC_KLEIN_NISHINA);
	//InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, Emin, Emax, Ephmin, Ephmax, photonDistribution, ComptonSolverType::ISOTROPIC_THOMSON);

	//comptonEvaluator2->outputDifferentialFlux("output3.dat");
	//comptonEvaluator->outputDifferentialFluxJones("output2.dat", photonDistribution, electrons);
	//return;

	double minEev = 0.3 * 1000 * 1.6E-12;
	double maxEev = 10 * 1000 * 1.6E-12;
	double kevFlux = comptonEvaluator1->evaluateTotalFluxInEnergyRange(minEev, maxEev, 20, source);
	double totalLuminosity = kevFlux * 4 * pi * distance * distance;
	FILE* outFile = fopen("SNRtoWindData.dat", "w");
	printf("total luminosity = %g erg/s \n", totalLuminosity);
	fprintf(outFile, "total luminosity = %g erg/s \n", totalLuminosity);
	printf("total flux = %g erg/s cm^2 \n", kevFlux);
	printf("99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	fprintf(outFile, "99 days F = 1.33+-0.76 10^-15 L = 3.4+-1.9 10^39\n");
	fclose(outFile);

	//initializing photon energy grid for output
	int Nnu = 200;
	double* E = new double[Nnu];
	double* F = new double[Nnu];

	double EphFinalmin = 0.0001 * kBoltzman * Tstar;
	double EphFinalmax = 10 * Emax + Emin;
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
		fprintf(output_ev_EFE1, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator1->evaluateFluxFromSource(E[i], source));
		fprintf(output_ev_EFE2, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator2->evaluateFluxFromSource(E[i], source));
		fprintf(output_ev_EFE3, "%g %g\n", E[i] / (1.6E-9), comptonEvaluator3->evaluateFluxFromSource(E[i], source));
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
	delete comptonEvaluator1;
	delete comptonEvaluator2;
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

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 2.0, 100*Energy, concentration);

	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance);
	//RadiationSource* source = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance, Udownstream);
	RadiationSource* source = new TabulatedSectoralSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, 0, concentrations, R, Rin,Rmin, 2*pi, distance, Udownstream);
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

	MassiveParticleIsotropicDistribution* electrons = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 2.0, 100 * Energy, concentration);

	//RadiationSource* source = new TabulatedSphericalLayerSource(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance);
	//RadiationSource* source = new TabulatedSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, concentrations, R, (1.0 - widthFraction) * R, distance, Udownstream);
	RadiationSource* source = new TabulatedSectoralSLSourceWithSynchCutoff(Nrho, Nz, Nphi, electrons, B, theta, 0, concentrations, R, Rin, Rmin, 2 * pi, distance, Udownstream);
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
	double R = 2E17;
	double B = 0.24;
	double f = 0.4;
	double electronConcentration = 200;
	double distance = 150 * 1000000 * parsec;

	double velocity = 0.7 * speed_of_light;
	double downstreamV = 0.25 * velocity;
	
	MassiveParticleTabulatedIsotropicDistribution* electronDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, "./examples_data/gamma1.5_combined_cutoff/Ee3.dat", "./examples_data/gamma1.5_combined_cutoff/Fs3.dat", electronConcentration, DistributionInputType::GAMMA_KIN_FGAMMA);
	
	int Nrho = 1;
	int Nz = 1000;
	int Nphi = 1;


	TabulatedDiskSourceWithSynchCutoff* source = new TabulatedDiskSourceWithSynchCutoff(Nrho, Nz, Nphi, electronDistribution, B, pi / 2, 0, electronConcentration, R, f * R, distance, downstreamV, velocity);

	int Ne = 1000;
	double Emin = me_c2;
	double Emax = me_c2 * 1E8;
	SynchrotronEvaluator* evaluator = new SynchrotronEvaluator(Ne, Emin, Emax, true, true);

	double kevFlux = evaluator->evaluateTotalFluxInEnergyRange(0.3 * keV, 10 * keV, 100, source);

	double mevFlux = evaluator->evaluateTotalFluxInEnergyRange(0.1 * MeV, 3 * MeV, 100, source);

	printf("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);
	printLog("keV flux = %g, luminosity = %g\n", kevFlux, kevFlux * 4 * pi * distance * distance);

	printf("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);
	printLog("MeV flux = %g, luminocity = %g\n", mevFlux, mevFlux * 4 * pi * distance * distance);

	evaluator->writeFluxFromSourceToFile("wideRangeSynch.dat", source, 1E8 * hplank, 100 * MeV, 1000);
}


int main() {
	resetLog();
	//evaluateSimpleSynchrotron();
	//evaluateComptonWithPowerLawDistribution();
	//fitCSS161010withPowerLawDistribition();
	//fitCSS161010withTabulatedDistributions();
	fitTimeDependentCSS161010();
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

	//evaluateFluxSNRtoWind();
	//evaluateComtonFromWind();
	//evaluateTychoProfile();
	//fitTychoProfile();
	//evaluateSynchrotronInWideRange();

	return 0;
}
