#include "stdio.h"
#include "math.h"
#include <omp.h>
#include <time.h>

#include "W50.h"

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

void evaluateW50bremsstrahlung() {
	double distance = (18000 / 3.26) * parsec;

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


	ThermalRectangularSource* source = RadiationSourceFactory::readThermalRectangularSourceFromFile(minX, maxX, minZ, maxZ, minY, maxY, Nx, Nz, Ny, distance, SourceInputGeometry::CYLINDRICAL, BFileName, concentrationFileName, temperatureFileName, 0.8 * pi / 2, 0, 0);

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
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 250, 0.1 * kBoltzman * 2.75, 10 * kBoltzman * 140, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", source, 1.6E-10, 1.6E3, 100);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50highenergy.dat", source, 1.6E-1, 1.6E3, 100);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", source, 1.6E-9, 50 * 1.6E-9, 100);

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
	sumEvaluator->writeEFEFromSourceToFile("./output/W50kev.dat", source, 1.6E-9, 50 * 1.6E-9, 200);

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
				concentrationArray[Nx - i - 1][j][k] = concentration[i + minIndex] * electronToProtonCorrection;
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
					Bper[Nx - i - 1] = Beff[i + minIndex] / sqrt(2.0);
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
	MassiveParticleIsotropicDistribution* distributionRight = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2 + 200, 0, 0));
	distributionRight->writeDistribution("./output/distributionRight.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionMiddle = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2, 0, 0));
	distributionMiddle->writeDistribution("./output/distributionMiddle.dat", 200, me_c2, 1E10 * me_c2);
	MassiveParticleIsotropicDistribution* distributionLeft = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(Nx / 2 - 200, 0, 0));
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
			double F = distribution->distributionNormalized(E) * concentrationArray[i][0][0];
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
				Btheta[i][j][k] = atan2(Bper[i], sqrt(Bpar[i] * Bpar[i] + Bper[i] * Bper[i]));
				//downstreamBtheta[i][j][k] = pi / 2;
				Bphi[i][j][k] = pi / 4;
				concentrationArray[i][j][k] = 1.0;
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
				concentrationArray[i][j][k] = concentration1 * electrons2->distributionNormalized(E0) * 0.02 * E0;
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


	int Ne = 10;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 140 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

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

	const char* fileName = "./examples_data/W50/B15FEB6/electrons.dat";
	const char* farFileName = "./examples_data/W50/B15FEB6/fardownstreamelectrons.dat";
	const char* farUpFileName = "./examples_data/W50/B15FEB6/farupstreamelectrons.dat";
	const char* protonsFileName = "./examples_data/W50/B15FEB6/protons.dat";




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
	electrons1->transformToThickRegime(photonTotalEnergyDensity + magneticEnergyDensity, norm);
	double norm2 = 0;
	double concentration2 = concentration1;
	fardownstreamDistribution->transformToThickRegime(photonTotalEnergyDensity + magneticEnergyDensity, norm2);
	double norm3 = 0;
	farupstreamDistribution->transformToThickRegime(photonTotalEnergyDensity + magneticEnergyDensity, norm3);
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
	//u = 10.3E8;
	double outOfJetFactor = 1.0;
	concentration1 *= outOfJetFactor * u * pi * size * size * electronToProtonCorrection * norm;
	concentration2 *= u * pi * size * size * electronToProtonCorrection * norm2;

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, downstreamXgrid, Ny, Nz, upstreamElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	//RectangularSourceInhomogenousDistribution* downstreamSource = new RectangularSourceInhomogenousDistribution(Nx, downstreamXgrid, Ny, Nz, electrons2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	RectangularSource* source = new RectangularSource(1, 1, 1, electrons1, B0, pi / 2, 0, concentration1, 0, 1, 0, 1, 0, 1, distance);
	RectangularSource* source2 = new RectangularSource(1, 1, 1, fardownstreamDistribution, B0, pi / 2, 0, concentration2, 0, 1, 0, 1, 0, 1, distance);


	int Ne = 2000;
	int Nmu = 100;
	int Nphi = 4;
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 2.75 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E11 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 2.75, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E11 * me_c2, false, false);
	RadiationEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E11 * me_c2, comptonEvaluator, synchrotronEvaluator, false, false);

	sumEvaluator->writeEFEFromSourceToFile("./output/W50thickCompton.dat", source, 1.6E-12, 1.6E4, 200);
	sumEvaluator->writeEFEFromSourceToFile("./output/W50thickCompton2.dat", source2, 1.6E-12, 1.6E4, 200);

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
	RectangularSourceWithSynchAndComptCutoffFromRight* source = new RectangularSourceWithSynchAndComptCutoffFromRight(Nx, xgrid, Ny, Nz, electrons2, B, Btheta, Bphi, concentrationArray, 0, size, 0, pi * size, distance, 0.25 * 0.2 * speed_of_light, 0.25 * 0.2 * speed_of_light, photonTotalEnergyDensity);
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

	double size = 0.5 * fabs(headMaxX);
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
		if (sqrt(Bpar[i] * Bpar[i] + 2 * Bper[i] * Bper[i]) < minField) {
			Bpar[i] = minField / sqrt(3.0);
			Bper[i] = minField / sqrt(3.0);
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
			if ((downstreamXgrid[j] >= coneMaxX) && (downstreamXgrid[j] <= coneMinX)) {
				fluxCone += flux;
			}
		}
		fprintf(outFile, "%g %g %g %g\n", currentE / 1.6E-12, currentE * (flux1 + flux2), currentE * fluxHead, currentE * fluxCone);
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

	const char* xfileName = "./examples_data/W50/B15FEB6/x_grid.dat";
	const char* BfileName = "./examples_data/W50/B15FEB6/Beff.dat";

	const char* distributionFileName = "./examples_data/W50/newdistribution/electrons_full.dat";
	const char* pfileName = "./examples_data/W50/B15FEB6/p_grid.dat";

	const char* fileName = "./examples_data/W50/B15FEB6/electrons.dat";
	const char* protonsFileName = "./examples_data/W50/B15FEB6/protons.dat";

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
	double sizeBrinkmann = 1.5 * size;

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
	double* Bpar = getUvarovBpar2(downstreamNx, downstreamXgrid, L0, 0.25);
	double* Bper = getUvarovBper2(downstreamNx, downstreamXgrid, L0, 0.25);
	//double* Bpar = getUvarovBpar2new(downstreamNx, downstreamXgrid, L0, 0.6);
	//double* Bper = getUvarovBper2new(downstreamNx, downstreamXgrid, L0, 0.4);
	double L1 = 3E19;
	//double* Bpar1 = getUvarovBpar2(downstreamNx, downstreamXgrid, L1, 0.125);
	//double* Bper1 = getUvarovBper2(downstreamNx, downstreamXgrid, L1, 0.125);
	double* Bpar1 = getUvarovBpar2new(downstreamNx, downstreamXgrid, L1, 0.5);
	double* Bper1 = getUvarovBper2new(downstreamNx, downstreamXgrid, L1, 0.5);

	/*for (int i = 0; i < downstreamNx; ++i) {
		Bpar[i] = Bpar[i] + Bpar1[i];
		Bper[i] = Bper[i] + Bper1[i];
	}*/

	double L0Brinkmann = 2.0E18;
	double* BparBrinkmann = getUvarovBpar2(downstreamNxBrinkmann, downstreamXgridBrinkmann, L0Brinkmann, 9.3 / 60.0);
	double* BperBrinkmann = getUvarovBper2(downstreamNxBrinkmann, downstreamXgridBrinkmann, L0Brinkmann, 9.3 / 60.0);

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

	double minField = 8.0E-6;
	double sinField = 0.0 * minField;

	int minFieldIndex = 0;
	for (int i = 1; i < downstreamNx; ++i) {
		if (sqrt(Bpar[downstreamNx - i - 1] * Bpar[downstreamNx - i - 1] + 2 * Bper[downstreamNx - i - 1] * Bper[downstreamNx - i - 1]) < minField) {
			minFieldIndex = downstreamNx - i;
			break;
		}
	}
	for (int i = 0; i < minFieldIndex; ++i) {
		Bpar[i] = Bpar[minFieldIndex];
		Bper[i] = Bper[minFieldIndex];
	}

	/*for (int i = 0; i < downstreamNx; ++i) {
		if (sqrt(Bpar[i] * Bpar[i] + 2 * Bper[i] * Bper[i]) < minField) {
			Bpar[i] = minField / sqrt(3.0);
			Bper[i] = minField / sqrt(3.0);
		}
	}*/

	double minFieldBrinkmann = 3.0E-6;
	double sinFieldBrinkmann = 0.0 * minFieldBrinkmann;

	int minFieldIndexBrinkmann = 0;
	for (int i = 1; i < downstreamNxBrinkmann; ++i) {
		if (sqrt(BparBrinkmann[downstreamNxBrinkmann - i - 1] * BparBrinkmann[downstreamNxBrinkmann - i - 1] + 2 * BperBrinkmann[downstreamNxBrinkmann - i - 1] * BperBrinkmann[downstreamNxBrinkmann - i - 1]) < minFieldBrinkmann) {
			minFieldIndexBrinkmann = downstreamNxBrinkmann - i;
			break;
		}
	}
	for (int i = 0; i < minFieldIndexBrinkmann; ++i) {
		BparBrinkmann[i] = BparBrinkmann[minFieldIndexBrinkmann];
		BperBrinkmann[i] = BperBrinkmann[minFieldIndexBrinkmann];
	}

	/*for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		if (sqrt(BparBrinkmann[i] * BparBrinkmann[i] + 2 * BperBrinkmann[i] * BperBrinkmann[i]) < minFieldBrinkmann) {
			BparBrinkmann[i] = minFieldBrinkmann / sqrt(3.0);
			BperBrinkmann[i] = minFieldBrinkmann / sqrt(3.0);
		}
	}*/

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

	FILE* BoutputFileBrinkmann = fopen("./output/BturbBrinkmann.dat", "w");

	for (int i = 0; i < downstreamNxBrinkmann; ++i) {
		//fprintf(Bfile, "%g %g %g\n", -downstreamXgrid[i]/L0 + 1.5, Bpar[i], Bper[i]);
		fprintf(BoutputFileBrinkmann, "%g %g %g\n", downstreamXgridBrinkmann[i], BparBrinkmann[i], BperBrinkmann[i]);
		//fprintf(BoutputFile, "%g %g %g\n", downstreamXgrid[i], 0.0, downstreamB[i][0][0]);
	}

	fclose(BoutputFileBrinkmann);


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
	RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 10.3E8, 10.3E8, photonEnergyDensity);
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


	int Ne = 100;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 200, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", downstreamSource, 1.6E-12, 1.6E4, 1000);


	double Ephmin = 1.6E-12;
	double Ephmax = 1.6E4;
	int Nph = 200;
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

		fprintf(outFile, "%g %g %g %g %g\n", currentE / 1.6E-12, currentE * (flux1 + flux2), currentE * fluxHead, currentE * fluxCone, currentE * fluxBrinkmann);
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
		p[i] = x * massProton / massElectron;

		distribution[i] = y * cube(massElectron * speed_of_light) / (pow(massProton * speed_of_light, 3) * x * x * x * x);

	}
	fclose(file);

	double norm = distribution[0] * p[0] * p[0] * (p[1] - p[0]);
	for (int i = 1; i < Nmomentum; ++i) {
		norm += distribution[i] * p[i] * p[i] * (p[i] - p[i - 1]);
	}
	norm *= 4 * pi;
	for (int i = 0; i < Nmomentum; ++i) {
		//substitution f = p^3 f
		distribution[i] *= p[i] * p[i] * p[i];
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
		p[i] = p[i] * massElectron * speed_of_light;
	}

	MassiveParticleDistribution**** electrons = new MassiveParticleDistribution * **[Nx];
	double*** concentrationArray = new double** [Nx];
	for (int i = 0; i < Nx; ++i) {
		electrons[i] = new MassiveParticleDistribution * *[Nz];
		concentrationArray[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			electrons[i][j] = new MassiveParticleDistribution * [Ny];
			concentrationArray[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				concentrationArray[i][j][k] = 4 * pi * F[i][j][k][0] * (p[1] - p[0]) / p[0];
				F[i][j][k][0] /= (p[0] * p[0] * p[0]);
				for (int l = 1; l < Nmomentum; ++l) {
					concentrationArray[i][j][k] += 4 * pi * F[i][j][k][l] * (p[l] - p[l - 1]) / p[l];
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
				fprintf(distributionFile2, "%g\n", concentrationArray[i][j][k] * F[i][j][k][Nmomentum - 1]);
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
	return pi * r * r * concentration * u;
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

	const char* concentrationFileName = "./examples_data/W50/density.dat";
	const char* BFileName = "./examples_data/W50/B.dat";

	FILE* concentrationFile = fopen(concentrationFileName, "r");
	int Nx, Ny, Nz;
	fscanf(concentrationFile, "%d %d %d", &Nz, &Nx, &Ny);
	Ny = Nz;
	double minX, maxX, minY, maxY, minZ, maxZ;
	fscanf(concentrationFile, "%lf %lf %lf", &minZ, &minX, &minY);
	fscanf(concentrationFile, "%lf %lf %lf", &maxZ, &maxX, &maxY);

	//minZ = -maxZ;
	//minY = 0;
	//minY = -maxZ;
	//maxY = maxZ;

	minY = minX / 2;
	maxY = maxX / 2;
	minZ = minY;
	maxZ = maxY;

	Nx = 160;
	Nz = 80;
	Ny = 80;

	double dx = (maxX - minX) / Nx;
	double dy = (maxY - minY) / Ny;
	double dz = (maxZ - minZ) / Nz;

	double time = 100000 * pi * 1E7;
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
				if (r < 0.5E19) {
					ambientConcentration[i][j][k] = 0.0;
				}
			}
		}
	}

	const char* protonsFileName = "./examples_data/W50/B15FEB6/protons.dat";

	MassiveParticleTabulatedIsotropicDistribution* protons;
	double concentration1;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileName, protons, concentration1);


	//int Ne = protons->getN();
	//double* energy = protons->getEnergyArray();
	//double* sourceDistribution = protons->getDistributionArray();

	//MassiveParticleDistribution**** distributions;
	//double*** concentration;

	//MassiveParticleDistribution**** distributions2;
	//double*** concentration2;

	//double* D = getDiffusionCoefficient(energy, Ne);

	//RadiationSourceFactory::createRectangularSourceArraysFromDiffusion(massProton, energy, sourceDistribution, Ne, D, Nx, Ny, Nz, minX, maxX, minY, maxY, minZ, maxZ, time, Nt, sourceCoordinates, sourcePower, distributions, concentration);
	//RadiationSourceFactory::createRectangularSourceArraysFromDiffusion(massProton, energy, sourceDistribution, Ne, D, Nx, Ny, Nz, minX, maxX, minY, maxY, minZ, maxZ, time, Nt, sourceCoordinates2, sourcePower, distributions2, concentration2);

	/*FILE* concentrationFile2 = fopen("concentration.dat", "w");
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
	fclose(distributionFile2);*/

	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	double size = fabs(headMaxX);

	double u = 0.26 * 0.15 * speed_of_light;
	double Nprotons = 2 * time * concentration1 * u * pi * size * size;
	double averageConcentration = Nprotons / ((maxX - minX) * (maxY - minY) * (maxZ - minZ));

	//hack - ambient and cr concentrations participate as multiplication so we use namb = 1 ncr = ncr*namb
	/*for (int i = 0; i < Nx; ++i) {
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
	}*/

	//hack swap ambient and acelerated concentration

	//RectangularSourceInhomogenousDistribution* source = new RectangularSourceInhomogenousDistribution(Nx, Ny, Nz, protons, B, Btheta, Bphi, ambientConcentration, minX, maxX, minY, maxY, minZ, maxZ, distance);
	RectangularSource* source = new RectangularSource(Nx, Ny, Nz, protons, B, Btheta, Bphi, ambientConcentration, minX, maxX, minY, maxY, minZ, maxZ, distance);

	PionDecayEvaluatorKelner* evaluator = new PionDecayEvaluatorKelner(200, massProton * speed_of_light2, 1E8 * massProton * speed_of_light2, averageConcentration);

	evaluator->writeEFEFromSourceToFile("W50pion.dat", source, 1.6E-12, 1.6E4, 200);

	printf("start writing image\n");
	evaluator->writeImageFromSourceToFile("W50pionImageGeV.dat", source, 1.6E-3, 1.6E-2, 20);
	evaluator->writeImageFromSourceToFile("W50pionImageTeV.dat", source, 1.6E0, 1.6E1, 20);
	evaluator->writeImageFromSourceToFile("W50pionImagePeV.dat", source, 1.6E3, 1.6E4, 20);
}