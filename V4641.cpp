#include "stdio.h"
#include "math.h"
#include <omp.h>
#include <time.h>

#include "V4641.h"

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

void evaluateV4641comptonAndSynchrotronAdvectionfunction() {
	double distance = (20200 / 3.26) * parsec;

	double* energy;
	double* downstreamXgrid;
	double* upstreamXgrid;
	double* concentration;
	double** distributions;

	int Nenergy;


	double secondToRadian = pi / (180 * 3600);
	double headMinSec = 0;
	double headMaxSec = 12 * 15;
	double coneMinSec = headMaxSec;
	double coneMaxSec = 26 * 15;

	double headMinX = -headMinSec * secondToRadian * distance;
	double headMaxX = -headMaxSec * secondToRadian * distance;
	double coneMinX = -coneMinSec * secondToRadian * distance;
	double coneMaxX = -coneMaxSec * secondToRadian * distance;

	const char* xfileName = "./examples_data/V4641/x_grid.dat";

	const char* pfileName = "./examples_data/V4641/p_grid.dat";

	const char* fileName = "./examples_data/V4641/electrons.dat";

	int Nx = 0;
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

	double* concentration1;
	double** distributions1;

	//double electronToProtonCorrection = 3E-7;

	MassiveParticleTabulatedIsotropicDistribution* frontElectrons;
	double concentration2;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, frontElectrons, concentration2);
	//frontElectrons = new MassiveParticleTabulatedIsotropicDistribution(new MassiveParticlePowerLawDistribution(massElectron, 2.0, me_c2), me_c2, 1600, 1000);


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



	double size = 0.5 * fabs(headMaxX);
	double B0 = 1E-6;

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


	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[i] = -downstreamXgrid[i];
	}



	/*for (int i = 0; i < downstreamNx; ++i) {
		if (sqrt(Bpar[i] * Bpar[i] + 2 * Bper[i] * Bper[i]) < minField) {
			Bpar[i] = minField / sqrt(3.0);
			Bper[i] = minField / sqrt(3.0);
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
				downstreamB[i][j][k] = B0;
				downstreamB1[i] = downstreamB[i][j][k];
				//downstreamB[i][j][k] = 2E-5;
				//par - x, per - y and z
				downstreamBtheta[i][j][k] = pi / 2;
				//downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphi[i][j][k] = 0;
				downstreamConcentrationArray[i][j][k] = 1.0;
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


	int Nediff = frontElectrons->getN();
	double* energyGrid = frontElectrons->getEnergyArray();
	double* frontDistribution = frontElectrons->getDistributionArray();
	//from SS433
	double electronToProtonCorrection = 1.24E-8;
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


	//MassiveParticleIsotropicDistribution* distribution0 = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
	//MassiveParticleTabulatedIsotropicDistribution* electrons2 = new MassiveParticleTabulatedIsotropicDistribution(distribution0, me_c2, 1E16 * 1.6E-12, 1000);
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
	RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 10.3E8, 10.3E8, photonEnergyDensity);

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


	FILE* outXfile = fopen("./output/x_grid.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(outXfile, "%g\n", downstreamXgrid[i]);
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
	fclose(outDistributionFile);

	FILE* concentrationFile = fopen("./output/concentration.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(concentrationFile, "%g %g\n", downstreamXgrid[i], downstreamConcentrationArray[i][0][0]);
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


	double Ephmin = 1.6E-18;
	double Ephmax = 1.6E4;
	int Nph = 400;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	FILE* outFile = fopen("./output/V4641synchandcompt.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		//omp_set_lock(&my_lock);
		printf("writeEFEFromSourceToFile iph = %d\n", i);
		printLog("writeEFEFromSourceToFile iph = %d\n", i);
		//omp_unset_lock(&my_lock);
		//double flux1 = sumEvaluator->evaluateFluxFromSource(currentE, downstreamSource);
		double flux1 = 0;
		int j;
#pragma omp parallel for private(j) shared(currentE, downstreamSource, downstreamXgrid, downstreamNx) reduction(+:flux1)
		for (j = 0; j < downstreamNx; ++j) {
			double flux = sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, j, 0);
			flux1 += flux;
		}

		fprintf(outFile, "%g %g %g %g %g\n", currentE / 1.6E-12, currentE * flux1, 0, 0, 0);
		currentE = currentE * factor;
	}
	fclose(outFile);


	return;
}

void evaluateV4641comptonThickRegime() {
	double distance = (20200 / 3.26) * parsec;

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
	double B0 = 1E-6;
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

	const char* fileName = "./examples_data/V4641/B1FEB6/electrons.dat";
	const char* farFileName = "./examples_data/V4641/B1FEB6/fardownstreamelectrons.dat";




	MassiveParticleTabulatedIsotropicDistribution* electrons1;
	double concentration1;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileName, electrons1, concentration1);

	//from SS433
	double electronToProtonCorrection = 1.24E-8;

	MassiveParticleTabulatedIsotropicDistribution* fardownstreamDistribution = new MassiveParticleTabulatedIsotropicDistribution(massElectron, farFileName, DistributionInputType::ENERGY_FE);
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

	sumEvaluator->writeEFEFromSourceToFile("./output/V4641thickCompton.dat", source, 1.6E-18, 1.6E4, 400);
	sumEvaluator->writeEFEFromSourceToFile("./output/V4641thickCompton2.dat", source2, 1.6E-18, 1.6E4, 400);

	return;
}


void evaluateV4641comptonAndSynchrotronAdvectionfunctionChangingB() {
	double distance = (20200 / 3.26) * parsec;

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

	//double size = 0.5 * fabs(headMaxX);
	double size = 1.5E19;

	const char* xfileName = "./examples_data/V4641/B14FEB6/x_grid.dat";
	const char* BfileName = "./examples_data/V4641/B14FEB6/Beff.dat";


	const char* fileName = "./examples_data/V4641/B14FEB6/electrons.dat";
	const char* protonsFileName = "./examples_data/V4641/B14FEB6/protons.dat";


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

	//double electronToProtonCorrection = concentration3 * frontProtons->getDistributionArray()[70] / (concentration2 * frontElectrons->getDistributionArray()[70]);

	int Ne = frontElectrons->getN();
	double* electronDistributionArray = frontElectrons->getDistributionArray();
	double* electronEnergy = frontElectrons->getEnergyArray();

	int leftBound = 0;
	double leftEnergy;
	for (int i = 0; i < Ne; ++i) {
		if (electronDistributionArray[i] > 0) {
			leftBound = i;
			leftEnergy = electronEnergy[i] - me_c2 + massProton * speed_of_light2;
			break;
		}
	}
	int rightBound = Ne - 1;
	double rightEnergy;
	for (int i = Ne - 1; i >= 0; --i) {
		if (electronDistributionArray[i] > 0) {
			rightBound = i - 1;
			rightEnergy = electronEnergy[i - 1] - me_c2 + massProton * speed_of_light2;
			break;
		}
	}

	double electronToProtonCorrection = frontProtons->evaluateDistributionInRange(200, leftEnergy, rightEnergy);


	double* Beff = new double[Nx];
	FILE* Bfile = fopen(BfileName, "r");
	for (int i = 0; i < Nx; ++i) {
		fscanf(Bfile, "%lf", &Beff[i]);
		//Beff[i] = 3E-6;
	}
	fclose(Bfile);

	int zeroIndex = 0;
	for (int i = 0; i < Nx; ++i) {
		if (xgrid1[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}

	double downstreamVelocity = 10.3E8;
	double timeSource = 12000 * 3.14E7;
	//double downstreamSize = 1E20;
	double downstreamSize = downstreamVelocity*timeSource;
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
	double* downstreamB1 = new double[downstreamNx];
	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[downstreamNx - i - 1] = xgrid1[i + zeroIndex];
		downstreamB1[downstreamNx - 1 - i] = Beff[i + zeroIndex];
	}

	
	
	double B0 = 1E-6;

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
	double* Bpar = getUvarovBpar2(downstreamNx, downstreamXgrid, L0, 14.0);
	double* Bper = getUvarovBper2(downstreamNx, downstreamXgrid, L0, 14.0);
	//double* Bpar = getUvarovBpar2new(downstreamNx, downstreamXgrid, L0, 0.6);
	//double* Bper = getUvarovBper2new(downstreamNx, downstreamXgrid, L0, 0.4);
	double L1 = 3E19;
	//double* Bpar1 = getUvarovBpar2(downstreamNx, downstreamXgrid, L1, 0.125);
	//double* Bper1 = getUvarovBper2(downstreamNx, downstreamXgrid, L1, 0.125);
	//double* Bpar1 = getUvarovBpar2new(downstreamNx, downstreamXgrid, L1, 0.125);
	//double* Bper1 = getUvarovBper2new(downstreamNx, downstreamXgrid, L1, 0.125);


	for (int i = 0; i < downstreamNx; ++i) {
		downstreamXgrid[i] = -downstreamXgrid[i];
	}


	double minField = 1.0E-6;
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
				/*downstreamB[i][j][k] = downstreamB1[i];
				downstreamBtheta[i][j][k] = pi / 2;
				downstreamBphi[i][j][k] = 0;
				downstreamConcentrationArray[i][j][k] = 1.0;*/

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

	FILE* BoutputFile = fopen("./output/Bturb.dat", "w");

	for (int i = 0; i < downstreamNx; ++i) {
		//fprintf(Bfile, "%g %g %g\n", -downstreamXgrid[i]/L0 + 1.5, Bpar[i], Bper[i]);
		fprintf(BoutputFile, "%g %g %g %g\n", downstreamXgrid[i], Bpar[i], Bper[i], downstreamB1[i]);
		//fprintf(BoutputFile, "%g %g %g\n", downstreamXgrid[i], 0.0, downstreamB[i][0][0]);
	}

	fclose(BoutputFile);



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


	/*int Nediff = frontElectrons->getN();
	double* energyGrid = frontElectrons->getEnergyArray();
	double* frontDistribution = frontElectrons->getDistributionArray();
	for (int i = 0; i < Nediff; ++i) {
		frontDistribution[i] = frontDistribution[i] * concentration3 * electronToProtonCorrection;
	}
	double norm1 = 4 * pi * MassiveParticleDistributionFactory::evaluateNorm(energyGrid, frontDistribution, Nediff);
	double** diffDistributions = NULL;*/
	double Uph[1];
	double Eph[1];
	Uph[0] = photonEnergyDensity;
	//Uph[1] = photonIRenergyDensity;
	Eph[0] = 2.8 * kBoltzman * 2.725;
	//Eph[1] = 2.8 * kBoltzman * 140;
	



	//MassiveParticleIsotropicDistribution* distribution0 = new MassiveParticlePowerLawCutoffDistribution(massElectron, 2.0, me_c2, 1.0, 1E15 * 1.6E-12);
	//MassiveParticleTabulatedIsotropicDistribution* electrons2 = new MassiveParticleTabulatedIsotropicDistribution(distribution0, me_c2, 1E16 * 1.6E-12, 1000);
	//upstreamElectrons = frontElectrons;
	//double E0 = 1.6E-1;
	//MassiveParticleIsotropicDistribution* upstreamElectrons = new MassiveParticleMonoenergeticDistribution(massElectron, E0, 0.01 * E0);

	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				downstreamConcentrationArray[i][j][k] = concentration3 * electronToProtonCorrection;
			}
		}
	}

	double timeOf = 3000 * 3.14E7;

	//TabulatedDiskSourceWithSynchAndComptCutoff* downstreamSource = new TabulatedDiskSourceWithSynchAndComptCutoff(Nrho, Nz, 1, upstreamElectrons, B0, pi / 2, 0, concentration, size, size, distance, 0.25 * 0.1 * speed_of_light, photonEnergyDensity);
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, 0.15 * 0.2 * speed_of_light, photonTotalEnergyDensity);
	
	//RectangularSourceWithSynchAndComptCutoffFromRight* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, downstreamVelocity, downstreamVelocity, photonEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRightFieldDecay* downstreamSource = new RectangularSourceWithSynchAndComptCutoffFromRightFieldDecay(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance, downstreamVelocity, downstreamVelocity, timeOf, 200*3.14E7, 1E-6, photonEnergyDensity);
	

	/*MassiveParticleDistribution**** distributions2 = new MassiveParticleDistribution * **[downstreamNx];


	for (int i = 0; i < downstreamNx; ++i) {
		distributions2[i] = new MassiveParticleDistribution * *[1];
		distributions2[i][0] = new MassiveParticleDistribution * [1];
		distributions2[i][0][0] =new MassiveParticleTabulatedIsotropicDistribution(* dynamic_cast<MassiveParticleTabulatedIsotropicDistribution*>(downstreamSource->getParticleDistribution(i, 0, 0)));
	}


	double pmin = 0.1 * massProton / massElectron;
	double pmax = 5E6 * massProton / massElectron;
	int Np = 100;
	double factorp = pow(pmax / pmin, 1.0 / (Np - 1.0));

	FILE* outDistributionFile = fopen("./output/pdf.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		double p = pmin;
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
	fclose(outDistributionFile);*/

	FILE* concentrationFile = fopen("./output/concentration.dat", "w");
	for (int i = 0; i < downstreamNx; ++i) {
		fprintf(concentrationFile, "%g %g\n", downstreamXgrid[i], downstreamConcentrationArray[i][0][0]);
	}
	fclose(concentrationFile);

	for (int i = 0; i < downstreamNx; ++i) {
		for (int j = 0; j < Nz; ++j) {
			for (int k = 0; k < Ny; ++k) {
				downstreamB[i][j][k] = B0;
			}
		}
	}

	/*for (int i = 0; i < downstreamNx; ++i) {
		if (downstreamXgrid[i] > -timeOf * downstreamVelocity) {
			for (int j = 0; j < Nz; ++j) {
				for (int k = 0; k < Ny; ++k) {
					downstreamConcentrationArray[i][j][k] = 0;
				}
			}
		}
	}*/

	//RectangularSourceInhomogenousDistribution* downstreamSource2 = new RectangularSourceInhomogenousDistribution(downstreamNx, downstreamXgrid, Ny, Nz, distributions2, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);
	//RectangularSource* downstreamSource2 = new RectangularSource(downstreamNx, downstreamXgrid, Ny, Nz, frontElectrons, downstreamB, downstreamBtheta, downstreamBphi, downstreamConcentrationArray, 0, size, 0, pi * size, distance);


	Ne = 100;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 200, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", downstreamSource, 1.6E-12, 1.6E4, 1000);


	double Ephmin = 1.6E-18;
	double Ephmax = 1.6E4;
	int Nph = 400;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	FILE* outFile = fopen("./output/V4641synchandcompt.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		//omp_set_lock(&my_lock);
		printf("writeEFEFromSourceToFile iph = %d\n", i);
		printLog("writeEFEFromSourceToFile iph = %d\n", i);
		//omp_unset_lock(&my_lock);
		//double flux1 = sumEvaluator->evaluateFluxFromSource(currentE, downstreamSource);
		double flux1 = 0;
		int j;
#pragma omp parallel for private(j) shared(currentE, downstreamSource, downstreamXgrid, downstreamNx) reduction(+:flux1)
		for (j = 0; j < downstreamNx; ++j) {
			double flux = sumEvaluator->evaluateFluxFromSourceAtPoint(currentE, downstreamSource, j, 0);
			flux1 += flux;
		}

		fprintf(outFile, "%g %g %g %g %g\n", currentE / 1.6E-12, currentE * flux1 , 0, 0, 0);
		currentE = currentE * factor;
	}
	fclose(outFile);


	return;
}

void evaluateV4641comptonAndSynchrotronWind()
{
	double distance = (20200 / 3.26) * parsec;

	double* energy;
	double* downstreamXgrid;
	double* upstreamXgrid;
	double* concentration;
	double** distributions;

	int Nenergy;
	int Nx;
	int Ny = 1;
	int Nz = 1;



	//4 for quasispherical?
	double sizeForward = 4E19;
	double sizeTermination = 4E19;

	double rForward = 1E20;
	double rTermination = 1E20;

	double B0 = 4E-6;

	//const char* BfileNameF = "./examples_data/V4641/ForwardWind/Beff.dat";


	const char* fileNameF = "./examples_data/V4641/ForwardWind/electrons.dat";
	const char* protonsFileNameF = "./examples_data/V4641/ForwardWind/protons.dat";
	const char* xfileNameF = "./examples_data/V4641/ForwardWind/x_grid.dat";

	//const char* BfileNameT = "./examples_data/V4641/TerminationWind/Beff.dat";


	const char* fileNameT = "./examples_data/V4641/TerminationWind/electrons.dat";
	const char* protonsFileNameT = "./examples_data/V4641/TerminationWind/protons.dat";
	const char* xfileNameT = "./examples_data/V4641/TerminationWind/x_grid.dat";
	

	int NxF = 0;
	FILE* xfile = fopen(xfileNameF, "r");
	while (!feof(xfile)) {
		double a;
		fscanf(xfile, "%lf", &a);
		NxF = NxF + 1;
	}
	fclose(xfile);
	NxF = NxF - 1;
	double* xgrid1F = new double[NxF];
	xfile = fopen(xfileNameF, "r");
	for (int i = 0; i < NxF; ++i) {
		fscanf(xfile, "%lf", &xgrid1F[i]);
	}
	fclose(xfile);

	int NxT = 0;
	xfile = fopen(xfileNameT, "r");
	while (!feof(xfile)) {
		double a;
		fscanf(xfile, "%lf", &a);
		NxT = NxT + 1;
	}
	fclose(xfile);
	NxT = NxT - 1;
	double* xgrid1T = new double[NxT];
	xfile = fopen(xfileNameT, "r");
	for (int i = 0; i < NxT; ++i) {
		fscanf(xfile, "%lf", &xgrid1T[i]);
	}
	fclose(xfile);

	int zeroIndex = 0;
	for (int i = 0; i < NxT; ++i) {
		if (xgrid1F[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}

	int maxIndex = NxF - 1;
	for (int i = 0; i < NxF; ++i) {
		if (xgrid1F[i] >= sizeForward) {
			maxIndex = i;
			break;
		}
	}

	int downstreamNxF = maxIndex + 1 - zeroIndex;

	double* downstreamXgridF = new double[downstreamNxF];
	double* downstreamB1F = new double[downstreamNxF];
	for (int i = 0; i < downstreamNxF; ++i) {
		downstreamXgridF[downstreamNxF - i - 1] = xgrid1F[i + zeroIndex];
		downstreamB1F[downstreamNxF - 1 - i] = B0;
	}

	zeroIndex = 0;
	for (int i = 0; i < NxT; ++i) {
		if (xgrid1T[i] >= 0) {
			zeroIndex = i;
			break;
		}
	}

	maxIndex = NxT - 1;
	for (int i = 0; i < NxT; ++i) {
		if (xgrid1T[i] >= sizeTermination) {
			maxIndex = i;
			break;
		}
	}

	int downstreamNxT = maxIndex + 1 - zeroIndex;

	double* downstreamXgridT = new double[downstreamNxT];
	double* downstreamB1T = new double[downstreamNxT];
	for (int i = 0; i < downstreamNxT; ++i) {
		downstreamXgridT[downstreamNxT - i - 1] = xgrid1T[i + zeroIndex];
		downstreamB1T[downstreamNxT - 1 - i] = B0;
	}

	for (int i = 0; i < downstreamNxF; ++i) {
		downstreamXgridF[i] = -downstreamXgridF[i];
	}

	for (int i = 0; i < downstreamNxT; ++i) {
		downstreamXgridT[i] = -downstreamXgridT[i];
	}
	

	MassiveParticleTabulatedIsotropicDistribution* frontElectronsF;
	double concentration2F;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileNameF, frontElectronsF, concentration2F);
	//frontElectrons = new MassiveParticleTabulatedIsotropicDistribution(new MassiveParticlePowerLawDistribution(massElectron, 2.0, me_c2), me_c2, 1600, 1000);
	MassiveParticleTabulatedIsotropicDistribution* frontProtonsF;
	double concentration3F;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileNameF, frontProtonsF, concentration3F);

	MassiveParticleTabulatedIsotropicDistribution* frontElectronsT;
	double concentration2T;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massElectron, fileNameT, frontElectronsT, concentration2T);
	//frontElectrons = new MassiveParticleTabulatedIsotropicDistribution(new MassiveParticlePowerLawDistribution(massElectron, 2.0, me_c2), me_c2, 1600, 1000);
	MassiveParticleTabulatedIsotropicDistribution* frontProtonsT;
	double concentration3T;
	MassiveParticleDistributionFactory::readTabulatedIsotropicDistributionFromMonteCarlo(massProton, protonsFileNameT, frontProtonsT, concentration3T);

	frontElectronsF->writeDistribution("forward.dat", 200, me_c2, 1E10 * me_c2);
	frontElectronsT->writeDistribution("backward.dat", 200, me_c2, 1E10 * me_c2);

	int Ne = frontElectronsF->getN();
	double* electronDistributionArrayF = frontElectronsF->getDistributionArray();
	double* electronEnergyF = frontElectronsF->getEnergyArray();

	int leftBoundF = 0;
	double leftEnergyF;
	for (int i = 0; i < Ne; ++i) {
		if (electronDistributionArrayF[i] > 0) {
			leftBoundF = i;
			leftEnergyF = electronEnergyF[i] - me_c2 + massProton * speed_of_light2;
			break;
		}
	}
	int rightBoundF = Ne - 1;
	double rightEnergyF;
	for (int i = Ne - 1; i >= 0; --i) {
		if (electronDistributionArrayF[i] > 0) {
			rightBoundF = i - 1;
			rightEnergyF = electronEnergyF[i - 1] - me_c2 + massProton * speed_of_light2;
			break;
		}
	}

	Ne = frontElectronsT->getN();
	double* electronDistributionArrayT = frontElectronsT->getDistributionArray();
	double* electronEnergyT = frontElectronsT->getEnergyArray();

	int leftBoundT = 0;
	double leftEnergyT;
	for (int i = 0; i < Ne; ++i) {
		if (electronDistributionArrayT[i] > 0) {
			leftBoundT = i;
			leftEnergyT = electronEnergyT[i] - me_c2 + massProton * speed_of_light2;
			break;
		}
	}
	int rightBoundT = Ne - 1;
	double rightEnergyT;
	for (int i = Ne - 1; i >= 0; --i) {
		if (electronDistributionArrayT[i] > 0) {
			rightBoundT = i - 1;
			rightEnergyT = electronEnergyT[i - 1] - me_c2 + massProton * speed_of_light2;
			break;
		}
	}

	//double electronToProtonCorrectionF = frontProtonsF->getDistributionArray()[70] / (frontElectronsF->getDistributionArray()[70]);
	double electronToProtonCorrectionF = frontProtonsF->evaluateDistributionInRange(200, leftEnergyF, rightEnergyF);
	//double electronToProtonCorrectionT = frontProtonsT->getDistributionArray()[70] / (frontElectronsT->getDistributionArray()[70]);
	double electronToProtonCorrectionT = frontProtonsT->evaluateDistributionInRange(200, leftEnergyT,rightEnergyT);

	double concentration0F = 4 * 2E-3*electronToProtonCorrectionF;
	double concentration0T = 2E-3 * electronToProtonCorrectionT;



	double*** downstreamBF = new double** [downstreamNxF];
	double*** downstreamBthetaF = new double** [downstreamNxF];
	double*** downstreamBphiF = new double** [downstreamNxF];
	double*** downstreamConcentrationArrayF = new double** [downstreamNxF];
	for (int i = 0; i < downstreamNxF; ++i) {
		downstreamBF[i] = new double* [Nz];
		downstreamBthetaF[i] = new double* [Nz];
		downstreamBphiF[i] = new double* [Nz];
		downstreamConcentrationArrayF[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			downstreamBF[i][j] = new double[Ny];
			downstreamBthetaF[i][j] = new double[Ny];
			downstreamBphiF[i][j] = new double[Ny];
			downstreamConcentrationArrayF[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				downstreamBF[i][j][k] = B0;
				downstreamBthetaF[i][j][k] = pi / 2;
				downstreamBphiF[i][j][k] = 0;
				downstreamConcentrationArrayF[i][j][k] = concentration0F;
			}
		}
	}

	double*** downstreamBT = new double** [downstreamNxT];
	double*** downstreamBthetaT = new double** [downstreamNxT];
	double*** downstreamBphiT = new double** [downstreamNxT];
	double*** downstreamConcentrationArrayT = new double** [downstreamNxT];
	for (int i = 0; i < downstreamNxT; ++i) {
		downstreamBT[i] = new double* [Nz];
		downstreamBthetaT[i] = new double* [Nz];
		downstreamBphiT[i] = new double* [Nz];
		downstreamConcentrationArrayT[i] = new double* [Nz];
		for (int j = 0; j < Nz; ++j) {
			downstreamBT[i][j] = new double[Ny];
			downstreamBthetaT[i][j] = new double[Ny];
			downstreamBphiT[i][j] = new double[Ny];
			downstreamConcentrationArrayT[i][j] = new double[Ny];
			for (int k = 0; k < Ny; ++k) {
				downstreamBT[i][j][k] = B0;
				downstreamBthetaT[i][j][k] = pi / 2;
				downstreamBphiT[i][j][k] = 0;
				downstreamConcentrationArrayT[i][j][k] = concentration0T;
			}
		}
	}

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

	double downstreamVelocityF = 100000000.0 / 4.0;
	double downstreamVelocityT = 300000000.0 / 4.0;

	RectangularSourceWithSynchAndComptCutoffFromRight* forwardSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNxF, downstreamXgridF, Ny, Nz, frontElectronsF, downstreamBF, downstreamBthetaF, downstreamBphiF, downstreamConcentrationArrayF, 0, 4*rForward, 0, pi * rForward, distance, downstreamVelocityF, downstreamVelocityF, photonEnergyDensity);
	RectangularSourceWithSynchAndComptCutoffFromRight* backwardSource = new RectangularSourceWithSynchAndComptCutoffFromRight(downstreamNxT, downstreamXgridT, Ny, Nz, frontElectronsT, downstreamBT, downstreamBthetaT, downstreamBphiT, downstreamConcentrationArrayT, 0, 4*rTermination, 0, pi * rTermination, distance, downstreamVelocityT, downstreamVelocityT, photonEnergyDensity);

	//DiskSource* forwardSource = new SimpleFlatSource(frontElectronsT, B0, pi / 2, 0, concentration0F, rForward, sizeForward, distance);
	//DiskSource* backwardSource = new SimpleFlatSource(frontElectronsT, B0, pi / 2, 0, concentration0T, rTermination, sizeTermination, distance);

	Ne = 100;
	int Nmu = 100;
	int Nphi = 4;
	//RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 2000, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photonsTotal, photonTotalConcentration, ComptonSolverType::ISOTROPIC_JONES);
	RadiationEvaluator* comptonEvaluator = new InverseComptonEvaluator(Ne, Nmu, Nphi, me_c2 * 500, 1E10 * me_c2, 200, 0.1 * kBoltzman * 2.75, 30 * kBoltzman * 20, photons, photonConcentration, ComptonSolverType::ISOTROPIC_JONES);

	//comptonEvaluator->writeEFEFromSourceToFile("W50compton.dat", downstreamSource, 1.6E-10, 1.6E3, 2000);

	RadiationEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, false);

	//synchrotronEvaluator->writeEFEFromSourceToFile("W50synchrotron.dat", downstreamSource, 1.6E-18, 1.6E-5, 2000);

	RadiationSumEvaluator* sumEvaluator = new RadiationSumEvaluator(Ne, me_c2 * 500, 1E10 * me_c2, comptonEvaluator, synchrotronEvaluator, false);

	//sumEvaluator->writeEFEFromSourceToFile("./output/W50synchandcompt.dat", downstreamSource, 1.6E-12, 1.6E4, 1000);


	double Ephmin = 1.6E-18;
	double Ephmax = 1.6E4;
	int Nph = 400;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nph - 1));
	double currentE = Ephmin;
	FILE* outFile = fopen("./output/V4641wind.dat", "w");
	for (int i = 0; i < Nph; ++i) {
		//omp_set_lock(&my_lock);
		printf("writeEFEFromSourceToFile iph = %d\n", i);
		printLog("writeEFEFromSourceToFile iph = %d\n", i);
		//omp_unset_lock(&my_lock);
		//double flux1 = sumEvaluator->evaluateFluxFromSource(currentE, downstreamSource);
		double flux1 = 0;
		int j;
#pragma omp parallel for private(j) shared(currentE, downstreamSource, downstreamXgrid, downstreamNx) reduction(+:flux1)
		double fluxF = sumEvaluator->evaluateFluxFromSource(currentE, forwardSource);
		double fluxT = sumEvaluator->evaluateFluxFromSource(currentE, backwardSource);

		fprintf(outFile, "%g %g %g %g %g\n", currentE / 1.6E-12, currentE * (fluxF + fluxT), currentE*fluxF, currentE*fluxT, 0);
		currentE = currentE * factor;
	}
	fclose(outFile);


	return;
}
