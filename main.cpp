#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "electronDistribution.h"
#include "photonDistribution.h"
#include "util.h"
#include "inverseCompton.h"
#include "radiationSource.h"
#include "synchrotron.h"


int main() {
	//CSS161010 99 days
	double electronConcentration = 150;
	double B = 0.6;
	double sinTheta = sin(pi / 6);
	double rmax = 1.3E17;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 3.08 * 1.0E24;
	double Emin = me_c2;
	double Emax = 10000 * me_c2;

	PhotonPlankDistribution* CMBradiation = PhotonPlankDistribution::getCMBradiation();
	PhotonMultiPlankDistribution* galacticRadiation = PhotonMultiPlankDistribution::getGalacticField();
	ElectronPowerLawDistribution* electrons = new ElectronPowerLawDistribution(3.5, Emin, electronConcentration);
	RadiationSource* source = new SimpleFlatSource(electrons, B, sinTheta, rmax, rmax, distance);
	RadiationSource* sphericalSource = new TabulatedSphericalLayerSource(20, 20, 4, electrons, B, 1.0, rmax, 0.5 * rmax, distance);
	double volume = source->getTotalVolume();

	InverseComptonEvaluator comptonEvaluator = InverseComptonEvaluator(200, 20, 20, Emin, Emax);
	SynchrotronEvaluator synchrotronEvaluator = SynchrotronEvaluator(200, Emin, Emax);

	int Ndistributions = 10;
	ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs",".dat",10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, 1.0, 0, rmax, 0.5 * rmax, distance);


	int Nnu = 200;
	double* E = new double[Nnu];
	double* I = new double[Nnu];

	double Ephmin = 0.0001 * kBoltzman * 2.7;
	double Ephmax = 2*10000*me_c2 + Emin;
	double factor = pow(Ephmax / Ephmin, 1.0 / (Nnu - 1));
	E[0] = Ephmin;
	I[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		I[i] = 0;
	}

	printLog("evaluating\n");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		//I[i] = comptonEvaluator.evaluateComptonIsotropicFluxFromSource(E[i], CMBradiation, source);
	}

	double* Nu = new double[Nnu];
	double* Is = new double[Nnu];

	double Numin = 1E8;
	double Numax = 1E18;
	factor = pow(Numax / Numin, 1.0 / (Nnu - 1));
	Nu[0] = Numin;
	Is[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		Nu[i] = Nu[i - 1] * factor;
		Is[i] = 0;
	}

	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		//Is[i] = evaluateSynchrotronFluxFromSource(source, Nu[i], Emin, Emax, 200);
		//Is[i] = synchrotronEvaluator.evaluateSynchrotronFluxFromSource(sphericalSource, Nu[i]);
		Is[i] = synchrotronEvaluator.evaluateSynchrotronFluxFromSource(angleDependentSource, Nu[i]);
	}

	printLog("output\n");

	FILE* output_ev_EFE = fopen("outputE.dat", "w");
	FILE* output_GHz_Jansky = fopen("outputNu.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i] / (1.6E-12), E[i] * E[i] * I[i]);
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * E[i] * I[i]);
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	FILE* output_synchr = fopen("outputSynch.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		fprintf(output_synchr, "%g %g\n", Nu[i]/1E9, Is[i]*1E26);
	}
	fclose(output_synchr);

	delete electrons;
	delete source;
	delete[] E;
	delete[] I;

	return 0;
}