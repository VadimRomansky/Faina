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

	srand(10);

	PhotonPlankDistribution* CMBradiation = PhotonPlankDistribution::getCMBradiation();
	PhotonMultiPlankDistribution* galacticRadiation = PhotonMultiPlankDistribution::getGalacticField();
	ElectronPowerLawDistribution* electrons = new ElectronPowerLawDistribution(3.5, Emin, electronConcentration);
	RadiationSource* source = new SimpleFlatSource(electrons, B, sinTheta, electronConcentration, rmax, rmax, distance);
	RadiationSource* sphericalSource = new TabulatedSphericalLayerSource(20, 20, 4, electrons, B, 1.0, electronConcentration, rmax, 0.5 * rmax, distance);
	double volume = source->getTotalVolume();

	InverseComptonEvaluator* comptonEvaluator = new InverseComptonEvaluator(200, 20, 20, Emin, Emax);
	SynchrotronEvaluator* synchrotronEvaluator = new SynchrotronEvaluator(200, Emin, Emax);

	int Ndistributions = 10;
	ElectronIsotropicDistribution** angleDependentDistributions = ElectronDistributionFactory::readTabulatedIsotropicDistributions("./input/Ee", "./input/Fs",".dat",10, ElectronInputType::GAMMA_KIN_FGAMMA, electronConcentration, 200);
	for (int i = 0; i < Ndistributions; ++i) {
		(dynamic_cast<ElectronTabulatedIsotropicDistribution*>(angleDependentDistributions[i]))->rescaleDistribution(sqrt(18));
	}
	AngleDependentElectronsSphericalSource* angleDependentSource = new AngleDependentElectronsSphericalSource(20, 20, 4, Ndistributions, angleDependentDistributions, B, 1.0, 0, electronConcentration, rmax, 0.5 * rmax, distance);

	const int Nparams = 4;
	double minParameters[Nparams] = {1E16, 0.01, 0.01, 0.1};
	double maxParameters[Nparams] = {2E17, 10, 1000, 1.0};
	double vector[Nparams] = { 1.3E17, 0.6, 150, 0.5 };
	for (int i = 0; i < Nparams; ++i) {
		vector[i] = vector[i] / maxParameters[i];
	}
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
	const int Nnu1 = 4;
	double nu1[Nnu1] = { 1.5E9, 3.0E9, 6.1E9, 0.97E9 };
	double observedFlux[Nnu1] = { 1.5, 4.3, 6.1, 4.2 };
	double observedError[Nnu1] = { 0.1, 0.2, 0.3, 0.2 };
	bool optPar[Nnu1] = { true, true, true, true };
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
	SynchrotronOptimizer* synchrotronOptimizer = new GradientDescentSynchrotronOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, 20, ErrorScale::LINEAR);
	int Npoints[Nparams] = { 5,5,5,5 };
	SynchrotronOptimizer* enumOptimizer = new EnumSynchrotronOptimizer(synchrotronEvaluator, minParameters, maxParameters, Nparams, ErrorScale::LINEAR, Npoints);
	enumOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, angleDependentSource);
	synchrotronOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, angleDependentSource);
	//synchrotronOptimizer->optimize(vector, optPar, nu1, observedFlux, observedError, Nnu1, source);
	//angleDependentSource->resetParameters(vector, maxParameters);
	source->resetParameters(vector, maxParameters);

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
		//I[i] = comptonEvaluator->evaluateComptonIsotropicFluxFromSource(E[i], CMBradiation, source);
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
		Is[i] = synchrotronEvaluator->evaluateSynchrotronFluxFromSource(angleDependentSource, Nu[i]);
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