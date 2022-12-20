#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "electronDistribution.h"
#include "photonDistribution.h"
#include "util.h"
#include "inverseCompton.h"


int main() {
	PhotonPlankDistribution* CMBradiation = PhotonPlankDistribution::getCMBradiation();
	PhotonMultiPlankDistribution* galacticRadiation = PhotonMultiPlankDistribution::getGalacticField();
	ElectronDistribution* electrons = new ElectronPowerLawDistribution(3.5, me_c2, 150);

	InverseComptonEvaluator evaluator = InverseComptonEvaluator(200, 20, 20, me_c2, 10000 * me_c2);

	int Nnu = 200;
	double* E = new double[Nnu];
	double* I = new double[Nnu];

	double Emin = 0.0001 * kBoltzman * 2.7;
	double Emax = 2*10000*me_c2 + Emin;
	double factor = pow(Emax / Emin, 1.0 / (Nnu - 1));
	E[0] = Emin;
	I[0] = 0;
	for (int i = 1; i < Nnu; ++i) {
		E[i] = E[i - 1] * factor;
		I[i] = 0;
	}

	double rmax = 1.3E17;
	double volume = 4 * pi * rmax * rmax * rmax / 3;

	//SN2009bb
	//const double distance = 40*3.08*1.0E24;
	//AT2018
	//const double distance = 60*3.08*1.0E24;
	//CSS161010
	const double distance = 150 * 3.08 * 1.0E24;

	printLog("evaluating\n");
	for (int i = 0; i < Nnu; ++i) {
		printf("%d\n", i);
		//I[i] = evaluateComptonLuminocity(E[i], 0, 0, CMBradiation, electrons, volume, distance, me_c2, 10000 * me_c2, 200, 20, 20);
		I[i] = evaluator.evaluateComptonLuminocity(E[i], 0, 0, CMBradiation, electrons, volume, distance);
	}

	printLog("output\n");

	FILE* output_ev_EFE = fopen("outputE.dat", "w");
	FILE* output_GHz_Jansky = fopen("outputNu.dat", "w");
	for (int i = 0; i < Nnu; ++i) {
		double nu = E[i] / hplank;
		fprintf(output_ev_EFE, "%g %g\n", E[i] / (1.6E-12), E[i] * E[i] * I[i] / sqr(distance));
		fprintf(output_GHz_Jansky, "%g %g\n", nu / 1E9, 1E26 * hplank * E[i] * I[i] / sqr(distance));
	}
	fclose(output_ev_EFE);
	fclose(output_GHz_Jansky);

	delete electrons;
	delete[] E;
	delete[] I;

	return 0;
}