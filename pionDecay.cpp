#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "pionDecay.h"

void PionDecayEvaluator::getCoefs(double& alpha, double& beta, double& gamma, double& lambda, const double& protonEnergy)
{
	if (protonEnergy < 1E9 * 1.6E-12) {
		lambda = 3.0;
		alpha = 1.0;
		beta = 3.29 - 0.2 * pow(protonEnergy / (massProton * speed_of_light2), -1.5);
		gamma = 0;
		return;
	}
	if (protonEnergy < 4E9 * 1.6E-12) {
		lambda = 3.0;
		alpha = 1.0;
		double q = (protonEnergy - 1E9 * 1.6E-12) / (massProton * speed_of_light2);
		double mu = 1.25 * pow(q, 1.25) * exp(-1.25 * q);
		beta = mu + 2.45;
		gamma = mu + 1.45;
		return;
	}
	if (protonEnergy < 20E9 * 1.6E-12) {
		lambda = 3.0;
		alpha = 1.0;
		double q = (protonEnergy - 1E9 * 1.6E-12) / (massProton * speed_of_light2);
		double mu = 1.25 * pow(q, 1.25) * exp(-1.25 * q);
		beta = 1.5*mu + 4.95;
		gamma = mu + 1.5;
		return;
	}
	if (protonEnergy < 50E9 * 1.6E-12) {
		lambda = 3.0;
		alpha = 0.5;
		beta = 4.2;
		gamma = 1;
		return;
	}
	if (protonEnergy < 100E9 * 1.6E-12) {
		//geant
		lambda = 3.0;
		alpha = 0.5;
		beta = 4.2;
		gamma = 1;
		return;
		//pythia
		/*lambda = 3.5;
		alpha = 0.5;
		beta = 4.0;
		gamma = 1;
		return;*/
	}
	//geant
	lambda = 3.0;
	alpha = 0.5;
	beta = 4.9;
	gamma = 1;
	//pythia
	/*lambda = 3.5;
	alpha = 0.5;
	beta = 4.0;
	gamma = 1;
	return; */
	//SIBYLL
	/*lambda = 3.55;
	alpha = 0.5;
	beta = 3.6;
	gamma = 1;*/
	//QGSJET
	/*lambda = 3.55;
	alpha = 0.5;
	beta = 4.5;
	gamma = 1;*/

}

void PionDecayEvaluator::getBcoefs(double& b1, double& b2, double& b3, const double& protonEnergy)
{
	if (protonEnergy < 5E9 * 1.6E-12) {
		b1 = 9.53;
		b2 = 0.52;
		b3 = 0.054;
		return;
	}
	if (protonEnergy < 50E9 * 1.6E-12) {
		b1 = 9.13;
		b2 = 0.35;
		b3 = 0.0097;
		return;
	}
	if (protonEnergy < 100E9 * 1.6E-12) {
		//geant
		/*b1 = 9.13;
		b2 = 0.35;
		b3 = 0.0097;*/
		//pythia
		b1 = 9.06;
		b2 = 0.3795;
		b3 = 0.01105;
		return;
	}
	//geant
	/*b1 = 9.13;
	b2 = 0.35;
	b3 = 0.0097;*/
	//pythia
	/*b1 = 9.06;
	b2 = 0.3795;
	b3 = 0.01105;*/
	//SIBYLL
	b1 = 10.77;
	b2 = 0.412;
	b3 = 0.01264;
	//QGSJET
	/*b1 = 13.16;
	b2 = 0.4419;
	b3 = 0.01439;*/
}

PionDecayEvaluator::PionDecayEvaluator(int Ne, double Emin, double Emax)
{
	my_Ne = Ne;
	my_Emin = Emin;
	my_Emax = Emax;

	double factor = pow(my_Emax / my_Emin, 1.0 / (my_Ne - 1));

	my_Ee = new double[my_Ne];
	my_Ee[0] = my_Emin;
	for (int i = 1; i < my_Ne; ++i) {
		my_Ee[i] = my_Ee[i - 1] * factor;
	}
}

PionDecayEvaluator::~PionDecayEvaluator()
{
	delete[] my_Ee;
}

double PionDecayEvaluator::sigmaInelastic(const double& energy)
{
	double thresholdEnergy = (2 * massPi0 + massPi0 * massPi0 / (2 * massProton)) * speed_of_light2;
	if (energy <= thresholdEnergy) {
		return 0;
	}
	double energyRatio = energy / thresholdEnergy;
	double logEnergy = log(energyRatio);
	double result = (30.7 - 0.9 * logEnergy + 0.18 * sqr(logEnergy)) * cube(1 - pow(energyRatio, -1.9)) * 1E-27;
	return result;
}

double PionDecayEvaluator::sigmaPion(const double& energy)
{
	double thresholdEnergy = (2 * massPi0 + massPi0 * massPi0 / (2 * massProton)) * speed_of_light2;
	if (energy <= thresholdEnergy) {
		return 0;
	}
	//depends on model
	double a1, a2, a3, a4, a5;
	//Geant 4
	double EtranGeant = (5E9) * 1.6E-12;
	a1 = 0.728;
	a2 = 0.596;
	a3 = 0.491;
	a4 = 0.2503;
	a5 = 0.117;
	//Pythia 8
	double EtranPythia = (50E9) * 1.6E-12;
	a1 = 0.652;
	a2 = 0.0016;
	a3 = 0.488;
	a4 = 0.1928;
	a5 = 0.482;
	//SIBYLL
	double EtranSIBYLL = (100E9) * 1.6E-12;
	a1 = 5.436;
	a2 = 0.254;
	a3 = 0.072;
	a4 = 0.075;
	a5 = 0.166;
	//QGSJET
	double EtranQGSJET = (100E9) * 1.6E-12;
	a1 = 0.908;
	a2 = 0.0009;
	a3 = 6.089;
	a4 = 0.176;
	a5 = 0.448;
	if (energy < thresholdEnergy) {
		return 0.0;
	}
	//from threshold to 2 GeV
	if (energy < (2E9) * 1.6E-12) {
		double Mres = 1.1883E9 * 1.6E-12;
		double Gres = 0.2264E9 * 1.6E-12;
		double sigma0 = 7.66E-30;
		double s = 2 * massProton * speed_of_light2 * (energy + 2 * massProton * speed_of_light2);
		double g = sqrt(Mres*Mres*(Mres*Mres + Gres*Gres));
		double K = sqrt(8) * Mres * Gres * g / (pi * sqrt(Mres * Mres + g));
		double fBW = massProton * K / (sqr(sqr(sqrt(s) - massProton * speed_of_light2) - Mres * Mres) + Mres * Mres * Gres * Gres);
		double eta = sqrt(sqr((s / speed_of_light4) - massPi0 * massPi0 - 4 * massProton * massProton) - 16 * massPi0 * massPi0 * massProton * massProton) / (2 * massPi0 * sqrt(s / speed_of_light4));
		double result = sigma0 * pow(eta, 1.95) * (1 + eta + eta * eta * eta * eta * eta) * pow(fBW, 1.86);
		return result + sigma2Pion(energy);
	}
	if (energy < (5E9) * 1.6E-12) {
		double Q = (energy - thresholdEnergy) / (massProton * speed_of_light2);
		double nmean = -0.006 + 0.237 * Q - 0.023 * Q * Q;
		return sigmaInelastic(energy) * nmean;
	}
	//from 5GeV to EtranGeant
	if (energy < EtranPythia) {
		//geant
		a1 = 0.728;
		a2 = 0.596;
		a3 = 0.491;
		a4 = 0.2503;
		a5 = 0.117;
		double dzeta = (energy - (3E9) * 1.6E-12) / (massProton * speed_of_light2);
		double nmean = a1*pow(dzeta, a4)*(1+exp(-a2*pow(dzeta, a5)))*(1-exp(-a3*pow(dzeta, 0.25)));
		return sigmaInelastic(energy) * nmean;
	}
	if (energy < EtranSIBYLL) {
		//pythia
		a1 = 0.652;
		a2 = 0.0016;
		a3 = 0.488;
		a4 = 0.1928;
		a5 = 0.482;
		double dzeta = (energy - (3E9) * 1.6E-12) / (massProton * speed_of_light2);
		double nmean = a1 * pow(dzeta, a4) * (1 + exp(-a2 * pow(dzeta, a5))) * (1 - exp(-a3 * pow(dzeta, 0.25)));
		return sigmaInelastic(energy) * nmean;
	}
	else {
		//SIBYLL
		a1 = 5.436;
		a2 = 0.254;
		a3 = 0.072;
		a4 = 0.075;
		a5 = 0.166;
		//QGSJET
		/*a1 = 0.908;
		a2 = 0.0009;
		a3 = 6.089;
		a4 = 0.176;
		a5 = 0.448;*/

		double dzeta = (energy - (3E9) * 1.6E-12) / (massProton * speed_of_light2);
		double nmean = a1 * pow(dzeta, a4) * (1 + exp(-a2 * pow(dzeta, a5))) * (1 - exp(-a3 * pow(dzeta, 0.25)));
		return sigmaInelastic(energy) * nmean;
	}
}

double PionDecayEvaluator::sigma2Pion(const double& energy)
{
	if (energy < (0.56E9) * 1.6E-12) {
		return 0.0;
	}
	if (energy < (2E9) * 1.6E-12) {
		return (5.7E-27) / (1 + exp(-9.3*((energy/(1.6E-3) - 1.4))));
	}
	return 0;
}

double PionDecayEvaluator::sigmaGamma(const double& photonEnergy, const double& protonEnergy)
{
	double thresholdEnergy = (2 * massPi0 + massPi0 * massPi0 / (2 * massProton)) * speed_of_light2;
	if (protonEnergy <= thresholdEnergy) {
		return 0;
	}

	double s = 2 * massProton * speed_of_light2 * (protonEnergy + 2 * massProton * speed_of_light2);
	double EpiCM = (s - 4 * massProton * massProton * speed_of_light4 + massPi0 * massPi0 * speed_of_light4) / (2 * sqrt(s));
	double gCM = (protonEnergy + 2 * massProton * speed_of_light2) / sqrt(s);
	double betaCM = sqrt(1 - 1.0 / (gCM * gCM));
	double PpiCM = sqrt(EpiCM * EpiCM - massPi0 * massPi0 * speed_of_light4);
	double EpiLabMax = gCM * (EpiCM + PpiCM * betaCM);
	double gLab = EpiLabMax / (massPi0 * speed_of_light2);
	double betaLab = sqrt(1 - 1.0 / (gLab * gLab));
	double EphMax = 0.5 * massPi0 * speed_of_light2 * gLab * (1 + betaLab);

	if (photonEnergy >= EphMax) {
		return 0;
	}

	double Yph = photonEnergy + massPi0 * massPi0 * speed_of_light4 / (4 * photonEnergy);
	double YphMax = EphMax + massPi0 * massPi0 * speed_of_light4 / (4 * EphMax);
	double Xph = (Yph - massPi0 * speed_of_light2) / (YphMax - massPi0 * speed_of_light2);

	double alpha, beta, gamma, lambda;

	getCoefs(alpha, beta, gamma, lambda, protonEnergy);

	double C = lambda * massPi0 * speed_of_light2 / YphMax;
	double F = pow(1 - pow(Xph, alpha), beta) / pow(1 + Xph / C, gamma);

	double Amax;
	double b0 = 5.9;
	if (protonEnergy < 1E9 * 1.6E-12) {
		Amax = b0 * sigmaPion(protonEnergy) / EpiLabMax;
	}
	else {
		double b1, b2, b3;
		double theta = protonEnergy / (massProton * speed_of_light2);
		getBcoefs(b1, b2, b3, protonEnergy);
		Amax = b1 * pow(theta, -b2) * exp(b3 * sqr(log(theta)) * sigmaPion(protonEnergy) / (massProton * speed_of_light2));
	}
	return Amax * F;
}

double PionDecayEvaluator::evaluatePionDecayLuminocityIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& ambientConcentration, const double& volume, const double& distance)
{
	double result = 0;
	double Yph = photonFinalEnergy + massPi0 * massPi0 * speed_of_light4 / (4 * photonFinalEnergy);
	
	for (int i = 0; i < my_Ne; ++i) {
		double protonEnergy = my_Ee[i];
		double protonGamma = protonEnergy / (massProton * speed_of_light2);
		double protonBeta = sqrt(1.0 - 1.0 / (protonGamma * protonGamma));
		double protonKineticEnergy = massProton * speed_of_light2 * (protonGamma - 1.0);

		double sigma = sigmaGamma(photonFinalEnergy, protonKineticEnergy);

		result += (speed_of_light * protonBeta/4*pi) * sigma * protonDistribution->distribution(protonEnergy) * ambientConcentration * volume / sqr(distance);
	}

	return result;
}

double PionDecayEvaluator::evaluatePionDecayIsotropicFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
				double ambientConcentration = source->getConcentration(irho, iz, iphi);
				result += evaluatePionDecayLuminocityIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), ambientConcentration, source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	return result;
}