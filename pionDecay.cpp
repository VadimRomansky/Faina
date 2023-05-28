#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "pionDecay.h"

PionDecayEvaluatorBase::PionDecayEvaluatorBase(int Ne, double Emin, double Emax, const double& ambientConcentration) : RadiationEvaluator(Ne, Emin, Emax){
    my_ambientConcentration = ambientConcentration;
}

PionDecayEvaluatorBase::~PionDecayEvaluatorBase(){

}

double PionDecayEvaluatorBase::sigmaInelastic(const double& energy)
{
    double thresholdEnergy = (2 * massPi0 + massPi0 * massPi0 / (2 * massProton)) * speed_of_light2;
    if (energy <= thresholdEnergy) {
        return 0;
    }
    double energyRatio = energy / thresholdEnergy;
    double logEnergy = log(energyRatio);
    double result = (30.7 - 0.96 * logEnergy + 0.18 * sqr(logEnergy)) * cube(1 - pow(energyRatio, -1.9)) * 1E-27;
    return result;
}



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
		b1 = 9.13;
		b2 = 0.35;
		b3 = 0.0097;
		//pythia
		/*b1 = 9.06;
		b2 = 0.3795;
		b3 = 0.01105;*/
		return;
	}
	//geant
	b1 = 9.13;
	b2 = 0.35;
	b3 = 0.0097;
	//pythia
	/*b1 = 9.06;
	b2 = 0.3795;
	b3 = 0.01105;*/
	//SIBYLL
	/*b1 = 10.77;
	b2 = 0.412;
	b3 = 0.01264;*/
	//QGSJET
	/*b1 = 13.16;
	b2 = 0.4419;
	b3 = 0.01439;*/
}

PionDecayEvaluator::PionDecayEvaluator(int Ne, double Emin, double Emax, const double& ambientConcentration) : PionDecayEvaluatorBase(Ne, Emin, Emax, ambientConcentration)
{
}

PionDecayEvaluator::~PionDecayEvaluator()
{
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
		double s = 2 * massProton * speed_of_light2 * (energy + 2 * massProton * speed_of_light2);//E2
		double g = sqrt(Mres*Mres*(Mres*Mres + Gres*Gres)); //E2
		double K = sqrt(8) * Mres * Gres * g / (pi * sqrt(Mres * Mres + g));//E3
		double fBW = massProton * speed_of_light2 * K / (sqr(sqr(sqrt(s) - massProton * speed_of_light2) - Mres * Mres) + Mres * Mres * Gres * Gres);//dimensionless
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
	//if (energy < EtranPythia) {
		//geant
		a1 = 0.728;
		a2 = 0.596;
		a3 = 0.491;
		a4 = 0.2503;
		a5 = 0.117;
		double dzeta = (energy - (3E9) * 1.6E-12) / (massProton * speed_of_light2);
		double nmean = a1*pow(dzeta, a4)*(1+exp(-a2*pow(dzeta, a5)))*(1-exp(-a3*pow(dzeta, 0.25)));
		return sigmaInelastic(energy) * nmean;
	//}
	if (energy < EtranSIBYLL) {
		//pythia
		a1 = 0.652;
		a2 = 0.0016;
		a3 = 0.488;
		a4 = 0.1928;
		a5 = 0.483;
		double dzeta = (energy - (3E9) * 1.6E-12) / (massProton * speed_of_light2);
		double nmean = a1 * pow(dzeta, a4) * (1 + exp(-a2 * pow(dzeta, a5))) * (1 - exp(-a3 * pow(dzeta, 0.25)));
		return sigmaInelastic(energy) * nmean;
	}
	else {
		//SIBYLL
		/*a1 = 5.436;
		a2 = 0.254;
		a3 = 0.072;
		a4 = 0.075;
		a5 = 0.166;*/
		//QGSJET
	    a1 = 0.908;
		a2 = 0.0009;
		a3 = 6.089;
		a4 = 0.176;
		a5 = 0.448;

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
		return (5.7E-27) / (1 + exp(-9.3*((energy/(1.6E-3)) - 1.4)));
	}
	return 0;
}

double PionDecayEvaluator::sigmaGamma(const double& photonEnergy, const double& protonEnergy)
{
	double thresholdEnergy = (2 * massPi0 + massPi0 * massPi0 / (2 * massProton)) * speed_of_light2;
	if (protonEnergy <= thresholdEnergy) {
		return 0;
	}

	double s = 2 * massProton * speed_of_light2 * (protonEnergy + 2 * massProton * speed_of_light2);//E2
	double EpiCM = (s - 4 * massProton * massProton * speed_of_light4 + massPi0 * massPi0 * speed_of_light4) / (2 * sqrt(s));//E
	double gCM = (protonEnergy + 2 * massProton * speed_of_light2) / sqrt(s);
	double betaCM = sqrt(1 - 1.0 / (gCM * gCM));
	double PpiCM = sqrt(EpiCM * EpiCM - massPi0 * massPi0 * speed_of_light4);
	double EpiLabMax = gCM * (EpiCM + PpiCM * betaCM);
	double gLab = EpiLabMax / (massPi0 * speed_of_light2);
	double betaLab = sqrt(1 - 1.0 / (gLab * gLab));
	double EphMax = 0.5 * massPi0 * speed_of_light2 * gLab * (1 + betaLab);

	double EphMin = 0.5 * massPi0 * speed_of_light2 * gLab * (1 - betaLab);

	if (photonEnergy >= EphMax) {
		return 0;
	}

	//todo what if energy less?
	/*if (photonEnergy <= 0.5 * massPi0 * speed_of_light2) {
		return 0;
	}*/
	double Yph = photonEnergy + massPi0 * massPi0 * speed_of_light4 / (4 * photonEnergy);
	double YphMax = EphMax + massPi0 * massPi0 * speed_of_light4 / (4 * EphMax);
	double Yphmin = EphMin + massPi0 * massPi0 * speed_of_light4 / (4 * EphMin);
	double Xph = (Yph - massPi0 * speed_of_light2) / (YphMax - massPi0 * speed_of_light2);
	if (Xph >= 1.0 || Xph <= 0) {
		return 0;
	}
	double alpha, beta, gamma, lambda;

	getCoefs(alpha, beta, gamma, lambda, protonEnergy);

	double C = lambda * massPi0 * speed_of_light2 / YphMax;
	double F = pow(1 - pow(Xph, alpha), beta) / pow(1 + Xph / C, gamma);

	double Amax;
	double b0 = 5.9;
	if (protonEnergy < (1E9) * 1.6E-12) {
		Amax = b0 * sigmaPion(protonEnergy) / EpiLabMax;
	}
	else {
		double b1, b2, b3;
		double theta = protonEnergy / (massProton * speed_of_light2);
		getBcoefs(b1, b2, b3, protonEnergy);
		Amax = b1 * pow(theta, -b2) * exp(b3 * sqr(log(theta))) * sigmaPion(protonEnergy) / (massProton * speed_of_light2);
	}
	return Amax * F;
}

void PionDecayEvaluator::resetParameters(const double *parameters, const double *normalizationUnits){
    //todo change ambient density?
}

double PionDecayEvaluator::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance)
{
	double result = 0;
	
	for (int i = 0; i < my_Ne; ++i) {
		double protonEnergy = my_Ee[i];
		double dprotonEnergy;
		if (i == 0) {
			dprotonEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			dprotonEnergy = my_Ee[i] - my_Ee[i - 1];
		}
		double protonGamma = protonEnergy / (massProton * speed_of_light2);
		double protonBeta = sqrt(1.0 - 1.0 / (protonGamma * protonGamma));
		double protonKineticEnergy = massProton * speed_of_light2 * (protonGamma - 1.0);

		double sigma = sigmaGamma(photonFinalEnergy, protonKineticEnergy);

        result += photonFinalEnergy*(speed_of_light * protonBeta/(4*pi)) * sigma * protonDistribution->distribution(protonEnergy) * my_ambientConcentration * volume * dprotonEnergy / sqr(distance);

		if (result != result) {
			printf("result = NaN in pion decay\n");
			printLog("result = NaN in pion decay\n");
			exit(0);
		}
	}

	return result;
}

double PionDecayEvaluator::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;
	int irho = 0;

	omp_init_lock(&my_lock);

#pragma omp parallel for private(irho) shared(photonFinalEnergy, source, Nrho, Nz, Nphi) reduction(+:result)

	for (irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
                result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	omp_destroy_lock(&my_lock);

	return result;
}

PionDecayEvaluatorKelner::PionDecayEvaluatorKelner(int Ne, double Emin, double Emax, const double& ambientConcentration) : PionDecayEvaluatorBase(Ne, Emin, Emax, ambientConcentration){

}

PionDecayEvaluatorKelner::~PionDecayEvaluatorKelner(){

}

void PionDecayEvaluatorKelner::resetParameters(const double *parameters, const double *normalizationUnits){
    //todo reset ambient concentration
}

double PionDecayEvaluatorKelner::functionKelner(const double& x, const double& protonEnergy)
{
    if (x >= 1.0) {
        return 0;
    }
    double L = log(protonEnergy / ((1E12) * 1.6E-12));
    double B = 1.30 + 0.14 * L + 0.011 * L * L;
    double beta = 1.0 / (1.79 + 0.11 * L + 0.008 * L * L);
    double k = 1.0 / (0.801 + 0.049 * L + 0.014 * L * L);
    double xbeta = pow(x, beta);
    double F = B * (log(x) / x) * pow((1 - xbeta) / (1 + k * xbeta * (1 - xbeta)), 4) * ((1 / log(x)) - 4 * beta * xbeta / (1 - xbeta) - 4 * k * beta * xbeta * (1 - 2 * xbeta) / (1 + k * xbeta * (1 - xbeta)));
    return F;
}

double PionDecayEvaluatorKelner::evaluateFluxFromIsotropicFunction(const double& photonFinalEnergy, MassiveParticleIsotropicDistribution* protonDistribution, const double& volume, const double& distance)
{
	double result = 0;

	for (int i = 0; i < my_Ne; ++i) {
		double protonEnergy = my_Ee[i];
		double dprotonEnergy;
		if (i == 0) {
			dprotonEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			dprotonEnergy = my_Ee[i] - my_Ee[i - 1];
		}
		double protonGamma = protonEnergy / (massProton * speed_of_light2);
		double protonBeta = sqrt(1.0 - 1.0 / (protonGamma * protonGamma));
		double protonKineticEnergy = massProton * speed_of_light2 * (protonGamma - 1.0);

		double sigma = (sigmaInelastic(protonKineticEnergy)/protonEnergy)*functionKelner(photonFinalEnergy/protonEnergy, protonEnergy);

        result += photonFinalEnergy*(speed_of_light * protonBeta / 4 * pi) * sigma * protonDistribution->distribution(protonEnergy) * my_ambientConcentration * volume * dprotonEnergy / sqr(distance);

		if (result != result) {
			printf("result = NaN in pion decay\n");
			printLog("result = NaN in pion decay\n");
			exit(0);
		}
	}

	return result;
}

double PionDecayEvaluatorKelner::evaluateFluxFromSource(const double& photonFinalEnergy, RadiationSource* source)
{
	int Nrho = source->getNrho();
	int Nz = source->getNz();
	int Nphi = source->getNphi();

	double result = 0;

	for (int irho = 0; irho < Nrho; ++irho) {
		for (int iz = 0; iz < Nz; ++iz) {
			for (int iphi = 0; iphi < Nphi; ++iphi) {
                result += evaluateFluxFromIsotropicFunction(photonFinalEnergy, source->getParticleDistribution(irho, iz, iphi), source->getVolume(irho, iz, iphi), source->getDistance());
			}
		}
	}

	return result;
}
