#include "stdio.h"
#include "math.h"
#include <omp.h>

#include "constants.h"
#include "util.h"
#include "coordinateTransform.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"

#include "bremsstrahlung.h"

//todo https://en.wikipedia.org/wiki/Kramers%27_opacity_law

BremsstrahlungThermalEvaluator::BremsstrahlungThermalEvaluator(bool absorption, bool doppler) : RadiationEvaluator(2, massElectron*speed_of_light2, 2*massElectron*speed_of_light2, absorption, doppler)
{
}

BremsstrahlungThermalEvaluator::~BremsstrahlungThermalEvaluator()
{
}

void BremsstrahlungThermalEvaluator::resetParameters(const double* parameters, const double* normalizationUnits)
{
}

double BremsstrahlungThermalEvaluator::evaluateGauntFactor(double eta, double theta)
{
	double gauntFactor = 1.0;

	if (eta >= 1.0) {
		if (theta >= 1.0) {
			gauntFactor = sqrt(3 / (pi * theta));
		}
		else {
			gauntFactor = (sqrt(3) / pi) * log(4 / (euler_mascheroni * theta));
		}
	}
	else {
		if (theta >= 1.0) {
			if (theta >= 1 / eta) {
				double dzeta = eta * theta;
				gauntFactor = sqrt(12 / dzeta);
			}
			else {
				gauntFactor = 1.0;
			}
		}
		else {
			if (theta < sqrt(eta)) {
				//where is 4? if in denominator, log can be < 0
				gauntFactor = (sqrt(3) / pi) * log(4 / (pow(euler_mascheroni, 2.5) * theta) * sqrt(eta));
			}
			else {
				gauntFactor = 1.0;
			}
		}
	}

	return gauntFactor;
}

double BremsstrahlungThermalEvaluator::evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	ThermalRectangularSource* thermalSource = dynamic_cast<ThermalRectangularSource*>(source);
	//MassiveParticleDistribution* electronDistribution = source->getParticleDistribution(ix1, iz, ix2);
	//MassiveParticleMaxwellDistribution* maxwellDistribution = dynamic_cast<MassiveParticleMaxwellDistribution*>(electronDistribution);
	//MassiveParticleMaxwellJuttnerDistribution* maxwellJuttnerDistribution = dynamic_cast<MassiveParticleMaxwellJuttnerDistribution*>(electronDistribution);
	double temperature;
	//if (maxwellDistribution != NULL) {
	//	temperature = maxwellDistribution->getTemperature();
	//}
	//else if (maxwellJuttnerDistribution != NULL) {
	//	temperature = maxwellJuttnerDistribution->getTemperature();
	//}
	//else {
	//	printf("primitive bremsstrahlung evaluator works only with maxwellian distribution\n");
	//	printLog("primitive bremsstrahlung evaluator works only with maxwellian distribution\n");
	//	exit(0);
	//}
	if (thermalSource == NULL) {
		printf("primitive bremsstrahlung evaluator works only with thermal sources\n");
		printLog("primitive bremsstrahlung evaluator works only with thermal sources\n");
		exit(0);
	}
	temperature = thermalSource->getTemperature(ix1, iz, ix2);
	double concentration = source->getConcentration(ix1, iz, ix2);
	double m = thermalSource->getParticleMass();
	double theta = photonFinalEnergy / (kBoltzman * temperature);
	double ritberg = m * electron_charge * electron_charge * electron_charge * electron_charge * 2 * pi * pi / (hplank * hplank);
	double eta = kBoltzman * temperature / ritberg;

	double gauntFactor = evaluateGauntFactor(eta, theta);

	double a = (32 * pi * pow(electron_charge, 6) / (3 * m * speed_of_light * speed_of_light2)) * sqrt(2 * pi / (3 * kBoltzman * m));
	double result = (1.0 / hplank) * (a / sqrt(temperature)) * concentration * concentration * exp(-theta) * gauntFactor;
	result *= 1 / (4 * pi);
	return result;
}

double BremsstrahlungThermalEvaluator::evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	ThermalRectangularSource* thermalSource = dynamic_cast<ThermalRectangularSource*>(source);
	//MassiveParticleDistribution* electronDistribution = source->getParticleDistribution(ix1, iz, ix2);
	//MassiveParticleMaxwellDistribution* maxwellDistribution = dynamic_cast<MassiveParticleMaxwellDistribution*>(electronDistribution);
	//MassiveParticleMaxwellJuttnerDistribution* maxwellJuttnerDistribution = dynamic_cast<MassiveParticleMaxwellJuttnerDistribution*>(electronDistribution);
	double temperature;
	//if (maxwellDistribution != NULL) {
	//	temperature = maxwellDistribution->getTemperature();
	//}
	//else if (maxwellJuttnerDistribution != NULL) {
	//	temperature = maxwellJuttnerDistribution->getTemperature();
	//}
	//else {
	//	printf("primitive bremsstrahlung evaluator works only with maxwellian distribution\n");
	//	printLog("primitive bremsstrahlung evaluator works only with maxwellian distribution\n");
	//	exit(0);
	//}
	if (thermalSource == NULL) {
		printf("primitive bremsstrahlung evaluator works only with thermal sources\n");
		printLog("primitive bremsstrahlung evaluator works only with thermal sources\n");
		exit(0);
	}

	temperature = thermalSource->getTemperature(ix1, iz, ix2);
	double concentration = source->getConcentration(ix1, iz, ix2);
	double m = thermalSource->getParticleMass();
	double theta = photonFinalEnergy / (kBoltzman * temperature);
	double ritberg = m * electron_charge * electron_charge * electron_charge * electron_charge * 2 * pi * pi / (hplank * hplank);
	double eta = kBoltzman * temperature / ritberg;

	double gauntFactor = evaluateGauntFactor(eta, theta);
	
	double photonFinalFrequency = photonFinalEnergy / hplank;

	double a = (4 * pow(electron_charge, 6.0) / (3 * m * hplank * speed_of_light)) * pow(2 * pi / (3 * kBoltzman * m), 0.5) * concentration * concentration* gauntFactor * (1.0 - exp(-photonFinalEnergy/(kBoltzman*temperature))) / (pow(temperature, 0.5) * cube(photonFinalFrequency));

	return a;
}

double BremsstrahlungEvaluator::evaluateSigma1(const double& gammaE, const double& epsilonG)
{
	if (epsilonG >= gammaE - 1.0) {
		return 0;
	}
	double sigma = (4*re2*alpha/epsilonG)*(1.0 + (1.0/3.0 - epsilonG/(gammaE))*(1.0 - epsilonG / (gammaE)))*(log(2*(gammaE)*(gammaE-epsilonG)/epsilonG)-0.5);

	if (sigma != sigma) {
		printf("sigma1 = NaN in bremsstrahlung evaluator\n");
		printLog("sigma1 = NaN in bremsstrahlung evaluator\n");
		exit(0);
	}

	if (sigma < 0) {
		printf("sigma1 < 0 in bremsstrahlung evaluator\n");
		printLog("sigma1 < 0 in bremsstrahlung evaluator\n");
		exit(0);
	}

	return sigma;
}

double BremsstrahlungEvaluator::evaluateSigma2(const double& gammaE, const double& epsilonG)
{
	if (epsilonG >= gammaE - 1) {
		return 0;
	}
	double coef = re2 * alpha / (3 * epsilonG);
	if (epsilonG < 0.5) {
		double sigma = coef * (16*(1.0-epsilonG+epsilonG*epsilonG)*log(gammaE/epsilonG) - 1.0/sqr(epsilonG) + 3.0/epsilonG - 4.0 + 4.0*epsilonG - 8.0*epsilonG*epsilonG - 2*(1 - 2*epsilonG)*log(1.0 - 2.0*epsilonG)*(1.0/(4.0*cube(epsilonG))-1/(2.0*sqr(epsilonG)+3.0/epsilonG-2.0+4.0*epsilonG)));
		if (sigma != sigma) {
			printf("sigma2 = NaN in bremsstrahlung evaluator\n");
			printLog("sigma2 = NaN in bremsstrahlung evaluator\n");
			exit(0);
		}

		if (sigma < 0) {
			printf("sigma2 < 0 in bremsstrahlung evaluator\n");
			printLog("sigma2 < 0 in bremsstrahlung evaluator\n");
			exit(0);
		}

		return sigma;
	}
	else {
		double sigma = coef * (2.0/epsilonG)*((4.0 - 1/epsilonG + 1/(4.0*sqr(epsilonG)))*log(2.0*gammaE)-2.0+2/epsilonG-5/(8.0*sqr(epsilonG)));
		if (sigma != sigma) {
			printf("sigma2 = NaN in bremsstrahlung evaluator\n");
			printLog("sigma2 = NaN in bremsstrahlung evaluator\n");
			exit(0);
		}

		if (sigma < 0) {
			printf("sigma2 < 0 in bremsstrahlung evaluator\n");
			printLog("sigma2 < 0 in bremsstrahlung evaluator\n");
			exit(0);
		}

		return sigma;
	}
}

double BremsstrahlungEvaluator::evaluateA(const double& gammaE, const double& epsilonG)
{
	if (epsilonG >= gammaE - 1) {
		return 0;
	}
	double A = 1.0 - (8.0/3.0)*pow(gammaE - 1.0, 1.0/5.0)/(gammaE + 1.0)*pow(epsilonG/gammaE, 1.0/3.0);

	if (A != A) {
		printf("BremsstrahlungEeEvaluator::evaluateA A = NaN gammaE = %g epsilonG = %g\n", gammaE, epsilonG);
		printLog("BremsstrahlungEeEvaluator::evaluateA A = NaN gammaE = %g epsilonG = %g\n", gammaE, epsilonG);
		exit(0);
	}

	if (A < 0) {
		printf("BremsstrahlungEeEvaluator::evaluateA A < 0 gammaE = %g epsilonG = %g\n", gammaE, epsilonG);
		printLog("BremsstrahlungEeEvaluator::evaluateA A < 0 gammaE = %g epsilonG = %g\n", gammaE, epsilonG);
		exit(0);
	}
	return A;
}

double BremsstrahlungEvaluator::evaluateSigmaNR(const double& gammaE, const double& epsilonG) {
	double x = 4 * epsilonG / (gammaE * gammaE - 1.0);
	if (x >= 1) {
		return 0;
	}
	if (x <= 0) {
		return 0;
	}

	double B = evaluateB(gammaE);
	double C = evaluateC(gammaE, x);

	double D = B * (17 - 3 * x * x / sqr(2 - x) - C);
	double E = (12 * (2 - x) - 7 * x * x / (2 - x) - 3 * x * x * x * x / cube(2 - x));

	double sigma = (4*re2 * alpha / (15.0 * epsilonG))*(D*sqrt(1-x)+E*log((1+sqrt(1-x))/sqrt(x)));
	if (sigma != sigma) {
		printf("sigmaNR = NaN in bremsstrahlung evaluator\n");
		printLog("sigmaNR = NaN in bremsstrahlung evaluator\n");
		exit(0);
	}

	//printf("sigma = %g\n", sigma);

	if (sigma < 0) {
		printf("sigmaNR = %g < 0 in bremsstrahlung evaluator\n", sigma);
		printLog("sigmaNR = %g < 0 in bremsstrahlung evaluator\n", sigma);
		//exit(0);
		sigma = 0;
	}

	return sigma;
}

double BremsstrahlungEvaluator::evaluateB(const double& gammaE) {
	return 1 + 0.5 * (gammaE * gammaE - 1.0);
}

double BremsstrahlungEvaluator::evaluateC(const double& gammaE, const double& x) {
	double betaE = sqrt(1.0 - 1.0 / (gammaE * gammaE));
	return 10 * x * gammaE * betaE * (2 + gammaE * betaE) / (1 + x * x * (gammaE * gammaE - 1.0));
}

double BremsstrahlungEvaluator::evaluateSigmaee(const double& gammaE, const double& epsilonG) {
	if (gammaE > 5.0) {
		return (evaluateSigma1(gammaE, epsilonG) + evaluateSigma2(gammaE, epsilonG)) * evaluateA(gammaE, epsilonG);
	}
	else {
		return evaluateSigmaNR(gammaE, epsilonG);
	}

}
double BremsstrahlungEvaluator::evaluateSigmape(const double& gammaE, const double& epsilonG) {
	//for ultrarelativistic limit
	if (gammaE > 1E8) {
		return evaluateSigma1(gammaE, epsilonG);
	}

	//Jauch Rorlich 15-95
	if (epsilonG > gammaE - 1.0) {
		return 0;
	}
	double betaE = sqrt(1.0 - 1.0 / (gammaE * gammaE));
	double gamma1 = gammaE - epsilonG;
	double beta1 = sqrt(1.0 - 1.0 / (gamma1 * gamma1));

	double l = (1.0 / (betaE * gammaE)) * log((1 + betaE) / (1 - betaE));
	double l1 = (1.0 / (beta1 * gamma1)) * log((1 + beta1) / (1 - beta1));
	double L = (2.0 / (gammaE * betaE * gamma1 * beta1)) * log((gammaE * gamma1 * (1 + betaE * beta1) - 1.0) / (epsilonG));

	double coef = (alpha * re2 / epsilonG) * beta1 * gamma1 / (betaE * gammaE);

	double sigma = coef * (4.0/3.0 - (2.0/sqr(betaE*beta1))*(gammaE/gamma1 + gamma1/gammaE - 2/(gammaE*gamma1)) 
		+ (l*gamma1/(gammaE*gammaE - 1.0) + l1*gammaE/(gamma1*gamma1 - 1.0) - l*l1)
		+ L*(8.0*gammaE*gamma1/3.0 + epsilonG*epsilonG*(1.0/sqr(betaE*beta1) + 1.0) + 0.5*epsilonG*(l*(1+gamma1/(betaE*betaE*gammaE)) - l1*(1 + gammaE/(beta1*beta1*gamma1)) + 2*epsilonG/(betaE*betaE*gammaE*beta1*beta1*gamma1))));

	if (sigma != sigma) {
		printf("sigma p-e = NaN in bremsstrahlung evaluator\n");
		printLog("sigma p-e = NaN in bremsstrahlung evaluator\n");
		exit(0);
	}

	if (sigma < 0) {
		printf("sigma p-e < 0 in bremsstrahlung evaluator\n");
		printLog("sigma p-e < 0 in bremsstrahlung evaluator\n");
		exit(0);
	}

	return sigma;
}

double BremsstrahlungEvaluator::evaluateSigma(const double& gammaE, const double& epsilonG) {
	double sigma = 0;
	if (my_ionNumber > 0) {
		sigma += evaluateSigmape(gammaE, epsilonG) * my_effectiveProtonConcentration;
	}

	sigma += evaluateSigmaee(gammaE, epsilonG);

	return sigma/me_c2;
}

//todo need ambient electron concentration?
BremsstrahlungEvaluator::BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, bool absorption, bool doppler) : RadiationEvaluator(Ne, Emin, Emax, absorption, doppler) {
	my_ionNumber = 0;
	my_ionConcentrations = NULL;
	my_ionConcentrations = NULL;
	my_effectiveProtonConcentration = 0;
}

BremsstrahlungEvaluator::BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, double protonsRelativeConcentration, bool absorption, bool doppler) : RadiationEvaluator(Ne, Emin, Emax, absorption, doppler) {
	my_ionNumber = 1;
	my_ionConcentrations = new double[my_ionNumber];
	my_ionConcentrations[0] = protonsRelativeConcentration;
	my_ionCharges = new int[my_ionNumber];
	my_ionCharges[0] = 1;
	my_effectiveProtonConcentration = protonsRelativeConcentration;
}


BremsstrahlungEvaluator::BremsstrahlungEvaluator(int Ne, const double& Emin, const double& Emax, int ionNumber, double* ionConcentrations, int* ionCharges, bool absorption, bool doppler) : RadiationEvaluator(Ne, Emin, Emax, absorption, doppler) {
	my_ionNumber = ionNumber;
	my_ionConcentrations = new double[my_ionNumber];
	my_ionCharges = new int[my_ionNumber];
	my_effectiveProtonConcentration = 0;
	for (int i = 0; i < my_ionNumber; ++i) {
		my_ionConcentrations[i] = ionConcentrations[i];
		my_ionCharges[i] = ionCharges[i];
		my_effectiveProtonConcentration += my_ionConcentrations[i] * my_ionCharges[i] * my_ionCharges[i];
	}
}

BremsstrahlungEvaluator::~BremsstrahlungEvaluator()
{
	if (my_ionNumber > 0) {
		delete[] my_ionConcentrations;
		delete[] my_ionCharges;
	}
}

void BremsstrahlungEvaluator::resetParameters(const double* parameters, const double* normalizationUnits)
{
}

void BremsstrahlungEvaluator::evaluateEmissivityAndAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source, double& I, double& A) {
	I = evaluateEmissivity(photonFinalEnergy, ix1, iz, ix2, source);

	MassiveParticleDistribution* electronDistribution = source->getParticleDistribution(ix1, iz, ix2);
	MassiveParticleMaxwellDistribution* maxwellDistribution = dynamic_cast<MassiveParticleMaxwellDistribution*>(electronDistribution);
	MassiveParticleMaxwellJuttnerDistribution* maxwellJuttnerDistribution = dynamic_cast<MassiveParticleMaxwellJuttnerDistribution*>(electronDistribution);
	double temperature;
	if (maxwellDistribution != NULL) {
		temperature = maxwellDistribution->getTemperature();
	}
	else if (maxwellJuttnerDistribution != NULL) {
		temperature = maxwellJuttnerDistribution->getTemperature();
	}
	else {
		printf("free-free absorption works only with maxwellian distribution\n");
		printLog("free-free absorption works only with maxwellian distribution\n");
		exit(0);
	}

	double photonFinalFrequency = photonFinalEnergy / hplank;
	double B = (2 * hplank * cube(photonFinalFrequency) / (speed_of_light2)) / (exp(photonFinalEnergy / (kBoltzman * temperature)) - 1.0);

	A = hplank * I / B;
}

double BremsstrahlungEvaluator::evaluateEmissivity(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	MassiveParticleIsotropicDistribution* electronDistribution = dynamic_cast<MassiveParticleIsotropicDistribution*>(source->getParticleDistribution(ix1, iz, ix2));
	if (electronDistribution == NULL) {
		printf("Bremsstrahlung e-e evaluator works only with isotropic electrons distribution\n");
		printLog("Bremsstrahlung e-e evaluator works only with isotropic electrons distribution\n");
		exit(0);
	}

	double result = 0;

	for (int i = 0; i < my_Ne; ++i) {
		double electronEnergy = my_Ee[i];
		double delectronEnergy;
		if (i == 0) {
			delectronEnergy = my_Ee[1] - my_Ee[0];
		}
		else {
			delectronEnergy = my_Ee[i] - my_Ee[i - 1];
		}
		double electronGamma = electronEnergy / (massElectron * speed_of_light2);
		double electronBeta = sqrt(1.0 - 1.0 / (electronGamma * electronGamma));
		double electronKineticEnergy = massElectron * speed_of_light2 * (electronGamma - 1.0);
		double epsilonG = photonFinalEnergy / (massElectron * speed_of_light2);
		double concentration = source->getConcentration(ix1, iz, ix2);

		double sigma = evaluateSigma(electronGamma, epsilonG);

		result += photonFinalEnergy * (speed_of_light * electronBeta) * sigma * concentration * concentration * (4 * pi * electronDistribution->distributionNormalized(electronEnergy)) * delectronEnergy / (4 * pi);

		if (result != result) {
			printf("result = NaN in bremsstrahlung\n");
			printLog("result = NaN in bremsstrahlung\n");
			exit(0);
		}
	}

	return result;
}

double BremsstrahlungEvaluator::evaluateAbsorption(const double& photonFinalEnergy, int ix1, int iz, int ix2, RadiationSource* source)
{
	double I = evaluateEmissivity(photonFinalEnergy, ix1, iz, ix2, source);

	MassiveParticleDistribution* electronDistribution = source->getParticleDistribution(ix1, iz, ix2);
	MassiveParticleMaxwellDistribution* maxwellDistribution = dynamic_cast<MassiveParticleMaxwellDistribution*>(electronDistribution);
	MassiveParticleMaxwellJuttnerDistribution* maxwellJuttnerDistribution = dynamic_cast<MassiveParticleMaxwellJuttnerDistribution*>(electronDistribution);
	double temperature;
	if (maxwellDistribution != NULL) {
		temperature = maxwellDistribution->getTemperature();
	}
	else if (maxwellJuttnerDistribution != NULL) {
		temperature = maxwellJuttnerDistribution->getTemperature();
	}
	else {
		printf("free-free absorption works only with maxwellian distribution\n");
		printLog("free-free absorption works only with maxwellian distribution\n");
		exit(0);
	}

	double photonFinalFrequency = photonFinalEnergy / hplank;
	double B = (2 * hplank * cube(photonFinalFrequency) / (speed_of_light2)) / exp(photonFinalEnergy / (kBoltzman * temperature) - 1);

	double a = hplank * I / B;

	return a;
}