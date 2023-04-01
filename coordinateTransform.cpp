#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "coordinateTransform.h"

void LorentzTransformationPhotonZ(const double& gamma, const double& Einit, const double& thetaInit, double& Eprime, double& thetaPrime) {
	//double gamma = 1.0 / sqrt(1.0 - beta * beta);
	double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
	double delta = relativisticDelta(gamma);
	double cosThetaInit = cos(thetaInit);
	double epsilon = versin(thetaInit);

	if (epsilon == 0) {
		//cosThetaPrime = 1.0;
		thetaPrime = 0;
		return;
	}
	if (delta == 0) {
		//cosThetaPrime = -1.0;
		thetaPrime = pi;
		return;
	}

	double factor = epsilon + delta - epsilon * delta;
	Eprime = gamma * factor * Einit;
	double versinThetaPrime = (2 * epsilon -epsilon*delta)/ factor;
	if (versinThetaPrime > 2.0 && versinThetaPrime < 2.0000001) {
		//printf("versin = %g > 2.0 reduced to 2.0\n", versinThetaPrime);
		//printLog("versin = %g > 2.0 reduced to 2.0\n", versinThetaPrime);
		versinThetaPrime = 2.0;
	}
	if (versinThetaPrime > 2.0) {
		printf("versin = %g > 2.0\n", versinThetaPrime);
		printLog("versin = %g > 2.0\n", versinThetaPrime);
		exit(0);
	}
	double sinThetaPrime = sqrt(2 * versinThetaPrime - versinThetaPrime * versinThetaPrime);
	//double cosThetaPrime = (delta - epsilon) / factor;
	//checkAndFixCosValue(cosThetaPrime);
	//thetaPrime = acos(cosThetaPrime);
	thetaPrime = asin(sinThetaPrime);
	if (versinThetaPrime > 1.0) {
		thetaPrime = pi - thetaPrime;
	}
}

void LorentzTransformationPhotonReverseZ(const double& gamma, const double& Einit, const double& thetaInit, double& Eprime, double& thetaPrime) {
	//double gamma = 1.0 / sqrt(1.0 - beta * beta);
	double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
	double delta = relativisticDelta(gamma);
	double cosThetaInit = cos(thetaInit);
	double epsilon = versin(pi - thetaInit);

	if (epsilon == 0) {
		//cosThetaPrime = -1.0;
		thetaPrime = pi;
		return;
	}
	if (delta == 0) {
		//cosThetaPrime = 1.0;
		thetaPrime = 0;
		return;
	}
	double factor = epsilon + delta - epsilon * delta;
	Eprime = gamma * factor * Einit;
	double versinThetaPrime = (2 * delta - epsilon*delta) / factor;
	if (versinThetaPrime > 2.0 && versinThetaPrime < 2.0000001) {
		//printf("versin = %g > 2.0 reduced to 2.0\n", versinThetaPrime);
		//printLog("versin = %g > 2.0 reduced to 2.0\n", versinThetaPrime);
		versinThetaPrime = 2.0;
	}
	if (versinThetaPrime > 2.0) {
		printf("versin = %g > 2.0\n", versinThetaPrime);
		printLog("versin = %g > 2.0\n", versinThetaPrime);
		exit(0);
	}
	double sinThetaPrime = sqrt(2 * versinThetaPrime - versinThetaPrime * versinThetaPrime);
	//double cosThetaPrime = (epsilon - delta) / factor;
	//checkAndFixCosValue(cosThetaPrime);
	//thetaPrime = acos(cosThetaPrime);
	thetaPrime = asin(sinThetaPrime);
	if (versinThetaPrime > 1.0) {
		thetaPrime = pi - thetaPrime;
	}
}

//transform from one spherical system to rotated. Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& thetar, const double& phir, const double& theta0, const double& phi0, double& theta1, double& phi1) {
	/*if (mur > 1.0) {
		printf("mu_r = %g > 1.0\n", mur);
		printLog("mu_r = %g > 1.0\n", mur);
		exit(0);
	}
	if (mur < -1.0) {
		printf("mu_r = %g < -1.0\n", mur);
		printLog("mu_r = %g < -1.0\n", mur);
		exit(0);
	}
	if (mu0 > 1.0) {
		printf("mu0 = %g > 1.0\n", mu0);
		printLog("mu0 = %g > 1.0\n", mu0);
		exit(0);
	}
	if (mu0 < -1.0) {
		printf("mu0 = %g < -1.0\n", mu0);
		printLog("mu0 = %g < -1.0\n", mu0);
		exit(0);
	}*/
	//double sinThetar = sqrt(1.0 - mur * mur);
	double sinThetar = sin(thetar);
	double mur = cos(thetar);
	//double sinTheta0 = sqrt(1.0 - mu0 * mu0);
	double sinTheta0 = sin(theta0);
	double mu0 = cos(theta0);
	double tempPhi = phi0 - phir;
	if (tempPhi < 0) {
		tempPhi = tempPhi + 2 * pi;
	}
	double mu1 = mu0 * mur + sinThetar * sinTheta0 * sin(tempPhi);

	if (mu1 > 1.0 && mu1 < 1.00001) {
		printf("mu1 = %g > 1.0 reduced to 1.0\n", mu1);
		printLog("mu1 = %g > 1.0 reduced to 1.0\n", mu1);
		mu1 = 1.0;
	}
	if (mu1 < -1.0 && mu1 > -1.00001) {
		printf("mu1 = %g < -1.0 reduced to =-1.0\n", mu1);
		printLog("mu1 = %g < -1.0 reduced to -1.0\n", mu1);
		mu1 = -1.0;
	}
	if (mu1 > 1.0) {
		printf("mu1 = %g > 1.0\n", mu1);
		printLog("mu1 = %g > 1.0\n", mu1);
		exit(0);
	}
	if (mu1 < -1.0) {
		printf("mu1 = %g < -1.0\n", mu1);
		printLog("mu1 = %g < -1.0\n", mu1);
		exit(0);
	}
	theta1 = acos(mu1);
	/*if (mu1 >= (1.0 - 1E-8) && mu0 >(1.0 - 1E-8) && mur > (1.0 - 1E-8)) {
		double Xidelta = 0.5 * (sinThetar * sinThetar + sinTheta0 * sinTheta0) - sinThetar * sinTheta0 * sin(tempPhi);
		theta1 = sqrt(2 * Xidelta);
	}*/

	/*if (mu1 <= (-1.0 + 1E-8) && sinThetar > 0 && sinTheta0 > 0) {
		double Xidelta = 0.5 * (sinThetar * sinThetar + sinTheta0 * sinTheta0) - sinThetar * sinTheta0 * sin(tempPhi);
		theta1 = sqrt(2 * Xidelta);
	}*/

	double sinTheta1 = sin(theta1);
	if (sinTheta1 == 0) {
		phi1 = 0;
		theta1 = 0;
		return;
	}
	double coshi = sinTheta0 * cos(tempPhi);
	double cosPhi1 = coshi / sinTheta1;
	if (cosPhi1 > 1.0 && cosPhi1 < 1.00001) {
		printf("cosPhi1 = %g > 1.0 reduced to 1.0\n", cosPhi1);
		printLog("cosPhi1 = %g > 1.0 reduced to 1.0\n", cosPhi1);
		cosPhi1 = 1.0;
	}
	if (cosPhi1 < -1.0 && cosPhi1 > -1.00001) {
		printf("cosPhi1 = %g < -1.0 reduced to -1.0\n", cosPhi1);
		printLog("cosPhi1 = %g < -1.0 reduced to -1.0\n", cosPhi1);
		cosPhi1 = -1.0;
	}
	if (cosPhi1 > 1.0) {
		printf("cosPhi1 = %g > 1.0\n", cosPhi1);
		printLog("cosPhi1 = %g > 1.0\n", cosPhi1);
		exit(0);
	}
	if (cosPhi1 < -1.0) {
		printf("cosPhi1 = %g < -1.0\n", cosPhi1);
		printLog("cosPhi1 = %g < -1.0\n", cosPhi1);
		exit(0);
	}
	double scalarMult = sin(tempPhi) * sinTheta0 * mur - mu0 * sinThetar;
	phi1 = acos(cosPhi1);
	if (scalarMult < 0) {
		phi1 = 2 * pi - phi1;
	}
}

void inverseRotationSphericalCoordinates(const double& thetar, const double& phir, const double& theta1, const double& phi1, double& theta0, double& phi0) {
	/*if (mur > 1.0) {
		printf("mu_r = %g > 1.0\n", mur);
		printLog("mu_r = %g > 1.0\n", mur);
		exit(0);
	}
	if (mur < -1.0) {
		printf("mu_r = %g < -1.0\n", mur);
		printLog("mu_r = %g < -1.0\n", mur);
		exit(0);
	}
	if (mu1 > 1.0) {
		printf("mu1 = %g > 1.0\n", mu1);
		printLog("mu1 = %g > 1.0\n", mu1);
		exit(0);
	}
	if (mu1 < -1.0) {
		printf("mu1 = %g < -1.0\n", mu1);
		printLog("mu1 = %g < -1.0\n", mu1);
		exit(0);
	}*/
	//double sinThetar = sqrt(1.0 - mur * mur);
	double sinThetar = sin(thetar);
	double mur = cos(thetar);
	//double sinTheta1 = sqrt(1.0 - mu1 * mu1);
	double sinTheta1 = sin(theta1);
	double mu1 = cos(theta1);


	double mu0 = mu1 * mur + sinThetar * sinTheta1 * cos(phi1 + pi / 2);
	if (mu0 > 1.0 && mu0 < 1.00001) {
		printf("mu0 = %g > 1.0 reduced to 1.0\n", mu0);
		printLog("mu0 = %g > 1.0 reduced to 1.0\n", mu0);
		mu0 = 1.0;
	}
	if (mu0 < -1.0 && mu0 > -1.00001) {
		printf("mu0 = %g < -1.0 reduced to =-1.0\n", mu0);
		printLog("mu0 = %g < -1.0 reduced to -1.0\n", mu0);
		mu0 = -1.0;
	}
	if (mu0 > 1.0) {
		printf("mu0 = %g > 1.0\n", mu0);
		printLog("mu0 = %g > 1.0\n", mu0);
		exit(0);
	}
	if (mu0 < -1.0) {
		printf("mu0 = %g < -1.0\n", mu0);
		printLog("mu0 = %g < -1.0\n", mu0);
		exit(0);
	}
	theta0 = acos(mu0);

	/*if (mu0 >= (1.0 - 1E-8) && mu1 >(1.0 - 1E-8) && mur > (1.0 - 1E-8)) {
		double Xidelta = 0.5 * (sinThetar * sinThetar + sinTheta1 * sinTheta1) - sinThetar * sinTheta1 * cos(phi1 + pi / 2);
		theta0 = sqrt(2 * Xidelta);
	}*/

	double sinTheta0 = sin(theta0);
	if (sinTheta0 == 0) {
		phi0 = 0;
		theta0 = 0;
		return;
	}

	double coshi = sinTheta1 * cos(phi1);

	double cosTempPhi = coshi / sinTheta0;
	if (cosTempPhi > 1.0 && cosTempPhi < 1.00001) {
		printf("cosTempPhi = %g > 1.0 reduced to 1.0\n", cosTempPhi);
		printLog("cosTempPhi = %g > 1.0 reduced to 1.0\n", cosTempPhi);
		cosTempPhi = 1.0;
	}
	if (cosTempPhi < -1.0 && cosTempPhi > -1.00001) {
		printf("cosTempPhi = %g < -1.0 reduced to -1.0\n", cosTempPhi);
		printLog("cosTempPhi = %g < -1.0 reduced to -1.0\n", cosTempPhi);
		cosTempPhi = -1.0;
	}
	if (cosTempPhi > 1.0) {
		printf("cosTempPhi = %g > 1.0\n", cosTempPhi);
		printLog("cosTempPhi = %g > 1.0\n", cosTempPhi);
		exit(0);
	}
	if (cosTempPhi < -1.0) {
		printf("cosTempPhi = %g < -1.0\n", cosTempPhi);
		printLog("cosTempPhi = %g < -1.0\n", cosTempPhi);
		exit(0);
	}
	double scalarMult = sin(phi1) * sinTheta1 * mur + mu1 * sinThetar;
	double tempPhi = acos(cosTempPhi);
	if (scalarMult < 0) {
		tempPhi = 2 * pi - tempPhi;
	}
	phi0 = tempPhi + phir;
	if (phi0 > 2 * pi) {
		phi0 = phi0 - 2 * pi;
	}
	if (phi0 < 0) {
		phi0 = phi0 + 2 * pi;
	}
}