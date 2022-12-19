#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "util.h"

#include "coordinateTransform.h"

void LorentzTransformationPhotonZ(const double& beta, const double& Einit, const double& cosThetaInit, double& Eprime, double& cosThetaPrime) {
	double gamma = 1.0 / sqrt(1.0 - beta * beta);

	Eprime = gamma*(1 - beta*cosThetaInit)*Einit;
	cosThetaPrime = (cosThetaInit - beta)/(1 - cosThetaInit*beta);
}

//transform from one spherical system to rotated. Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& mur, const double& phir, const double& mu0, const double& phi0, double& mu1, double& phi1) {
	if (mur > 1.0) {
		printf("mu_r = %g > 1.0\n", mur);
		exit(0);
	}
	if (mur < -1.0) {
		printf("mu_r = %g < -1.0\n", mur);
		exit(0);
	}
	if (mu0 > 1.0) {
		printf("mu0 = %g > 1.0\n", mu0);
		exit(0);
	}
	if (mu0 < -1.0) {
		printf("mu0 = %g < -1.0\n", mu0);
		exit(0);
	}
	double sinThetar = sqrt(1.0 - mur * mur);
	double sinTheta0 = sqrt(1.0 - mu0 * mu0);
	double tempPhi = phi0 - phir;
	if (tempPhi < 0) {
		tempPhi = tempPhi + 2 * pi;
	}
	mu1 = mu0 * mur + sinThetar * sinTheta0 * sin(tempPhi);

	if (mu1 > 1.0) {
		printf("mu1 = %g > 1.0\n", mu1);
		exit(0);
	}
	if (mu1 < -1.0) {
		printf("mu1 = %g < -1.0\n", mu1);
		exit(0);
	}

	double sinTheta1 = sqrt(1.0 - mu1 * mu1);
	if (sinTheta1 == 0) {
		phi1 = 0;
		return;
	}
	double coshi = sinTheta0 * cos(tempPhi);
	double cosPhi1 = coshi / sinTheta1;
	if (cosPhi1 > 1.0) {
		printf("cosPhi1 = %g > 1.0\n", cosPhi1);
		exit(0);
	}
	if (cosPhi1 < -1.0) {
		printf("cosPhi1 = %g < -1.0\n", cosPhi1);
		exit(0);
	}
	double scalarMult = sin(tempPhi) * sinTheta0 * mur - mu0 * sinThetar;
	phi1 = acos(cosPhi1);
	if (scalarMult < 0) {
		phi1 = 2 * pi - phi1;
	}
}

void inverseRotationSphericalCoordinates(const double& mur, const double& phir, const double& mu1, const double& phi1, double& mu0, double& phi0) {
	if (mur > 1.0) {
		printf("mu_r = %g > 1.0\n", mur);
		exit(0);
	}
	if (mur < -1.0) {
		printf("mu_r = %g < -1.0\n", mur);
		exit(0);
	}
	if (mu1 > 1.0) {
		printf("mu1 = %g > 1.0\n", mu1);
		exit(0);
	}
	if (mu1 < -1.0) {
		printf("mu1 = %g < -1.0\n", mu1);
		exit(0);
	}
	double sinThetar = sqrt(1.0 - mur * mur);
	double sinTheta1 = sqrt(1.0 - mu1 * mu1);

	mu0 = mu1 * mur + sinThetar * sinTheta1 * cos(phi1 + pi / 2);
	if (mu0 > 1.0) {
		printf("mu0 = %g > 1.0\n", mu0);
		exit(0);
	}
	if (mu0 < -1.0) {
		printf("mu0 = %g < -1.0\n", mu0);
		exit(0);
	}
	double sinTheta0 = sqrt(1.0 - mu0 * mu0);
	if (sinTheta0 == 0) {
		phi0 = 0;
		return;
	}

	double coshi = sinTheta1 * cos(phi1);

	double cosTempPhi = coshi / sinTheta0;
	if (cosTempPhi > 1.0) {
		printf("cosTempPhi = %g > 1.0\n", cosTempPhi);
		exit(0);
	}
	if (cosTempPhi < -1.0) {
		printf("cosTempPhi = %g < -1.0\n", cosTempPhi);
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