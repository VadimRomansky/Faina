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
		Eprime = gamma * delta * Einit;
		return;
	}
	if (delta == 0) {
		//cosThetaPrime = -1.0;
		thetaPrime = pi;
		Eprime = gamma *epsilon * Einit;
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

void LorentzTransformationPhotonReverseZalpha(const double& gamma, const double& Einit, const double& alphaInit, double& Eprime, double& alphaPrime) {
	//double gamma = 1.0 / sqrt(1.0 - beta * beta);
	double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
	double delta = relativisticDelta(gamma);
	//double cosThetaInit = cos(thetaInit);
	double epsilon = versin(alphaInit);

	if (epsilon == 0) {
		//cosThetaPrime = -1.0;
		//thetaPrime = pi;
		alphaPrime = 0;
		return;
	}
	if (delta == 0) {
		//cosThetaPrime = 1.0;
		//thetaPrime = 0;
		alphaPrime = pi;
		return;
	}
	double factor = epsilon + delta - epsilon * delta;
	Eprime = gamma * factor * Einit;
	double versinThetaPrime = (2 * delta - epsilon * delta) / factor;
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
	alphaPrime = asin(sinThetaPrime);
	if (versinThetaPrime > 1.0) {
		alphaPrime = pi - alphaPrime;
	}
}

//transform from one spherical system to rotated. Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& thetar, const double& phir, const double& theta0, const double& phi0, double& theta1, double& phi1) {
	double sinThetar = sin(thetar);
	double mur = cos(thetar);
	double versinr = versin(thetar);
	double versinr_c = versin(pi - thetar);
	double sinTheta0 = sin(theta0);
	double mu0 = cos(theta0);
	double versin0 = versin(theta0);
	double versin0_c = versin(pi - theta0);
	double versin1;
	double versin1_c;
	double sinTheta1;

	double tempPhi = phi0 - phir;
	if (tempPhi < 0) {
		tempPhi = tempPhi + 2 * pi;
	}

	//double mu1 = mu0 * mur + sinThetar * sinTheta0 * sin(tempPhi);

	if ((thetar <= pi / 2) || (theta0 <= pi / 2)) {
		versin1 = versinr + versin0 - versinr * versin0 - sinThetar * sinTheta0 * sin(tempPhi);

		checkAndFixVersin(versin1);
		sinTheta1 = sqrt(2 * versin1 - versin1 * versin1);
		theta1 = asin(sinTheta1);
		if (versin1 > 1) {
			theta1 = pi - theta1;
		}

		if (versin1 == 0) {
			phi1 = 0;
			theta1 = 0;
			return;
		}

		if (versin1 == 2) {
			phi1 = 0;
			theta1 = pi;
			return;
		}
	}
	else {
		versin1_c = 2-versinr_c - versin0_c + versinr_c * versin0_c + sinThetar * sinTheta0 * sin(tempPhi);

		checkAndFixVersin(versin1_c);
		sinTheta1 = sqrt(2 * versin1_c - versin1_c * versin1_c);
		theta1 = pi - asin(sinTheta1);
		if (versin1_c > 1) {
			theta1 = pi - theta1;
		}

		if (versin1_c == 0) {
			phi1 = 0;
			theta1 = pi;
			return;
		}

		if (versin1_c == 2) {
			phi1 = 0;
			theta1 = 0;
			return;
		}
	}

	double cosTempPhi = cos(tempPhi);
	double cosPhi1 = (sinTheta0 / sinTheta1)*cosTempPhi;
	checkAndFixCosValue(cosPhi1);
	double scalarMult = sin(tempPhi) * sinTheta0 * mur - mu0 * sinThetar;
	phi1 = acos(cosPhi1);
	if (scalarMult < 0) {
		phi1 = 2 * pi - phi1;
	}
}

void inverseRotationSphericalCoordinates(const double& thetar, const double& phir, const double& theta1, const double& phi1, double& theta0, double& phi0) {
	//todo
	/*if ((thetar > pi / 2)) {
		double tempThetar = pi - thetar;
		double tempTheta1 = pi - theta1;
		inverseRotationSphericalCoordinates(tempThetar, phi1, tempTheta1, phi1, theta0, phi0);
		theta0 = pi - theta0;
		return;
	}*/
	

	double sinThetar = sin(thetar);
	double mur = cos(thetar);
	double versinr = versin(thetar);
	double versinr_c = versin(pi - thetar);
	double sinTheta1 = sin(theta1);
	double mu1 = cos(theta1);
	double versin1 = versin(theta1);
	double versin1_c = versin(pi - theta1);
	double versin0;
	double versin0_c;
	double sinTheta0;


	//double mu0 = mu1 * mur + sinThetar * sinTheta1 * cos(phi1 + pi / 2);


	if ((thetar <= pi / 2) || (theta1 <= pi / 2)) {
		versin0 = versinr + versin1 - versinr * versin1 - sinThetar * sinTheta1 * cos(phi1 + pi / 2);


		checkAndFixVersin(versin0);
		sinTheta0 = sqrt(2 * versin0 - versin0 * versin0);
		theta0 = asin(sinTheta0);
		if (versin0 > 1) {
			theta0 = pi - theta0;
		}

		if (versin0 == 0) {
			phi0 = 0;
			theta0 = 0;
			return;
		}

		if (versin0 == 2) {
			phi0 = 0;
			theta0 = pi;
			return;
		}
	}
	else {
		versin0_c = 2-versinr_c - versin1_c + versinr_c * versin1_c + sinThetar * sinTheta1 * cos(phi1 + pi / 2);

		checkAndFixVersin(versin0_c);
		sinTheta0 = sqrt(2 * versin0_c - versin0_c * versin0_c);
		theta0 = pi - asin(sinTheta0);
		if (versin0_c > 1) {
			theta0 = pi - theta0;
		}

		if (versin0_c == 0) {
			phi0 = 0;
			theta0 = pi;
			return;
		}

		if (versin0_c == 2) {
			phi0 = 0;
			theta0 = 0;
			return;
		}
	}

	double cosPhi1 = cos(phi1);
	double cosTempPhi = (sinTheta1 / sinTheta0) * cosPhi1;
	checkAndFixCosValue(cosTempPhi);
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