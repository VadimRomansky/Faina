#ifndef inverse_compton_h
#define inverse_compton_h

#include "electronDistribution.h"
#include "photonDistribution.h"

//correct name?
double evaluateComptonLuminocity(const double& photonFinalEnergy, const double& photonFinalTheta, const double& photonFinalPhi, PhotonDistribution* photonDistribution, ElectronDistribution* electronDistribution, const double& volume, const double& distance, const double& Emin, const double& Emax, const int Ne, const int Nmu, const int Nphi);

#endif