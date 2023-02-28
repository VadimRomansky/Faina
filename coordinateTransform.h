#ifndef coordinate_transform_h
#define coordinate_transform_h

void LorentzTransformationPhotonZ(const double& gamma, const double& Einit, const double& cosThetaInit, double& Eprime, double& cosThetaPrime);
void LorentzTransformationPhotonReverseZ(const double& gamma, const double& Einit, const double& cosThetaInit, double& Eprime, double& cosThetaPrime);

// transform from one spherical system to rotated.Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& mur, const double& phir, const double& mu0, const double& phi0, double& mu1, double& phi1);
void inverseRotationSphericalCoordinates(const double& mur, const double& phir, const double& mu1, const double& phi1, double& mu0, double& phi0);

#endif