#ifndef coordinate_transform_h
#define coordinate_transform_h

void LorentzTransformationPhotonZ(const double& gamma, const double& Einit, const double& thetaInit, double& Eprime, double& thetaPrime);
void LorentzTransformationPhotonReverseZ(const double& gamma, const double& Einit, const double& thetaInit, double& Eprime, double& thetaPrime);
void LorentzTransformationPhotonReverseZalpha(const double& gamma, const double& Einit, const double& alphaInit, double& Eprime, double& alphaPrime);

// transform from one spherical system to rotated.Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& mur, const double& phir, const double& mu0, const double& phi0, double& mu1, double& phi1);
void inverseRotationSphericalCoordinates(const double& mur, const double& phir, const double& mu1, const double& phi1, double& mu0, double& phi0);

#endif