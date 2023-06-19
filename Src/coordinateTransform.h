#ifndef coordinate_transform_h
#define coordinate_transform_h

void LorentzTransformationPhotonZ(const double& gamma, const double& Einit, const double& thetaInit, double& Eprime, double& thetaPrime);
void LorentzTransformationPhotonReverseZ(const double& gamma, const double& Einit, const double& thetaInit, double& Eprime, double& thetaPrime);
void LorentzTransformationPhotonReverseZalpha(const double& gamma, const double& Einit, const double& alphaInit, double& Eprime, double& alphaPrime);

// transform from one spherical system to rotated.Rotation on phir around z, and then on mur around x' 
void rotationSphericalCoordinates(const double& thetar, const double& phir, const double& theta0, const double& phi0, double& theta1, double& phi1);
void inverseRotationSphericalCoordinates(const double& thetar, const double& phir, const double& theta1, const double& phi1, double& theta0, double& phi0);

#endif