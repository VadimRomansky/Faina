#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "electronDistribution.h"
#include "photonDistribution.h"
#include "util.h"


int main() {
	PhotonPlankDistribution* CMBradiation = PhotonPlankDistribution::getCMBradiation();
	PhotonMultiPlankDistribution* galacticRadiation = PhotonMultiPlankDistribution::getGalacticField();
	ElectronDistribution* electrons = new ElectronPowerLawDistribution(3.5, me_c2, 150);
	double nu = 2;
	double x = 0.1;
	printf("K(%g, %g) = %g", nu, x, McDonaldFunction(nu, x));
	return 0;
}