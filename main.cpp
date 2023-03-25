#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "massiveParticleDistribution.h"
#include "photonDistribution.h"
#include "util.h"
#include "inverseCompton.h"
#include "radiationSource.h"
#include "synchrotron.h"
#include "optimization.h"
#include "pionDecay.h"
#include "bremsstrahlung.h"
#include "examples.h"


int main() {
	//evaluateSimpleSynchrotron();
	//evaluateComtonWithPowerLawDistribution();
	//fitCSS161010withPowerLawDistribition();
	//fitCSS161010withTabulatedDistributions();
	fitTimeDependentCSS161010();
	//evaluatePionDecayWithPowerLawDistribution();
	//evaluateBremsstrahlung();
	//compareComptonSynchrotron();
	return 0;
}
