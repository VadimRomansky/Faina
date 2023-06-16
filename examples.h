#ifndef EXAMPLES_H
#define EXAMPLES_H

//example 0 evaluation simple synchrotron
void evaluateSimpleSynchrotron();

// example 1. Evaluating inverse compton flux of powerlaw distributed electrons on CMB radiation
void evaluateComtonWithPowerLawDistribution();

//example 2. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with simple flat disk source and powerlaw distribution
void fitCSS161010withPowerLawDistribition();

//example 3. Fitting observed synchrotron radio fluxes from CSS1601010 at one time moment with electron distributions read from files
void fitCSS161010withTabulatedDistributions();

//example 4. Fitting observed synchrotron radio fluxes from CSS161010 at 3 time moments
void fitTimeDependentCSS161010();

// example 5. Evaluating pion decay gamma flux of powerlaw distributed protons in cygnus cocoon
void evaluatePionDecay();

// example 6. Evaluating bremsstrahlung radiation from hot gas
void evaluateBremsstrahlung();

//example 7 compare compton and synchrotron
void compareComptonSynchrotron();

//example 8 evaluation synchrotron Image
void evaluateSynchrotronImage();



#endif
