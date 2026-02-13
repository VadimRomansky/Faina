#pragma once

void evaluateW50bremsstrahlung();

void evaluateW50synchrotron();

void evaluateW50comptonAndSynchrotron();

void evaluateW50comptonAndSynchrotron2();

void evaluateW50comptonAndSynchrotronMCfunctionUpstream();

void evaluateW50comptonAndSynchrotronAdvectionfunction();

void evaluateW50comptonThickRegime();

void evaluateW50comptonAdvectionBigSource();

void evaluateW50comptonAndSynchrotronMCwithoutupstream();

void evaluateW50comptonAndSynchrotronAdvectionfunctionWithUpstream();

void evaluateW50comptonAndSynchrotronAdvectionfunctionWithBrinkmann();

void evaluateW50comptonDiffusion();

void sourceCoordinates(const double& time, double& x, double& y, double& z);

void sourceCoordinates2(const double& time, double& x, double& y, double& z);

double sourcePower(const double& time);

double* getDiffusionCoefficient(double* energy, const int Ne);

void evaluateW50pion();




