#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <omp.h>

double uniformDistribution();
double sqr(const double& a);
double cube(const double& a);
double min(const double& a, const double& b);
double min3(const double& a, const double& b, const double& c);
double min4(const double& a, const double& b, const double& c, const double& d);
double max(const double& a, const double& b);
std::string convertIntToString(int a);
void printLog (const char *fmt, ...);
void resetLog();
double McDonaldFunction(double index, double x);
double versin(const double& x);
double relativisticDelta(const double& gamma);
void checkAndFixCosValue(double& mu);
void checkAndFixVersin(double& v);

double*** create3dArray(const int N1, const int N2, const int N3, const double& value = 0);
void delete3dArray(double*** a, const int N1, const int N2, const int N3);
void write3dArrayToFile(double*** a, const int N1, const int N2, const int N3, const char* fileName);
int readRadiationFromFile(double*& E, double*& F, double*& error, const char* fileName);
#endif
