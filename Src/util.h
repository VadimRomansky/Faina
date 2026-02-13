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
bool checkNaNorInfinity(const double& v);

double*** create3dArray(const int N1, const int N2, const int N3, const double& value = 0);
void delete3dArray(double*** a, const int N1, const int N2, const int N3);
double**** create4dArray(const int N1, const int N2, const int N3, const int N4, const double& value = 0);
void delete4dArray(double**** a, const int N1, const int N2, const int N3, const int N4);
void write3dArrayToFile(double*** a, const int N1, const int N2, const int N3, const char* fileName);
int readRadiationFromFile(double*& E, double*& F, double*& error, const char* fileName);

//from numerical recepies
double bessj0(double x);
double bessj1(double x);
double bessj(int n, double x);
double chebev(double a, double b, double c[], int m, double x);
void beschb(double x, double* gam1, double* gam2, double* gampl, double* gammi);
void bessik(double x, double xnu, double* ri, double* rk, double* rip, double* rkp);

double** inverseMatrix(double** matrix, int N);

double* getUvarovBpar(int Nx, double minX, double maxX, double L0);

double* getUvarovBper(int Nx, double minX, double maxX, double L0);

double* getUvarovBpar2(int Nx, double* xgrid, double L0, double factor);

double* getUvarovBper2(int Nx, double* xgrid, double L0, double factor);

double* getUvarovBpar2new(int Nx, double* xgrid, double L0, double factor);

double* getUvarovBper2new(int Nx, double* xgrid, double L0, double factor);
#endif
