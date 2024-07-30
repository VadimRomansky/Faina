#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include <omp.h>

#include "constants.h"

#include "util.h"

double uniformDistribution() {
    return (rand() % randomParameter + 0.5) / randomParameter;
}

double sqr(const double& a){
    return a*a;
}

double cube(const double& a){
    return a*a*a;
}

double min(const double& a, const double& b){
    if(a < b) {
        return a;
    } else {
        return b;
    }
}

double min3(const double& a, const double& b, const double& c){
    if(a < b) {
        return min(a, c);
    } else {
        return min(b, c);
    }
}

double min4(const double& a, const double& b, const double& c, const double& d){
    if(a < b) {
        if(c < d){
            return min(a, c);
        } else {
            return min(a, d);
        }
    } else {
        if(c < d){
            return min(b, c);
        } else {
            return min(b, d);
        }
    }
}

double max(const double& a, const double& b){
    if(a > b) {
        return a;
    } else {
        return b;
    }
}

std::string convertIntToString(int a) {
    if (a == 0) {
        std::string result = "0";
        return result;
    }
    if (a > 0) {
        std::string result = "";
        while (a > 0) {
            int last = a % 10;
            a = a / 10;
            char c = last + '0';
            result = c + result;
        }
        return result;
    }
    a = -a;
    std::string result = "-";
    return result + convertIntToString(a);
}


void printLog (const char *fmt, ...)
{
    FILE* logFile = fopen("log.dat", "a");
    char buffer[256];
    va_list args;
    va_start (args, fmt);
    vsprintf (buffer,fmt, args);
    fprintf(logFile, buffer);
    fclose(logFile);
}

void resetLog() {
    FILE* logFile = fopen("log.dat", "w");
    fclose(logFile);
}

double McDonaldFunction(double index, double x) {
    //todo approximation with small x!!!
    if (x < 0) {
        printf("x in McDonald < 0\n");
        printLog("x in McDonald < 0\n");
        exit(0);
    }

    if (x > 2 * index * index && x > 10) {
        double result = sqrt(pi / (2 * x)) * exp(-x) *
                (1 + ((4 * index * index - 1) / (8 * x)) +
                 ((4 * index * index - 1) * (4 * index * index - 9) / (128 * x * x)) +
                 ((4 * index * index - 1) * (4 * index * index - 9) * (4 * index * index - 25) / (6 * cube(8 * x))));
        return result;
    }
    if (x < sqrt(index + 1) / 100) {
        if (index == 0) {
            double eulerMaskeroni = 0.5772156649;
            return -log(x / 2) - eulerMaskeroni;
        } else {
            return (tgamma(index) / 2) * pow(2 / x, index);
        }
    }
    double dt;
    double t;
    double prevT = 0;
    double result = 0;
    double maxT;
    double maxLevel = 10;

    /*double discr = index * index - 2 * x * x + 2 * x * maxLevel;
    if (discr < 0) {
        maxT = log(maxLevel * maxLevel + sqrt(maxLevel * maxLevel - 1));
    }
    else {
        maxT = (index + sqrt(discr)) / (x);
    }*/
    maxT = max((1 - index / x + sqrt(sqr(1 - index / x) + 4 * (maxLevel - 1))), 5*log(index/x));
    ;
    double relation = cosh(maxT) - index / (x * maxT);

    /*while (x * cosh(maxT) - index * maxT < 10) {
        maxT = maxT + dt;
    }*/
    int Npoints = 10000;
    dt = maxT / Npoints;

    t = dt;
    prevT = 0;
    int i = 0;
    while (i < Npoints) {
        double middleT = 0.5 * (t + prevT);
        double dresult = exp(-x * cosh(middleT)) * cosh(index * middleT) * dt;
        result += dresult;
        if (dresult < result / 1E14) break;
        prevT = t;
        t += dt;
        ++i;
    }

    return result;
}

//1 - cos
double versin(const double& x) {
    if (fabs(x) > pi / 1000) {
        return 1.0 - cos(x);
    }
    double result = 0;
    double x2 = x * x;
    double temp = x2;
    int factorial = 2.0;
    int i = 2;
    double df = 1.0;
    while(fabs(df) > x2*1E-16){
        df = temp / factorial;
        result = result + df;
        temp = temp * x2;
        factorial = -factorial * (2 * i - 1) * (2 * i);
        ++i;
        if (i > 1000000) {
            printf("versin didn't converge, x = %g\n", x);
            printLog("versin didn't converge, x = %g\n", x);
            exit(0);
        }
    }
    return result;
}

//1 - beta
double relativisticDelta(const double& gamma) {
    if (gamma < 1000) {
        return 1.0 - sqrt(1.0 - 1.0 / (gamma * gamma));
    }
    double result = 0;
    double gamma2 = gamma * gamma;
    double temp = gamma2;
    double factor = 0.5;
    int i = 1;
    double df = 1.0;

    while(df > (1E-16)/gamma2) {
        df = factor / temp;
        result = result + df;
        temp = temp * gamma2;
        factor = factor * (2 * i - 1) / (2.0 * (i + 1));
        ++i;
        if (i > 1000000) {
            printf("relativisticDelta didn't converge, x = %g\n", gamma);
            printLog("relativisticDelta didn't converge, x = %g\n", gamma);
            exit(0);
        }
    }

    return result;
}

void checkAndFixCosValue(double& mu) {
    if (mu > 1.0 && mu < 1.00001) {
        //printf("mu = %g > 1.0 reduced to 1.0\n", mu);
        //printLog("mu = %g > 1.0 reduced to 1.0\n", mu);
        mu = 1.0;
    }
    if (mu < -1.0 && mu > -1.00001) {
       // printf("mu = %g < -1.0 reduced to =-1.0\n", mu);
       // printLog("mu = %g < -1.0 reduced to -1.0\n", mu);
        mu = -1.0;
    }
    if (mu > 1.0) {
        printf("mu = %g > 1.0\n", mu);
        printLog("mu = %g > 1.0\n", mu);
        exit(0);
    }
    if (mu < -1.0) {
        printf("mu = %g < -1.0\n", mu);
        printLog("mu = %g < -1.0\n", mu);
        exit(0);
    }
}

void checkAndFixVersin(double& v) {
    if (v > 2.0 && v < 2.0000001) {
        //printf("versin = %g > 2.0 reduced to 2.0\n", versinThetaPrime);
        //printLog("versin = %g > 2.0 reduced to 2.0\n", versinThetaPrime);
        v = 2.0;
    }
    if (v > 2.0) {
        printf("versin = %g > 2.0\n", v);
        printLog("versin = %g > 2.0\n", v);
        exit(0);
    }
    if (v < 0 && v > -5E-10) {
        v = 0;
    }
    if (v < 0) {
        printf("versin = %g < 0\n", v);
        printLog("versin = %g < 0\n", v);
        exit(0);
    }
}

double*** create3dArray(const int N1, const int N2, const int N3, const double& value)
{
    double*** a = new double** [N1];
    for (int i = 0; i < N1; ++i) {
        a[i] = new double*[N2];
        for (int j = 0; j < N2; ++j) {
            a[i][j] = new double[N3];
            for (int k = 0; k < N3; ++k) {
                a[i][j][k] = value;
            }
        }
    }
    return a;
}

void delete3dArray(double*** a, const int N1, const int N2, const int N3)
{
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            delete[] a[i][j];
        }
        delete[] a[i];
    }
    delete[] a;
}

void write3dArrayToFile(double*** a, const int N1, const int N2, const int N3, const char* fileName)
{
    FILE* outFile = fopen(fileName, "w");

    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            for (int k = 0; k < N3; ++k) {
                fprintf(outFile, "%g\n", a[i][j][k]);
            }
        }
    }

    fclose(outFile);
}

int readRadiationFromFile(double*& E, double*& F, double*& error, const char* fileName) {
    FILE* file = fopen(fileName, "r");
    double year, month, day, time, freq, band, flux, err;
    int ch;
    int N = 0;
    while (EOF != (ch = getc(file))){
        if ('\n' == ch)
            ++N;
    }
    fclose(file);
    if (N == 0) {
        printf("empty input file %s\n", fileName);
        printLog("empty input file %s\n", fileName);
        return N;
    }
    file = fopen(fileName, "r");
    E = new double[N];
    F = new double[N];
    error = new double[N];
    for (int i = 0; i < N; ++i) {
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf", &year, &month, &day, &time, &freq, &band, &flux, &err);
        /*fscanf(file, "%lf", &year);
        fscanf(file, "%lf", &month);
        fscanf(file, "%lf", &day);
        fscanf(file, "%lf", &time);
        fscanf(file, "%lf", &freq);
        fscanf(file, "%lf", &band);
        fscanf(file, "%lf", &flux);
        fscanf(file, "%lf", &err);*/
        E[i] = freq;
        F[i] = flux;
        error[i] = err;
    }
    fclose(file);
    return N;
}

bool checkNaNorInfinity(const double& v) {
    if (v != v) {
        return true;
    }

    if (0 * v != 0 * v) {
        return true;
    }

    return false;
}


//from numerical recepies
double bessj0(double x)
//Returns the Bessel function J0(x) for any real x.
{
    double ax, z;
    double xx, y, ans, ans1, ans2; //Accumulate polynomials in double precision.
    if ((ax = fabs(x)) < 8.0) { //Direct rational function fit.
        y = x * x;
        ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
            + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
        ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
            + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
        ans = ans1 / ans2;
    }
    else { //Fitting function (6.5.9).
        z = 8.0 / ax;
        y = z * z;
        xx = ax - 0.785398164;
        ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
            + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
            + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                - y * 0.934945152e-7)));
        ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
    }
    return ans;
}

double bessj1(double x)
//Returns the Bessel function J1(x) for any real x.
{
    double ax, z;
    double xx, y, ans, ans1, ans2; //Accumulate polynomials in double precision.
    if ((ax = fabs(x)) < 8.0) {
        //Direct rational approximation.
        y = x * x;
        ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
            + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
        ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
            + y * (99447.43394 + y * (376.9991397 + y * 1.0))));
        ans = ans1 / ans2;
    }
    else {
        //Fitting function(6.5.9).
        z = 8.0 / ax;
        y = z * z;
        xx = ax - 2.356194491;
        ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
            + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
        ans2 = 0.04687499995 + y * (-0.2002690873e-3
            + y * (0.8449199096e-5 + y * (-0.88228987e-6
                + y * 0.105787412e-6)));
        ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
        if (x < 0.0) ans = -ans;
    }
    return ans;
}
