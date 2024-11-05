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

#define ACC 40.0 //Make larger to increase accuracy.
#define BIGNO 1.0e10
#define BIGNI 1.0e-10
double bessj(int n, double x)
//Returns the Bessel function Jn(x) for any real xand n ≥ 2.
{
    int j, jsum, m;
    double ax, bj, bjm, bjp, sum, tox, ans;
    if (n < 2) {
        printf("Index n less than 2 in bessj");
        exit(0);
    }
    ax = fabs(x);
    if (ax == 0.0)
        return 0.0;
    else if (ax > (float) n) {
        //Upwards recurrence from J0and J1.
        tox = 2.0 / ax;
        bjm = bessj0(ax);
        bj = bessj1(ax);
        for (j = 1; j < n; j++) {
            bjp = j * tox * bj - bjm;
            bjm = bj;
            bj = bjp;
        }
        ans = bj;
    }
    else {
        //Downwards recurrence from an even m here computed
        tox = 2.0 / ax;
        m = 2 * ((n + (int)sqrt(ACC * n)) / 2);
        jsum = 0; //jsum will alternate between 0 and 1; when it is1, we accumulate in sum the even terms in(5.5.16).
        bjp = ans = sum = 0.0;
        bj = 1.0;
        for (j = m; j > 0; j--) {
            //The downward recurrence.
            bjm = j * tox * bj - bjp;
            bjp = bj;
            bj = bjm;
            if (fabs(bj) > BIGNO) {
                //Renormalize to prevent overflows.
                bj *= BIGNI;
                bjp *= BIGNI;
                ans *= BIGNI;
                sum *= BIGNI;
            }
            if (jsum) sum += bj; //Accumulate the sum.
                jsum = !jsum; //Change 0 to 1 or vice versa.
                if (j == n) ans = bjp; //Save the unnormalized answer.
        }
        sum = 2.0 * sum - bj; //Compute(5.5.16)
        ans /= sum;//and use it to normalize the answer.
    }
    return x < 0.0 && (n & 1) ? -ans : ans;
}

double chebev(double a, double b, double c[], int m, double x)
/*Chebyshev evaluation : All arguments are input.c[0..m - 1] is an array of Chebyshev coeffi -
cients, the first m elements of c output from chebft(which must have been called with the
    same aand b).The Chebyshev polynomial ?m - 1
    k = 0 ckTk(y) ? c0 / 2 is evaluated at a point
    y = [x ?(b + a) / 2] / [(b ? a) / 2], and the result is returned as the function value.*/
{
    //void nrerror(char error_text[]);
    double d = 0.0, dd = 0.0, sv, y, y2;
    int j;
    if ((x - a) * (x - b) > 0.0) { 
        printf("x not in range in routine chebev"); 
        exit(0);
    }
    y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a)); //Change of variable.
        for (j = m - 1; j >= 1; j--) {
            //Clenshaw’s recurrence.
            sv = d;
            d = y2 * d - dd + c[j];
            dd = sv;
        }
    return y * d - dd + 0.5 * c[0]; //Last step is different.
}

#define NUSE1 5
#define NUSE2 5
void beschb(double x, double* gam1, double* gam2, double* gampl, double* gammi)
/*Evaluates  and  by Chebyshev expansion for |x| ? 1/2. Also returns 1/(1 + x) and
1/(1 ? x). If converting to double precision, set NUSE1 = 7, NUSE2 = 8.*/
{
    //float chebev(float a, float b, float c[], int m, float x);
    double xx;
    static double c1[] = {
        -1.142022680371168e0,6.5165112670737e-3,
        3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
        3.67795e-11,-1.356e-13 };
    static double c2[] = {
        1.843740587300905e0,-7.68528408447867e-2,
        1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
        2.423096e-10,-1.702e-13,-1.49e-15 };
    xx = 8.0 * x * x - 1.0; //Multiply x by 2 to make range be ?1 to 1, and then apply transformation for evaluating even Chebyshev series.
    *gam1 = chebev(-1.0, 1.0, c1, NUSE1, xx);
    *gam2 = chebev(-1.0, 1.0, c2, NUSE2, xx);
    *gampl = *gam2 - x * (*gam1);
    *gammi = *gam2 + x * (*gam1);
}


#define EPS 1.0e-10
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define PI 3.141592653589793
void bessik(double x, double xnu, double* ri, double* rk, double* rip, double* rkp)
/*Returns the modified Bessel functions ri = I , rk = K and their derivatives rip = I ,
rkp = K, for positive x and for xnu. The relative accuracy is within one or two
significant digits of EPS. FPMIN is a number close to the machine’s smallest floating-point
number. All internal arithmetic is in double precision. To convert the entire routine to double
precision, change the float declarations above to double and decrease EPS to 10?16. Also
convert the function beschb.*/
{
    //void beschb(double x, double* gam1, double* gam2, double* gampl, double* gammi);
    //void nrerror(char error_text[]);
    int i, l, nl;
    double a, a1, b, c, d, del, del1, delh, dels, e, f, fact, fact2, ff, gam1, gam2,
        gammi, gampl, h, p, pimu, q, q1, q2, qnew, ril, ril1, rimu, rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1, x2, xi, xi2, xmu, xmu2;
    if (x <= 0.0 || xnu < 0.0) {
        printf("bad arguments in bessik\n");
        exit(0);
    }
    nl = (int)(xnu + 0.5); //nl is the number of downward recurrences of the I’s and upwardrecurrences of K’s. xmulies be-tween ?1/2 and 1/2.
    xmu = xnu - nl;
    xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    h = xnu * xi; //Evaluate CF1 by modified Lentz’smethod (5.2).
    if (h < FPMIN) h = FPMIN;
    b = xi2 * xnu;
    d = 0.0;
    c = h;
    for (i = 1; i <= MAXIT; i++) {
        b += xi2;
        d = 1.0 / (b + d); //Denominators cannot be zero here,so no need for special precautions.
        c = b + 1.0 / c;
        del = c * d;
        h = del * h;
        if (fabs(del - 1.0) < EPS) break;
    }
    if (i > MAXIT) {
        printf("x too large in bessik; try asymptotic expansion"); 
        exit(0);
    }
    ril = FPMIN; //Initialize I and I for downward recurrence.
    ripl = h * ril;
    ril1 = ril; //Store values for later rescaling.
    rip1 = ripl;
    fact = xnu * xi;
    for (l = nl; l >= 1; l--) {
        ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }
    f = ripl / ril; //Now have unnormalized I and I.
    if (x < XMIN) { //Use series.
        x2 = 0.5 * x;
        pimu = PI * xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
        d = -log(x2);
        e = xmu * d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e) / e);
        beschb(xmu, &gam1, &gam2, &gampl, &gammi); //Chebyshev evaluation of 1 and 2.
        ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        sum = ff;
        e = exp(e);
        p = 0.5 * e / gampl;
        q = 0.5 / (e * gammi);
        c = 1.0;
        d = x2 * x2;
        sum1 = p;
        for (i = 1; i <= MAXIT; i++) {
            ff = (i * ff + p + q) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q /= (i + xmu);
            del = c * ff;
            sum += del;
            del1 = c * (p - i * ff);
            sum1 += del1;
            if (fabs(del) < fabs(sum) * EPS) break;
        }
        if (i > MAXIT) {
            printf("bessk series failed to converge");
            exit(0);
        }
        rkmu = sum;
        rk1 = sum1 * xi2;
    }
    else { //Evaluate CF2 by Steed’s algorithm (5.2), which is OK because there can be no zero denominators.
        b = 2.0 * (1.0 + x);
        d = 1.0 / b;
        h = delh = d;
        q1 = 0.0; //Initializations for recurrence (6.7.35).
        q2 = 1.0;
        a1 = 0.25 - xmu2;
        q = c = a1; //First term in equation (6.7.34).
        a = -a1;
        s = 1.0 + q * delh;
        for (i = 2; i <= MAXIT; i++) {
            a -= 2 * (i - 1);
            c = -a * c / i;
            qnew = (q1 - b * q2) / a;
            q1 = q2;
            q2 = qnew;
            q += c * qnew;
            b += 2.0;
            d = 1.0 / (b + a * d);
            delh = (b * d - 1.0) * delh;
            h += delh;
            dels = q * delh;
            s += dels;
            if (fabs(dels / s) < EPS) break;
            //Need only test convergence of sum since CF2 itself converges more quickly.
        }
        if (i > MAXIT) {
            printf("bessik: failure to converge in cf2\n");
            exit(0);
        }
        h = a1 * h;
        rkmu = sqrt(PI / (2.0 * x)) * exp(-x) / s; //Omit the factor exp(?x) to scaleall the returned functions by exp(x) for x ? XMIN.
        rk1 = rkmu * (xmu + x + 0.5 - h) * xi;
    }
    rkmup = xmu * xi * rkmu - rk1;
    rimu = xi / (f * rkmu - rkmup); //Get I? from Wronskian.
    *ri = (rimu * ril1) / ril; //Scale original I and I .
    *rip = (rimu * rip1) / ril;
    for (i = 1; i <= nl; i++) { //Upward recurrence of K .
        rktemp = (xmu + i) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    *rk = rkmu;
    *rkp = xnu * xi * rkmu - rk1;
}

double** inverseMatrix(double** matrix, int N)
{
    double** temp = new double* [N];
    double** rightPart = new double* [N];
    for (int i = 0; i < N; ++i) {
        temp[i] = new double[N];
        rightPart[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            temp[i][j] = matrix[i][j];
            rightPart[i][j] = 0;
            if (i == j) {
                rightPart[i][j] = 1;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int k = i; k < N; ++k) {
            double divisor = temp[k][i];
            for (int j = 0; j < N; ++j) {
                temp[k][j] /= divisor;
                rightPart[k][j] /= divisor;
            }
        }

        for (int k = 0; k < N; ++k) {
            if (k != i) {
                double mult = temp[k][i];
                //temp[k][i] = 0;
                //rightPart[k][i] -= rightPart[i][i] * mult;
                for (int j = 0; j < N; ++j) {
                    temp[k][j] -= temp[i][j] * mult;
                    rightPart[k][j] -= rightPart[i][j] * mult;
                }
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        delete[] temp[i];
    }
    delete[] temp;

    return rightPart;
}
