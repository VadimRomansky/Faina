#include <cstdlib>
#include <math.h>
#include "stdio.h"
#include <time.h>
#include <cmath>
#include <string>
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

void printLog(const char* s){
	FILE* logFile = fopen("log.dat", "a");
	printf(s);
	fprintf(logFile, s);
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