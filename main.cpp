#include "stdio.h"
#include "math.h"

#include "constants.h"
#include "electronDistribution.h"
#include "util.h"


int main() {
	double nu = 2;
	double x = 0.1;
	printf("K(%g, %g) = %g", nu, x, McDonaldFunction(nu, x));
	return 0;
}