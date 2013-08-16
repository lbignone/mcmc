#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Define uniform prior*/
double mcmc_uniform(double x, double x_min, double x_max)
{
return 1.0/(x_max - x_min);
}

/* Define Jeffreys's prior */
double mcmc_jeffreys(double x, double x_min, double x_max)
{
return 1.0/(x*log(x_max/x_min));
}
