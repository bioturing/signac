#include <cmath>
#include <limits>
#include <stdexcept>
#include <boost/math/special_functions/gamma.hpp>
#include "chisq.h"

#define MAX_LOOP 200
#define ACCURACY_EPS 1e-20
#define M_SQRTPI 1.77245385090551602729816748334
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03

static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static double P[] = {
      2.97886562639399288862,
      7.409740605964741794425,
      6.1602098531096305440906,
      5.019049726784267463450058,
      1.275366644729965952479585264,
      0.5641895835477550741253201704
  };
  static double Q[] = {
      3.3690752069827527677,
      9.608965327192787870698,
      17.08144074746600431571095,
      12.0489519278551290360340491,
      9.396034016235054150430579648,
      2.260528520767326969591866945,
      1.0
  };
  double num=0.0, den=0.0;
  int i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}


inline
static double log_erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e = log(e) - x*x;
  return e;
}


double gsl_sf_log_erfc(double x)
{
  /* CHECK_POINTER(result) */

  if(x*x < 10.0*GSL_ROOT6_DBL_EPSILON) {
    const double y = x / M_SQRTPI;
    /* series for -1/2 Log[Erfc[Sqrt[Pi] y]] */
    const double c3 = (4.0 - M_PI)/3.0;
    const double c4 = 2.0*(1.0 - M_PI/3.0);
    const double c5 = -0.001829764677455021;  /* (96.0 - 40.0*M_PI + 3.0*M_PI*M_PI)/30.0  */
    const double c6 =  0.02629651521057465;   /* 2.0*(120.0 - 60.0*M_PI + 7.0*M_PI*M_PI)/45.0 */
    const double c7 = -0.01621575378835404;
    const double c8 =  0.00125993961762116;
    const double c9 =  0.00556964649138;
    const double c10 = -0.0045563339802;
    const double c11 =  0.0009461589032;
    const double c12 =  0.0013200243174;
    const double c13 = -0.00142906;
    const double c14 =  0.00048204;
    double series = c8 + y*(c9 + y*(c10 + y*(c11 + y*(c12 + y*(c13 + c14*y)))));
    series = y*(1.0 + y*(1.0 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*(c7 + y*series)))))));
    return -2.0 * series;
    //result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
  else if(x > 8.0) {
    return log_erfc8(x);
    //result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    return  log(std::erfc(x));
    //result->err  = fabs(result_erfc.err / result_erfc.val);
    //result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
}


double log1pexp(double x) {
    if (x <= -37)
        return exp(x);
    else if (x <= 18)
        return log1p(exp(x));
    else if (x <= 33.3)
        return x + exp(-x);
    else
        return x;
}

double log_chisqr2(int Dof, double x)
{
    return std::log(boost::math::gamma_q(Dof/2.0, x/2));
}


double log_chisqr(int Dof, double x)
{
    x /= 2;

    if (Dof == 1)
        return gsl_sf_log_erfc(sqrt(x));

    double f = -x;
    double i = 1;

    if (Dof & 1) {
        f -= log(M_PI * x)/2;
        i = 0.5;
    }

    double lx = log(x);

    double lvalue = lx - log(i);
    double lsum = lvalue;

    while (++i < Dof/2) {
        lvalue += lx - log(i);
        lsum += log1pexp(lvalue-lsum);
    }

    lsum += f;

    if (Dof & 1)
        lsum += log1pexp(gsl_sf_log_erfc(sqrt(x)) - lsum);
    
    return lsum;
}
