#include <cmath>
#include <limits>
#include <stdexcept>
#include "chisq.h"

#define M_SQRTPI 1.7724538509055160272981674833411451827975494561224
#define M_LN_SQRTPI 0.57236494292470008707171367567652935582364740645766

double log_erfc(double x)
{
    double x2 = x * x;

    const int order = 11;

    if (x < 0.4) {
        //continue fraction
        double tmp = 0;
        int sign = 2 * (order & 1) - 1;

        for (int i = order; i > 0; --i) {
            tmp = (2 * i * x2)/(2*i + 1 + sign * tmp);
            sign = -sign;
        }

        return log1p(-2/M_SQRTPI * exp(-x2) * x / (1-tmp));
    } else if (x > 12) {

        double tmp = 0;

        for (int i = order; i > 0; --i)
            tmp = (0.5 * i)/(x + tmp);

        return - M_LN_SQRTPI - x2 - log(x + tmp);
    } else {
        return  log(std::erfc(x));
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

double log_chisqr(int Dof, double x)
{
    if (Dof == 0)
        return 0;

    if (x <= 0)
        return 0;

    x /= 2;

    if (Dof == 1)
        return log_erfc(sqrt(x));

    if (Dof == 2)
        return -x;

    double f = -x;
    double i = 1;

    if (Dof & 1) {
        f -= log(M_PI * x)/2;
        i = 0.5;
    }

    double lx = log(x);

    double lvalue = lx - log(i);
    double lsum = (Dof & 1) ? lvalue : log1p(x/i);

    while (++i < Dof/2) {
        lvalue += lx - log(i);
        lsum += log1pexp(lvalue-lsum);
    }

    lsum += f;

    if (Dof & 1)
        lsum += log1pexp(log_erfc(sqrt(x)) - lsum);
    
    return lsum;
}
