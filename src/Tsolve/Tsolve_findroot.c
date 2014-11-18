#define TSOLVE_PRIVATE_DEFS
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "../Headers/Tsolve.h"

double tsolve_findroot_E(struct Tsolve * theTsolve){
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    gsl_function_fdf FDF;
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);

    FDF.f = &tsolve_E_eq;
    FDF.df = &tsolve_E_eq_deriv;
    FDF.fdf = &tsolve_E_eq_fdf;
    FDF.params = theTsolve;

    double Temp0; //store previous iteration value
    double Temp = theTsolve->guess; // set initial guess here
    gsl_root_fdfsolver_set (s, &FDF, Temp);

    do
    {
        iter++;
        status = gsl_root_fdfsolver_iterate (s);
        Temp0 = Temp;
        Temp = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta (Temp, Temp0, 0, 1e-12);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fdfsolver_free (s);

    return(Temp);
}

double tsolve_findroot_P(struct Tsolve * theTsolve){
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    gsl_function_fdf FDF;
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);

    FDF.f = &tsolve_P_eq;
    FDF.df = &tsolve_P_eq_deriv;
    FDF.fdf = &tsolve_P_eq_fdf;
    FDF.params = theTsolve;

    double Temp0; //store previous iteration value
    double Temp = theTsolve->guess; // set initial guess here
    gsl_root_fdfsolver_set (s, &FDF, Temp);

    do
    {
        iter++;
        status = gsl_root_fdfsolver_iterate (s);
        Temp0 = Temp;
        Temp = gsl_root_fdfsolver_root (s);
        status = gsl_root_test_delta (Temp, Temp0, 0, 1e-12);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fdfsolver_free (s);

    return(Temp);
}
