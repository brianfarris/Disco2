#ifndef TSOLVE_H
#define TSOLVE_H
struct Tsolve;

#ifdef TSOLVE_PRIVATE_DEFS
struct Tsolve{
    double guess;
    double rhoe_or_Ptot;    
    double Sigma;
    double Gamma;
    double xi_r;
    double xi_g;
    double Omega_eff;
};
#endif

//create and destroy
struct Tsolve * tsolve_create(double, double, double,double, double, double);
void tsolve_destroy(struct Tsolve *);

//other
double tsolve_E_eq (double , void *);
double tsolve_E_eq_deriv (double, void *);
void tsolve_E_eq_fdf (double, void *, double *, double *);
double tsolve_P_eq (double , void *);
double tsolve_P_eq_deriv (double, void *);
void tsolve_P_eq_fdf (double, void *, double *, double *);
double tsolve_findroot_E(struct Tsolve *);
double tsolve_findroot_P(struct Tsolve *);

#endif
