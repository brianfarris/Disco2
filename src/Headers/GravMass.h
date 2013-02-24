#ifndef PLANET_H
#define PLANET_H
struct GravMass;

#ifdef PLANET_PRIVATE_DEFS
struct GravMass{
   double r;
   double phi;
   double M;
   double E;
   double L;
   double vr;
   double RK_r;
   double RK_phi;
   double RK_M;
   double RK_E;
   double RK_L;
   double RK_vr;
   double Fr;
   double Fp;
};
#endif

struct GravMass *gravMass_create(int);
void gravMass_destroy(struct GravMass *);
void gravMass_initialize(struct GravMass *);
void gravMass_clean_pi(struct GravMass *);
void gravMass_copy(struct GravMass * );
double gravMass_r(struct GravMass * ,int);
double gravMass_phi(struct GravMass * ,int);
double gravMass_M(struct GravMass * ,int);
#endif 
