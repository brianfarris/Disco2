#ifndef PLANET_H
#define PLANET_H
struct GravMass;
struct Grid;

#ifdef PLANET_PRIVATE_DEFS
struct GravMass{
   double r;
   double phi;
   double M;
   double omega;
   //   double E;
   //   double L;
   //   double vr;
   double RK_r;
   double RK_phi;
   double RK_M;
   double RK_omega;
   //   double RK_E;
   //   double RK_L;
   //   double RK_vr;
   //   double Fr;
   //   double Fp;

};
#endif

//create and destroy
struct GravMass *gravMass_create(int);
void gravMass_destroy(struct GravMass *);
//initialization
void gravMass_initialize(struct GravMass *);
//access data
double gravMass_r(struct GravMass * ,int);
double gravMass_phi(struct GravMass * ,int);
double gravMass_M(struct GravMass * ,int);
//miscellaneous
void gravMass_clean_pi(struct GravMass *,struct Grid *);
void gravMass_copy(struct GravMass *,struct Grid *);
void gravMass_move(struct GravMass *,double);
void gravMass_update_RK( struct GravMass * ,struct Grid * , double );

#endif 
