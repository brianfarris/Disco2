#ifndef PLANET_H
#define PLANET_H
struct GravMass;
struct Sim;
struct Cell; // need Cell for live bnary update_RK()

#ifdef PLANET_PRIVATE_DEFS
struct GravMass{
   double M;
   double r;
   double phi;
   double omega;
   double E;
   double L;
   double vr;
   double Fr;
   double Fp;
   double RK_r;
   double RK_phi;
   double RK_M;
   double RK_omega;
   double RK_E;
   double RK_L;
   double RK_vr;
};
#endif

//create and destroy
struct GravMass *gravMass_create(int);
void gravMass_destroy(struct GravMass *);
//initialization
void gravMass_init_none(struct GravMass *,struct Sim *);
void gravMass_init_single(struct GravMass *,struct Sim *);
void gravMass_init_binary(struct GravMass *,struct Sim *);
void gravMass_init_livebinary(struct GravMass *,struct Sim *);
void (*gravMass_init_ptr(struct Sim * ))(struct GravMass *,struct Sim *);
//void gravMass_set_chkpt(struct GravMass * ,int ,double ,double ,double, double ); //w/o live binary params
void gravMass_set_chkpt(struct GravMass * ,int ,double ,double ,double, double, double, double, double, double, double );
//access data
double gravMass_r(struct GravMass * ,int);
double gravMass_phi(struct GravMass * ,int);
double gravMass_M(struct GravMass * ,int);
double gravMass_omega(struct GravMass * ,int );
  //For live binary -DD
double gravMass_L(struct GravMass * ,int);
double gravMass_E(struct GravMass * ,int);
double gravMass_vr(struct GravMass * ,int );
double gravMass_Fr(struct GravMass * ,int);
double gravMass_Fp(struct GravMass * ,int );
//set Data
void GravMass_set_Fr( struct GravMass *, int, double);
void GravMass_set_Fp( struct GravMass *, int, double);
void GravMass_set_omega( struct GravMass *, int, double);
//miscellaneous
void gravMass_clean_pi(struct GravMass *,struct Sim *);
void gravMass_copy(struct GravMass *,struct Sim *);
void gravMass_adv_anly( struct GravMass * , double, double );
//void gravMass_move(struct Sim *, struct GravMass *,double, double RK); //added RK as last argument
void gravMass_move(struct Sim *, struct GravMass *,double); //without RK as last argument
//void gravMass_move(struct GravMass *,double);
void gravMass_update_RK( struct Cell *** , struct GravMass * , struct Sim * , double, double );
//For live binary -DD
// OFF CM BINARY
void gravMass_beta(struct Sim *, struct GravMass *, double *, double *);
void gravMass_gam(struct Sim *, struct GravMass *, double *, double *);
#endif 
