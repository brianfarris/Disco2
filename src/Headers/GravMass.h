#ifndef PLANET_H
#define PLANET_H
struct GravMass;
struct Sim;

#ifdef PLANET_PRIVATE_DEFS
struct GravMass{
  double OrbShrinkTscale;
  double OrbShrinkT0;
  double r;
  double phi;
  double M;
  double omega;
  double E;
  double L;
  double Ltot;
  double vr;
  double vp;
  double Fr;
  double Fp;
  double RK_r;
  double RK_phi;
  double RK_M;
  double RK_omega;
  double RK_E;
  double RK_L;
  double RK_Ltot;
  double RK_vr;
  double Mdot;
  double Macc;
  double total_torque;
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
//void gravMass_set_chkpt(struct GravMass * ,int ,double ,double ,double, double );
void gravMass_set_chkpt(struct GravMass * ,int ,double ,double ,double, double, double, double, double, double, double, double ); //DD
//access data
double gravMass_r(struct GravMass * ,int);
double gravMass_phi(struct GravMass * ,int);
double gravMass_M(struct GravMass * ,int);
double gravMass_omega(struct GravMass * ,int );
double gravMass_Mdot(struct GravMass *,int);
double gravMass_Macc(struct GravMass *,int); 
double gravMass_total_torque(struct GravMass *,int);
//For live binary -DD
double gravMass_L(struct GravMass * ,int);
double gravMass_Ltot(struct GravMass * ,int);
double gravMass_E(struct GravMass * ,int);
double gravMass_vr(struct GravMass * ,int );
double gravMass_Fr(struct GravMass * ,int);
double gravMass_Fp(struct GravMass * ,int );
//set Data
//void GravMass_set_Fr( struct GravMass *, int, double);
//void GravMass_set_Fp( struct GravMass *, int, double);
void GravMass_set_omega( struct GravMass *, int, double);
void GravMass_set_Ltot( struct GravMass *, int, double);
//void GravMass_set_L( struct GravMass *, int, double);
//miscellaneous
void gravMass_clean_pi(struct GravMass *,struct Sim *);
void gravMass_copy(struct GravMass *,struct Sim *);
void gravMass_move(struct Sim *,struct GravMass *,double,double);
void gravMass_update_RK( struct GravMass * ,struct Sim * , double, double ); //DD added double for dt livebin
void gravMass_set_Mdot( struct GravMass * , double, int );
void gravMass_set_Macc( struct GravMass * , double, int );
void gravMass_set_total_torque( struct GravMass * , double , int );
//For moving Bin
void gravMass_adv_anly(struct GravMass * , double, double ); 
void gravMass_adv_arb( struct Sim *, struct GravMass * , double, double, double );
#endif 
