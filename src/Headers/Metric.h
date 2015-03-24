#ifndef METRIC_H
#define METRIC_H

struct Metric;
struct Sim;

#ifdef METRIC_PRIVATE_DEFS
struct Metric
{
    double x[4];
    double g_dd[10];
    double g_uu[10];
    double *dg_dd;
    double *dg_uu;
    double da[4];
    double db[12];
    int killing[4];
    int length_dg;
    struct Sim *theSim;
};
#endif

//create and destroy
struct Metric* (*metric_create)(double, double, double, double, struct Sim *);
struct Metric* metric_create_exact(double t, double r, double p, double z, 
                                struct Sim *theSim);
struct Metric* metric_create_adm(double t, double r, double p, double z, 
                                struct Sim *theSim);
void (*metric_create_der)(struct Metric *g);
void metric_create_der_exact(struct Metric *g);
void metric_create_der_adm(struct Metric *g);
void metric_destroy(struct Metric *g);

//Initialize
void metric_init_background(struct Sim *);
void metric_init_metric(struct Sim *);

//access
double metric_lapse(struct Metric *g);
double metric_shift_u(struct Metric *g, int i);
double metric_shift_d(struct Metric *g, int i);
double metric_gamma_dd(struct Metric *g, int i, int j);
double metric_gamma_uu(struct Metric *g, int i, int j);
double metric_g_dd(struct Metric *g, int i, int j);
double metric_g_uu(struct Metric *g, int i, int j);
double metric_sqrtgamma(struct Metric *g);
double metric_sqrtg(struct Metric *g);
double metric_dg_dd(struct Metric *g, int k, int i, int j);
double metric_dg_uu(struct Metric *g, int k, int i, int j);
double metric_dlapse(struct Metric *g, int k);
double metric_dshift_u(struct Metric *g, int k, int i);
int metric_killcoord(struct Metric *g, int i);

//routines
double metric_square3_u(struct Metric *g, double *vec);
double metric_square3_d(struct Metric *g, double *vec);
double metric_square4_u(struct Metric *g, double *vec);
double metric_square4_d(struct Metric *g, double *vec);
double metric_dot3_u(struct Metric *g, double *a, double *b);
double metric_dot3_d(struct Metric *g, double *a, double *b);
double metric_dot4_u(struct Metric *g, double *a, double *b);
double metric_dot4_d(struct Metric *g, double *a, double *b);
double metric_conn(struct Metric *g, int tau, int mu, int nu);
void metric_shear_uu(struct Metric *g, double *v, double *dv, double *shear, struct Sim *theSim);

//Frames
double (*metric_frame_U_u)(struct Metric *, int, struct Sim *);
double (*metric_frame_dU_du)(struct Metric *, int, int, struct Sim *);

double metric_frame_U_u_euler(struct Metric *g, int mu, struct Sim *theSim);
double metric_frame_dU_du_euler(struct Metric *g, int mu, int nu, struct Sim *theSim);
double metric_frame_U_u_kep(struct Metric *g, int mu, struct Sim *theSim);
double metric_frame_dU_du_kep(struct Metric *g, int mu, int nu, struct Sim *theSim);
double metric_frame_U_u_acc(struct Metric *g, int mu, struct Sim *theSim);
double metric_frame_dU_du_acc(struct Metric *g, int mu, int nu, struct Sim *theSim);

//Exact Metrics
double (*metric_g_dd_exact)(int i, int j, double t, double r, double p, double z, struct Sim *theSim);
double (*metric_g_uu_exact)(int i, int j, double t, double r, double p, double z, struct Sim *theSim);
double (*metric_dg_dd_exact)(int k, int i, int j, double t, double r, double p, double z, struct Sim *theSim);
double (*metric_dg_uu_exact)(int k, int i, int j, double t, double r, double p, double z, struct Sim *theSim);
void (*metric_killing_exact)(int *k);
double (*metric_horizon)(struct Sim *);

//SR
double metric_g_dd_exact_sr(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_g_uu_exact_sr(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_dd_exact_sr(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_uu_exact_sr(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
void metric_killing_exact_sr(int *k);
double metric_horizon_exact_sr(struct Sim *);

//Schwarzschild - Schwarzschild coords
double metric_g_dd_exact_schw_sc(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_g_uu_exact_schw_sc(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_dd_exact_schw_sc(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_uu_exact_schw_sc(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
void metric_killing_exact_schw_sc(int *k);
double metric_horizon_exact_schw_sc(struct Sim *);

//Schwarzschild - Kerr-Schild coords
double metric_g_dd_exact_schw_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_g_uu_exact_schw_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_dd_exact_schw_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_uu_exact_schw_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
void metric_killing_exact_schw_ks(int *k);
double metric_horizon_exact_schw_ks(struct Sim *);

//SR Cartesian coords
double metric_g_dd_exact_sr_cart(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_g_uu_exact_sr_cart(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_dd_exact_sr_cart(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_uu_exact_sr_cart(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
void metric_killing_exact_sr_cart(int *k);
double metric_horizon_exact_sr_cart(struct Sim *);

//Kerr - Kerr-Schild coords
double metric_g_dd_exact_kerr_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_g_uu_exact_kerr_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_dd_exact_kerr_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
double metric_dg_uu_exact_kerr_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim);
void metric_killing_exact_kerr_ks(int *k);
double metric_horizon_exact_kerr_ks(struct Sim *);

//ADM Metrics
double (*metric_lapse_adm)(double, double, double, double, struct Sim *);
double (*metric_shift_adm_exact)(int, double, double, double, double, 
                                struct Sim *);
double (*metric_spatial_adm)(int, int, double, double, double, double, 
                            struct Sim *);
double (*metric_ispatial_adm)(int, int, double, double, double, double, 
                            struct Sim *);
double (*metric_dlapse_adm)(int, double, double, double, double, struct Sim *);
double (*metric_dshift_adm_exact)(int, int, double, double, double, double, 
                                struct Sim *);
double (*metric_dspatial_adm)(int, int, int, double, double, double, double, 
                            struct Sim *);
double metric_shift_adm(int i, double t, double r, double p, double z, 
                        struct Sim *theSim);
double metric_dshift_adm(int mu, int i, double t, double r, double p, 
                        double z, struct Sim *theSim);

// SR
double metric_lapse_adm_sr(double t, double r, double p, double z, 
                            struct Sim *theSim);
double metric_shift_adm_sr(int i, double t, double r, double p, double z, 
                                struct Sim *theSim);
double metric_spatial_adm_sr(int i, int j, double t, double r, double p, 
                            double z, struct Sim *theSim);
double metric_ispatial_adm_sr(int i, int j, double t, double r, double p, 
                            double z, struct Sim *theSim);
double metric_dlapse_adm_sr(int mu, double t, double r, double p, double z, 
                            struct Sim *theSim);
double metric_dshift_adm_sr(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim);
double metric_dspatial_adm_sr(int mu, int i, int j, double t, double r, 
                            double p, double z, struct Sim *theSim);
void metric_killing_adm_sr(int *k);
double metric_horizon_adm_sr(struct Sim *theSim);

// Schwarzschild KS
double metric_lapse_adm_schw_ks(double t, double r, double p, double z, 
                            struct Sim *theSim);
double metric_shift_adm_schw_ks(int i, double t, double r, double p, double z, 
                                struct Sim *theSim);
double metric_spatial_adm_schw_ks(int i, int j, double t, double r, double p, 
                            double z, struct Sim *theSim);
double metric_ispatial_adm_schw_ks(int i, int j, double t, double r, double p, 
                            double z, struct Sim *theSim);
double metric_dlapse_adm_schw_ks(int mu, double t, double r, double p, double z, 
                            struct Sim *theSim);
double metric_dshift_adm_schw_ks(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim);
double metric_dspatial_adm_schw_ks(int mu, int i, int j, double t, double r, 
                            double p, double z, struct Sim *theSim);
void metric_killing_adm_schw_ks(int *k);
double metric_horizon_adm_schw_ks(struct Sim *theSim);

//Boosts

double (*metric_shift_adm_boost)(int, double, double, double, double, 
                                struct Sim *);
double (*metric_dshift_adm_boost)(int, int, double, double, double, double, 
                                struct Sim *);
void (*metric_killing_boost)(int *);

double metric_boost_none(int i, double t, double r, double p, double z, 
                                struct Sim *theSim);
double metric_dboost_none(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim);
void metric_killing_boost_none(int *k);
double metric_boost_rigid(int i, double t, double r, double p, double z, 
                                struct Sim *theSim);
double metric_dboost_rigid(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim);
void metric_killing_boost_rigid(int *k);
double metric_boost_bin(int i, double t, double r, double p, double z, 
                                struct Sim *theSim);
double metric_dboost_bin(int mu, int i, double t, double r, double p, 
                            double z, struct Sim *theSim);
void metric_killing_boost_bin(int *k);
#endif
