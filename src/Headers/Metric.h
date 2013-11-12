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
    int killing[4];
    int length_dg;
};
#endif

//create and destroy
struct Metric* metric_create(double t, double r, double p, double z);
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

//Exact Metrics
double (*metric_g_dd_exact)(int i, int j, double t, double r, double p, double z);
double (*metric_g_uu_exact)(int i, int j, double t, double r, double p, double z);
double (*metric_dg_dd_exact)(int k, int i, int j, double t, double r, double p, double z);
double (*metric_dg_uu_exact)(int k, int i, int j, double t, double r, double p, double z);
void (*metric_killing_exact)(int *k);

//SR
double metric_g_dd_exact_sr(int mu, int nu, double t, double r, double p, double z);
double metric_g_uu_exact_sr(int mu, int nu, double t, double r, double p, double z);
double metric_dg_dd_exact_sr(int k, int mu, int nu, double t, double r, double p, double z);
double metric_dg_uu_exact_sr(int k, int mu, int nu, double t, double r, double p, double z);
void metric_killing_exact_sr(int *k);
#endif
