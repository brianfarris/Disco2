#ifndef METRIC_H
#define METRIC_H

struct Metric;
struct Sim;

#ifdef METRIC_PRIVATE_DEFS
struct Metric
{
    double lapse;
    double shift[3]; //Contravariant shift \beta^i
    double gamma_dd[6]; //Spatial Metric tensor \gamma_{ij}
    double gamma_uu[6]; //Inverse spatial metric \gamma^{ij}
};
#endif

//create and destroy
struct Metric* metric_create(struct Sim *, double t, double r, double p, double z);
void metric_destroy(struct Metric *g);

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

#endif
