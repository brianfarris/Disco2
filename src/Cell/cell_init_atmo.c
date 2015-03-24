#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

//ATMO: A Tenuous atMOsphere

void cell_init_atmo_normal(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho0, T0, rho, T, vr, vp;
    double R0, q;

    rho0 = sim_InitPar1(theSim);
    T0 = sim_InitPar2(theSim);
    R0 = sim_InitPar3(theSim);

    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, b[3];
    a = metric_lapse(g);
    b[0] = metric_shift_u(g,0);
    b[1] = metric_shift_u(g,1);
    b[2] = metric_shift_u(g,2);
    metric_destroy(g);

    rho = rho0;
    T = T0;
    vr = -b[0];
    vp = -b[1];
    q = r < R0 ? 1.0 : 0.0;
/*
    if(sim_Background(theSim) != GRDISC)
    {
        double tp[5];
        tp[RHO] = rho0;
        tp[TTT] = T0;
        tp[URR] = vr;
        tp[UPP] = vp;
        tp[UZZ] = -b[2];
        double eps = eos_eps(tp, theSim);
        double P = eos_ppp(tp, theSim);
        double H = sqrt(r*r*r*P/(M*(rho+rho*eps+P))) * a;
        rho = rho0 * H;
        T = P * H;
    }
*/
    c->prim[RHO] = rho;
    c->prim[TTT] = T;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = -b[2];

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}

void cell_init_atmo_nozzle(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho0, T0, rho, T, vr, vp;
    double rho1, T1, v1, impactpar, r1, q;

    rho0 = sim_InitPar1(theSim);
    T0 = sim_InitPar2(theSim);
    rho1 = sim_InitPar3(theSim);
    T1 = sim_InitPar4(theSim);
    r1 = sim_InitPar5(theSim);
    impactpar = sim_InitPar6(theSim);

    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, b[3];
    a = metric_lapse(g);
    b[0] = metric_shift_u(g,0);
    b[1] = metric_shift_u(g,1);
    b[2] = metric_shift_u(g,2);
    metric_destroy(g);

    rho = rho0;
    T = T0;
    vr = -b[0];
    vp = -b[1];
    q = 0.0;

    double x = r*cos(phi);
    double y = r*sin(phi);
    double sina = impactpar/r1;
    double yp = fabs(y + sina/sqrt(1-sina*sina)*(x - r1)) / sqrt(1.0+sina*sina/(1-sina*sina));

    if(yp < 20*M && r > 0.8*r1 && (phi < 0.5*M_PI-asin(sina) || phi > 1.5*M_PI-asin(sina)))
    {
        double v = sqrt(2*M/r1); //Newtonian parabolic orbit.
        double scale = exp(-yp*yp/(2*(3*M)*(3*M))) * 0.5*(atan((r-r1-20*M)/(20*M))+1.0);
        vr += sina * v * scale;
        vp += sqrt(1-sina*sina) * v/r * scale;
        rho += rho1 * scale;
        T += T1 * scale;
        q = 1.0;
    }
/*
    if(sim_Background(theSim) != GRDISC)
    {
        double tp[5];
        tp[RHO] = rho0;
        tp[TTT] = T0;
        tp[URR] = vr;
        tp[UPP] = vp;
        tp[UZZ] = -b[2];
        double eps = eos_eps(tp, theSim);
        double P = eos_ppp(tp, theSim);
        double H = sqrt(r*r*r*P/(M*(rho+rho*eps+P))) * a;
        rho = rho0 * H;
        T = P * H;
    }
*/
    c->prim[RHO] = rho;
    c->prim[TTT] = T;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = -b[2];

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}

void cell_init_atmo_disc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int i;
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double alpha = sim_AlphaVisc(theSim);

    double rho0, T0, rho1, T1, rho, T, vr, vp;
    double R0, q;

    rho0 = sim_InitPar1(theSim);
    T0 = sim_InitPar2(theSim);
    R0 = sim_InitPar3(theSim);
    rho1 = sim_InitPar4(theSim);
    T1 = sim_InitPar5(theSim);

    struct Metric *g = metric_create(time_global, r, phi, z, theSim);
    double a, b[3];
    a = metric_lapse(g);
    b[0] = metric_shift_u(g,0);
    b[1] = metric_shift_u(g,1);
    b[2] = metric_shift_u(g,2);
    metric_destroy(g);

    rho = rho0;
    T = T0;
    vr = -b[0];
    vp = -b[1];
    q = 0.0;

    if(r < R0 && r > 6*M)
    {
        rho = rho1;
        T = T1;
        vr = (-2*M/r) / (1+2*M/r);
        vp = sqrt(M/(r*r*(r+2*M)));
        q = 1.0;
    }
/*
    if(sim_Background(theSim) != GRDISC)
    {
        double tp[5];
        tp[RHO] = rho0;
        tp[TTT] = T0;
        tp[URR] = vr;
        tp[UPP] = vp;
        tp[UZZ] = -b[2];
        double eps = eos_eps(tp, theSim);
        double P = eos_ppp(tp, theSim);
        double H = sqrt(r*r*r*P/(M*(rho+rho*eps+P))) * a;
        rho = rho0 * H;
        T = P * H;
    }
*/
    c->prim[RHO] = rho;
    c->prim[TTT] = T;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = -b[2];

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q;
}
void cell_init_atmo_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int opt = sim_InitPar0(theSim);

    if(opt == 0)
        cell_init_atmo_normal(c, r, phi, z, theSim);
    else if(opt == 1)
        cell_init_atmo_nozzle(c, r, phi, z, theSim);
    else if(opt == 2)
        cell_init_atmo_disc(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_atmo given bad option.\n");
}

void cell_single_init_atmo(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_atmo_calc(theCell, r, t, z, theSim);
}

void cell_init_atmo(struct Cell ***theCells,struct Sim *theSim,
                    struct MPIsetup *theMPIsetup)
{
    int i, j, k;
    for (k = 0; k < sim_N(theSim,Z_DIR); k++)
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double z = 0.5*(zm+zp);

        for (i = 0; i < sim_N(theSim,R_DIR); i++)
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = 0.5*(rm+rp);

            for (j = 0; j < sim_N_p(theSim,i); j++)
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;

                cell_init_atmo_calc(&(theCells[k][i][j]), r, t, z, theSim);

                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n",
                          i,j,k,r,t,z,theCells[k][i][j].prim[RHO], theCells[k][i][j].prim[URR],
                          theCells[k][i][j].prim[UPP], theCells[k][i][j].prim[UZZ],
                          theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
