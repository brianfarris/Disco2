#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/Face.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

void riemann_visc_flux_gr(struct Riemann *theRiemann, struct Sim *theSim)
{
    int i,j,dir;
    double a, sqrtg, hn, cs2, rhoh, visc, height, r, u0;
    double v[3], dv[12], b[3], shear[16], U[4], Ud[4];
    
    double alpha = sim_AlphaVisc(theSim);
    double GAMMA = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    int NUMQ = sim_NUM_Q(theSim);
    double *prim = (double *)malloc(NUMQ * sizeof(double));
    double *F = (double *)malloc(NUMQ * sizeof(double));
    struct Metric *g = metric_create(time_global, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    r = theRiemann->pos[R_DIR];

    for(i=0; i<NUMQ; i++)
        prim[i] = 0.5*(theRiemann->primL[i] + theRiemann->primR[i]);

    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    for(i=0; i<3; i++)
        if(theRiemann->n[i] == 1)
            dir = i;

    dv[0] = 0.0;
    dv[1] = 0.0;
    dv[2] = 0.0;
    dv[3] = 0.5*(cell_gradr(theRiemann->cL, URR) + cell_gradr(theRiemann->cR, URR));
    dv[4] = 0.5*(cell_gradr(theRiemann->cL, UPP) + cell_gradr(theRiemann->cR, UPP));
    dv[5] = 0.5*(cell_gradr(theRiemann->cL, UZZ) + cell_gradr(theRiemann->cR, UZZ));
    dv[6] = 0.5*(cell_gradp(theRiemann->cL, URR) + cell_gradp(theRiemann->cR, URR));
    dv[7] = 0.5*(cell_gradp(theRiemann->cL, UPP) + cell_gradp(theRiemann->cR, UPP));
    dv[8] = 0.5*(cell_gradp(theRiemann->cL, UZZ) + cell_gradp(theRiemann->cR, UZZ));
    dv[9] = 0.5*(cell_gradz(theRiemann->cL, URR) + cell_gradz(theRiemann->cR, URR));
    dv[10] = 0.5*(cell_gradz(theRiemann->cL, UPP) + cell_gradz(theRiemann->cR, UPP));
    dv[11] = 0.5*(cell_gradz(theRiemann->cL, UZZ) + cell_gradz(theRiemann->cR, UZZ));

    if(dir == 0)
    {
        double dphiR = theRiemann->pos[P_DIR] - (cell_tiph(theRiemann->cR) - 0.5*cell_dphi(theRiemann->cR));
        double dphiL = theRiemann->pos[P_DIR] - (cell_tiph(theRiemann->cL) - 0.5*cell_dphi(theRiemann->cL));
        while(dphiR<-M_PI) dphiR += 2*M_PI;
        while(dphiR> M_PI) dphiR -= 2*M_PI;
        while(dphiL<-M_PI) dphiL += 2*M_PI;
        while(dphiL> M_PI) dphiL -= 2*M_PI;
        double idr = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);

        dv[3] = idr * (cell_prim(theRiemann->cR,URR)+cell_gradp(theRiemann->cR,URR)*dphiR - cell_prim(theRiemann->cL,URR)-cell_gradp(theRiemann->cL,URR)*dphiL);
        dv[4] = idr * (cell_prim(theRiemann->cR,UPP)+cell_gradp(theRiemann->cR,UPP)*dphiR - cell_prim(theRiemann->cL,UPP)-cell_gradp(theRiemann->cL,UPP)*dphiL);
        dv[5] = idr * (cell_prim(theRiemann->cR,UZZ)+cell_gradp(theRiemann->cR,UZZ)*dphiR - cell_prim(theRiemann->cL,UZZ)-cell_gradp(theRiemann->cL,UZZ)*dphiL);
    }
    else if (dir == 1)
    {
        double dphi = cell_tiph(theRiemann->cR)-0.5*cell_dphi(theRiemann->cR) - cell_tiph(theRiemann->cL)+0.5*cell_dphi(theRiemann->cL);
        while (dphi < -M_PI) dphi += 2*M_PI;
        while (dphi >  M_PI) dphi -= 2*M_PI;
        double idp = 1.0/dphi;

        dv[6] = idp * (cell_prim(theRiemann->cR,URR) - cell_prim(theRiemann->cL,URR));
        dv[7] = idp * (cell_prim(theRiemann->cR,UPP) - cell_prim(theRiemann->cL,UPP));
        dv[8] = idp * (cell_prim(theRiemann->cR,UZZ) - cell_prim(theRiemann->cL,UZZ)); 
    }
    //TODO: FIX THIS.  Doesn't account for misaligned zones in z faces.
    else
    {
        double idz = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[9] = idz * (cell_prim(theRiemann->cR,URR) - cell_prim(theRiemann->cL,URR));
        dv[10] = idz * (cell_prim(theRiemann->cR,UPP) - cell_prim(theRiemann->cL,UPP));
        dv[11] = idz * (cell_prim(theRiemann->cR,UZZ) - cell_prim(theRiemann->cL,UZZ));
    }
    
    if(PRINTTOOMUCH && 0)
    {
        FILE *gradfile = fopen("grad_face.out","a");
        fprintf(gradfile, "%d, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g\n", dir, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], v[0], v[1], v[2], dv[3], dv[4], dv[5], dv[6], dv[7], dv[8], dv[9], dv[10], dv[11]);
        fclose(gradfile);
    }
    
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g)/r;
    u0 = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
    for(i=0; i<4; i++)
    {
        Ud[i] = 0.0;
        U[i] = metric_frame_U_u(g,i,theSim);
    }
    for(i=0; i<4; i++)
    {
        Ud[0] += metric_g_dd(g,0,i)*U[i];
        Ud[1] += metric_g_dd(g,1,i)*U[i];
        Ud[2] += metric_g_dd(g,2,i)*U[i];
        Ud[3] += metric_g_dd(g,3,i)*U[i];
    }

    metric_shear_uu(g, v, dv, shear, theSim);
    if(sim_N(theSim,Z_DIR)==1)
    {
        for(i=0; i<4; i++)
        {
            shear[4*i+3] = 0.0;
            shear[4*3+i] = 0.0;
        }
    }
    
    //TODO: CHECK THIS!  Probably not relativistic and/or isothermal.
    
    rhoh = prim[RHO] + GAMMA*prim[PPP]/(GAMMA-1.0);
    cs2 = GAMMA*prim[PPP] / rhoh;
    height = sqrt(prim[PPP]*r*r*r/(rhoh*M)) / u0;
    if(alpha > 0)
        visc = alpha * sqrt(cs2) * height * prim[RHO];
    else
        visc = -alpha * prim[RHO];

    if(dir == 1)
        hn = theRiemann->pos[R_DIR];
    else
        hn = 1.0;

    for(i=0; i<NUMQ; i++)
        F[i] = 0.0;

    for(i=0; i<4; i++)
    {
        F[SRR] += metric_g_dd(g,i,1)*shear[4*(dir+1)+i];
        F[LLL] += metric_g_dd(g,i,2)*shear[4*(dir+1)+i];
        F[SZZ] += metric_g_dd(g,i,3)*shear[4*(dir+1)+i];
        F[TAU] += Ud[i] * shear[4*(dir+1)+i];
    }
    F[SRR] *= -a*sqrtg*hn * visc;
    F[LLL] *= -a*sqrtg*hn * visc;
    F[SZZ] *= -a*sqrtg*hn * visc;
    F[TAU] *=  a*sqrtg*hn * visc;
    //F[TAU] = -a*a*sqrtg*hn * visc * shear[dir+1];

    if(PRINTTOOMUCH)
    {
        printf("v: %.12g %.12g %.12g\n", v[0],v[1],v[2]);
        for(i=0; i<4; i++)
        {
            printf("d%dv ",i);
            for(j=0; j<3; j++)
                printf("%.12g ", dv[3*i+j]);
            printf("\n");
        }
        for(i=0; i<4; i++)
        {
            for(j=0; j<4; j++)
                printf("%.12g ", shear[4*i+j]);
            printf("\n");
        }
        printf("prim: %.12g, %.12g, %.12g, %.12g, %.12g,\n", prim[RHO], prim[URR], prim[UPP], prim[UZZ], prim[PPP]);
        printf("nonvisc flux: %.12g, %.12g, %.12g, %.12g,\n", theRiemann->F[SRR], theRiemann->F[LLL], theRiemann->F[SZZ], theRiemann->F[TAU]);
        printf("   visc flux: %.12g, %.12g, %.12g, %.12g,\n", F[SRR], F[LLL], F[SZZ], F[TAU]);
        printf("nu: %.12g\n", visc);
        printf("g: ");
        for(i=0; i<4; i++)
        {
            for(j=0; j<4; j++)
                printf("%.12g ", metric_g_dd(g,i,j));
            printf("\n");
        }

    }

    theRiemann->F[SRR] += F[SRR];
    theRiemann->F[LLL] += F[LLL];
    theRiemann->F[SZZ] += F[SZZ];
    theRiemann->F[TAU] += F[TAU];

    free(prim);
    free(F);
    metric_destroy(g);
}

void riemann_visc_flux_LR_gr(struct Riemann *theRiemann, struct Sim *theSim, int state, double *F)
{
    int i,j,dir;
    double a, sqrtg, hn, cs2, rhoh, visc, height, r, u0;
    double v[3], dv[12], b[3], shear[16];
    
    double alpha = sim_AlphaVisc(theSim);
    double GAMMA = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    int NUMQ = sim_NUM_Q(theSim);
    double *prim = (double *)malloc(NUMQ * sizeof(double));
    struct Metric *g = metric_create(time_global, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    r = theRiemann->pos[R_DIR];
    for(i=0; i<NUMQ; i++)
    {
        if(state == LEFT)
            prim[i] = theRiemann->primL[i];
        else
            prim[i] = theRiemann->primR[i];
    }

    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    for(i=0; i<3; i++)
        if(theRiemann->n[i] == 1)
            dir = i;

    //TODO: FIX THIS.  Doesn't account for misaligned zones in r/z faces.
    struct Cell *mycell;
    if(state == LEFT)
        mycell = theRiemann->cL;
    else
        mycell = theRiemann->cR;

    dv[0] = 0.0;
    dv[1] = 0.0;
    dv[2] = 0.0;
    dv[3] = cell_gradr(mycell, URR);
    dv[4] = cell_gradr(mycell, UPP);
    dv[5] = cell_gradr(mycell, UZZ);
    dv[6] = cell_gradp(mycell, URR);
    dv[7] = cell_gradp(mycell, UPP);
    dv[8] = cell_gradp(mycell, UZZ);
    dv[9] = cell_gradz(mycell, URR);
    dv[10] = cell_gradz(mycell, UPP);
    dv[11] = cell_gradz(mycell, UZZ);

    /*
    if(dir == 0)
    {
        double idr = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[3] = idr * (theRiemann->cR->prim[URR] - theRiemann->cL->prim[URR]);
        dv[4] = idr * (theRiemann->cR->prim[UPP] - theRiemann->cL->prim[UPP]);
        dv[5] = idr * (theRiemann->cR->prim[UZZ] - theRiemann->cL->prim[UZZ]);
    }
    else if (dir == 1)
    {
        double idp = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[6] = idp * (theRiemann->cR->prim[URR] - theRiemann->cL->prim[URR]);
        dv[7] = idp * (theRiemann->cR->prim[UPP] - theRiemann->cL->prim[UPP]);
        dv[8] = idp * (theRiemann->cR->prim[UZZ] - theRiemann->cL->prim[UZZ]);
    }
    else
    {
        double idz = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[9] = idz * (theRiemann->cR->prim[URR] - theRiemann->cL->prim[URR]);
        dv[10] = idz * (theRiemann->cR->prim[UPP] - theRiemann->cL->prim[UPP]);
        dv[11] = idz * (theRiemann->cR->prim[UZZ] - theRiemann->cL->prim[UZZ]);
    }
    */
    
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g)/r;
    u0 = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));

    metric_shear_uu(g, v, dv, shear, theSim);
    if(PRINTTOOMUCH)
    {
        printf("v: %.12g %.12g %.12g\n", v[0],v[1],v[2]);
        for(i=0; i<4; i++)
        {
            printf("d%dv ",i);
            for(j=0; j<3; j++)
                printf("%.12g ", dv[3*i+j]);
            printf("\n");
        }
        for(i=0; i<4; i++)
        {
            for(j=0; j<4; j++)
                printf("%.12g ", shear[4*i+j]);
            printf("\n");
        }
    }

    //TODO: CHECK THIS!  Probably not relativistic and/or isothermal.
    rhoh = prim[RHO] + GAMMA*prim[PPP]/(GAMMA-1.0);
    cs2 = GAMMA*prim[PPP] / rhoh;
    height = sqrt(prim[PPP]*r*r*r/(rhoh*M)) / u0;
    if(alpha > 0)
        visc = alpha * sqrt(cs2) * height * prim[RHO];
    else
        visc = -alpha * prim[RHO];

    if(dir == 1)
        hn = theRiemann->pos[R_DIR];
    else
        hn = 1.0;

    for(i=0; i<NUMQ; i++)
        F[i] = 0.0;

    for(i=0; i<4; i++)
    {
        F[SRR] += metric_g_dd(g,i,1)*shear[4*(dir+1)+i];
        F[LLL] += metric_g_dd(g,i,2)*shear[4*(dir+1)+i];
        F[SZZ] += metric_g_dd(g,i,3)*shear[4*(dir+1)+i];
    }
    F[SRR] *= -a*sqrtg*hn * visc;
    F[LLL] *= -a*sqrtg*hn * visc;
    F[SZZ] *= -a*sqrtg*hn * visc;
    F[TAU] = -a*a*sqrtg*hn * visc * shear[dir+1];

    free(prim);
    metric_destroy(g);
}

