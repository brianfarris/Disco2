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

void riemann_visc_flux(struct Riemann *theRiemann, struct Sim *theSim)
{
    int i,j,k,dir;
    double a, sqrtg, du0, hn, cs2, rhoh, visc, height, r;
    double v[3], dv[12], b[3], u[4], du[16], shear[16];
    
    double alpha = sim_AlphaVisc(theSim);
    double GAMMA = sim_GAMMALAW(theSim);
    double M = sim_GravM(theSim);
    int NUMQ = sim_NUM_Q(theSim);
    double *prim = (double *)malloc(NUMQ * sizeof(double));
    double *F = (double *)malloc(NUMQ * sizeof(double));
    struct Metric *g = metric_create(time_global, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    r = theRiemann->pos[R_DIR];

    //TODO: account for z-shearing too.
    /*
    if(dir == 0)
    {
        double dphiR = theRiemann->pos[P_DIR] - (cell_tiph(theRiemann->cR) - 0.5*cell_dphi(theRiemann->cR));
        double dphiL = theRiemann->pos[P_DIR] - (cell_tiph(theRiemann->cL) - 0.5*cell_dphi(theRiemann->cL));
        while(dphiR<-M_PI) dphiR += 2*M_PI;
        while(dphiR> M_PI) dphiR -= 2*M_PI;
        while(dphiL<-M_PI) dphiL += 2*M_PI;
        while(dphiL> M_PI) dphiL -= 2*M_PI;

        for(i=0; i<NUMQ; i++)
            prim[i] = 0.5*(cell_prim(theRiemann->cL,i)+cell_gradp(theRiemann->cL,i)*dphiL + cell_prim(theRiemann->cR,i)+cell_gradp(theRiemann->cR,i)*dphiR);
    }
    else
    */
        for(i=0; i<NUMQ; i++)
            prim[i] = 0.5*(theRiemann->primL[i] + theRiemann->primR[i]);

    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    for(i=0; i<3; i++)
        if(theRiemann->n[i] == 1)
            dir = i;

    //TODO: FIX THIS.  Doesn't account for misaligned zones in r/z faces.
    
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
/*
        double PLM = sim_PLM(theSim);
        if(fabs(dv[3]) > PLM*fabs(cell_gradr(theRiemann->cR,URR)))
            dv[3] = PLM*cell_gradr(theRiemann->cR,URR);
        if(fabs(dv[3]) > PLM*fabs(cell_gradr(theRiemann->cL,URR)))
            dv[3] = PLM*cell_gradr(theRiemann->cL,URR);
        if(cell_gradr(theRiemann->cL,URR)*cell_gradr(theRiemann->cR,URR) < 0.0)
            dv[3] = 0.0;
        if(fabs(dv[4]) > PLM*fabs(cell_gradr(theRiemann->cR,UPP)))
            dv[4] = PLM*cell_gradr(theRiemann->cR,UPP);
        if(fabs(dv[4]) > PLM*fabs(cell_gradr(theRiemann->cL,UPP)))
            dv[4] = PLM*cell_gradr(theRiemann->cL,UPP);
        if(cell_gradr(theRiemann->cL,UPP)*cell_gradr(theRiemann->cR,UPP) < 0.0)
            dv[4] = 0.0;
        if(fabs(dv[5]) > PLM*fabs(cell_gradr(theRiemann->cR,UZZ)))
            dv[5] = PLM*cell_gradr(theRiemann->cR,UZZ);
        if(fabs(dv[5]) > PLM*fabs(cell_gradr(theRiemann->cL,UZZ)))
            dv[5] = PLM*cell_gradr(theRiemann->cL,UZZ);
        if(cell_gradr(theRiemann->cL,UZZ)*cell_gradr(theRiemann->cR,UZZ) < 0.0)
            dv[5] = 0.0;
  */
    }
    else if (dir == 1)
    {
        double dphi = cell_tiph(theRiemann->cR)-0.5*cell_dphi(theRiemann->cR) - cell_tiph(theRiemann->cL)+0.5*cell_dphi(theRiemann->cL);
        while (dphi < -M_PI) dphi += 2*M_PI;
        while (dphi >  M_PI) dphi -= 2*M_PI;
        double idp = 1.0/dphi;
        //double idp = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[6] = idp * (cell_prim(theRiemann->cR,URR) - cell_prim(theRiemann->cL,URR));
        dv[7] = idp * (cell_prim(theRiemann->cR,UPP) - cell_prim(theRiemann->cL,UPP));
        dv[8] = idp * (cell_prim(theRiemann->cR,UZZ) - cell_prim(theRiemann->cL,UZZ));
    /*    
        double PLM = sim_PLM(theSim);
        if(fabs(dv[6]) > PLM*fabs(cell_gradp(theRiemann->cR,URR)))
            dv[6] = PLM*cell_gradp(theRiemann->cR,URR);
        if(fabs(dv[6]) > PLM*fabs(cell_gradp(theRiemann->cL,URR)))
            dv[6] = PLM*cell_gradp(theRiemann->cL,URR);
        if(cell_gradp(theRiemann->cL,URR)*cell_gradp(theRiemann->cR,URR) < 0.0)
            dv[6] = 0.0;
        if(fabs(dv[7]) > PLM*fabs(cell_gradp(theRiemann->cR,UPP)))
            dv[7] = PLM*cell_gradp(theRiemann->cR,UPP);
        if(fabs(dv[7]) > PLM*fabs(cell_gradp(theRiemann->cL,UPP)))
            dv[7] = PLM*cell_gradp(theRiemann->cL,UPP);
        if(cell_gradp(theRiemann->cL,UPP)*cell_gradp(theRiemann->cR,UPP) < 0.0)
            dv[7] = 0.0;
        if(fabs(dv[8]) > PLM*fabs(cell_gradp(theRiemann->cR,UZZ)))
            dv[8] = PLM*cell_gradp(theRiemann->cR,UZZ);
        if(fabs(dv[8]) > PLM*fabs(cell_gradp(theRiemann->cL,UZZ)))
            dv[8] = PLM*cell_gradp(theRiemann->cL,UZZ);
        if(cell_gradp(theRiemann->cL,UZZ)*cell_gradp(theRiemann->cR,UZZ) < 0.0)
            dv[8] = 0.0;
    */
    }
    else
    {
        double idz = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[9] = idz * (cell_prim(theRiemann->cR,URR) - cell_prim(theRiemann->cL,URR));
        dv[10] = idz * (cell_prim(theRiemann->cR,UPP) - cell_prim(theRiemann->cL,UPP));
        dv[11] = idz * (cell_prim(theRiemann->cR,UZZ) - cell_prim(theRiemann->cL,UZZ));
      /*  
        double PLM = sim_PLM(theSim);
        if(fabs(dv[9]) > PLM*fabs(cell_gradz(theRiemann->cR,URR)))
            dv[9] = PLM*cell_gradz(theRiemann->cR,URR);
        if(fabs(dv[9]) > PLM*fabs(cell_gradz(theRiemann->cL,URR)))
            dv[9] = PLM*cell_gradz(theRiemann->cL,URR);
        if(cell_gradz(theRiemann->cL,URR)*cell_gradz(theRiemann->cR,URR) < 0.0)
            dv[9] = 0.0;
        if(fabs(dv[10]) > PLM*fabs(cell_gradz(theRiemann->cR,UPP)))
            dv[10] = PLM*cell_gradz(theRiemann->cR,UPP);
        if(fabs(dv[10]) > PLM*fabs(cell_gradz(theRiemann->cL,UPP)))
            dv[10] = PLM*cell_gradz(theRiemann->cL,UPP);
        if(cell_gradz(theRiemann->cL,UPP)*cell_gradz(theRiemann->cR,UPP) < 0.0)
            dv[10] = 0.0;
        if(fabs(dv[11]) > PLM*fabs(cell_gradz(theRiemann->cR,UZZ)))
            dv[11] = PLM*cell_gradz(theRiemann->cR,UZZ);
        if(fabs(dv[11]) > PLM*fabs(cell_gradz(theRiemann->cL,UZZ)))
            dv[11] = PLM*cell_gradz(theRiemann->cL,UZZ);
        if(cell_gradz(theRiemann->cL,UZZ)*cell_gradz(theRiemann->cR,UZZ) < 0.0)
            dv[11] = 0.0;
    */
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

    metric_shear_uu(g, v, dv, shear, theSim);
    
    //TODO: CHECK THIS!  Probably not relativistic and/or isothermal.
    rhoh = prim[RHO] + GAMMA*prim[PPP]/(GAMMA-1.0);
    cs2 = GAMMA*prim[PPP] / rhoh;
    height = 2 * sqrt(prim[PPP]*r*r*r*(1-3*M/r)/(rhoh*M));
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

void riemann_visc_flux_LR(struct Riemann *theRiemann, struct Sim *theSim, int state, double *F)
{
    int i,j,k,dir;
    double a, sqrtg, du0, hn, cs2, rhoh, visc, height, r;
    double v[3], dv[12], b[3], u[4], du[16], shear[16];
    
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
    height = 2*sqrt(prim[PPP]*r*r*r*(1-3*M/r)/(rhoh*M));
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

