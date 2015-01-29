#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

// this routine is only called by riemann_set_vel.
// It is used to find various L/R quantities. 
void LR_speed_grdisc(double *prim, double r, int *n, double *p_vn, 
                    double *p_cf2, double *Fm, double *p_mn, 
                    struct Sim *theSim)
{
    double vn = prim[URR]*n[0] + prim[UPP]*n[1] + prim[UZZ]*n[2];
    double cf2 = eos_cs2(prim, theSim);

    *p_vn = vn;
    *p_cf2 = cf2;

    //TODO: These are only used for HLLC and are WRONG.
    *p_mn = 0;
    Fm[0] = 0;
    Fm[1] = 0;
    Fm[2] = 0;
}

// Find velocities needed for the Riemann problem
void riemann_set_vel_grdisc(struct Riemann *theRiemann, struct Sim *theSim, 
                            double r, double GAMMALAW)
{
    int i, dir;
    double Sl, Sr, SlL, SrL, SlR, SrR;
    double a, b[3], bn, sig, gam, w, v[3], dv;
    struct Metric *g;

    g = metric_create(time_global, theRiemann->pos[R_DIR], 
                theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    //Get needed Metric quantities.
    dir = -1;
    a = metric_lapse(g);
    for(i=0; i<3; i++)
    {
        b[i] = metric_shift_u(g,i);
        if(theRiemann->n[i] != 0)
            dir = i;
    }
    bn = b[dir];
    gam = metric_gamma_uu(g, dir, dir);

    //Leftside wave speeds
    double vnL, cf2L, mnL;
    double FmL[3];
    LR_speed_grdisc(theRiemann->primL, r, theRiemann->n, &vnL, &cf2L, FmL, 
                    &mnL, theSim);
    v[0] = theRiemann->primL[URR];
    v[1] = theRiemann->primL[UPP];
    v[2] = theRiemann->primL[UZZ];
    w = a / sqrt(-metric_g_dd(g,0,0) - 2.0*metric_dot3_u(g,b,v) 
                - metric_square3_u(g,v));
    sig = cf2L / (w*w*(1.0-cf2L));
    dv = sqrt(sig*(1.0+sig)*a*a*gam - sig*(vnL+bn)*(vnL+bn));

    SlL = (vnL - sig*bn - dv) / (1.0 + sig);
    SrL = (vnL - sig*bn + dv) / (1.0 + sig);

    //Rightside wave speeds
    double vnR, cf2R, mnR;
    double FmR[3];
    LR_speed_grdisc(theRiemann->primR, r, theRiemann->n, &vnR, &cf2R, FmR, 
                    &mnR, theSim);
    v[0] = theRiemann->primR[URR];
    v[1] = theRiemann->primR[UPP];
    v[2] = theRiemann->primR[UZZ];
    w = a / sqrt(-metric_g_dd(g,0,0) - 2.0*metric_dot3_u(g,b,v) 
                - metric_square3_u(g,v));
    sig = cf2R / (w*w*(1.0-cf2R));
    dv = sqrt(sig*(1.0+sig)*a*a*gam - sig*(vnR+bn)*(vnR+bn));

    SlR = (vnR - sig*bn - dv) / (1.0 + sig);
    SrR = (vnR - sig*bn + dv) / (1.0 + sig);

    //Assign largest wave speeds.
    Sl = (SlL < SlR) ? SlL : SlR;
    Sr = (SrL > SrR) ? SrL : SrR;

    //Fluxes are evaluated in orthonormal basis of coordinates.
    if(dir == PDIRECTION)
    {
        Sl *= r;
        Sr *= r;
    }

    theRiemann->Sl = Sl;
    theRiemann->Sr = Sr;
    //TODO: This is only used for HLLC and is WRONG.
    theRiemann->Ss = (Sl+Sr)/2;

    metric_destroy(g);
}

void riemann_set_flux_grdisc(struct Riemann *theRiemann, struct Sim *theSim, double GAMMALAW, int SetState)
{
    double r = theRiemann->pos[R_DIR];
    double *prim;
    double *F;

    if(SetState==LEFT)
    {
        prim = theRiemann->primL;
        F = theRiemann->FL;
    }
    else if(SetState==RIGHT)
    {
        prim = theRiemann->primR;
        F = theRiemann->FR;
    }
    else
    {
        printf("ERROR: riemann_set_flux given unrecognized state.\n");
        exit(0);
    }

    int i,j;
    struct Metric *g;
    double a, b[3], sqrtg, U[4];
    double u0, u_d[4];
    double rho, v[3];
    double Pp, eps, rhoh, vn, bn, hn, Un, M, H;

    //Get Hydro primitives
    rho = prim[RHO];
    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    //Get needed metric values
    //TODO: Make Metric obj part of Riemann.
    g = metric_create(time_global, theRiemann->pos[R_DIR], 
                    theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g,i);
    sqrtg = metric_sqrtgamma(g) / r;
    for(i=0; i<4; i++)
        U[i] = metric_frame_U_u(g, i, theSim);

    //Calculate 4-velocity
    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v)
                    - metric_square3_u(g,v));
    u_d[0] = metric_g_dd(g,0,0)*u0 + metric_dot3_u(g,b,v)*u0;
    for(i=1; i<4; i++)
    {
        u_d[i] = 0.0;
        for(j=0; j<3; j++)
            u_d[i] += metric_gamma_dd(g,i-1,j) * (v[j] + b[j]);
        u_d[i] *= u0;
    }

    if(-metric_g_dd(g,0,0)-2*metric_dot3_u(g,b,v)-metric_square3_u(g,v) < 0)
            printf("ERROR: Velocity too high in flux. r=%.12g, vr=%.12g, vp=%.12g, vz=%.12g\n", r, v[0], v[1], v[2]);

    //Calculate beta & v normal to face.
    vn = v[0]*theRiemann->n[0] + v[1]*theRiemann->n[1] + v[2]*theRiemann->n[2];
    bn = b[0]*theRiemann->n[0] + b[1]*theRiemann->n[1] + b[2]*theRiemann->n[2];
    Un = U[1]*theRiemann->n[0] + U[2]*theRiemann->n[1] + U[3]*theRiemann->n[2];
    if(theRiemann->n[1] == 1)
        hn = r;
    else
        hn = 1.0;

    //Thermal things
    Pp = eos_ppp(prim, theSim);
    eps = eos_eps(prim, theSim);
    rhoh = rho*(1.0 + eps) + Pp;
    M = sim_GravM(theSim);
    H = sqrt(r*r*r*Pp / (rhoh*M)) / u0;

    //The Fluxes.
    F[DDD] = hn*sqrtg*a*u0 * rho*vn * H;
    F[SRR] = hn*sqrtg*a*(u0*rhoh*u_d[1]*vn + Pp*theRiemann->n[0])*H;
    F[LLL] = hn*sqrtg*a*(u0*rhoh*u_d[2]*vn + Pp*theRiemann->n[1])*H;
    F[SZZ] = hn*sqrtg*a*(u0*rhoh*u_d[3]*vn + Pp*theRiemann->n[2])*H;
    F[TAU] = hn*sqrtg*a*(-u0*rhoh*(U[0]*u_d[0]+U[1]*u_d[1]+U[2]*u_d[2]
                        +U[3]*u_d[3])*vn - u0*rho*vn - Un*Pp) * H;

    //Passive Fluxes.
    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        F[i] = prim[i] * F[DDD];

    metric_destroy(g);
}

void riemann_visc_flux_grdisc(struct Riemann *theRiemann, struct Sim *theSim)
{
    int i,j,dir;
    double a, sqrtg, hn, cs2, Pp, eps, rhoh, visc, r, u0, H;
    double v[3], dv[12], b[3], shear[16], U[4], Ud[4];
    
    double alpha = sim_AlphaVisc(theSim);
    double M = sim_GravM(theSim);
    int NUMQ = sim_NUM_Q(theSim);
    double *prim = (double *)malloc(NUMQ * sizeof(double));
    double *F = (double *)malloc(NUMQ * sizeof(double));
    struct Metric *g = metric_create(time_global, theRiemann->pos[R_DIR], 
            theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], theSim);

    r = theRiemann->pos[R_DIR];

    //Face centered prims averaged from L/R.
    for(i=0; i<NUMQ; i++)
        prim[i] = 0.5*(theRiemann->primL[i] + theRiemann->primR[i]);

    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    dir = -1;
    for(i=0; i<3; i++)
        if(theRiemann->n[i] == 1)
            dir = i;
    if(dir == -1)
        printf("ERROR: n[i] != 1 in riemann_visc_flux_grdisc\n");

    //v-derivatives normal to face are averages of L/R gradients.
    // dv/dt = 0 in viscous terms... for now.
    dv[0] = 0.0;
    dv[1] = 0.0;
    dv[2] = 0.0;
    // dv/dr
    dv[3] = 0.5*(cell_gradr(theRiemann->cL, URR)
               + cell_gradr(theRiemann->cR, URR));
    dv[4] = 0.5*(cell_gradr(theRiemann->cL, UPP) 
               + cell_gradr(theRiemann->cR, UPP));
    dv[5] = 0.5*(cell_gradr(theRiemann->cL, UZZ) 
               + cell_gradr(theRiemann->cR, UZZ));
    // dv/dphi
    dv[6] = 0.5*(cell_gradp(theRiemann->cL, URR) 
               + cell_gradp(theRiemann->cR, URR));
    dv[7] = 0.5*(cell_gradp(theRiemann->cL, UPP) 
               + cell_gradp(theRiemann->cR, UPP));
    dv[8] = 0.5*(cell_gradp(theRiemann->cL, UZZ) 
               + cell_gradp(theRiemann->cR, UZZ));
    // dv/dz
    dv[9] = 0.5*(cell_gradz(theRiemann->cL, URR) 
               + cell_gradz(theRiemann->cR, URR));
    dv[10] = 0.5*(cell_gradz(theRiemann->cL, UPP) 
                + cell_gradz(theRiemann->cR, UPP));
    dv[11] = 0.5*(cell_gradz(theRiemann->cL, UZZ) 
                + cell_gradz(theRiemann->cR, UZZ));

    //Derivatives across face are calculated directly using 2nd order stencil.
    if(dir == RDIRECTION)
    {
        double dphiR = theRiemann->pos[P_DIR] - (cell_tiph(theRiemann->cR)
                                            -0.5*cell_dphi(theRiemann->cR));
        double dphiL = theRiemann->pos[P_DIR] - (cell_tiph(theRiemann->cL)
                                            -0.5*cell_dphi(theRiemann->cL));
        while(dphiR<-M_PI) dphiR += 2*M_PI;
        while(dphiR> M_PI) dphiR -= 2*M_PI;
        while(dphiL<-M_PI) dphiL += 2*M_PI;
        while(dphiL> M_PI) dphiL -= 2*M_PI;
        double idr = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);


        dv[3] = idr * (cell_prim(theRiemann->cR,URR)
                        + cell_gradp(theRiemann->cR,URR)*dphiR
                    - cell_prim(theRiemann->cL,URR)
                        - cell_gradp(theRiemann->cL,URR)*dphiL);
        dv[4] = idr * (cell_prim(theRiemann->cR,UPP)
                        + cell_gradp(theRiemann->cR,UPP)*dphiR 
                    - cell_prim(theRiemann->cL,UPP)
                        - cell_gradp(theRiemann->cL,UPP)*dphiL);
        dv[5] = idr * (cell_prim(theRiemann->cR,UZZ) 
                        + cell_gradp(theRiemann->cR,UZZ)*dphiR 
                    - cell_prim(theRiemann->cL,UZZ)
                        - cell_gradp(theRiemann->cL,UZZ)*dphiL);
    }
    else if (dir == PDIRECTION)
    {
        double dphi = cell_tiph(theRiemann->cR)-0.5*cell_dphi(theRiemann->cR) 
                    - cell_tiph(theRiemann->cL)+0.5*cell_dphi(theRiemann->cL);
        while (dphi < -M_PI) dphi += 2*M_PI;
        while (dphi >  M_PI) dphi -= 2*M_PI;
        double idp = 1.0/dphi;

        dv[6] = idp * (cell_prim(theRiemann->cR,URR)
                        - cell_prim(theRiemann->cL,URR));
        dv[7] = idp * (cell_prim(theRiemann->cR,UPP) 
                        - cell_prim(theRiemann->cL,UPP));
        dv[8] = idp * (cell_prim(theRiemann->cR,UZZ) 
                        - cell_prim(theRiemann->cL,UZZ)); 
    }
    //TODO: FIX THIS.  Doesn't account for misaligned zones in z faces.
    else
    {
        double idz = 1.0/(theRiemann->x_cell_R - theRiemann->x_cell_L);
        dv[9] = idz * (cell_prim(theRiemann->cR,URR)
                        - cell_prim(theRiemann->cL,URR));
        dv[10] = idz * (cell_prim(theRiemann->cR,UPP)
                        - cell_prim(theRiemann->cL,UPP));
        dv[11] = idz * (cell_prim(theRiemann->cR,UZZ)
                        - cell_prim(theRiemann->cL,UZZ));
    }
    
    if(PRINTTOOMUCH && 0)
    {
        FILE *gradfile = fopen("grad_face.out","a");
        fprintf(gradfile, "%d, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g, %.12g\n", dir, theRiemann->pos[R_DIR], theRiemann->pos[P_DIR], theRiemann->pos[Z_DIR], v[0], v[1], v[2], dv[3], dv[4], dv[5], dv[6], dv[7], dv[8], dv[9], dv[10], dv[11]);
        fclose(gradfile);
    }
   
    //Metric quantities
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g) / r;
    u0 = 1.0/sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v)
                - metric_square3_u(g,v));
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

    //Shear tensor!
    metric_shear_uu(g, v, dv, shear, theSim);

    //No z-shear if we're in 2D.
    if(sim_N(theSim,Z_DIR)==1)
    {
        for(i=0; i<4; i++)
        {
            shear[4*i+3] = 0.0;
            shear[4*3+i] = 0.0;
        }
    }
   
    //Thermal quantities
    Pp = eos_ppp(prim, theSim);
    eps = eos_eps(prim, theSim);
    rhoh = prim[RHO]*(1.0 + eps) + Pp;
    cs2 = eos_cs2(prim, theSim);
    H = sqrt(r*r*r*Pp / (rhoh*M)) / u0;
    if(alpha > 0)
        visc = alpha * sqrt(cs2) * H * prim[RHO];
    else
        visc = -alpha * prim[RHO];

    if(dir == PDIRECTION)
        hn = r;
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
    F[SRR] *= -a*sqrtg*hn * visc * H;
    F[LLL] *= -a*sqrtg*hn * visc * H;
    F[SZZ] *= -a*sqrtg*hn * visc * H;
    F[TAU] *=  a*sqrtg*hn * visc * H;

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

