#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

void cell_prim2cons_grdisc(double *prim, double *cons, double *pos, 
                            double dV, struct Sim *theSim)
{
    struct Metric *g;
    double a, b[3], sqrtg;
    double u0, u_d[4], U[4];
    double rho, v[3];
    double Pp, eps, rhoh, H, M;
    double r;
    int i,j;

    r = pos[R_DIR];

    //Get hydro primitives
    rho = prim[RHO];
    v[0] = prim[URR];
    v[1] = prim[UPP];
    v[2] = prim[UZZ];

    //Get ADM metric values
    g = metric_create(time_global, pos[R_DIR], pos[P_DIR], pos[Z_DIR], theSim);
    a = metric_lapse(g);
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g, i);
    sqrtg = metric_sqrtgamma(g) / r;
    for(i=0; i<4; i++)
        U[i] = metric_frame_U_u(g, i, theSim);
    M = sim_GravM(theSim);

    //Calculate 4-velocity
    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v) - metric_square3_u(g,v));
    u_d[0] = metric_g_dd(g,0,0)*u0 + metric_dot3_u(g,b,v)*u0;
    for(i=1; i<4; i++)
    {
        u_d[i] = 0.0;
        for(j=0; j<3; j++)
            u_d[i] += metric_gamma_dd(g,i-1,j) * (v[j] + b[j]);
        u_d[i] *= u0;
    }

    Pp = eos_ppp(prim, theSim);
    eps = eos_eps(prim, theSim);
    rhoh = rho + rho*eps + Pp;

    H = sqrt(r*r*r*Pp / (rhoh*M)) / u0;

    u0 = 1.0 / sqrt(-metric_g_dd(g,0,0) - 2*metric_dot3_u(g,b,v)
                    - metric_square3_u(g,v));

    cons[DDD] = a*sqrtg*u0 * rho * H * dV;
    cons[SRR] = a*sqrtg*rhoh*u0 * u_d[1] * H * dV;
    cons[LLL] = a*sqrtg*rhoh*u0 * u_d[2] * H * dV;
    cons[SZZ] = a*sqrtg*rhoh*u0 * u_d[3] * H * dV;
    cons[TAU] = a*sqrtg*(-rhoh*u0*(U[0]*u_d[0]+U[1]*u_d[1]+U[2]*u_d[2]
                        +U[3]*u_d[3]) - rho*u0 - U[0]*Pp) * H * dV;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        cons[i] = prim[i] * cons[DDD];

    metric_destroy(g);
}

void cell_cons2prim_grdisc(double *cons, double *prim, double *pos, double dV,
                            struct Sim *theSim)
{
    int mu, nu, i;
    double D, S[3], tau;
    double rho, T, h, P, eps, l, dpdr, dpdt, dedr, dedt, dhdr, dhdt;
    double w, dwdr, dwdt, H, dHdr, dHdt;
    double tmp_prim[5], res;
    double r, a, b[3], sqrtg, U[4], Ud[4], ud[4], M;
    double S2, SU, A, B, C;
    double f1, f2, df1dr, df1dt, df2dr, df2dt, detj;
    struct Metric *g;
    double tol = 1.0e-10;
    int maxIter = 100;

    r = pos[R_DIR];

    D = cons[RHO] / dV;
    tau = cons[TAU] / dV;
    S[0] = cons[SRR] / dV;
    S[1] = cons[LLL] / dV;
    S[2] = cons[SZZ] / dV;

    g = metric_create(time_global, pos[R_DIR], pos[P_DIR], pos[Z_DIR], theSim);

    //Get some metric quantities
    a = metric_lapse(g);
    sqrtg = metric_sqrtgamma(g) / r;
    for(i=0; i<3; i++)
        b[i] = metric_shift_u(g,i);
    for(mu=0; mu<4; mu++)
        U[mu] = metric_frame_U_u(g, mu, theSim);
    for(mu=0; mu<4; mu++)
    {
        Ud[mu] = 0;
        for(nu=0; nu<4; nu++)
            Ud[mu] += metric_g_dd(g,mu,nu) * U[nu];
    }
    M = sim_GravM(theSim);
    
    //Useful invariants
    S2 = metric_square3_d(g, S);
    SU = metric_dot3_d(g, S, &(Ud[1]));

    //Constants for NR scheme
    A = (tau + SU + D) / (a*U[0]*D);
    B = D / sqrtg;
    C = S2/(D*D);

    //Initial Guess: prim.
    for(i=0; i<5; i++)
        tmp_prim[i] = prim[i];
    rho = prim[RHO];
    T = prim[TTT];

    double rho1,rho2,rho3,rho4,T1,T2,T3,T4;
    rho1 = -1;
    rho2 = -1;
    rho3 = -1;
    rho4 = -1;
    T1 = -1;
    T2 = -1;
    T3 = -1;
    T4 = -1;

    i = 0;
    //2D Newton-Raphson to find rho and T
    do
    {
        rho4 = rho3;
        rho3 = rho2;
        rho2 = rho1;
        rho1 = rho;
        T4 = T3;
        T3 = T2;
        T2 = T1;
        T1 = T;

        //Update from previous guess.
        tmp_prim[RHO] = rho;
        tmp_prim[TTT] = T;

        //Calculate for this iteration.
        P = eos_ppp(tmp_prim, theSim);
        eps = eos_eps(tmp_prim, theSim);
        dedr = eos_depsdrho(tmp_prim, theSim);
        dedt = eos_depsdttt(tmp_prim, theSim);
        dpdr = eos_dpppdrho(tmp_prim, theSim);
        dpdt = eos_dpppdttt(tmp_prim, theSim);
        h = 1.0 + eps + P/rho;
        dhdr = dedr + (dpdr*rho-P)/(rho*rho);
        dhdt = dedt + dpdt/rho;

        l = S[2] / (D*h);
        w = sqrt(1.0+C/(h*h));
        dwdr = -C / (h*h*h * w) * dhdr;
        dwdt = -C / (h*h*h * w) * dhdt;

        H = sqrt(r*r*r*P / (rho*h*M)) * (a/w);
        dHdr = (-dwdr/w + 0.5*dpdr/P - 0.5/rho - 0.5*dhdr/h) * H; 
        dHdt = (-dwdt/w + 0.5*dpdt/P - 0.5*dhdt/h) * H; 
        
        //f-values
        f1 = (h*w*w - P/rho) / w - A;
        f2 = w*rho*H - B;
        printf("    %.12lg (%.12lg %.12lg %.12lg) %.12lg (%.12lg %.12lg)\n", 
                f1, h*w,-P/(rho*w),-A,f2,w*rho*H,-B);

        //jacobian
        df1dr = ((dhdr*w*w+2*h*w*dwdr-(dpdr*rho-P)/(rho*rho))*w
                    - (h*w*w-P/rho)*dwdr) / (w*w);
        df1dt = ((dhdt*w*w+2*h*w*dwdt-dpdt/rho)*w-(h*w*w-P/rho)*dwdt) / (w*w);
        df2dr = dwdr*rho*H + w*H + w*rho*dHdr;
        df2dt = dwdt*rho*H + w*rho*dHdt;
        detj = df1dr*df2dt - df1dt*df2dr;

        //NR update for rho and T
        rho += -( df2dt*f1 - df1dt*f2) / detj;
        T   += -(-df2dr*f1 + df1dr*f2) / detj;

        if(rho <= 0.0)
        {
            T = 0.5*tmp_prim[RHO]/(tmp_prim[RHO]-rho) * (T-tmp_prim[TTT]) + tmp_prim[TTT];
            rho = 0.5*tmp_prim[RHO];
        }
        else if(T <= 0.0)
        {
            rho = 0.5*tmp_prim[TTT]/(tmp_prim[TTT]-T) * (rho-tmp_prim[RHO]) + tmp_prim[RHO];
            T = 0.5*tmp_prim[TTT];
        }
        
        if(rho==rho2 && T==T2)
        {
            printf("2-Cycle: taking average\n");
            rho = 0.5*(rho1+rho2);
            T = 0.5*(T1+T2);
        }
        else if(rho==rho3 && T==T3)
        {
            printf("3-Cycle: taking average\n");
            rho = (rho1+rho2+rho3)/3.0;
            T = (T1+T2+T3)/3.0;
        }
        else if(rho==rho4 && T==T4)
        {
            printf("4-Cycle: taking average\n");
            rho = 0.25*(rho1+rho2+rho3+rho4);
            T = 0.25*(T1+T2+T3+T4);
        }

        // 2-norm of change in rho & T.
        res = sqrt((rho-tmp_prim[RHO])*(rho-tmp_prim[RHO])/(rho*rho)
                + (T-tmp_prim[TTT])*(T-tmp_prim[TTT])/(T*T));
        i++;

        printf("%d rho: %.20lg T: %.20lg res: %.20lg\n", i, rho, T, res);
    }
    while(res > tol && i < maxIter);

    if(i == maxIter)
    {
        printf("ERROR: NR failed to converge after %d iterations.  Res = %.12lg\n",
                maxIter, res);
        printf("r: %.12lg rho: %.12lg T: %.12lg\n", r, rho, T);
        printf("D: %.12lg tau: %.12lg SR: %.12lg SP: %.12lg SZ: %.12lg\n",
                D, tau, S[0], S[1], S[2]);
        printf("S2: %.12lg SU: %.12lg\n", S2, SU);
        printf("A: %.20lg B: %.20lg C: %.20lg\n", A, B, C);

        double tc[5];
        cell_prim2cons_grdisc(prim, tc, pos, dV, theSim);
        printf("rho: %.12lg T: %.12lg vr: %.12lg vp: %.12lg vz: %.12lg\n", 
                prim[RHO], prim[TTT], prim[URR], prim[UPP], prim[UZZ]);
        printf("D: %.12lg tau: %.12lg Sr: %.12lg Sp: %.12lg Sz: %.12lg\n", 
                tc[RHO]/dV, tc[TTT]/dV, tc[URR]/dV, tc[UPP]/dV, tc[UZZ]/dV);
    }
    else
    {
        printf("NR converged after %d iterations.  Res = %.12lg\n", i, res);
        printf("r: %.12lg rho: %.12lg T: %.12lg\n", r, rho, T);
        printf("D: %.12lg tau: %.12lg SR: %.12lg SP: %.12lg SZ: %.12lg\n",
                D, tau, S[0], S[1], S[2]);
        printf("S2: %.12lg SU: %.12lg\n", S2, SU);
        printf("A: %.20lg B: %.20lg C: %.20lg\n", A, B, C);

        double tc[5];
        cell_prim2cons_grdisc(prim, tc, pos, dV, theSim);
        printf("rho: %.12lg T: %.12lg vr: %.12lg vp: %.12lg vz: %.12lg\n", 
                prim[RHO], prim[TTT], prim[URR], prim[UPP], prim[UZZ]);
        printf("D: %.12lg tau: %.12lg Sr: %.12lg Sp: %.12lg Sz: %.12lg\n", 
                tc[RHO]/dV, tc[TTT]/dV, tc[URR]/dV, tc[UPP]/dV, tc[UZZ]/dV);
    }

    //Prim recovery
    if(rho < sim_RHO_FLOOR(theSim))
        rho = sim_RHO_FLOOR(theSim);

    //TODO: enforce CS FLOOR/CAP
   
    prim[RHO] = rho;
    prim[TTT] = T;

    P = eos_ppp(prim, theSim);
    eps = eos_eps(prim, theSim);
    h = 1.0 + eps + P/rho;

    //Recover v^i
    w = sqrt(1.0+C/(h*h));
    for(i=0; i<3; i++)
        ud[i+1] = S[i] / (D*h);
    ud[0] = (w/a - metric_g_uu(g,0,1)*ud[1] - metric_g_uu(g,0,2)*ud[2]
            - metric_g_uu(g,0,3)*ud[3]) / metric_g_uu(g,0,0);

    prim[URR] = 0;
    prim[UPP] = 0;
    prim[UZZ] = 0;
    for(mu = 0; mu < 4; mu++)
    {
        prim[URR] += metric_g_uu(g,1,mu)*ud[mu];
        prim[UPP] += metric_g_uu(g,2,mu)*ud[mu];
        prim[UZZ] += metric_g_uu(g,3,mu)*ud[mu];
    }
    prim[URR] /= w/a;
    prim[UPP] /= w/a;
    prim[UZZ] /= w/a;

    for(i = sim_NUM_C(theSim); i < sim_NUM_Q(theSim); i++)
        prim[i] = cons[i]/cons[DDD];

    metric_destroy(g);
}

