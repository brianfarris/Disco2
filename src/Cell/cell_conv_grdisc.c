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

//Local Functions
static int cons2prim_prep(double *cons, double *prim, double *pos, double dV, 
                        struct Metric *g, struct Sim *theSim);
static int cons2prim_solve(double *cons, double *prim, double *pos, double dV, 
                        struct Metric *g, struct Sim *theSim);

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

    if(sim_Metric(theSim) == KERR_KS)
    {
        double A = M*sim_GravA(theSim);

        H = r*r * sqrt(Pp / (rhoh*(u_d[2]*u_d[2]-A*A*(u_d[0]*u_d[0]-1.0))));
    }
    else
        H = sqrt(r*r*r*Pp / (rhoh*M)) / u0;

    cons[DDD] = a*sqrtg * u0 * rho * H * dV;
    cons[SRR] = a*sqrtg * rhoh*u0 * u_d[1] * H * dV;
    cons[LLL] = a*sqrtg * rhoh*u0 * u_d[2] * H * dV;
    cons[SZZ] = a*sqrtg * rhoh*u0 * u_d[3] * H * dV;
    cons[TAU] = a*sqrtg * (-rhoh*u0*(U[0]*u_d[0]+U[1]*u_d[1]+U[2]*u_d[2]
                            +U[3]*u_d[3]) - rho*u0 - U[0]*Pp) * H * dV;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        cons[i] = prim[i] * cons[DDD];

    metric_destroy(g);
}

void cell_cons2prim_grdisc(double *cons, double *prim, double *pos, double dV,
                            struct Sim *theSim)
{
    int err;
    struct Metric *g;
    g = metric_create(time_global, pos[R_DIR], pos[P_DIR], pos[Z_DIR], theSim);

    err = cons2prim_prep(cons, prim, pos, dV, g, theSim);
    err = cons2prim_solve(cons, prim, pos, dV, g, theSim);

    if(err)
    {
        printf("cons:  %.12lg %.12lg %.12lg %.12lg %.12lg\n", 
                cons[DDD], cons[TAU], cons[SRR], cons[LLL], cons[SZZ]);
        printf("prim:  %.12lg %.12lg %.12lg %.12lg %.12lg\n", 
                prim[RHO], prim[TTT], prim[URR], prim[UPP], prim[UZZ]);
        if(prim[RHO] < sim_RHO_FLOOR(theSim))
        {
            prim[RHO] = sim_RHO_FLOOR(theSim);
            printf("Floored RHO.\n.");
        }
        if(prim[TTT] < sim_CS_FLOOR(theSim))
        {
            prim[TTT] = sim_CS_FLOOR(theSim);
            printf("Floored TTT.\n.");
        }
        printf("prim2: %.12lg %.12lg %.12lg %.12lg %.12lg\n", 
                prim[RHO], prim[TTT], prim[URR], prim[UPP], prim[UZZ]);
        
        cell_prim2cons_grdisc(prim, cons, pos, dV, theSim);
        
        printf("cons2: %.12lg %.12lg %.12lg %.12lg %.12lg\n", 
                cons[DDD], cons[TAU], cons[SRR], cons[LLL], cons[SZZ]);
    }

    metric_destroy(g);
}

static int cons2prim_prep(double *cons, double *prim, double *pos, double dV, 
                        struct Metric *g, struct Sim *theSim)
{
    // Enforces the inequalities E^2 > rho^2 + S^2 and E > rho. This *should*
    // guarantee a successful prim recovery in cons2prim_solve.  This 
    // inequality is not always true, depends on the EoS.
    // TODO: Check inequality w/ Gas,Rad,Deg EoSes.

    int err = 0;
    int q;
    int NUMQ = sim_NUM_Q(theSim);

    double D, tau, S[3], S2, SU, U[4], Ud[3], E, W;
    int i,j;

    D = cons[DDD]/dV;
    tau = cons[TAU]/dV;
    S[0] = cons[SRR]/dV;
    S[1] = cons[LLL]/dV;
    S[2] = cons[SZZ]/dV;

    S2 = metric_square3_d(g, S);
    for(i=0; i<4; i++)
        U[i] = metric_frame_U_u(g, i, theSim);
    for(i=0; i<3; i++)
    {
        Ud[i] = 0;
        for(j=0; j<4; j++)
            Ud[i] += metric_g_dd(g,i+1,j)*U[j];
    }

    SU = metric_dot3_d(g, S, Ud);
    W = metric_lapse(g) * U[0];

    E = (tau + SU + D)/W;

    if(S2 + D*D >= E*E || E < D)
    {
    //    if(pos[R_DIR] >= metric_horizon(theSim))
            printf("2 fast 2 cold. r = %.12lg\n", pos[R_DIR]);
        err = 1;
        E = 1.01 * D * (sqrt(1.0+S2/(D*D)));
        cons[TAU] = (W*E - SU - D) * dV;
    }

    int nan = 0;
    for(q=0; q<NUMQ; q++)
        if(cons[q] != cons[q])
            nan = 1;
    if(nan)
    {
        printf("ERROR: NaN in cons2prim. r=%.12g\n", pos[R_DIR]);
        printf("   DDD=%.12g TAU=%.12g SRR=%.12g LLL=%.12g SZZ=%.12g\n",
                 cons[DDD], cons[TAU], cons[SRR], cons[LLL], cons[SZZ]);
        err = 1;
    }
    
    return err;
}

static int cons2prim_solve(double *cons, double *prim, double *pos, double dV, 
                        struct Metric *g, struct Sim *theSim)
{
    int mu, nu, i;
    double D, S[3], tau;
    double rho, T, h, P, eps, dpdr, dpdt, dedr, dedt, dhdr, dhdt;
    double w, dwdr, dwdt, H, dHdr, dHdt;
    double tmp_prim[5], res;
    double r, a, b[3], sqrtg, U[4], Ud[4], ud[4], M;
    double S2, SU, A, B, C;
    double f1, f2, df1dr, df1dt, df2dr, df2dt, detj;
    int err = 0;
    double tol = 1.0e-10;
    int maxIter = 100;

    r = pos[R_DIR];

    D = cons[RHO] / dV;
    tau = cons[TAU] / dV;
    S[0] = cons[SRR] / dV;
    S[1] = cons[LLL] / dV;
    S[2] = cons[SZZ] / dV;

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

        w = sqrt(1.0+C/(h*h));
        dwdr = -C / (h*h*h * w) * dhdr;
        dwdt = -C / (h*h*h * w) * dhdt;

        if(sim_Metric(theSim) == KERR_KS)
        {
            double AA, l2, dl2dr, dl2dt, e2, de2dr, de2dt, guu;
            AA = M*sim_GravA(theSim);
            l2 = S[1]*S[1] / (D*D*h*h);
            dl2dr = -2*l2/h*dhdr;
            dl2dt = -2*l2/h*dhdt;
            guu = (metric_g_uu(g,1,1)*S[0]*S[0] + metric_g_uu(g,2,2)*S[1]*S[1]
                    + metric_g_uu(g,3,3)*S[2]*S[2]) / (D*D*h*h);
            e2 = (-1.0 - guu) / metric_g_uu(g,0,0);
            de2dr = 2*guu / ( h * metric_g_uu(g,0,0)) * dhdr;
            de2dt = 2*guu / ( h * metric_g_uu(g,0,0)) * dhdt;

            H = r*r*sqrt(P/(rho*h*(l2 - AA*AA*(e2-1.0))));
            dHdr = (-0.5*(dl2dr-AA*AA*de2dr)/(l2-AA*AA*(e2-1.0))
                    + 0.5*dpdr/P - 0.5/rho - 0.5*dhdr/h ) * H;
            dHdt = (-0.5*(dl2dt-AA*AA*de2dt)/(l2-AA*AA*(e2-1.0))
                    + 0.5*dpdt/P - 0.5*dhdt/h ) * H;
        }
        else
        {
            H = sqrt(r*r*r*P / (rho*h*M)) * (a/w);
            dHdr = (-dwdr/w + 0.5*dpdr/P - 0.5/rho - 0.5*dhdr/h) * H; 
            dHdt = (-dwdt/w + 0.5*dpdt/P - 0.5*dhdt/h) * H; 
        }
        
        //f-values
        //f1 = (h*w*w - P/rho) / w - A;
        f1 = (1.0+eps)*w + P*C/(h*h*rho*w) - A;
        f2 = w*rho*H - B;
        if(PRINTTOOMUCH)
            printf("    %.12lg (%.12lg %.12lg %.12lg) %.12lg (%.12lg %.12lg)\n", 
                f1, h*w,-P/(rho*w),-A,f2,w*rho*H,-B);

        //jacobian
        //df1dr = ((dhdr*w*w+2*h*w*dwdr-(dpdr*rho-P)/(rho*rho))*w
        //            - (h*w*w-P/rho)*dwdr) / (w*w);
        //df1dt = ((dhdt*w*w+2*h*w*dwdt-dpdt/rho)*w-(h*w*w-P/rho)*dwdt) / (w*w);
        df1dr = dedr*w + (1+eps)*dwdr + C * (dpdr*h*rho*w - P*(2*dhdr*rho*w
                    + h*w + h*rho*dwdr)) / (h*h*h*rho*rho*w*w);
        df1dt = dedt*w + (1+eps)*dwdt + C * (dpdt*h*w - P*(2*dhdt*w + h*dwdt))
                                            / (h*h*h*rho*w*w);
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
        if(T <= 0.0)
        {
            rho = 0.5*tmp_prim[TTT]/(tmp_prim[TTT]-T) * (rho-tmp_prim[RHO]) + tmp_prim[RHO];
            T = 0.5*tmp_prim[TTT];
        }
        
        
        if(rho==rho2 && T==T2)
        {
            //if(PRINTTOOMUCH)
                printf("2-Cycle: taking average\n");
            rho = 0.5*(rho1+rho2);
            T = 0.5*(T1+T2);
        }
        else if(rho==rho3 && T==T3)
        {
            //if(PRINTTOOMUCH)
                printf("3-Cycle: taking average\n");
            rho = (rho1+rho2+rho3)/3.0;
            T = (T1+T2+T3)/3.0;
        }
        else if(rho==rho4 && T==T4)
        {
            //if(PRINTTOOMUCH)
                printf("4-Cycle: taking average\n");
            rho = 0.25*(rho1+rho2+rho3+rho4);
            T = 0.25*(T1+T2+T3+T4);
        }

        // 2-norm of change in rho & T.
        res = sqrt((rho-tmp_prim[RHO])*(rho-tmp_prim[RHO])/(rho*rho)
                + (T-tmp_prim[TTT])*(T-tmp_prim[TTT])/(T*T));
        i++;

        if(PRINTTOOMUCH)
            printf("%d rho: %.20lg T: %.20lg res: %.20lg\n", i, rho, T, res);
    }
    while(res > tol && i < maxIter);

    if(i == maxIter)
    {
        err = 1;

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
    else if(PRINTTOOMUCH)
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

    return err;
}

