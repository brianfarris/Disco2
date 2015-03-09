#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/Tsolve.h"
#include "../Headers/header.h"

void cell_prim2cons( double * prim , double * cons , double r , double phi, double dV ,struct Sim * theSim, struct GravMass * theGravMasses){
    double GAMMALAW = sim_GAMMALAW(theSim);

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[URR];
    double vp  = prim[UPP]*r;
    double vz  = prim[UZZ];
    double v2  = vr*vr + vp*vp + vz*vz;

    double xi_g = 0.1;
    double xi_r_times_xi_g_pow_8 = 1.e-1;
    double xi_r =  xi_r_times_xi_g_pow_8*pow(xi_g,-8);

    int NumGravMass = sim_NumGravMass(theSim);
    double Omega_eff = gravMass_Omega_eff(r,phi,NumGravMass,theGravMasses); 
    struct Tsolve * params = tsolve_create(Pp,rho,GAMMALAW,xi_r,Omega_eff,0.5);
    double Temp = tsolve_findroot_P(params);
    tsolve_destroy(params);
    double Pgas = rho*Temp;
    double Prad = Pp - Pgas;
    double rhoe = Pgas/(GAMMALAW-1) + Prad*3.;

    //double rhoe = Pp/(GAMMALAW - 1.);


    cons[DDD] = rho*dV;
    cons[SRR] = rho*vr*dV;
    cons[LLL] = r*rho*vp*dV; // note that the conserved variable is angular momentum, not linear momentum
    cons[SZZ] = rho*vz*dV;
    cons[TAU] = ( .5*rho*v2 + rhoe )*dV;

    int q;
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
        cons[q] = prim[q]*cons[DDD];
    }

}

void cell_calc_cons( struct Cell *** theCells,struct Sim *theSim , struct GravMass * theGravMasses){
    double GAMMALAW = sim_GAMMALAW(theSim);

    int i,j,k;
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double dz = zp-zm;
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
            double rp = sim_FacePos(theSim,i,R_DIR);
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double r = .5*(rp+rm);
            for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
                struct Cell * c = &(theCells[k][i][j]);
                double phi = c->tiph-.5*c->dphi;                                
                double dV = .5*(rp*rp - rm*rm)*c->dphi*dz;
                cell_prim2cons( c->prim , c->cons , r ,phi, dV,theSim,theGravMasses);
            }
        }    
    }
}

void cell_cons2prim( double * cons , double * prim , double r , double phi, double dV ,struct Sim * theSim, struct GravMass * theGravMasses){
    double CS_FLOOR = sim_CS_FLOOR(theSim);
    double CS_CAP = sim_CS_CAP(theSim);
    double RHO_FLOOR = sim_RHO_FLOOR(theSim);
    double VEL_CAP = sim_VEL_CAP(theSim);
    double GAMMALAW = sim_GAMMALAW(theSim);
    double rho = cons[DDD]/dV;
    if( rho < RHO_FLOOR ) {
        rho = RHO_FLOOR;
    }
    double Sr  = cons[SRR]/dV;
    double Sp  = cons[LLL]/dV/r;
    double Sz  = cons[SZZ]/dV;
    double E   = cons[TAU]/dV;

    double vr = Sr/rho;
    double vp = Sp/rho;
    double vz = Sz/rho;
    double KE = .5*( Sr*vr + Sp*vp + Sz*vz );
    double v_magnitude = sqrt(2.*KE/rho);

    double rhoe = E - KE;

    double xi_g = 0.1;
    double xi_r_times_xi_g_pow_8 = 1.e-1;
    double xi_r =  xi_r_times_xi_g_pow_8*pow(xi_g,-8);

    double Tguess;
    if (prim[PGAS]>1.e-12){
       Tguess = 100*prim[PGAS]/prim[RHO];
    } else{
        Tguess = 0.5;
    }
    int NumGravMass = sim_NumGravMass(theSim);
    double Omega_eff = gravMass_Omega_eff(r,phi,NumGravMass,theGravMasses); 
    struct Tsolve * params = tsolve_create(rhoe,rho,GAMMALAW,xi_r,Omega_eff,Tguess);
    double Temp = tsolve_findroot_E(params);
    tsolve_destroy(params);
    double Pgas = rho*Temp;
    double Ptot = rhoe/3. + (GAMMALAW-4./3.)/(GAMMALAW-1.)*Pgas;
    double Prad = Ptot - Pgas; 
    double radfrac = Prad/Ptot;
    double gasfrac = Pgas/Ptot;

    double Pp = Ptot;
    //double Pp = (GAMMALAW-1.)*rhoe;

    if( Pp < CS_FLOOR*CS_FLOOR*rho/GAMMALAW ) {
        Pp = CS_FLOOR*CS_FLOOR*rho/GAMMALAW;
        Prad = CS_FLOOR*CS_FLOOR*rho/GAMMALAW*radfrac;
        Pgas = CS_FLOOR*CS_FLOOR*rho/GAMMALAW*gasfrac;
    }
    if ( Pp > CS_CAP*CS_CAP*rho/GAMMALAW) {
        Pp = CS_CAP*CS_CAP*rho/GAMMALAW;
        Prad = CS_CAP*CS_CAP*rho/GAMMALAW*radfrac;
        Pgas = CS_CAP*CS_CAP*rho/GAMMALAW*gasfrac;
    }
    if ( v_magnitude>VEL_CAP ){
        vr = vr*VEL_CAP/v_magnitude;
        vp = vp*VEL_CAP/v_magnitude;
        vz = vz*VEL_CAP/v_magnitude;
    }
    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = vr;
    prim[UPP] = vp/r;
    prim[UZZ] = vz;
    prim[PRAD] = Prad;
    prim[PGAS] = Pgas;

    int q;
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
        prim[q] = cons[q]/cons[DDD];
    }

}

void cell_calc_prim( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses){
    int i,j,k;
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double dz = zp-zm;
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = .5*(rp+rm);
            for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
                struct Cell * c = cell_single(theCells,i,j,k);
                double phi = c->tiph-.5*c->dphi;                
                double dV = .5*(rp*rp-rm*rm)*c->dphi*dz;
                cell_cons2prim( c->cons , c->prim , r , phi, dV ,theSim, theGravMasses);
                //printf("%e %e %e %e\n",r,c->P_components[RAD],c->P_components[GAS],c->P_components[RAD]+c->P_components[GAS]);
            }
        }
    }
}



