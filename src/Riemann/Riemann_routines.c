#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/GravMass.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/Face.h"
#include "../Headers/header.h"

// this routine is only called by riemann_set_vel.
// It is used to find sound speed and normal velocity. 
void LR_speed(double *prim,double r,int * n,double GAMMALAW,double * p_vn,double * p_cf2){
    double P   = prim[PPP];
    double rho = prim[RHO];
    double vr  =   prim[URR];
    double vp  = r*prim[UPP];
    double vz  =   prim[UZZ];
    double vn  = vr*n[0] + vp*n[1] + vz*n[2];
    double cf2  = GAMMALAW*(P/rho);

    //normal velocity
    *p_vn = vn;
    //sound speed squared
    *p_cf2 = cf2;
}

// Find velocities needed for the Riemann problem
void riemann_set_vel(struct Riemann * theRiemann,struct Sim * theSim,struct GravMass * theGravMasses, double r,double GAMMALAW){
    double Sl, Sr, Ss;

    double vnL,cf2L;
    // get sound speed squared, normal velocity on L face
    LR_speed(theRiemann->primL,r,theRiemann->n,GAMMALAW,&vnL,&cf2L);
    Sl = vnL - sqrt( cf2L );
    Sr = vnL + sqrt( cf2L );

    // get sound speed squared, normal velocity on L face
    double vnR,cf2R;
    LR_speed(theRiemann->primR,r,theRiemann->n,GAMMALAW,&vnR,&cf2R);

    //See Eqns 10.47 and 10.48 of Toro
    if( Sl > vnR - sqrt( cf2R ) ) Sl = vnR - sqrt( cf2R );
    if( Sr < vnR + sqrt( cf2R ) ) Sr = vnR + sqrt( cf2R );

    // See Eq 10.37 of Toro
    Ss = (theRiemann->primR[RHO]*vnR*(Sr-vnR)-theRiemann->primL[RHO]*vnL*(Sl-vnL)+theRiemann->primL[PPP]-theRiemann->primR[PPP])
        /(theRiemann->primR[RHO]*(Sr-vnR)-theRiemann->primL[RHO]*(Sl-vnL));

    theRiemann->Sl = Sl;
    theRiemann->Sr = Sr;
    theRiemann->Ss = Ss;

}

// Which state of the riemann problem are we in?
// See Eq 10.26 of Toro
void riemann_set_state(struct Riemann * theRiemann,double w ){
    if (w < theRiemann->Sl){
        theRiemann->state=LEFT;
    }else if( w > theRiemann->Sr ){
        theRiemann->state=RIGHT;
    }else{
        if( w < theRiemann->Ss ){
            theRiemann->state=LEFTSTAR;
        }else{
            theRiemann->state=RIGHTSTAR;
        }
    }
}

// 
void riemann_set_star_hll(struct Riemann * theRiemann,struct Sim * theSim){
    double Sr =  theRiemann->Sr;
    double Sl = theRiemann->Sl;
    int q;
    // don't be confused by the fact that we refer to these values as "star". We just don't want to store more stuff in memory than we need to.
    for( q=0 ; q<sim_NUM_Q(theSim) ; ++q ){
        //See Eq 10.13 of Toro
        theRiemann->Ustar[q] = ( Sr*theRiemann->UR[q] - Sl*theRiemann->UL[q] + theRiemann->FL[q] - theRiemann->FR[q] )/( Sr - Sl );
        //See Eq 10.20 of Toro
        theRiemann->Fstar[q] = ( Sr*theRiemann->FL[q] - Sl*theRiemann->FR[q] - Sr*Sl*( theRiemann->UL[q] - theRiemann->UR[q] ) )/( Sr - Sl );
    }
}

void riemann_set_star_hllc(struct Riemann * theRiemann,struct Sim * theSim,double GAMMALAW){
    double r = theRiemann->r;
    double *prim;
    double Sk;
    double *Uk;
    double *Fk;
    // Set K values according to which state we are in
    if (theRiemann->state==LEFTSTAR){
        prim = theRiemann->primL;
        Sk = theRiemann->Sl;
        Uk = theRiemann->UL;
        Fk = theRiemann->FL;
    }else{
        prim = theRiemann->primR;
        Sk = theRiemann->Sr;
        Uk = theRiemann->UR;
        Fk = theRiemann->FR;
    }
    double Ss=theRiemann->Ss;

    // See Eq 10.39 of Toro
    double rho = prim[RHO];
    double vr  = prim[URR];
    double vp  = prim[UPP]*r;
    double vz  = prim[UZZ];
    double Pp  = prim[PPP];
    double v2 = vr*vr+vp*vp+vz*vz;
    double vn = vr*theRiemann->n[0] + vp*theRiemann->n[1] + vz*theRiemann->n[2];
    double rhoe = Pp/(GAMMALAW-1.);
    double D  = rho;
    double mr = rho*vr;
    double mp = rho*vp;
    double mz = rho*vz;
    double E_hydro  = .5*rho*v2 + rhoe;
    double Ps  = rho*( Sk - vn )*( Ss - vn ) + Pp;
    double Dstar = ( Sk - vn )*D/( Sk - Ss );
    double Msn = Dstar*Ss;
    double Msr   = ( Sk - vn )*mr / ( Sk - Ss );
    double Msp   = ( Sk - vn )*mp / ( Sk - Ss );
    double Msz   = ( Sk - vn )*mz / ( Sk - Ss );
    double Estar = ( ( Sk - vn )*E_hydro + Ps*Ss - Pp*vn ) / ( Sk - Ss );

    double mn  = Msr*theRiemann->n[0] + Msp*theRiemann->n[1] + Msz*theRiemann->n[2];

    Msr += theRiemann->n[0]*( Msn - mn );
    Msp += theRiemann->n[1]*( Msn - mn );
    Msz += theRiemann->n[2]*( Msn - mn );

    theRiemann->Ustar[DDD] = Dstar;
    theRiemann->Ustar[SRR] = Msr;
    theRiemann->Ustar[LLL] = r*Msp;
    theRiemann->Ustar[SZZ] = Msz;
    theRiemann->Ustar[TAU] = Estar;

    int q;
    
    //Now set Fstar
    // See Eq 10.38 of Toro
    for (q=0;q<sim_NUM_Q(theSim);++q){
        theRiemann->Fstar[q] = Fk[q] + Sk*( theRiemann->Ustar[q] - Uk[q] ) ;
    }

    // passive scalar stuff
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
        theRiemann->Ustar[q] = prim[q]*theRiemann->Ustar[DDD];
    }
}

// Here we compute the Fluxes on the L or R face, depending on what SetState is
void riemann_set_flux(struct Riemann * theRiemann, struct Sim * theSim,double GAMMALAW,int SetState){
    double r = theRiemann->r;
    double *prim;
    double *F;
    if (SetState==LEFT){
        prim = theRiemann->primL;
        F = theRiemann->FL;
    }else if (SetState==RIGHT){
        prim = theRiemann->primR;
        F = theRiemann->FR;
    } else{
        printf("ERROR\n");
        exit(0);
    }
    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double vr  = prim[URR];
    double vp  = prim[UPP]*r;
    double vz  = prim[UZZ];
    double vn = vr*theRiemann->n[0] + vp*theRiemann->n[1] + vz*theRiemann->n[2];
    double rhoe = Pp/(GAMMALAW-1.);
    double v2 = vr*vr + vp*vp + vz*vz;
    F[DDD] = rho*vn;
    F[SRR] =     rho*vr*vn + Pp*theRiemann->n[0] ;
    F[LLL] = r*( rho*vp*vn + Pp*theRiemann->n[1] );
    F[SZZ] =     rho*vz*vn + Pp*theRiemann->n[2] ;
    F[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn ;

    //passive scalar stuff
    int q;
    for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
        F[q] = prim[q]*F[DDD];
    }

}

void riemann_visc_flux(struct Riemann * theRiemann,struct Sim * theSim,struct GravMass * theGravMasses ){
    int NUM_Q = sim_NUM_Q(theSim);

    double r = theRiemann->r;
    double tiph = theRiemann->cm;

    double Mtotal = 1.0;
    double sep = 1.0;
    double M0 = gravMass_M(theGravMasses,0);
    double M1 = gravMass_M(theGravMasses,1);

    double r0 = gravMass_r(theGravMasses,0);
    double r1 = gravMass_r(theGravMasses,1);

    double a = r0 + r1;

    double phi_bh0 = gravMass_phi(theGravMasses,0);
    double phi_bh1 = gravMass_phi(theGravMasses,1);

    double xbh0 = r0*cos(phi_bh0);
    double ybh0 = r0*sin(phi_bh0);

    double xbh1 = r1*cos(phi_bh1);
    double ybh1 = r1*sin(phi_bh1);

    double xpos = r*cos(tiph);
    double ypos = r*sin(tiph);

    double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
    double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));

    double * VFlux = malloc(NUM_Q*sizeof(double));
    double * AvgPrim = malloc(NUM_Q*sizeof(double));
    double * Grad_r_prim = malloc(NUM_Q*sizeof(double));
    double * Grad_ph_prim = malloc(NUM_Q*sizeof(double));
    int q;
    for (q=0;q<NUM_Q;++q){
        AvgPrim[q] = .5*(theRiemann->primL[q]+theRiemann->primR[q]);
        //AvgPrim[q] = 0.5*(cell_prim(theRiemann->cL,q)+cell_prim(theRiemann->cR,q));
        Grad_ph_prim[q] = .5*(cell_gradp(theRiemann->cL,q)+cell_gradp(theRiemann->cR,q))/r;    
        Grad_r_prim[q] = .5*(cell_grad(theRiemann->cL,q)+cell_grad(theRiemann->cR,q));

        if (1==1){
            if (theRiemann->n[0]==1){
                double pL = cell_tiph(theRiemann->cL) - .5*cell_dphi(theRiemann->cL);
                double pR = cell_tiph(theRiemann->cR) - .5*cell_dphi(theRiemann->cR);   

                double dpL =  theRiemann->cm - pL;
                double dpR = -theRiemann->cm + pR;
                while( dpL >  PHIMAX/2. ) dpL -= PHIMAX;
                while( dpL < -PHIMAX/2. ) dpL += PHIMAX;
                while( dpR >  PHIMAX/2. ) dpR -= PHIMAX;
                while( dpR < -PHIMAX/2. ) dpR += PHIMAX;

                double WL = cell_prim(theRiemann->cL,q) + dpL*cell_gradp(theRiemann->cL,q);
                double WR = cell_prim(theRiemann->cR,q) - dpR*cell_gradp(theRiemann->cR,q);

                double deltaL = theRiemann->r-theRiemann->r_cell_L;
                double deltaR = theRiemann->r_cell_R - theRiemann->r;
                AvgPrim[q] = 0.5*(WL+WR);
                Grad_r_prim[q] = (WR-WL)/(deltaR+deltaL);
            }
            if (theRiemann->n[1]==1){
                AvgPrim[q] = 0.5*(cell_prim(theRiemann->cR,q) + cell_prim(theRiemann->cL,q));
                Grad_ph_prim[q] = (cell_prim(theRiemann->cR,q) - cell_prim(theRiemann->cL,q))/ (cell_tiph(theRiemann->cR)-cell_tiph(theRiemann->cL))/r;
            }
        }

        VFlux[q] = 0.0;
    }

    double nu;
    double alpha =  sim_EXPLICIT_VISCOSITY(theSim);
    double Gamma =  sim_GAMMALAW(theSim);
    if (sim_VISC_CONST(theSim)==1){
        nu = alpha;
    } else{
        if (sim_InitialDataType(theSim)==SHEAR){
            double HoR = 0.1;
            //nu = alpha*HoR*HoR*pow(fabs((r*cos(tiph))),1.5);
            nu = alpha*Gamma*AvgPrim[PPP]/AvgPrim[RHO]*pow(fabs(r*cos(tiph)),2.0);    
            if (r*cos(tiph)>20.) nu=0.0;
        } else{
            //nu = alpha*Gamma*AvgPrim[PPP]/AvgPrim[RHO]*pow(r,1.5);
            //nu = alpha*Gamma*AvgPrim[PPP]/AvgPrim[RHO]*(sqrt(M0)+sqrt(M1))/(sqrt(M0)*pow(dist_bh0,-1.5)+sqrt(M1)*pow(dist_bh1,-1.5));
            double eps = sim_G_EPS(theSim);
            nu = alpha*AvgPrim[PPP]/AvgPrim[RHO]/sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
        }
    }

    double rho = AvgPrim[RHO];
    double vr  = AvgPrim[URR];
    double om  = AvgPrim[UPP];
    double vz  = AvgPrim[UZZ];

    if (fabs(vz)>1.e-12){
        printf("you should not be using viscosity and having a nonzero component of velocity in z direction. Viscosity is only for 2d planar motion right now\n");
        exit(0);
    }

    double Gr_vr = Grad_r_prim[URR];
    double Gr_om = Grad_r_prim[UPP];
    double Gr_vz = Grad_r_prim[UZZ];

    double Gp_vr = Grad_ph_prim[URR];
    double Gp_om = Grad_ph_prim[UPP];
    double Gp_vz = Grad_ph_prim[UZZ];

    double om_cell = sim_rOm_a(theSim,r,a)/r;
    double rdr_om_cell = sim_rdrOm_a(theSim,r,a);


    //r direction
    if (theRiemann->n[0] ==1){
        VFlux[SRR] = -nu*rho*(r*Gr_vr - vr - r*r*Gp_om);
        VFlux[LLL] = -nu*rho*(r*(r*Gr_om+rdr_om_cell) + r*Gp_vr);
        VFlux[SZZ] = 0.0; //deal with this later
        VFlux[TAU] = - nu * rho * ( vr*( Gr_vr - vr/r - r*Gp_om ) + r*om*( Gp_vr + (rdr_om_cell + r*Gr_om) ));
    }
    //phi direction
    if (theRiemann->n[1] ==1){
        VFlux[SRR] = -nu*rho*(r*(r*Gr_om+rdr_om_cell) + r*Gp_vr);
        VFlux[LLL] = -nu*rho*(-r*Gr_vr + vr + r*r*Gp_om );
        VFlux[SZZ] = 0.0; //deal with this later
        VFlux[TAU] = -nu * rho * ( -r*om*( Gr_vr - vr/r - r*Gp_om ) + vr*(Gp_vr + (rdr_om_cell + r*Gr_om ) ));
    }

    for (q=0;q<NUM_Q;++q){
        theRiemann->Fvisc[q] = VFlux[q];
    }

    free(VFlux);
    free(Grad_r_prim);
    free(Grad_ph_prim);
    free(AvgPrim);
}


void riemann_AddFlux(struct Riemann * theRiemann, struct Sim *theSim,struct GravMass *theGravMasses,double dt ){
    int NUM_Q = sim_NUM_Q(theSim);
    double GAMMALAW = sim_GAMMALAW(theSim);

    riemann_set_vel(theRiemann,theSim,theGravMasses,theRiemann->r,GAMMALAW);

    double Bk_face,Psi_face;
    double w;
    if (theRiemann->n[PDIRECTION]){
        w = cell_wiph(theRiemann->cL);
    } else{
        w = 0.0;
    }

    // which state of the riemann problem are we in?
    riemann_set_state(theRiemann,w);


    if (theRiemann->state==LEFT){
        riemann_set_flux( theRiemann , theSim, GAMMALAW,LEFT);//in this case, we only need FL
        cell_prim2cons( theRiemann->primL , theRiemann->UL , theRiemann->r , 1.0 ,theSim);
        int q;
        for (q=0;q<sim_NUM_Q(theSim) ; ++q ){
            theRiemann->F[q] = theRiemann->FL[q] - w*theRiemann->UL[q];// w is only nonzero when we are in phi direction
        }
    } else if (theRiemann->state==RIGHT){
        riemann_set_flux( theRiemann , theSim, GAMMALAW,RIGHT);//in this case, we only need FR
        cell_prim2cons( theRiemann->primR , theRiemann->UR , theRiemann->r , 1.0 ,theSim);
        int q;
        for (q=0;q<sim_NUM_Q(theSim) ; ++q ){
            theRiemann->F[q] = theRiemann->FR[q] - w*theRiemann->UR[q];// w is only nonzero when we are in phi direction
        }
    } else{
        if (sim_Riemann(theSim)==HLL){
            riemann_set_flux( theRiemann , theSim, GAMMALAW,LEFT);  //we need both
            riemann_set_flux( theRiemann , theSim, GAMMALAW,RIGHT);    
            cell_prim2cons( theRiemann->primL , theRiemann->UL , theRiemann->r , 1.0 ,theSim);//we need both
            cell_prim2cons( theRiemann->primR , theRiemann->UR , theRiemann->r , 1.0 ,theSim);
            riemann_set_star_hll(theRiemann,theSim);// get Ustar and Fstar
        } else if (sim_Riemann(theSim)==HLLC){
            if (theRiemann->state==LEFTSTAR){
                riemann_set_flux( theRiemann , theSim, GAMMALAW,LEFT);//in this case, we only need FL
                cell_prim2cons( theRiemann->primL , theRiemann->UL , theRiemann->r , 1.0 ,theSim);
                riemann_set_star_hllc(theRiemann,theSim,GAMMALAW);// get Ustar and Fstar
            } else if (theRiemann->state==RIGHTSTAR){
                riemann_set_flux( theRiemann , theSim, GAMMALAW,RIGHT);//in this case, we only need FR      
                cell_prim2cons( theRiemann->primR , theRiemann->UR , theRiemann->r , 1.0 ,theSim);
                riemann_set_star_hllc(theRiemann,theSim,GAMMALAW);// get Ustar and Fstar
            }
        } else{
            printf("ERROR\n");
            exit(0);
        }
        int q;
        for (q=0;q<sim_NUM_Q(theSim) ; ++q ){
            theRiemann->F[q] = theRiemann->Fstar[q] - w*theRiemann->Ustar[q];// w is only nonzero when we are in phi direction
        }
    }


    // compute viscous flux terms
    if (sim_EXPLICIT_VISCOSITY(theSim)>0.0){
        riemann_visc_flux(theRiemann,theSim,theGravMasses );  
    }

    // Now actually update the conserved quantities using to the fluxes we have computed
    int q;
    for( q=0 ; q<NUM_Q ; ++q ){
        cell_add_cons(theRiemann->cL,q,-dt*theRiemann->dA*theRiemann->F[q]);
        cell_add_cons(theRiemann->cR,q,dt*theRiemann->dA*theRiemann->F[q]);
    }

    if (sim_EXPLICIT_VISCOSITY(theSim)>0.){
        cell_add_cons(theRiemann->cL,SRR,-dt*theRiemann->dA*theRiemann->Fvisc[SRR]/theRiemann->r_cell_L);
        cell_add_cons(theRiemann->cR,SRR, dt*theRiemann->dA*theRiemann->Fvisc[SRR]/theRiemann->r_cell_R);

        cell_add_cons(theRiemann->cL,LLL,-dt*theRiemann->dA*theRiemann->Fvisc[LLL]);
        cell_add_cons(theRiemann->cR,LLL, dt*theRiemann->dA*theRiemann->Fvisc[LLL]);

        //viscous heating
        cell_add_cons(theRiemann->cL,TAU,-dt*theRiemann->dA*theRiemann->Fvisc[TAU]);
        cell_add_cons(theRiemann->cR,TAU,+dt*theRiemann->dA*theRiemann->Fvisc[TAU]);
    }
}



