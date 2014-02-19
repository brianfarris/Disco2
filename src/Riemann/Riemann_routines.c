#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/Face.h"
#include "../Headers/header.h"


// ********************************************************************************************
// WE REALLY SHOULD IMPROVE THE COMMENTING OF ALL OF THESE ROUTINES. 
// THERE IS SOME COMPLICATED STUFF HERE. 
// LETS CHOOSE A REFERENCE SUCH AS TORO AND IDENTIFY LINES OF CODE WITH EQUATIONS IN THE BOOK.
// ********************************************************************************************

// Which state of the riemann problem are we in?
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

void riemann_set_star_hll(struct Riemann * theRiemann,struct Sim * theSim){
  double Sr =  theRiemann->Sr;
  double Sl = theRiemann->Sl;
  int q;
  for( q=0 ; q<sim_NUM_Q(theSim) ; ++q ){
    theRiemann->Ustar[q] = ( Sr*theRiemann->UR[q] - Sl*theRiemann->UL[q] + theRiemann->FL[q] - theRiemann->FR[q] )/( Sr - Sl );
    theRiemann->Fstar[q] = ( Sr*theRiemann->FL[q] - Sl*theRiemann->FR[q] - Sr*Sl*( theRiemann->UL[q] - theRiemann->UR[q] ) )/( Sr - Sl );
  }
}


void riemann_set_star_hllc(struct Riemann * theRiemann,struct Sim * theSim,double GAMMALAW){
  double r = theRiemann->pos[R_DIR];
  double *prim;
  double Sk;
  double *Uk;
  double *Fk;
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
  for( q=sim_NUM_C(theSim) ; q<sim_NUM_Q(theSim) ; ++q ){
    theRiemann->Ustar[q] = prim[q]*theRiemann->Ustar[DDD];
  }

  //Now set Fstar
  for (q=0;q<sim_NUM_Q(theSim);++q){
    theRiemann->Fstar[q] = Fk[q] + Sk*( theRiemann->Ustar[q] - Uk[q] ) ;
  }
}

void riemann_setup_rz(struct Riemann * theRiemann,struct Face * theFaces,struct Sim * theSim,int FaceNumber,int direction){
  theRiemann->n[direction]=1; // set
  int NUM_Q = sim_NUM_Q(theSim);
  double deltaL = face_deltaL(theFaces,FaceNumber);
  double deltaR = face_deltaR(theFaces,FaceNumber);
  theRiemann->cL = face_L_pointer(theFaces,FaceNumber);
  theRiemann->cR = face_R_pointer(theFaces,FaceNumber);
  double pL = cell_tiph(theRiemann->cL) - .5*cell_dphi(theRiemann->cL);
  double pR = cell_tiph(theRiemann->cR) - .5*cell_dphi(theRiemann->cR);   
  double dpL =  face_cm(theFaces,FaceNumber) - pL;
  double dpR = -face_cm(theFaces,FaceNumber) + pR;
  while( dpL >  M_PI ) dpL -= 2.*M_PI;
  while( dpL < -M_PI ) dpL += 2.*M_PI;
  while( dpR >  M_PI ) dpR -= 2.*M_PI;
  while( dpR < -M_PI ) dpR += 2.*M_PI;
  dpL = dpL;
  dpR = dpR;
  theRiemann->pos[R_DIR] = face_r(theFaces,FaceNumber);
  theRiemann->pos[Z_DIR] = face_z(theFaces,FaceNumber);
  theRiemann->pos[P_DIR] = face_phi(theFaces,FaceNumber);
  theRiemann->dA = face_dA(theFaces,FaceNumber);
  theRiemann->cm = face_cm(theFaces,FaceNumber);

  if (direction==0){
    theRiemann->r_cell_L = theRiemann->pos[R_DIR]-deltaL;
    theRiemann->r_cell_R = theRiemann->pos[R_DIR]+deltaR;
  } else{
    theRiemann->r_cell_L = theRiemann->pos[R_DIR];
    theRiemann->r_cell_R = theRiemann->pos[R_DIR];
  }

  int q;
  for (q=0;q<NUM_Q;++q){
    theRiemann->primL[q] = cell_prim(theRiemann->cL,q) + cell_grad(theRiemann->cL,q)*deltaL + cell_gradp(theRiemann->cL,q)*dpL;
    theRiemann->primR[q] = cell_prim(theRiemann->cR,q) - cell_grad(theRiemann->cR,q)*deltaR - cell_gradp(theRiemann->cR,q)*dpR;
  }
}

void riemann_setup_p(struct Riemann * theRiemann,struct Cell *** theCells,struct Sim * theSim,int i,int j_low,int k,int direction){
  theRiemann->n[direction]=1; // set
  int NUM_Q = sim_NUM_Q(theSim);

  int j_hi;
  if (j_low == sim_N_p(theSim,i)-1){
    j_hi = 0;
  } else{
    j_hi = j_low+1;
  }
  theRiemann->cL = cell_single(theCells,i,j_low,k);
  theRiemann->cR = cell_single(theCells,i,j_hi ,k);
  double dpL = cell_dphi(theRiemann->cL);
  double dpR = cell_dphi(theRiemann->cR);
  double zm = sim_FacePos(theSim,k-1,Z_DIR);
  double zp = sim_FacePos(theSim,k,Z_DIR);
  double dz = zp-zm;
  double z = 0.5*(zp+zm);
  double rm = sim_FacePos(theSim,i-1,R_DIR);
  double rp = sim_FacePos(theSim,i,R_DIR);
  double dr = rp-rm;
  double r = .5*(rp+rm);
  theRiemann->dA = dr*dz;
  theRiemann->pos[R_DIR] = r; 
  theRiemann->pos[Z_DIR] = z; 
  theRiemann->pos[P_DIR] = cell_tiph(theRiemann->cL); 
  theRiemann->cm = cell_tiph(theRiemann->cL);

  theRiemann->r_cell_L = r;
  theRiemann->r_cell_R = r;

  int q;
  for (q=0;q<NUM_Q;++q){
    theRiemann->primL[q] = cell_prim(theRiemann->cL,q) + 0.5*cell_gradp(theRiemann->cL,q)*dpL;
    theRiemann->primR[q] = cell_prim(theRiemann->cR,q) - 0.5*cell_gradp(theRiemann->cR,q)*dpR;
  }

}

void riemann_AddFlux(struct Riemann * theRiemann, struct Sim *theSim,double dt ){
  int NUM_Q = sim_NUM_Q(theSim);
  double GAMMALAW = sim_GAMMALAW(theSim);

  riemann_set_vel(theRiemann,theSim,theRiemann->pos[R_DIR],GAMMALAW);

  double w;
  if (theRiemann->n[PDIRECTION]){
    if( sim_MOVE_CELLS(theSim) == C_WRIEMANN ) cell_add_wiph(theRiemann->cL,theRiemann->Ss);
    w = cell_wiph(theRiemann->cL);
  } else{
    w = 0.0;
  }
  // which state of the riemann problem are we in?
  riemann_set_state(theRiemann,w);

  if (theRiemann->state==LEFT){
    riemann_set_flux( theRiemann , theSim, GAMMALAW,LEFT);//in this case, we only need FL
    cell_prim2cons( theRiemann->primL , theRiemann->UL , theRiemann->pos , 1.0 ,theSim);
    int q;
    for (q=0;q<sim_NUM_Q(theSim) ; ++q ){
      theRiemann->F[q] = theRiemann->FL[q] - w*theRiemann->UL[q];// w is only nonzero when we are in phi direction
    }
  } else if (theRiemann->state==RIGHT){
    riemann_set_flux( theRiemann , theSim, GAMMALAW,RIGHT);//in this case, we only need FR
    cell_prim2cons( theRiemann->primR , theRiemann->UR , theRiemann->pos , 1.0 ,theSim);
    int q;
    for (q=0;q<sim_NUM_Q(theSim) ; ++q ){
      theRiemann->F[q] = theRiemann->FR[q] - w*theRiemann->UR[q];// w is only nonzero when we are in phi direction
    }
  } else{
    if (sim_Riemann(theSim)==HLL){
      riemann_set_flux( theRiemann , theSim, GAMMALAW,LEFT);  //we need both
      riemann_set_flux( theRiemann , theSim, GAMMALAW,RIGHT);    
      cell_prim2cons( theRiemann->primL , theRiemann->UL , theRiemann->pos , 1.0 ,theSim);//we need both
      cell_prim2cons( theRiemann->primR , theRiemann->UR , theRiemann->pos , 1.0 ,theSim);
      riemann_set_star_hll(theRiemann,theSim);// get Ustar and Fstar
    } else if (sim_Riemann(theSim)==HLLC){
      if (theRiemann->state==LEFTSTAR){
        riemann_set_flux( theRiemann , theSim, GAMMALAW,LEFT);//in this case, we only need FL
        cell_prim2cons( theRiemann->primL , theRiemann->UL , theRiemann->pos , 1.0 ,theSim);
        riemann_set_star_hllc(theRiemann,theSim,GAMMALAW);// get Ustar and Fstar
      } else if (theRiemann->state==RIGHTSTAR){
        riemann_set_flux( theRiemann , theSim, GAMMALAW,RIGHT);//in this case, we only need FR      
        cell_prim2cons( theRiemann->primR , theRiemann->UR , theRiemann->pos , 1.0 ,theSim);
        riemann_set_star_hllc(theRiemann,theSim,GAMMALAW);// get Ustar and Fstar
      }
    } else{
      printf("ERROR: Unknown state given to riemann_add_flux.\n");
      exit(0);
    }
    int q;
    for (q=0;q<sim_NUM_Q(theSim) ; ++q ){
      theRiemann->F[q] = theRiemann->Fstar[q] - w*theRiemann->Ustar[q];// w is only nonzero when we are in phi direction
    }
  }

  int q;
  for( q=0 ; q<NUM_Q ; ++q ){
    cell_add_cons(theRiemann->cL,q,-dt*theRiemann->dA*theRiemann->F[q]);
    cell_add_cons(theRiemann->cR,q,dt*theRiemann->dA*theRiemann->F[q]);

  }
}



