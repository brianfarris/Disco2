#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Grid.h"
#include "../Headers/Cell.h"
#include "../Headers/Face.h"
#include "../Headers/header.h"

//try to come up with a better name for this routine
void LR_speed(double *prim,double *n,double r,double GAMMALAW,double * p_vn,double * p_cf2,double *Fm,double * p_mn){
  double P   = prim[PPP];
  double rho = prim[RHO];
  double vr  =   prim[URR];
  double vp  = r*prim[UPP];
  double vz  =   prim[UZZ];
  double vn  = vr*n[0] + vp*n[1] + vz*n[2];
  double cf2  = GAMMALAW*(P/rho);
  double mr = rho*vr;
  double mp = rho*vp;
  double mz = rho*vz;
  double mn = mr*n[0]+mp*n[1]+mz*n[2];

  Fm[0] = rho*vr*vn + P*n[0];
  Fm[1] = rho*vp*vn + P*n[1];
  Fm[2] = rho*vz*vn + P*n[2];

  *p_vn = vn;
  *p_cf2 = cf2;
  *p_mn = mn;
}

void LR_speed_mhd(double *prim,double *n,double r,double ch, double * p_cf2,double *F,double *Fm,double * p_Bn,double * p_B2){
  double rho  =   prim[RHO];
  double vr  =   prim[URR];
  double vp  = r*prim[UPP];
  double vz  =   prim[UZZ];
  double vn  = vr*n[0] + vp*n[1] + vz*n[2];
 
  double Br = prim[BRR];
  double Bp = prim[BPP];
  double Bz = prim[BZZ];
  double psi = prim[PSI];

  double Bn =  Br*n[0] + Bp*n[1] + Bz*n[2];
  double B2 = (Br*Br  + Bp*Bp  + Bz*Bz);
  double Fps = pow(ch,2.)*Bn;

  F[0] = vn*Br - Bn*vr + psi*n[0];
  F[1] = vn*Bp - Bn*vp + psi*n[1];
  F[2] = vn*Bz - Bn*vz + psi*n[2];

  Fm[0] +=  .5*B2*n[0] - Br*Bn;
  Fm[1] +=  .5*B2*n[1] - Bp*Bn;
  Fm[2] +=  .5*B2*n[2] - Bz*Bn;

  *p_cf2 = .5*( *p_cf2 + B2/rho + sqrt(fabs(  (*p_cf2+B2/rho)*(*p_cf2+B2/rho) - 4.0*(*p_cf2)*Bn*Bn/rho )) );

  *p_Bn = Bn;
  *p_B2 = B2;
}

void riemann_set_vel(struct Riemann * theRiemann,struct Grid * theGrid, double *n,double r,double *Bpack,double GAMMALAW,double DIVB_CH){
  double L_Mins, L_Plus, L_Star;

  double vnL,cf21,mnL,BnL,B2L;
  double FL[3], FmL[3];
  LR_speed(theRiemann->primL,n,r,GAMMALAW,&vnL,&cf21,FmL,&mnL);
  if (grid_runtype(theGrid)==MHD){
    LR_speed_mhd(theRiemann->primL,n,r,DIVB_CH,&cf21,FL,FmL,&BnL,&B2L);
  }
  L_Mins = vnL - sqrt( cf21 );
  L_Plus = vnL + sqrt( cf21 );
//  printf("1 L_Mins: %e L_Plus: %e\n",L_Mins,L_Plus);

  double vnR,cf22,mnR,BnR,B2R;
  double FR[3],FmR[3];
  LR_speed(theRiemann->primR,n,r,GAMMALAW,&vnR,&cf22,FmR,&mnR);
  if (grid_runtype(theGrid)==MHD){
    LR_speed_mhd(theRiemann->primR,n,r,DIVB_CH,&cf22,FR,FmR,&BnR,&B2R);
  }
  if( L_Mins > vnR - sqrt( cf22 ) ) L_Mins = vnR - sqrt( cf22 );
  if( L_Plus < vnR + sqrt( cf22 ) ) L_Plus = vnR + sqrt( cf22 );
//  printf("2 L_Mins: %e L_Plus: %e\n",L_Mins,L_Plus);
  
  if (grid_runtype(theGrid)==MHD){
    if( L_Mins > -DIVB_CH ) L_Mins = -DIVB_CH;
    if( L_Plus <  DIVB_CH ) L_Plus =  DIVB_CH;
  }

  double aL = L_Plus;
  double aR = -L_Mins;

  double mr = ( aR*theRiemann->primL[RHO]*theRiemann->primL[URR] + aL*theRiemann->primR[RHO]*theRiemann->primR[URR] + FmL[0] - FmR[0] )/( aL + aR );
  double mp = ( aR*theRiemann->primL[RHO]*theRiemann->primL[UPP]*r + aL*theRiemann->primR[RHO]*theRiemann->primR[UPP]*r + FmL[1] - FmR[1] )/( aL + aR );
  double mz = ( aR*theRiemann->primL[RHO]*theRiemann->primL[UZZ] + aL*theRiemann->primR[RHO]*theRiemann->primR[UZZ] + FmL[2] - FmR[2] )/( aL + aR );
  double rho = ( aR*theRiemann->primL[RHO] + aL*theRiemann->primR[RHO] + mnL - mnR )/( aL + aR );

  L_Star = (theRiemann->primR[RHO]*vnR*(L_Plus-vnR)-theRiemann->primL[RHO]*vnL*(L_Mins-vnL)+theRiemann->primL[PPP]-theRiemann->primR[PPP])
    /(theRiemann->primR[RHO]*(L_Plus-vnR)-theRiemann->primL[RHO]*(L_Mins-vnL));
//  printf("3a L_Star: %e\n",L_Star);

  if (grid_runtype(theGrid)==MHD){  
    double Br = ( aR*theRiemann->primL[BRR] + aL*theRiemann->primR[BRR] + FL[0] - FR[0] )/( aL + aR );
    double Bp = ( aR*theRiemann->primL[BPP] + aL*theRiemann->primR[BPP] + FL[1] - FR[1] )/( aL + aR );
    double Bz = ( aR*theRiemann->primL[BZZ] + aL*theRiemann->primR[BZZ] + FL[2] - FR[2] )/( aL + aR );
    double psi = ( aR*theRiemann->primL[PSI] + aL*theRiemann->primR[PSI] + pow(DIVB_CH,2.)*BnL - pow(DIVB_CH,2.)*BnR )/( aL + aR );
    Bpack[0] = Br*n[0] + Bp*n[1] + Bz*n[2]; // Bn
    Bpack[1] = Br;
    Bpack[2] = Bp;
    Bpack[3] = Bz;
    Bpack[4] = (mr*Br + mp*Bp + mz*Bz)/rho; // v dot B
    Bpack[5] = psi;

    L_Star += ((.5*B2L-BnL*BnL)-.5*B2R+BnR*BnR)
      /(theRiemann->primR[RHO]*(L_Plus-vnR)-theRiemann->primL[RHO]*(L_Mins-vnL));
// printf("blah: %e\n",((.5*B2L-BnL*BnL)-.5*B2R+BnR*BnR)
 //     /(theRiemann->primR[RHO]*(L_Plus-vnR)-theRiemann->primL[RHO]*(L_Mins-vnL)));
 //printf("3b L_Star: %e\n",L_Star);
  }
//  printf("mr: %e, mp: %e, mz: %e, rho: %e\n",mr,mp,mz,rho);
/*
  int q;
  for (q=0;q<6;++q){
    printf("Bpack[%d]: %e\n",q,Bpack[q]);
  }
  */
  theRiemann->Sl = L_Mins;
  theRiemann->Sr = L_Plus;
  theRiemann->Ss = L_Star;

}

void riemann_set_state(struct Riemann * theRiemann,int w ){
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


void riemann_set_Ustar(struct Riemann * theRiemann,struct Grid * theGrid, double *n,double r,double *Bpack,double GAMMALAW){
  double *prim;
  double Sk;
  if (theRiemann->state==LEFTSTAR){
    prim = theRiemann->primL;
    Sk = theRiemann->Sl;
  }else{
    prim = theRiemann->primR;
    Sk = theRiemann->Sr;
  }
  double Ss=theRiemann->Ss;

  double rho = prim[RHO];
  double vr  = prim[URR];
  double vp  = prim[UPP]*r;
  double vz  = prim[UZZ];
  double Pp  = prim[PPP];
  double v2 = vr*vr+vp*vp+vz*vz;
  double vn = vr*n[0] + vp*n[1] + vz*n[2];
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

  if (grid_runtype(theGrid)==MHD){
    double Bsn = Bpack[0];
    double Bsr = Bpack[1];
    double Bsp = Bpack[2];
    double Bsz = Bpack[3];
    double vBs = Bpack[4];
    double psi = Bpack[5];

    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double B2 = Br*Br+Bp*Bp+Bz*Bz;
    double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];
    double vB = vr*Br   + vp*Bp   + vz*Bz;
    double Bs2 = Bsr*Bsr+Bsp*Bsp+Bsz*Bsz;
    //Ps    +=  (.5*B2-Bn*Bn) - .5*Bs2 + Bsn*Bsn;
    double Ps_mag = (.5*B2-Bn*Bn) - .5*Bs2 + Bsn*Bsn;
    Msr   += ( Br*Bn - Bsr*Bsn ) / ( Sk - Ss );
    Msp   += ( Bp*Bn - Bsp*Bsn ) / ( Sk - Ss );
    Msz   += ( Bz*Bn - Bsz*Bsn ) / ( Sk - Ss );
    Estar += ( ( Sk - vn )*.5*B2 + (Ps_mag+.5*Bs2)*Ss - .5*B2*vn - vBs*Bsn + vB*Bn ) / ( Sk - Ss );
    /*
    if (fabs(r-2.875)<0.00001){
      printf("Estar: %e\n",Estar);
      printf("Sk - vn: %e\n",Sk-vn);
      printf("Sk - Ss: %e\n",Sk-Ss);
      printf("Ps: %e\n",Ps+Ps_mag);
      printf("Pp: %e\n",Pp);
      printf("vBs: %e\n",vBs);
      printf("Bn: %e\n",Bn);
      printf("Bsn: %e\n",Bsn);
      printf("vn: %e\n",vn);
      printf("B2: %e\n",B2);
      printf("Bs2: %e\n",Bs2);
      printf("rhoe: %e\n",rhoe);
      printf("vB: %e\n",vB);
      printf("rho: %e\n",rho);
      printf("v2: %e\n",v2);
  exit(0); 
    }
    */

    theRiemann->Ustar[BRR] = Bsr/r;
    theRiemann->Ustar[BPP] = Bsp/r;
    theRiemann->Ustar[BZZ] = Bsz;
    theRiemann->Ustar[PSI] = psi;
  }
  double mn  = Msr*n[0] + Msp*n[1] + Msz*n[2];

  Msr += n[0]*( Msn - mn );
  Msp += n[1]*( Msn - mn );
  Msz += n[2]*( Msn - mn );

  theRiemann->Ustar[DDD] = Dstar;
  theRiemann->Ustar[SRR] = Msr;
  theRiemann->Ustar[LLL] = r*Msp;
  theRiemann->Ustar[SZZ] = Msz;
  theRiemann->Ustar[TAU] = Estar;

}

void riemann_set_flux(struct Riemann * theRiemann, struct Grid * theGrid, double r , double * n ,double GAMMALAW,double DIVB_CH){
  double *prim;
  if ((theRiemann->state==LEFT)||(theRiemann->state==LEFTSTAR)){
    prim = theRiemann->primL;
  }else{
    prim = theRiemann->primR;
  }

  double rho = prim[RHO];
  double Pp  = prim[PPP];
  double vr  = prim[URR];
  double vp  = prim[UPP]*r;
  double vz  = prim[UZZ];
  double vn = vr*n[0] + vp*n[1] + vz*n[2];
  double rhoe = Pp/(GAMMALAW-1.);
  double v2 = vr*vr + vp*vp + vz*vz;
  theRiemann->F[DDD] = rho*vn;
  theRiemann->F[SRR] =     rho*vr*vn + Pp*n[0] ;
  theRiemann->F[LLL] = r*( rho*vp*vn + Pp*n[1] );
  theRiemann->F[SZZ] =     rho*vz*vn + Pp*n[2] ;
  theRiemann->F[TAU] = ( .5*rho*v2 + rhoe + Pp )*vn ;

  if (grid_runtype(theGrid)==MHD){ 
    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];
    double vB = vr*Br + vp*Bp + vz*Bz;
    double B2 = Br*Br + Bp*Bp + Bz*Bz;

    theRiemann->F[SRR] +=     .5*B2*n[0] - Br*Bn;
    theRiemann->F[LLL] += r*( .5*B2*n[1] - Bp*Bn );
    theRiemann->F[SZZ] +=     .5*B2*n[2] - Bz*Bn;
    theRiemann->F[TAU] += B2*vn - vB*Bn;
    double psi = prim[PSI];
    theRiemann->F[BRR] =(Br*vn - vr*Bn + psi*n[0])/r;
    theRiemann->F[BPP] =(Bp*vn - vp*Bn + psi*n[1])/r;
    theRiemann->F[BZZ] = Bz*vn - vz*Bn + psi*n[2];
    theRiemann->F[PSI] = pow(DIVB_CH,2.)*Bn;
  }
}

void riemann_addto_flux_general(struct Riemann * theRiemann,double w,int NUM_Q){
  int q;
  for (q=0;q<NUM_Q;++q){
    if ((theRiemann->state==LEFT)||(theRiemann->state==RIGHT)){
      theRiemann->F[q] -= w*theRiemann->Uk[q];
    }else if(theRiemann->state==LEFTSTAR){
      theRiemann->F[q] += theRiemann->Sl*( theRiemann->Ustar[q] - theRiemann->Uk[q] ) - w*theRiemann->Ustar[q];
    }else{
      theRiemann->F[q] += theRiemann->Sr*( theRiemann->Ustar[q] - theRiemann->Uk[q] ) - w*theRiemann->Ustar[q]; 
    }
  }
}

void riemann_setup_rz(struct Riemann * theRiemann,struct Face * theFaces,struct Grid * theGrid,int n){
  int NUM_Q = grid_NUM_Q(theGrid);
  double deltaL = face_deltaL(face_pointer(theFaces,n));
  double deltaR = face_deltaR(face_pointer(theFaces,n));
  theRiemann->cL = face_L_pointer(theFaces,n);
  theRiemann->cR = face_R_pointer(theFaces,n);
  double pL = cell_tiph(theRiemann->cL) - .5*cell_dphi(theRiemann->cL);
  double pR = cell_tiph(theRiemann->cR) - .5*cell_dphi(theRiemann->cR);   
  double dpL =  face_cm(face_pointer(theFaces,n)) - pL;
  double dpR = -face_cm(face_pointer(theFaces,n)) + pR;
  while( dpL >  M_PI ) dpL -= 2.*M_PI;
  while( dpL < -M_PI ) dpL += 2.*M_PI;
  while( dpR >  M_PI ) dpR -= 2.*M_PI;
  while( dpR < -M_PI ) dpR += 2.*M_PI;
  dpL = dpL;
  dpR = dpR;
  theRiemann->r = face_r(face_pointer(theFaces,n));
  theRiemann->dA = face_dA(face_pointer(theFaces,n));

  int q;
  for (q=0;q<NUM_Q;++q){
    theRiemann->primL[q] = cell_prims(theRiemann->cL)[q] + cell_grad(theRiemann->cL)[q]*deltaL + cell_gradp(theRiemann->cL)[q]*dpL;
    theRiemann->primR[q] = cell_prims(theRiemann->cR)[q] + cell_grad(theRiemann->cR)[q]*deltaR + cell_gradp(theRiemann->cR)[q]*dpR;
  }
}

void riemann_setup_p(struct Riemann * theRiemann,struct Cell *** theCells,struct Grid * theGrid,int i,int j_low,int j_hi,int k){
  int NUM_Q = grid_NUM_Q(theGrid);
  theRiemann->cL = cell_single(theCells,i,j_low,k);
  theRiemann->cR = cell_single(theCells,i,j_hi ,k);
  double dpL = cell_dphi(theRiemann->cL);
  double dpR = cell_dphi(theRiemann->cR);
  double zm = grid_z_faces(theGrid,k-1);
  double zp = grid_z_faces(theGrid,k);
  double dz = zp-zm;
  double rm = grid_r_faces(theGrid,i-1);
  double rp = grid_r_faces(theGrid,i);
  double dr = rp-rm;
  double r = .5*(rp+rm);
  theRiemann->dA = dr*dz;
  theRiemann->r = r; 
  int q;
  for (q=0;q<NUM_Q;++q){
    theRiemann->primL[q] = cell_prims(theRiemann->cL)[q] + cell_gradp(theRiemann->cL)[q]*dpL;
    theRiemann->primR[q] = cell_prims(theRiemann->cR)[q] + cell_gradp(theRiemann->cR)[q]*dpR;
  }

}


//void riemann_hllc(struct Riemann * theRiemann, struct Cell * cL , struct Cell * cR, struct Grid *theGrid, double dA,double dt, double r,double deltaL,double deltaR, double dpL, double dpR , int direction ){
void riemann_hllc(struct Riemann * theRiemann, struct Grid *theGrid,double dt, int direction ){
  int NUM_Q = grid_NUM_Q(theGrid);
  double GAMMALAW = grid_GAMMALAW(theGrid);
  double DIVB_CH = grid_DIVB_CH(theGrid);

  /*
     int q;
     for (q=0;q<NUM_Q;++q){
     theRiemann->primL[q] = cell_prims(cL)[q] + cell_grad(cL)[q]*deltaL + cell_gradp(cL)[q]*dpL;
     theRiemann->primR[q] = cell_prims(cR)[q] + cell_grad(cR)[q]*deltaR + cell_gradp(cR)[q]*dpR;
     }
     */
  //initialize
  double n[3];int i;
  for (i=0;i<3;++i){
    n[i]=0;
  }
  //set
  n[direction]=1.0;

  double Bpack[6];
  riemann_set_vel(theRiemann,theGrid,n,theRiemann->r,Bpack,GAMMALAW,DIVB_CH);

  //printf("Sl: %e, Sr: %e, Ss: %e\n",theRiemann->Sl,theRiemann->Sr,theRiemann->Ss);

  double Bk_face;
  if (direction==0){
    Bk_face = 0.5*(theRiemann->primL[BRR]+theRiemann->primR[BRR]);  
  } else if (direction==1){
    Bk_face = 0.5*(theRiemann->primL[BPP]+theRiemann->primR[BPP]);
  } else if (direction==2){
    Bk_face = 0.5*(theRiemann->primL[BZZ]+theRiemann->primR[BZZ]);
  }
  double Psi_face = 0.5*(theRiemann->primL[PSI]+theRiemann->primR[PSI]);

  double w;
  if (direction==1){
    if( grid_MOVE_CELLS(theGrid) == C_WRIEMANN ) cell_add_wiph(theRiemann->cL,theRiemann->Ss);
    w = cell_wiph(theRiemann->cL);
  } else{
    w = 0.0;
  }

  riemann_set_state(theRiemann,w);

  riemann_set_flux( theRiemann , theGrid, theRiemann->r , n,GAMMALAW,DIVB_CH );
  if (theRiemann->state==LEFTSTAR){
    cell_prim2cons( theRiemann->primL , theRiemann->Uk , theRiemann->r , 1.0 ,GAMMALAW,grid_runtype(theGrid));
    riemann_set_Ustar(theRiemann,theGrid,n,theRiemann->r,Bpack,GAMMALAW);
  
  } else if(theRiemann->state==RIGHTSTAR){
    cell_prim2cons( theRiemann->primR , theRiemann->Uk , theRiemann->r , 1.0 ,GAMMALAW,grid_runtype(theGrid));
    riemann_set_Ustar(theRiemann,theGrid,n,theRiemann->r,Bpack,GAMMALAW);
  }
  riemann_addto_flux_general(theRiemann,w,grid_NUM_Q(theGrid));
/*
 if ((fabs(theRiemann->r-2.6875)<0.00001)&&(direction==2)){
      printf("r: %e, Sl: %e, Sr: %e, Ss: %e\n",theRiemann->r,theRiemann->Sl,theRiemann->Sr,theRiemann->Ss);
      printf("Flux[TAU]: %e, Uk[TAU]: %e\n",theRiemann->F[TAU],theRiemann->Uk[TAU]);
      int q;
      for (q=0;q<NUM_Q;++q){
        printf("primL[%d]: %e\n",q,theRiemann->primL[q]);
        printf("primR[%d]: %e\n",q,theRiemann->primR[q]);
      }
exit(0);
    }
    */
  //  if ((theRiemann->primL[BZZ]>0.00001)&&(direction==0)){
  int q;
  for( q=0 ; q<NUM_Q ; ++q ){
    cell_add_cons(theRiemann->cL,q,-dt*theRiemann->dA*theRiemann->F[q]);
    cell_add_cons(theRiemann->cR,q,dt*theRiemann->dA*theRiemann->F[q]);
  }

  cell_add_divB(theRiemann->cL,theRiemann->dA*Bk_face);
  cell_add_divB(theRiemann->cR,-theRiemann->dA*Bk_face);
  cell_add_GradPsi(theRiemann->cL,direction,Psi_face);
  cell_add_GradPsi(theRiemann->cR,direction,-Psi_face);

}



