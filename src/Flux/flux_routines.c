#define FLUX_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Flux.h"
#include "../Headers/header.h"

void flux_set_primL(struct Flux * theFlux,int q,double input){
  theFlux->primL[q] = input;
}
void flux_set_primR(struct Flux * theFlux,int q,double input){
  theFlux->primR[q] = input;
}
void flux_set_vel(struct Flux * theFlux,double *n,double r,double *Bpack,double GAMMALAW,double DIVB_CH){
  double wp = 0.0*n[1];
  double ch = DIVB_CH;

  double L_Mins, L_Plus, L_Star;

  double * prim1 = theFlux->primL;
  double * prim2 = theFlux->primR;

  double P1   = prim1[PPP];
  double rho1 = prim1[RHO];

  double vr1  =   prim1[URR];
  double vp1  = r*prim1[UPP];
  double vz1  =   prim1[UZZ];

  double vn1  = vr1*n[0] + vp1*n[1] + vz1*n[2];

  double cs1  = sqrt( GAMMALAW*(P1/rho1) );

  double Br1 = prim1[BRR];
  double Bp1 = prim1[BPP];
  double Bz1 = prim1[BZZ];

  double Bn1 =  Br1*n[0] + Bp1*n[1] + Bz1*n[2];
  double B21 = (Br1*Br1  + Bp1*Bp1  + Bz1*Bz1);
  double b21 = B21/rho1;
  double psiL = prim1[PSI];
  double FpsL = wp*psiL + pow(DIVB_CH,2.)*Bn1;

  double FrL = vn1*Br1 - Bn1*vr1 + psiL*n[0];
  double FpL = vn1*Bp1 - Bn1*vp1 + psiL*n[1];
  double FzL = vn1*Bz1 - Bn1*vz1 + psiL*n[2];

  double mrL = rho1*vr1;
  double mpL = rho1*vp1;
  double mzL = rho1*vz1;

  double FmrL = rho1*vr1*vn1 + (P1+.5*B21)*n[0] - Br1*Bn1;
  double FmpL = rho1*vp1*vn1 + (P1+.5*B21)*n[1] - Bp1*Bn1;
  double FmzL = rho1*vz1*vn1 + (P1+.5*B21)*n[2] - Bz1*Bn1;

  double cf21 = .5*( cs1*cs1 + b21 + sqrt(fabs(  (cs1*cs1+b21)*(cs1*cs1+b21) - 4.0*cs1*cs1*Bn1*Bn1/rho1 )) );

  L_Mins = vn1 - sqrt( cf21 );
  L_Plus = vn1 + sqrt( cf21 );

  double P2   = prim2[PPP];
  double rho2 = prim2[RHO];

  double vr2  =   prim2[URR];
  double vp2  = r*prim2[UPP];
  double vz2  =   prim2[UZZ];

  double vn2  = vr2*n[0] + vp2*n[1] + vz2*n[2];

  double cs2  = sqrt( GAMMALAW*(P2/rho2) );

  double Br2 = prim2[BRR];
  double Bp2 = prim2[BPP];
  double Bz2 = prim2[BZZ];

  double Bn2 =  Br2*n[0] + Bp2*n[1] + Bz2*n[2];
  double B22 = (Br2*Br2  + Bp2*Bp2  + Bz2*Bz2);
  double b22 = B22/rho2;
  double psiR = prim2[PSI];
  double FpsR = wp*psiR + pow(DIVB_CH,2.)*Bn2;

  double FrR = vn2*Br2 - Bn2*vr2 + psiR*n[0];
  double FpR = vn2*Bp2 - Bn2*vp2 + psiR*n[1];
  double FzR = vn2*Bz2 - Bn2*vz2 + psiR*n[2];

  double mrR = rho2*vr2;
  double mpR = rho2*vp2;
  double mzR = rho2*vz2;

  double FmrR = rho2*vr2*vn2 + (P2+.5*B22)*n[0] - Br2*Bn2;
  double FmpR = rho2*vp2*vn2 + (P2+.5*B22)*n[1] - Bp2*Bn2;
  double FmzR = rho2*vz2*vn2 + (P2+.5*B22)*n[2] - Bz2*Bn2;

  double cf22 = .5*( cs2*cs2 + b22 + sqrt(fabs(  (cs2*cs2+b22)*(cs2*cs2+b22) - 4.0*cs2*cs2*Bn2*Bn2/rho2 )) );

  if( L_Mins > vn2 - sqrt( cf22 ) ) L_Mins = vn2 - sqrt( cf22 );
  if( L_Plus < vn2 + sqrt( cf22 ) ) L_Plus = vn2 + sqrt( cf22 );
  if( L_Mins > -ch+wp ) L_Mins = -ch+wp;
  if( L_Plus <  ch+wp ) L_Plus =  ch+wp;

  double aL = L_Plus;
  double aR = -L_Mins;

  double Br = ( aR*Br1 + aL*Br2 + FrL - FrR )/( aL + aR );
  double Bp = ( aR*Bp1 + aL*Bp2 + FpL - FpR )/( aL + aR );
  double Bz = ( aR*Bz1 + aL*Bz2 + FzL - FzR )/( aL + aR );
  double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];

  double mr = ( aR*mrL + aL*mrR + FmrL - FmrR )/( aL + aR );
  double mp = ( aR*mpL + aL*mpR + FmpL - FmpR )/( aL + aR );
  double mz = ( aR*mzL + aL*mzR + FmzL - FmzR )/( aL + aR );

  double mnL = mrL*n[0]+mpL*n[1]+mzL*n[2];
  double mnR = mrR*n[0]+mpR*n[1]+mzR*n[2];
  double rho = ( aR*rho1 + aL*rho2 + mnL - mnR )/( aL + aR );
  double psi = ( aR*psiL + aL*psiR + FpsL - FpsR )/( aL + aR );

  L_Star = ( rho2*vn2*(L_Plus-vn2) - rho1*vn1*(L_Mins-vn1) + (P1+.5*B21-Bn1*Bn1) - (P2+.5*B22-Bn2*Bn2) )/( rho2*(L_Plus-vn2) - rho1*(L_Mins-vn1) );

  double vr = mr/rho;
  double vp = mp/rho;
  double vz = mz/rho;
  double vdotB = vr*Br + vp*Bp + vz*Bz;

  Bpack[0] = Bn;
  Bpack[1] = Br;
  Bpack[2] = Bp;
  Bpack[3] = Bz;
  Bpack[4] = vdotB;
  Bpack[5] = psi;

  theFlux->Sl = L_Mins;
  theFlux->Sr = L_Plus;
  theFlux->Ss = L_Star;

}

flux_set_Ustar(struct Flux * theFlux,double *n,double r,double *Bpack,double GAMMALAW,double w){
  double *prim;
  double Sk;
  if (w<(theFlux->Ss)){
    prim = theFlux->primL;
    Sk = theFlux->Sl;
  }else{
    prim = theFlux->primR;
    Sk = theFlux->Sr;
  }
  double Ss=theFlux->Ss;
  
  double Bsn = Bpack[0];
  double Bsr = Bpack[1];
  double Bsp = Bpack[2];
  double Bsz = Bpack[3];
  double vBs = Bpack[4];
  double psi = Bpack[5];

  double rho = prim[RHO];
  double vr  = prim[URR];
  double vp  = prim[UPP]*r;
  double vz  = prim[UZZ];
  double Pp  = prim[PPP];

  double Br  = prim[BRR];
  double Bp  = prim[BPP];
  double Bz  = prim[BZZ];

  double v2 = vr*vr+vp*vp+vz*vz;
  double B2 = Br*Br+Bp*Bp+Bz*Bz;

  double vn = vr*n[0] + vp*n[1] + vz*n[2];
  double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];
  double vB = vr*Br   + vp*Bp   + vz*Bz;

  double rhoe = Pp/(GAMMALAW-1.);

  double D  = rho;
  double mr = rho*vr;
  double mp = rho*vp;
  double mz = rho*vz;
  double E  = .5*rho*v2 + rhoe + .5*B2;

  double Bs2 = Bsr*Bsr+Bsp*Bsp+Bsz*Bsz;
  double Ps  = rho*( Sk - vn )*( Ss - vn ) + (Pp+.5*B2-Bn*Bn) - .5*Bs2 + Bsn*Bsn;

  double Dstar = ( Sk - vn )*D/( Sk - Ss );
  double Msr   = ( ( Sk - vn )*mr + ( Br*Bn - Bsr*Bsn ) ) / ( Sk - Ss );
  double Msp   = ( ( Sk - vn )*mp + ( Bp*Bn - Bsp*Bsn ) ) / ( Sk - Ss );
  double Msz   = ( ( Sk - vn )*mz + ( Bz*Bn - Bsz*Bsn ) ) / ( Sk - Ss );
  double Estar = ( ( Sk - vn )*E + (Ps+.5*Bs2)*Ss - (Pp+.5*B2)*vn - vBs*Bsn + vB*Bn ) / ( Sk - Ss );

  double Msn = Dstar*Ss;
  double mn  = Msr*n[0] + Msp*n[1] + Msz*n[2];

  Msr += n[0]*( Msn - mn );
  Msp += n[1]*( Msn - mn );
  Msz += n[2]*( Msn - mn );

  theFlux->Ustar[DDD] = Dstar;
  theFlux->Ustar[SRR] = Msr;
  theFlux->Ustar[LLL] = r*Msp;
  theFlux->Ustar[SZZ] = Msz;
  theFlux->Ustar[TAU] = Estar;

  theFlux->Ustar[BRR] = Bsr/r;
  theFlux->Ustar[BPP] = Bsp/r;
  theFlux->Ustar[BZZ] = Bsz;

  theFlux->Ustar[PSI] = psi;


}

void flux_set_flux(struct Flux * theFlux, double r , double * n ,double GAMMALAW,double DIVB_CH, double w){
  double *prim;
  if (w<(theFlux->Ss)){
    prim = theFlux->primL;
  }else{
    prim = theFlux->primR;
  }

  double rho = prim[RHO];
  double Pp  = prim[PPP];
  double vr  = prim[URR];
  double vp  = prim[UPP]*r;
  double vz  = prim[UZZ];

  double Br  = prim[BRR];
  double Bp  = prim[BPP];
  double Bz  = prim[BZZ];

  double vn = vr*n[0] + vp*n[1] + vz*n[2];
  double Bn = Br*n[0] + Bp*n[1] + Bz*n[2];
  double vB = vr*Br + vp*Bp + vz*Bz;

  double rhoe = Pp/(GAMMALAW-1.);
  double v2 = vr*vr + vp*vp + vz*vz;
  double B2 = Br*Br + Bp*Bp + Bz*Bz;

  theFlux->F[DDD] = rho*vn;
  theFlux->F[SRR] =     rho*vr*vn + (Pp+.5*B2)*n[0] - Br*Bn;
  theFlux->F[LLL] = r*( rho*vp*vn + (Pp+.5*B2)*n[1] - Bp*Bn );
  theFlux->F[SZZ] =     rho*vz*vn + (Pp+.5*B2)*n[2] - Bz*Bn;
  theFlux->F[TAU] = ( .5*rho*v2 + rhoe + Pp + B2 )*vn - vB*Bn;

  double psi = prim[PSI];
  theFlux->F[BRR] =(Br*vn - vr*Bn + psi*n[0])/r;
  theFlux->F[BPP] =(Bp*vn - vp*Bn + psi*n[1])/r;
  theFlux->F[BZZ] = Bz*vn - vz*Bn + psi*n[2];

  double wp = 0.0;
  theFlux->F[PSI] = wp*psi*n[1] + pow(DIVB_CH,2.)*Bn;

  /*
     int q;
     for( q=NUM_C ; q<NUM_Q ; ++q ){
     flux[q] = prim[q]*flux[DDD];
     }
     */

}

void flux_addto_flux(struct Flux * theFlux,double w,int NUM_Q){
  double Sk;
  if (w<(theFlux->Ss)){
    Sk = theFlux->Sl;
  }else{
    Sk = theFlux->Sr;
  }
  int q;
  for( q=0 ; q<NUM_Q ; ++q ){
    theFlux->F[q] += Sk*( theFlux->Ustar[q] - theFlux->Uk[q] ) - w*theFlux->Ustar[q];
  }

}



void flux_set_Uk(struct Flux *theFlux,double * Uk_in){
  theFlux->Uk = Uk_in; 
}

double * flux_primL(struct Flux * theFlux){
  return(theFlux->primL);
}
double * flux_primR(struct Flux * theFlux){
  return(theFlux->primR);
}