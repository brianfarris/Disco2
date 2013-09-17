#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"

void cell_EMF(struct Cell *** theCells, struct Sim *theSim,double dt ){
  int NUM_Q = sim_NUM_Q(theSim);
  double GAMMALAW = sim_GAMMALAW(theSim);

  int i,j,k;
  for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double dz = zp-zm;
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double dr = rp-rm;
      double r = 0.5*(rm+rp);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){

        struct Cell * theCell = cell_single(theCells,i,j,k);

        double dp = cell_dphi(theCell);

        double * prim_minsr = malloc(NUM_Q*sizeof(double));
        double * prim_plusr = malloc(NUM_Q*sizeof(double));
        double * prim_minsp = malloc(NUM_Q*sizeof(double));
        double * prim_plusp = malloc(NUM_Q*sizeof(double));
        double * prim_minsz = malloc(NUM_Q*sizeof(double));
        double * prim_plusz = malloc(NUM_Q*sizeof(double));

        int q;
        for (q=0;q<NUM_Q;++q){
          prim_minsp[q] = theCell->prim[q] - 0.5*cell_gradp(theCell,q)*dp;
          prim_plusp[q] = theCell->prim[q] + 0.5*cell_gradp(theCell,q)*dp;
          prim_minsr[q] = theCell->prim[q] - 0.5*cell_gradr(theCell,q)*dr;
          prim_plusr[q] = theCell->prim[q] + 0.5*cell_gradr(theCell,q)*dr;
          prim_minsz[q] = theCell->prim[q] - 0.5*cell_gradz(theCell,q)*dz;
          prim_plusz[q] = theCell->prim[q] + 0.5*cell_gradz(theCell,q)*dz;
        }

        double FpBz_plusp = r*prim_plusp[UPP]*prim_plusp[BZZ]-prim_plusp[UZZ]*prim_plusp[BPP];
        double FpBz_minsp = r*prim_minsp[UPP]*prim_minsp[BZZ]-prim_minsp[UZZ]*prim_minsp[BPP];
        double FzBp_plusz = prim_plusz[UZZ]*prim_plusz[BPP]-r*prim_plusz[UPP]*prim_plusz[BZZ]; 
        double FzBp_minsz = prim_minsz[UZZ]*prim_minsz[BPP]-r*prim_minsz[UPP]*prim_minsz[BZZ];

        double FrBz_plusr = prim_plusr[URR]*prim_plusr[BZZ]-prim_plusr[UZZ]*prim_plusr[BRR];
        double FrBz_minsr = prim_minsr[URR]*prim_minsr[BZZ]-prim_minsr[UZZ]*prim_minsr[BRR];
        double FzBr_plusz = prim_plusz[UZZ]*prim_plusz[BRR]-prim_plusz[URR]*prim_plusz[BZZ];
        double FzBr_minsz = prim_minsz[UZZ]*prim_minsz[BRR]-prim_minsz[URR]*prim_minsz[BZZ];

        double FrBp_plusr = prim_plusr[URR]*prim_plusr[BPP]-rp*prim_plusr[UPP]*prim_plusr[BRR];
        double FrBp_minsr = prim_minsr[URR]*prim_minsr[BPP]-rm*prim_minsr[UPP]*prim_minsr[BRR];
        double FpBr_plusp = rp*prim_plusp[UPP]*prim_plusp[BRR]-prim_plusp[URR]*prim_plusp[BPP];
        double FpBr_minsp = rm*prim_minsp[UPP]*prim_minsp[BRR]-prim_minsp[URR]*prim_minsp[BPP];

        double Er = 0.25*(-FpBz_plusp - FpBz_minsp + FzBp_plusz + FzBp_minsz);
        double Ep = 0.25*( FrBz_plusr + FrBz_minsr - FzBr_plusz - FzBr_minsz);
        double Ez = 0.25*(-FrBp_plusr - FrBp_minsr + FpBr_plusp + FpBr_minsp);

        //printf("prim_minsp[BRR]: %e\n",prim_minsp[BRR]);
        //printf("cell_gradp(theCell,BRR): %e\n",cell_gradp(theCell,BRR));
        //printf("theCell->prim[BRR]: %e\n",theCell->prim[BRR]);
        double By = theCell->prim[BRR]*sin(cell_tiph(theCell)) +  theCell->prim[BPP]*cos(cell_tiph(theCell));
        double vx = theCell->prim[URR]*cos(cell_tiph(theCell)) -  theCell->prim[UPP]*r*sin(cell_tiph(theCell));
        //printf("By: %e, vx: %e, -20y: %e, Ez: %e\n",By,vx,-20.*r*sin(cell_tiph(theCell)+0.5*cell_dphi(theCell)),Ez);
        //printf("By: %e, Ez: %e, -vx: %e\n",By,Ez,20*r*sin(cell_tiph(theCell)+0.5*cell_dphi(theCell)));
        //printf("FrPp_plusr: %e, FrBp_minsr: %e, FpBr_plusp: %e, FpBr_minsp: %e, Ez: %e, -vx: %e\n",FrBp_plusr,FrBp_minsr, FpBr_plusp, FpBr_minsp,Ez,20*r*sin(cell_tiph(theCell)+0.5*cell_dphi(theCell)));
        free(prim_minsr);
        free(prim_plusr);
        free(prim_minsp);
        free(prim_plusp);
        free(prim_minsz);
        free(prim_plusz);

        theCell->prim[ARR] += Er * dt;
        theCell->prim[APP] += Ep * dt;
        theCell->prim[AZZ] += Ez * dt;

      }
    }
  }
}

