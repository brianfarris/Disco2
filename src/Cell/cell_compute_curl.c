#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/header.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


void cell_compute_curl( struct Cell *** theCells ,struct Sim * theSim ){
  double Mtotal = 1.0;
  double sep = 1.0;
  
  int imin = sim_Nghost_min(theSim,R_DIR);
  int kmin = sim_Nghost_min(theSim,Z_DIR);
  int imax = sim_N(theSim,R_DIR) - sim_Nghost_max(theSim,R_DIR);
  int kmax = sim_N(theSim,Z_DIR) - sim_Nghost_max(theSim,Z_DIR);

  int i,j,k;
  for( k=kmin ; k<kmax ; ++k ){
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    for( i=imin ; i<imax ; ++i ){
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      double rp = sim_FacePos(theSim,i,R_DIR);
      double r = 0.5*(rm+rp);
      double * phi_plus = malloc((sim_N_p(theSim,i+1)+4)*sizeof(double));
      double * Ar_plus = malloc((sim_N_p(theSim,i+1)+4)*sizeof(double));
      double * Ap_plus = malloc((sim_N_p(theSim,i+1)+4)*sizeof(double));
      double * Az_plus = malloc((sim_N_p(theSim,i+1)+4)*sizeof(double));
      double * phi_mins = malloc((sim_N_p(theSim,i-1)+4)*sizeof(double));
      double * Ar_mins = malloc((sim_N_p(theSim,i-1)+4)*sizeof(double));
      double * Ap_mins = malloc((sim_N_p(theSim,i-1)+4)*sizeof(double));
      double * Az_mins = malloc((sim_N_p(theSim,i-1)+4)*sizeof(double));
      double * phi_interp = malloc(sim_N_p(theSim,i)*sizeof(double));
      double * Ar_plus_interp = malloc(sim_N_p(theSim,i)*sizeof(double));
      double * Ap_plus_interp = malloc(sim_N_p(theSim,i)*sizeof(double));
      double * Az_plus_interp = malloc(sim_N_p(theSim,i)*sizeof(double));
      double * Ar_mins_interp = malloc(sim_N_p(theSim,i)*sizeof(double));
      double * Ap_mins_interp = malloc(sim_N_p(theSim,i)*sizeof(double));
      double * Az_mins_interp = malloc(sim_N_p(theSim,i)*sizeof(double));

      int next_to_last_plus = sim_N_p(theSim,i+1)-2;
      int last_plus = sim_N_p(theSim,i+1)-1;
      int next_to_last_mins = sim_N_p(theSim,i-1)-2;
      int last_mins = sim_N_p(theSim,i-1)-1;

      phi_plus[0] = theCells[k][i+1][next_to_last_plus].tiph - 0.5*theCells[k][i+1][next_to_last_plus].dphi;
      phi_plus[1] = theCells[k][i+1][last_plus].tiph - 0.5*theCells[k][i+1][next_to_last_plus].dphi;
      Ap_plus[0] = theCells[k][i+1][next_to_last_plus].prim[APP];
      Ap_plus[1] = theCells[k][i+1][last_plus].prim[APP];
      Az_plus[0] = theCells[k][i+1][next_to_last_plus].prim[AZZ];
      Az_plus[1] = theCells[k][i+1][last_plus].prim[AZZ];
      for( j=0 ; j<sim_N_p(theSim,i+1) ; ++j ){
        phi_plus[j+2] = theCells[k][i+1][j].tiph - 0.5*theCells[k][i+1][j].dphi;
        Ap_plus[j+2] = theCells[k][i+1][j].prim[APP];
        Az_plus[j+2] = theCells[k][i+1][j].prim[AZZ];
      }
      phi_plus[sim_N_p(theSim,i+1)+2] = theCells[k][i+1][0].tiph - 0.5*theCells[k][i+1][0].dphi;
      phi_plus[sim_N_p(theSim,i+1)+3] = theCells[k][i+1][1].tiph - 0.5*theCells[k][i+1][1].dphi;
      Ap_plus[sim_N_p(theSim,i+1)+2] = theCells[k][i+1][0].prim[APP];
      Ap_plus[sim_N_p(theSim,i+1)+3] = theCells[k][i+1][1].prim[APP];
      Az_plus[sim_N_p(theSim,i+1)+2] = theCells[k][i+1][0].prim[AZZ];
      Az_plus[sim_N_p(theSim,i+1)+3] = theCells[k][i+1][1].prim[AZZ];

      for ( j=1 ; j<sim_N_p(theSim,i+1)+4 ; ++j ){
        while (phi_plus[j]<phi_plus[j-1]) phi_plus[j] += 2.*M_PI;
      }

      phi_mins[0] = theCells[k][i-1][next_to_last_mins].tiph - 0.5*theCells[k][i-1][next_to_last_mins].dphi;
      phi_mins[1] = theCells[k][i-1][last_mins].tiph - 0.5*theCells[k][i-1][next_to_last_mins].dphi;
      Ap_mins[0] = theCells[k][i-1][next_to_last_mins].prim[APP];
      Ap_mins[1] = theCells[k][i-1][last_mins].prim[APP];
      Az_mins[0] = theCells[k][i-1][next_to_last_mins].prim[AZZ];
      Az_mins[1] = theCells[k][i-1][last_mins].prim[AZZ];
      for( j=0 ; j<sim_N_p(theSim,i-1) ; ++j ){
        phi_mins[j+2] = theCells[k][i-1][j].tiph - 0.5*theCells[k][i-1][j].dphi;
        Ap_mins[j+2] = theCells[k][i-1][j].prim[APP];
        Az_mins[j+2] = theCells[k][i-1][j].prim[AZZ];
      }
      phi_mins[sim_N_p(theSim,i-1)+2] = theCells[k][i-1][0].tiph - 0.5*theCells[k][i-1][0].dphi;
      phi_mins[sim_N_p(theSim,i-1)+3] = theCells[k][i-1][1].tiph - 0.5*theCells[k][i-1][1].dphi;
      Ap_mins[sim_N_p(theSim,i-1)+2] = theCells[k][i-1][0].prim[APP];
      Ap_mins[sim_N_p(theSim,i-1)+3] = theCells[k][i-1][1].prim[APP];
      Az_mins[sim_N_p(theSim,i-1)+2] = theCells[k][i-1][0].prim[AZZ];
      Az_mins[sim_N_p(theSim,i-1)+3] = theCells[k][i-1][1].prim[AZZ];

      for ( j=1 ; j<sim_N_p(theSim,i-1)+4 ; ++j ){
        while (phi_mins[j]<phi_mins[j-1]) phi_mins[j] += 2.*M_PI;
      }

      gsl_interp_accel *acc_plus = gsl_interp_accel_alloc ();
      gsl_spline *spline_Ap_plus = gsl_spline_alloc (gsl_interp_cspline_periodic, sim_N_p(theSim,i+1)+4);
      gsl_spline *spline_Az_plus = gsl_spline_alloc (gsl_interp_cspline_periodic, sim_N_p(theSim,i+1)+4);
      gsl_spline_init (spline_Ap_plus, phi_plus, Ap_plus, sim_N_p(theSim,i+1)+4);
      gsl_spline_init (spline_Az_plus, phi_plus, Az_plus, sim_N_p(theSim,i+1)+4);

      gsl_interp_accel *acc_mins = gsl_interp_accel_alloc ();
      gsl_spline *spline_Ap_mins = gsl_spline_alloc (gsl_interp_cspline_periodic, sim_N_p(theSim,i-1)+4);
      gsl_spline *spline_Az_mins = gsl_spline_alloc (gsl_interp_cspline_periodic, sim_N_p(theSim,i-1)+4);
      gsl_spline_init (spline_Ap_mins, phi_mins, Ap_mins, sim_N_p(theSim,i-1)+4);
      gsl_spline_init (spline_Az_mins, phi_mins, Az_mins, sim_N_p(theSim,i-1)+4);

      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double two_dphi = 2.*theCells[k][i][j].dphi;
        double two_dr = 2.*(rp-rm);
        double two_dz = 2.*(zp-zm);
        struct Cell * c = &(theCells[k][i][j]);
        double phi = c->tiph - 0.5*c->dphi;
        while (phi<phi_plus[0]) phi += 2.*M_PI;
        double Ap_interp_r_plus = gsl_spline_eval(spline_Ap_plus,phi,acc_plus);
        double Az_interp_r_plus = gsl_spline_eval(spline_Az_plus,phi,acc_plus);
        phi = c->tiph - 0.5*c->dphi;
        while (phi<phi_mins[0]) phi += 2.*M_PI;
        double Ap_interp_r_mins = gsl_spline_eval(spline_Ap_mins,phi,acc_mins); 
        double Az_interp_r_mins = gsl_spline_eval(spline_Az_mins,phi,acc_mins);
        double dAp_dr = (Ap_interp_r_plus - Ap_interp_r_mins)/two_dr;
        //printf("Ap_interp_r_plus: %e, Ap_interp_r_mins: %e, dAp_dr: %e\n",Ap_interp_r_plus,Ap_interp_r_mins,dAp_dr);
        double dAz_dr = (Az_interp_r_plus - Az_interp_r_mins)/two_dr;
        //printf("Az_interp_r_plus: %e, Az_interp_r_mins: %e, dAz_dr: %e\n",Az_interp_r_plus,Az_interp_r_mins,dAz_dr);

        double Ar_interp_p_mins;
        double Ar_interp_p_plus;
        double Az_interp_p_mins;
        double Az_interp_p_plus;
        if (j==0){
          Ar_interp_p_mins = theCells[k][i][sim_N_p(theSim,i)-1].prim[ARR];
          Az_interp_p_mins = theCells[k][i][sim_N_p(theSim,i)-1].prim[AZZ];
        } else{
          Ar_interp_p_mins = theCells[k][i][j-1].prim[ARR];
          Az_interp_p_mins = theCells[k][i][j-1].prim[AZZ];
        }          
        if (j==sim_N_p(theSim,i)-1){
          Ar_interp_p_plus = theCells[k][i][0].prim[ARR];
          Az_interp_p_plus = theCells[k][i][0].prim[AZZ];
        } else{
          Ar_interp_p_plus = theCells[k][i][j+1].prim[ARR];
          Az_interp_p_plus = theCells[k][i][j+1].prim[AZZ];
        }
        double dAr_dp = (Ar_interp_p_plus - Ar_interp_p_mins)/two_dphi;
        double dAz_dp = (Az_interp_p_plus - Az_interp_p_mins)/two_dphi;
        double dAr_dz;
        double dAp_dz;
        if (sim_N(theSim,Z_DIR)>1){
          double Ar_interp_z_mins = theCells[k-1][i][j].prim[ARR];
          double Ap_interp_z_mins = theCells[k-1][i][j].prim[APP];
          double Ar_interp_z_plus = theCells[k+1][i][j].prim[ARR];
          double Ap_interp_z_plus = theCells[k+1][i][j].prim[APP];
          dAr_dz = (Ar_interp_z_plus - Ar_interp_z_mins)/two_dz;
          dAp_dz = (Ap_interp_z_plus - Ap_interp_z_mins)/two_dz;
        } else{
          dAr_dz = 0.0;
          dAp_dz = 0.0;
        }

//        theCells[k][i][j].prim[BRR] = dAz_dp/r - dAp_dz;
//        theCells[k][i][j].prim[BPP] = dAr_dz - dAz_dr;
//        theCells[k][i][j].prim[BZZ] = dAp_dr - dAr_dp/r + theCells[k][i][j].prim[APP];

        theCells[k][i][j].prim[BRR] = dAz_dp/r - dAp_dz;
        theCells[k][i][j].prim[BPP] = dAr_dz - dAz_dr;
        theCells[k][i][j].prim[BZZ] = dAp_dr - dAr_dp/r + theCells[k][i][j].prim[APP]/r;

        //printf("i: %d, j: %d, k: %d, dAz_dp: %e, dAp_dz: %e, dAr_dz: %e, dAz_dr: %e, dAp_dr: %e, dAr_dp: %e, A_phi: %e\n",i,j,k,dAz_dp,dAp_dz,dAr_dz,dAz_dr,dAp_dr,dAr_dp,theCells[k][i][j].prim[APP]);
        //printf("i: %d, j: %d, k: %d, r: %e, Bz: %e\n",i,j,k,r,theCells[k][i][j].prim[BZZ]);
       // theCells[k][i][j].cons[BRR] = theCells[k][i][j].prim[BRR]*dV/r;
       // theCells[k][i][j].cons[BPP] = theCells[k][i][j].prim[BPP]*dV/r;
       // theCells[k][i][j].cons[BZZ] = theCells[k][i][j].prim[BZZ]*dV;
        

      }
      gsl_spline_free (spline_Ap_mins);
      gsl_spline_free (spline_Az_mins);
      gsl_spline_free (spline_Ap_plus);
      gsl_spline_free (spline_Az_plus);
      gsl_interp_accel_free (acc_mins);
      gsl_interp_accel_free (acc_plus);

    }
  } 
}
