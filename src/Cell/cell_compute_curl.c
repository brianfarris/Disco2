#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


void set_plus_mins(double * phi_pm, double * Ap_pm, double * Az_pm, struct Cell *** theCells, struct Sim * theSim,int next_to_last_pm, int last_pm,int i,int k,int shift){
  phi_pm[0] = theCells[k][i+shift][next_to_last_pm].tiph - 0.5*theCells[k][i+shift][next_to_last_pm].dphi;
  phi_pm[1] = theCells[k][i+shift][last_pm].tiph - 0.5*theCells[k][i+shift][next_to_last_pm].dphi;
  Ap_pm[0] = theCells[k][i+shift][next_to_last_pm].prim[APP];
  Ap_pm[1] = theCells[k][i+shift][last_pm].prim[APP];
  Az_pm[0] = theCells[k][i+shift][next_to_last_pm].prim[AZZ];
  Az_pm[1] = theCells[k][i+shift][last_pm].prim[AZZ];
  int j;
  for( j=0 ; j<sim_N_p(theSim,i+shift) ; ++j ){
    phi_pm[j+2] = theCells[k][i+shift][j].tiph - 0.5*theCells[k][i+shift][j].dphi;
    Ap_pm[j+2] = theCells[k][i+shift][j].prim[APP];
    Az_pm[j+2] = theCells[k][i+shift][j].prim[AZZ];
  }
  phi_pm[sim_N_p(theSim,i+shift)+2] = theCells[k][i+shift][0].tiph - 0.5*theCells[k][i+shift][0].dphi;
  phi_pm[sim_N_p(theSim,i+shift)+3] = theCells[k][i+shift][1].tiph - 0.5*theCells[k][i+shift][1].dphi;
  Ap_pm[sim_N_p(theSim,i+shift)+2] = theCells[k][i+shift][0].prim[APP];
  Ap_pm[sim_N_p(theSim,i+shift)+3] = theCells[k][i+shift][1].prim[APP];
  Az_pm[sim_N_p(theSim,i+shift)+2] = theCells[k][i+shift][0].prim[AZZ];
  Az_pm[sim_N_p(theSim,i+shift)+3] = theCells[k][i+shift][1].prim[AZZ];

  for ( j=1 ; j<sim_N_p(theSim,i+shift)+4 ; ++j ){
    while (phi_pm[j]<phi_pm[j-1]) phi_pm[j] += 2.*M_PI;
    //printf("j: %d, phi_pm: %e, Ap_pm: %e\n",j,phi_pm[j],Ap_pm[j]);
  }

}

void cell_compute_curl( struct Cell *** theCells ,struct Sim * theSim , struct MPIsetup * theMPIsetup){
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
      double rmm = sim_FacePos(theSim,i-2,R_DIR);
      double rpp = sim_FacePos(theSim,i+1,R_DIR);
      double r_plus = 0.5*(rp+rpp);
      double r_mins = 0.5*(rm+rmm);
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

      if (1==0){
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
      } else{
        set_plus_mins(phi_plus, Ap_plus, Az_plus, theCells, theSim,next_to_last_plus, last_plus,i,k,1);
      }

      if (1==0){
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
      } else{
        set_plus_mins(phi_mins, Ap_mins, Az_mins, theCells, theSim,next_to_last_mins, last_mins,i,k,-1);
      }

      gsl_interp *interp_Ap_plus;
      gsl_interp *interp_Az_plus;
      gsl_interp *interp_Ap_mins;
      gsl_interp *interp_Az_mins;

      int interp_method=LINEAR;

      if (interp_method==LINEAR){
        interp_Ap_plus = gsl_interp_alloc(gsl_interp_linear,sim_N_p(theSim,i+1)+4);
        interp_Az_plus = gsl_interp_alloc(gsl_interp_linear,sim_N_p(theSim,i+1)+4);
        interp_Ap_mins = gsl_interp_alloc(gsl_interp_linear,sim_N_p(theSim,i-1)+4);
        interp_Az_mins = gsl_interp_alloc(gsl_interp_linear,sim_N_p(theSim,i-1)+4);
      } else if (interp_method==SPLINE){
        interp_Ap_plus = gsl_interp_alloc(gsl_interp_cspline,sim_N_p(theSim,i+1)+4);
        interp_Az_plus = gsl_interp_alloc(gsl_interp_cspline,sim_N_p(theSim,i+1)+4);
        interp_Ap_mins = gsl_interp_alloc(gsl_interp_cspline,sim_N_p(theSim,i-1)+4);
        interp_Az_mins = gsl_interp_alloc(gsl_interp_cspline,sim_N_p(theSim,i-1)+4);
      } else{
        printf("invalid interpolation method\n");
        exit(1);
      }

      gsl_interp_accel *acc_plus = gsl_interp_accel_alloc ();
      gsl_interp_accel *acc_mins = gsl_interp_accel_alloc ();

      gsl_interp_init (interp_Ap_plus, phi_plus, Ap_plus, sim_N_p(theSim,i+1)+4);
      gsl_interp_init (interp_Az_plus, phi_plus, Az_plus, sim_N_p(theSim,i+1)+4);
      gsl_interp_init (interp_Ap_mins, phi_mins, Ap_mins, sim_N_p(theSim,i-1)+4);
      gsl_interp_init (interp_Az_mins, phi_mins, Az_mins, sim_N_p(theSim,i-1)+4);

      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double two_dphi = 2.*theCells[k][i][j].dphi;
        double dr = rp-rm;
        double two_dr = 2.*(rp-rm);
        double dz = zp-zm;
        double two_dz = 2.*(zp-zm);
        struct Cell * c = &(theCells[k][i][j]);
        double phi = c->tiph - 0.5*c->dphi;
        while (phi<phi_plus[0]) phi += 2.*M_PI;
        double Ap_interp_r_plus = gsl_interp_eval(interp_Ap_plus,phi_plus, Ap_plus,phi,acc_plus);
        double Az_interp_r_plus = gsl_interp_eval(interp_Az_plus,phi_plus, Az_plus,phi,acc_plus);
        phi = c->tiph - 0.5*c->dphi;
        while (phi<phi_mins[0]) phi += 2.*M_PI;
        double Ap_interp_r_mins = gsl_interp_eval(interp_Ap_mins,phi_mins, Ap_mins,phi,acc_mins); 
        double Az_interp_r_mins = gsl_interp_eval(interp_Az_mins,phi_mins, Az_mins,phi,acc_mins);
      
        double one_o_r_drAp_dr;
        double dAz_dr;
        if (i==imin && mpisetup_check_rin_bndry(theMPIsetup)){
          one_o_r_drAp_dr = (Ap_interp_r_plus*r_plus - c->prim[APP]*r)/dr/r;
          dAz_dr = (Az_interp_r_plus - c->prim[AZZ])/dr;
        } else if (i==imax-1 && mpisetup_check_rout_bndry(theMPIsetup)){
          one_o_r_drAp_dr = (c->prim[APP]*r - Ap_interp_r_mins*r_mins)/dr/r;
          dAz_dr = (c->prim[AZZ] - Az_interp_r_mins)/dr;
        } else{
          one_o_r_drAp_dr = (Ap_interp_r_plus*r_plus - Ap_interp_r_mins*r_mins)/two_dr/r;
          dAz_dr = (Az_interp_r_plus - Az_interp_r_mins)/two_dr;
        }
        //get phi derivatives of Ar and Az
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

        //get z derivatives of Ar and Ap. 
        //Assume that cells at different heights with same j lie on top of each other
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

        theCells[k][i][j].prim[BRR] = dAz_dp/r - dAp_dz;
        theCells[k][i][j].prim[BPP] = dAr_dz - dAz_dr;
        theCells[k][i][j].prim[BZZ] = one_o_r_drAp_dr - dAr_dp/r;// + theCells[k][i][j].prim[APP]/r;

      }
      gsl_interp_free (interp_Ap_mins);
      gsl_interp_free (interp_Az_mins);
      gsl_interp_free (interp_Ap_plus);
      gsl_interp_free (interp_Az_plus);
      gsl_interp_accel_free (acc_mins);
      gsl_interp_accel_free (acc_plus);
      //if(k==4){
      //  printf("%e %e\n",r,theCells[k][i][0].prim[BZZ]);
      //}
    }
  } 
}
