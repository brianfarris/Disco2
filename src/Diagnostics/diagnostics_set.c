#define DIAGNOSTICS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/GravMass.h"
#include "../Headers/Metric.h"
#include "../Headers/Diagnostics.h"
#include "../Headers/TimeStep.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

void diagnostics_set(struct Diagnostics * theDiagnostics,struct Cell *** theCells,struct Sim * theSim,struct TimeStep * theTimeStep,struct MPIsetup * theMPIsetup,struct GravMass * theGravMasses){
  if (timestep_get_t(theTimeStep)>=diagnostics_tdiag_measure(theDiagnostics)){
    int num_r_points_global = sim_N_global(theSim,R_DIR);

    int NUM_SCAL = theDiagnostics->NUM_DIAG;
    int NUM_VEC = theDiagnostics->NUM_DIAG+1;
    int NUM_EQ = theDiagnostics->NUM_DIAG+2;

    int i,j,k,n;

    double * EquatDiag_temp = malloc(sizeof(double) * theDiagnostics->N_eq_cells*NUM_EQ);
    double * EquatDiag_reduce = malloc(sizeof(double) * theDiagnostics->N_eq_cells*NUM_EQ);
    double * VectorDiag_temp = malloc(sizeof(double) * num_r_points_global*NUM_VEC);
    double * VectorDiag_reduce = malloc(sizeof(double) * num_r_points_global*NUM_VEC);
    double * ScalarDiag_temp = malloc(sizeof(double)*NUM_SCAL);
    double * ScalarDiag_reduce = malloc(sizeof(double)*NUM_SCAL);

    double dtout = timestep_get_t(theTimeStep)-theDiagnostics->toutprev;

    int imin = sim_Nghost_min(theSim,R_DIR);
    int imax = sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
    int kmin = sim_Nghost_min(theSim,Z_DIR);
    int kmax = sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR);

    for (n=0;n<NUM_SCAL;++n){
      ScalarDiag_temp[n]=0.0;
      ScalarDiag_reduce[n]=0.0;
    }
    for (i=0;i<num_r_points_global;++i){
      for (n=0;n<NUM_VEC;++n){
        VectorDiag_temp[i*NUM_VEC+n]=0.0;
      }
    }

    int position=0;
    for (i=0;i<num_r_points_global;++i){
      for(j = 0; j < theDiagnostics->N_p_global[i]; j++){
        for (n=0;n<NUM_EQ;++n){
          EquatDiag_temp[position]=0.0;
          ++position;
        }
      }
    }

    position=0;
    for (k=kmin;k<kmax;++k){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double z = 0.5*(zm+zp);
      double dz = zp-zm;
      for (i=imin;i<imax;++i){
        double rp = sim_FacePos(theSim,i,R_DIR);
        double rm = sim_FacePos(theSim,i-1,R_DIR);
        double r = 0.5*(rm+rp);
        for (j=0;j<sim_N_p(theSim,i);++j){
          double phi = cell_tiph(cell_single(theCells,i,j,k));
          double rho = cell_prim(cell_single(theCells,i,j,k),RHO);
          double press = cell_prim(cell_single(theCells,i,j,k),PPP);
          double vr = cell_prim(cell_single(theCells,i,j,k),URR);
          double vp = cell_prim(cell_single(theCells,i,j,k),UPP)*r;
          double vz = cell_prim(cell_single(theCells,i,j,k),UZZ);

          double passive_scalar;
          int NUMQ = sim_NUM_Q(theSim);
          int NUMC = sim_NUM_C(theSim);
          if(NUMQ > NUMC)
            passive_scalar = cell_prim(cell_single(theCells,i,j,k),NUMC);
          else
            passive_scalar = 0.0;

          if ((fabs(zp)<0.0000001)||(fabs(z)<0.0000001)){
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+0] = r;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+1] = phi;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+2] = rho;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+3] = vr;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+4] = vp;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+5] = press;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+6] = passive_scalar;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+7] = rho*vr;
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+8] = rho*vp;            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+9] = rho*vr*cos(phi);            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+10] = rho*vr*sin(phi);            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+11] = rho*vr*cos(2.*phi);            
            EquatDiag_temp[(theDiagnostics->offset_eq+position)*NUM_EQ+12] = rho*vr*sin(2.*phi);            
            ++position;
          }
          // divide by number of phi cells to get phi average, mult by dz because we are doing a z integration;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+0] += (r/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+1] += (rho/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+2] += (vr/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+3] += (vp/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+4] += (press/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+5] += (passive_scalar/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+6] += (rho*vr/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+7] += (rho*vp/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+8] += (rho*vr*cos(phi)/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+9] += (rho*vr*sin(phi)/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+10] += (rho*vr*cos(2.*phi)/sim_N_p(theSim,i)*dz) ;
          VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+11] += (rho*vr*sin(2.*phi)/sim_N_p(theSim,i)*dz) ;

          if(sim_Background(theSim) != NEWTON)
          {
              struct Metric *g = metric_create(time_global, r, phi, z, theSim);
              double tmp_prim[5];
              tmp_prim[RHO] = rho;
              tmp_prim[PPP] = press;
              tmp_prim[URR] = vr;
              tmp_prim[UPP] = vp/r;
              tmp_prim[UZZ] = vz;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+12] += (eos_cs2(tmp_prim,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+13] += (eos_eps(tmp_prim,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+14] += (eos_ppp(tmp_prim,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+15] += 
                  (metric_frame_U_u(g,0,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+16] += 
                  (metric_frame_U_u(g,1,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+17] += 
                  (metric_frame_U_u(g,2,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+18] += 
                  (metric_frame_dU_du(g,1,0,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+19] += 
                  (metric_frame_dU_du(g,1,1,theSim)/sim_N_p(theSim,i)*dz) ;
              VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+20] += 
                  (metric_frame_dU_du(g,1,2,theSim)/sim_N_p(theSim,i)*dz) ;

              metric_destroy(g);
          }
          // the above are just placeholders. Put the real diagnostics you want here, then adjust NUM_DIAG accordingly.
        }
      }
    }

    for (i=imin;i<imax;++i){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      for (n=0;n<NUM_SCAL;++n){
        // mult by delta r^2 because we are doing an r integration
        ScalarDiag_temp[n] += VectorDiag_temp[(sim_N0(theSim,R_DIR)+i-imin)*NUM_VEC+n+1] * (rp*rp-rm*rm); 
      }
    }

    MPI_Allreduce( ScalarDiag_temp,ScalarDiag_reduce , NUM_SCAL, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( VectorDiag_temp,VectorDiag_reduce , num_r_points_global*NUM_VEC, MPI_DOUBLE, MPI_SUM, sim_comm);
    MPI_Allreduce( EquatDiag_temp,EquatDiag_reduce ,theDiagnostics->N_eq_cells*NUM_EQ, MPI_DOUBLE, MPI_SUM, sim_comm);

    double RMIN = sim_MIN(theSim,R_DIR);
    double RMAX = sim_MAX(theSim,R_DIR);
    double ZMIN = sim_MIN(theSim,Z_DIR);
    double ZMAX = sim_MAX(theSim,Z_DIR);

    for (n=0;n<NUM_SCAL;++n){
      ScalarDiag_reduce[n] *= dtout/((ZMAX-ZMIN)*(RMAX*RMAX-RMIN*RMIN));
    }

    for (i=0;i<num_r_points_global;++i){
      for (n=0;n<NUM_VEC;++n){
        VectorDiag_reduce[i*NUM_VEC+n] *= dtout/(ZMAX-ZMIN);
      }
    }

    //We are doing time averaged diagnostics, so mult by delta t and add it
    //We will divide by the total delta next time we save to disk;
    for (i=0;i<num_r_points_global;++i){
      for (n=0;n<NUM_VEC;++n){
        theDiagnostics->VectorDiag[i][n] += VectorDiag_reduce[i*NUM_VEC+n] ;
      }
    }
    for (n=0;n<NUM_SCAL;++n){
      theDiagnostics->ScalarDiag[n] += ScalarDiag_reduce[n] ;
    }

    position=0;
    for (i=0;i<num_r_points_global;++i){
      for(j = 0; j < theDiagnostics->N_p_global[i]; j++){
        for (n=0;n<NUM_EQ;++n){
          theDiagnostics->EquatDiag[position][n] = EquatDiag_reduce[position*NUM_EQ+n];
        }
        ++position;
      }
    }

    //update output time;
    theDiagnostics->toutprev = timestep_get_t(theTimeStep);

    free(ScalarDiag_temp);
    free(ScalarDiag_reduce);
    free(VectorDiag_temp);
    free(VectorDiag_reduce);
    free(EquatDiag_temp);
    free(EquatDiag_reduce);

    theDiagnostics->tdiag_measure += theDiagnostics->dtdiag_measure;
  } 
}






