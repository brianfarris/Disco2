#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"


void cell_boundary_outflow_r( struct Cell *** theCells , struct Face * theFaces ,struct Sim * theSim,struct MPIsetup * theMPIsetup, struct TimeStep * theTimeStep ){
  int Nf = timestep_nri(theTimeStep,sim_N(theSim,R_DIR)-1);
  int NUM_Q = sim_NUM_Q(theSim);
  int n1 = timestep_nri(theTimeStep,sim_N(theSim,R_DIR)-2);

  int n,q;
  int i,j,k;
  double r;

  if( mpisetup_check_rin_bndry(theMPIsetup) ){

    for( i=0 ; i>=0 ; --i ){
      r=sim_FacePos(theSim,i,R_DIR);
      for( n=timestep_nri(theTimeStep,i) ; n<timestep_nri(theTimeStep,i+1) ; ++n ){
        for( q=0 ; q<NUM_Q ; ++q ){
          face_L_pointer(theFaces,n)->prim[q] = 0.0;
        }
      } 
      for( n=timestep_nri(theTimeStep,i) ; n<timestep_nri(theTimeStep,i+1) ; ++n ){
        struct Cell * cL = face_L_pointer(theFaces,n);
        struct Cell * cR = face_R_pointer(theFaces,n);
        for( q=0 ; q<NUM_Q ; ++q ){
          cL->prim[q] += cR->prim[q]*face_dA(theFaces,n);
        }
      }

      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double dz = zp-zm;
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          double dA = dz*r*theCells[k][i][j].dphi;
          for( q=0 ; q<NUM_Q ; ++q ){
            theCells[k][i][j].prim[q] /= dA;
          }
          if( theCells[k][i][j].prim[URR] > 0.0 ) theCells[k][i][j].prim[URR] = 0.0;
        }
      }
    }
  }

  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( n=n1 ; n<Nf ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_R_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=n1 ; n<Nf ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(theFaces,n);
      }
    }
    r = sim_FacePos(theSim,sim_N(theSim,R_DIR)-2,R_DIR);

    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double dz = zp-zm;
      for( j=0 ; j<sim_N_p(theSim,sim_N(theSim,R_DIR)-1) ; ++j ){
        double dA = dz*r*theCells[k][sim_N(theSim,R_DIR)-1][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][sim_N(theSim,R_DIR)-1][j].prim[q] /= dA;
        }
        theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] *= -1.;
      }
    }

  }

}

void cell_boundary_outflow_z( struct Cell *** theCells , struct Face * theFaces, struct Sim * theSim,struct MPIsetup * theMPIsetup,struct TimeStep * theTimeStep ){
  int NUM_Q = sim_NUM_Q(theSim);

  int j,i;
  int Nf = timestep_nzk(theTimeStep,sim_N(theSim,Z_DIR)-1);
  int n0 = timestep_nzk(theTimeStep,1);
  int n1 = timestep_nzk(theTimeStep,sim_N(theSim,Z_DIR)-2);

  int n,q;

  if(mpisetup_check_zbot_bndry(theMPIsetup)){
    for( n=0 ; n<n0 ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_L_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=0 ; n<n0 ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cL->prim[q] += cR->prim[q]*face_dA(theFaces,n);
      }
    }

    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double dA = .5*(rp*rp-rm*rm)*theCells[0][i][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[0][i][j].prim[q] /= dA;
        }
        theCells[0][i][j].prim[UZZ] *= -1.0;
      }
    }
  }

  if(mpisetup_check_ztop_bndry(theMPIsetup)){
    for( n=n1 ; n<Nf ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_R_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=n1 ; n<Nf ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(theFaces,n);
      }
    }


    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double dA = .5*(rp*rp-rm*rm)*theCells[sim_N(theSim,Z_DIR)-1][i][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[sim_N(theSim,Z_DIR)-1][i][j].prim[q] /= dA;
        }
        theCells[sim_N(theSim,Z_DIR)-1][i][j].prim[UZZ] *= -1.0;
      }
    }
  }

}

void cell_boundary_fixed_r( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup, 
    void (*cell_single_init_ptr)(struct Cell ***,struct Sim *,int,int,int) ){
  int i,j,k;
  if(mpisetup_check_rin_bndry(theMPIsetup)){
    for( i=0 ; i<sim_Nghost_min(theSim,R_DIR) ; ++i ){
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*cell_single_init_ptr)(theCells,theSim,i,j,k);
        }
      }
    }
  }

  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( i=sim_N(theSim,R_DIR)-1 ; i>sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR)-1 ; --i ){
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*cell_single_init_ptr)(theCells,theSim,i,j,k);
        }
      }
    }
  }

}

void cell_boundary_fixed_z( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup,
    void (*cell_single_init_ptr)(struct Cell ***,struct Sim *,int,int,int)  ){
  int i,j,k;
  if(mpisetup_check_zbot_bndry(theMPIsetup)){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( k=0 ; k<sim_Nghost_min(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*cell_single_init_ptr)(theCells,theSim,i,j,k);
        }
      }
    }
  }

  if( mpisetup_check_ztop_bndry(theMPIsetup) ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( k=sim_N(theSim,Z_DIR)-1 ; k>sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR)-1 ; --k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*cell_single_init_ptr)(theCells,theSim,i,j,k);
        }
      }
    }
  }

}


