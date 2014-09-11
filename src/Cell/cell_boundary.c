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
  int Nf = timestep_n(theTimeStep,sim_N(theSim,R_DIR)-1,R_DIR);
  int NUM_Q = sim_NUM_Q(theSim);
  int n1 = timestep_n(theTimeStep,sim_N(theSim,R_DIR)-2,R_DIR);
  int n2 = timestep_n(theTimeStep,sim_N(theSim,R_DIR)-3,R_DIR);

  int n,q;
  int i,j,k;
  double r_face,r_face_m1,r_face_p1,r_cell,r_cellR;

  if( mpisetup_check_rin_bndry(theMPIsetup) ){
    if (sim_NoInnerBC(theSim)!=1){ // if the global inner radius is set negative, we don't apply an inner BC
      for( i=0 ; i>=0 ; --i ){
        r_face=sim_FacePos(theSim,i,R_DIR);
        r_face_m1=sim_FacePos(theSim,i-1,R_DIR);
        r_face_p1=sim_FacePos(theSim,i+1,R_DIR);
        r_cell = 0.5*(r_face+r_face_m1);
        r_cellR = 0.5*(r_face+r_face_p1);
        for( n=timestep_n(theTimeStep,i,R_DIR) ; n<timestep_n(theTimeStep,i+1,R_DIR) ; ++n ){
          for( q=0 ; q<NUM_Q ; ++q ){
            face_L_pointer(theFaces,n)->prim[q] = 0.0;
          }
        } 
        for( n=timestep_n(theTimeStep,i,R_DIR) ; n<timestep_n(theTimeStep,i+1,R_DIR) ; ++n ){
          struct Cell * cL = face_L_pointer(theFaces,n);
          struct Cell * cR = face_R_pointer(theFaces,n);
          if (SS_BCS==1){
            for( q=4 ; q<NUM_Q ; ++q ){
              cL->prim[q] += (cR->prim[q]-EXTRAP_BC*cR->grad[q]*(sim_FacePos(theSim,i+1,R_DIR)-sim_FacePos(theSim,i-1,R_DIR))/2.)*face_dA(theFaces,n);
            }
            cL->prim[RHO] += (cR->prim[RHO]/pow(r_cellR,-3./5.)*pow(r_cell,-3./5.))*face_dA(theFaces,n);
            cL->prim[PPP] += (cR->prim[PPP]/pow(r_cellR,-3./2.)*pow(r_cell,-3./2.))*face_dA(theFaces,n);
            cL->prim[URR] += (cR->prim[URR]-EXTRAP_BC*cR->grad[URR]*(sim_FacePos(theSim,i+1,R_DIR)-sim_FacePos(theSim,i-1,R_DIR))/2.)*face_dA(theFaces,n);
            //cL->prim[UPP] += (cR->prim[UPP]-EXTRAP_BC*cR->grad[UPP]*(sim_FacePos(theSim,i+1,R_DIR)-sim_FacePos(theSim,i-1,R_DIR))/2.)*face_dA(theFaces,n);
            cL->prim[UPP] += (cR->prim[UPP]/pow(r_cellR,-7./5.)*pow(r_cell,-7./5.))*face_dA(theFaces,n);
            //cL->prim[UPP] += ((cR->prim[UPP]-pow(r_cellR,-1.5)*exp(-pow(r_cellR/0.5,1.5)))/pow(r_cellR,-7./5.)*pow(r_cell,-7./5.) + pow(r_cell,-1.5)*exp(-pow(r_cell/0.5,1.5)))*face_dA(theFaces,n);
          } else{
            for( q=0 ; q<NUM_Q ; ++q ){
              cL->prim[q] += (cR->prim[q]-EXTRAP_BC*cR->grad[q]*(sim_FacePos(theSim,i+1,R_DIR)-sim_FacePos(theSim,i-1,R_DIR))/2.)*face_dA(theFaces,n);
            }
          }
        }
        for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
          double zp = sim_FacePos(theSim,k,Z_DIR);
          double zm = sim_FacePos(theSim,k-1,Z_DIR);
          double dz = zp-zm;
          for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
            double dA = dz*r_face*theCells[k][i][j].dphi;
            for( q=0 ; q<NUM_Q ; ++q ){
              theCells[k][i][j].prim[q] /= dA;
            }
            if (fabs(r_face-0.6)<0.00001){
              //printf("r_face: %e, wiph left: %e, theCells[k][i][j].prim[URR]: %e, theCells[k][i+1][j].prim[URR]: %e, theCells[k][i+1][j].grad[URR]: %e\n",r_face,theCells[k][i][j].wiph,theCells[k][i][j].prim[URR],theCells[k][i+1][j].prim[URR],theCells[k][i+1][j].grad[URR]);
              //exit(1);
            }
            if( theCells[k][i][j].prim[URR] > 0.0 && diode==1) theCells[k][i][j].prim[URR] = 0.0;
            if (KEP_BNDRY==1){
              theCells[k][i][j].prim[UPP] = 0.0;// pow(r_cell,-1.5);          
              //printf("hello rin. r_cell: %e\n",r_cell);
            }
          }
        }
      }
    }
    //exit(1);
  }


  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( n=n2 ; n<n1 ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_R_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=n2 ; n<n1 ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(theFaces,n);
      }
    }
    //r_face = sim_FacePos(theSim,sim_N(theSim,R_DIR)-2,R_DIR);
    //r_face_p1 = sim_FacePos(theSim,sim_N(theSim,R_DIR)-1,R_DIR);
    //r_cell = 0.5*(r_face+r_face_p1);
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double dz = zp-zm;
      for( j=0 ; j<sim_N_p(theSim,sim_N(theSim,R_DIR)-2) ; ++j ){
        r_face = sim_FacePos(theSim,sim_N(theSim,R_DIR)-3,R_DIR);
        r_face_p1 = sim_FacePos(theSim,sim_N(theSim,R_DIR)-2,R_DIR);
        r_cell = 0.5*(r_face+r_face_p1); 
        double dA = dz*r_face*theCells[k][sim_N(theSim,R_DIR)-2][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][sim_N(theSim,R_DIR)-2][j].prim[q] /= dA;
        }
        if( theCells[k][sim_N(theSim,R_DIR)-2][j].prim[URR] < 0.0 && diode==1 ) theCells[k][sim_N(theSim,R_DIR)-2][j].prim[URR] = 0.0;
        if (KEP_BNDRY==1){
          theCells[k][sim_N(theSim,R_DIR)-2][j].prim[UPP] = 0.0;// pow(r_cell,-1.5);
        }
      }
    }
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


    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double dz = zp-zm;
      for( j=0 ; j<sim_N_p(theSim,sim_N(theSim,R_DIR)-1) ; ++j ){
        r_face = sim_FacePos(theSim,sim_N(theSim,R_DIR)-2,R_DIR);
        r_face_p1 = sim_FacePos(theSim,sim_N(theSim,R_DIR)-1,R_DIR);
        r_cell = 0.5*(r_face+r_face_p1); 
        double dA = dz*r_face*theCells[k][sim_N(theSim,R_DIR)-1][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][sim_N(theSim,R_DIR)-1][j].prim[q] /= dA;
        }
        //theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] *= -1.; 
        if( theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] < 0.0 && diode==1 ) theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] = 0.0;
        if (KEP_BNDRY==1){
          theCells[k][sim_N(theSim,R_DIR)-1][j].prim[UPP] = 0.0;// pow(r_cell,-1.5);
        }
      }

    }
  }
}

void cell_boundary_outflow_z( struct Cell *** theCells , struct Face * theFaces, struct Sim * theSim,struct MPIsetup * theMPIsetup,struct TimeStep * theTimeStep ){
  int NUM_Q = sim_NUM_Q(theSim);

  int j,i;
  int Nf = timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-1,Z_DIR);
  int n0 = timestep_n(theTimeStep,1,Z_DIR);
  int n1 = timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-2,Z_DIR);

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
    void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int) ){
  int i,j,k;
  if (sim_NoInnerBC(theSim)!=1){ 
    if(mpisetup_check_rin_bndry(theMPIsetup)){
      for( i=0 ; i<sim_Nghost_min(theSim,R_DIR) ; ++i ){
        for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
          for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
            (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
            //printf("theCells[0][%d][%d].prim[PPP]: %e\n",i,j,theCells[0][i][j].prim[PPP]);
          }
        }
      }
    }
  }
  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( i=sim_N(theSim,R_DIR)-1 ; i>sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR)-1 ; --i ){
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
        }
      }
    }
  }

}

void cell_boundary_fixed_z( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup,
    void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int)  ){
  int i,j,k;
  if(mpisetup_check_zbot_bndry(theMPIsetup)){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( k=0 ; k<sim_Nghost_min(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
        }
      }
    }
  }

  if( mpisetup_check_ztop_bndry(theMPIsetup) ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( k=sim_N(theSim,Z_DIR)-1 ; k>sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR)-1 ; --k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
        }
      }
    }
  }

}


