#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"


void cell_boundary_outflow_r( struct Cell *** theCells , struct Face * theFaces ,struct Grid * theGrid,struct MPIsetup * theMPIsetup, int * nri ){
  int Nf = nri[grid_N_r(theGrid)-1];
  int NUM_Q = grid_NUM_Q(theGrid);
  int n1 = nri[grid_N_r(theGrid)-2];

  int n,q;
  int i,j,k;
  double r;

  if( mpisetup_check_rin_bndry(theMPIsetup) ){

    for( i=0 ; i>=0 ; --i ){
      r=grid_r_faces(theGrid,i);
      for( n=nri[i] ; n<nri[i+1] ; ++n ){
        for( q=0 ; q<NUM_Q ; ++q ){
          face_L_pointer(theFaces,n)->prim[q] = 0.0;
        }
      } 
      for( n=nri[i] ; n<nri[i+1] ; ++n ){
        struct Face * f = face_pointer(theFaces,n);//&(theFaces[n]);
        struct Cell * cL = face_L_pointer(theFaces,n);
        struct Cell * cR = face_R_pointer(theFaces,n);
        for( q=0 ; q<NUM_Q ; ++q ){
          cL->prim[q] += cR->prim[q]*face_dA(f);
        }
      }

      for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
        double zp = grid_z_faces(theGrid,k);
        double zm = grid_z_faces(theGrid,k-1);
        double dz = zp-zm;
        for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
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
      struct Face * f = face_pointer(theFaces,n);
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(f);
      }
    }
    r = grid_r_faces(theGrid,grid_N_r(theGrid)-2);

    for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
      double zp = grid_z_faces(theGrid,k);
      double zm = grid_z_faces(theGrid,k-1);
      double dz = zp-zm;
      for( j=0 ; j<grid_N_p(theGrid,grid_N_r(theGrid)-1) ; ++j ){
        double dA = dz*r*theCells[k][grid_N_r(theGrid)-1][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][grid_N_r(theGrid)-1][j].prim[q] /= dA;
        }
        theCells[k][grid_N_r(theGrid)-1][j].prim[URR] *= -1.;
      }
    }

  }

}

void cell_boundary_outflow_z( struct Cell *** theCells , struct Face * theFaces, struct Grid * theGrid,struct MPIsetup * theMPIsetup,int * nzk ){
  int NUM_Q = grid_NUM_Q(theGrid);

  int j,i;
  int Nf = nzk[grid_N_z(theGrid)-1];
  int n0 = nzk[1];
  int n1 = nzk[grid_N_z(theGrid)-2];

  int n,q;

  if(mpisetup_check_zbot_bndry(theMPIsetup)){
    for( n=0 ; n<n0 ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_L_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=0 ; n<n0 ; ++n ){
      struct Face * f = face_pointer(theFaces,n);
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cL->prim[q] += cR->prim[q]*face_dA(f);
      }
    }

    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
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
      struct Face * f = face_pointer(theFaces,n);
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(f);
      }
    }


    for( i=0 ; i<grid_N_r(theGrid) ; ++i ){
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        double dA = .5*(rp*rp-rm*rm)*theCells[grid_N_z(theGrid)-1][i][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[grid_N_z(theGrid)-1][i][j].prim[q] /= dA;
        }
        theCells[grid_N_z(theGrid)-1][i][j].prim[UZZ] *= -1.0;
      }
    }
  }

}

void cell_boundary_fixed_r( struct Cell *** theCells, struct Grid *theGrid,struct MPIsetup * theMPIsetup ){
  int i,j,k;

  if(mpisetup_check_rin_bndry(theMPIsetup)){
    for( i=0 ; i<grid_Nghost_rmin(theGrid) ; ++i ){
      for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
        for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
          cell_single_init_shear(theCells,theGrid,i,j,k);
        }
      }
    }
  }

  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( i=grid_N_r(theGrid)-1 ; i>grid_N_r(theGrid)-grid_Nghost_rmax(theGrid)-1 ; --i ){
      for( k=0 ; k<grid_N_z(theGrid) ; ++k ){
        for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
          cell_single_init_shear(theCells,theGrid,i,j,k);
        }
      }
    }
  }

}


