#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdIO.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"


void cell_boundary_outflow_r( struct Cell *** theCells , struct Face * theFaces ,struct Grid * theGrid, int * nri ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);
  int Nf = nri[N_r_withghost-1];
  //   int n0 = nri[1];
  int n1 = nri[N_r_withghost-2];

  int n,q;
  int i,j,k;
  double r;

  if( dim_MyProc[0] == 0 ){

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

      for( k=0 ; k<N_z_withghost ; ++k ){
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

  if( dim_MyProc[0] == dim_NumProcs[0]-1 ){
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
    r = grid_r_faces(theGrid,N_r_withghost-2);

    for( k=0 ; k<N_z_withghost ; ++k ){
      double zp = grid_z_faces(theGrid,k);
      double zm = grid_z_faces(theGrid,k-1);
      double dz = zp-zm;
      for( j=0 ; j<grid_N_p(theGrid,N_r_withghost-1) ; ++j ){
        double dA = dz*r*theCells[k][N_r_withghost-1][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][N_r_withghost-1][j].prim[q] /= dA;
        }
        theCells[k][N_r_withghost-1][j].prim[URR] *= -1.;
      }
    }

  }

}

void cell_boundary_outflow_z( struct Cell *** theCells , struct Face * theFaces, struct Grid * theGrid,int * nzk ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int j,i;
  int Nf = nzk[N_z_withghost-1];
  int n0 = nzk[1];
  int n1 = nzk[N_z_withghost-2];

  int n,q;

  if( dim_MyProc[1] == 0 ){
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

    for( i=0 ; i<N_r_withghost ; ++i ){
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

  if( dim_MyProc[1] == dim_NumProcs[1]-1 ){
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


    for( i=0 ; i<N_r_withghost ; ++i ){
      double rp = grid_r_faces(theGrid,i);
      double rm = grid_r_faces(theGrid,i-1);
      for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
        double dA = .5*(rp*rp-rm*rm)*theCells[N_z_withghost-1][i][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[N_z_withghost-1][i][j].prim[q] /= dA;
        }
        theCells[N_z_withghost-1][i][j].prim[UZZ] *= -1.0;
      }
    }
  }

}

void cell_boundary_fixed_r( struct Cell *** theCells, struct Grid *theGrid ){
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i,j,k;

  if( dim_MyProc[0] == 0 ){
    for( i=0 ; i<grid_Nghost_rmin(theGrid) ; ++i ){
      for( k=0 ; k<N_z_withghost ; ++k ){
        for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
          cell_single_init(theCells,theGrid,i,j,k);
        }
      }
    }
  }

  if( dim_MyProc[0] == dim_NumProcs[0]-1 ){
    for( i=N_r_withghost-1 ; i>N_r_withghost-grid_Nghost_rmax(theGrid)-1 ; --i ){
      for( k=0 ; k<N_z_withghost ; ++k ){
        for( j=0 ; j<grid_N_p(theGrid,i) ; ++j ){
          cell_single_init(theCells,theGrid,i,j,k);
        }
      }
    }
  }

}


