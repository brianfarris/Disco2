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

void cell_single_init(struct Cell ***theCells, struct Grid *theGrid,int i,int j,int k){
  double DISK_MACH = 10.;
  double GAMMALAW = grid_GAMMALAW(theGrid);
  double rm =  grid_r_faces(theGrid,i-1);
  double rp =  grid_r_faces(theGrid,i);
  double r = .5*(rm+rp);
  double omega = 1./pow(r,1.5);
  double cs = 1.0/DISK_MACH*omega*r;
  double rho = 1.0;
  double Pp = cs*cs*rho/GAMMALAW;

  theCells[k][i][j].prim[RHO] = 1.0;
  theCells[k][i][j].prim[PPP] = Pp;
  theCells[k][i][j].prim[URR] = 0.0;
  theCells[k][i][j].prim[UPP] = omega;
  theCells[k][i][j].prim[UZZ] = 0.0;
  theCells[k][i][j].prim[BRR] = 0.0;
  theCells[k][i][j].prim[BPP] = 0.0;
  theCells[k][i][j].prim[BZZ] = 0.0;
  theCells[k][i][j].prim[PSI] = 0.0;
  theCells[k][i][j].wiph = omega*r;
  theCells[k][i][j].divB = 0.0;
  theCells[k][i][j].GradPsi[0] = 0.0;
  theCells[k][i][j].GradPsi[1] = 0.0;
  theCells[k][i][j].GradPsi[2] = 0.0;
}

void cell_init(struct Cell ***theCells,struct Grid *theGrid,struct MPIsetup * theMPIsetup) {
  int N_r_withghost = grid_N_r(theGrid)+grid_Nghost_rmin(theGrid)+grid_Nghost_rmax(theGrid);
  int N_z_withghost = grid_N_z(theGrid)+grid_Nghost_zmin(theGrid)+grid_Nghost_zmax(theGrid);

  int i, j,k;

  double DISK_MACH = 10.;
  double GAMMALAW = grid_GAMMALAW(theGrid);
  
  srand(666+mpisetup_MyProc(theMPIsetup));
  double rho0 = 1.0;
  for (k = 0; k < N_z_withghost; k++) {
    for (i = 0; i < N_r_withghost; i++) {
      double rm = grid_r_faces(theGrid,i-1);
      double rp = grid_r_faces(theGrid,i);
      double r = 0.5*(rm+rp);
      double omega = 1./pow(r,1.5);
      double cs = 1.0/DISK_MACH*omega*r;
      double rho = rho0;
      double Pp = cs*cs*rho/GAMMALAW;

      double Bz;
      if ((r>2.)&&(r<3.)){
        double n=4.0;
        Bz = 0.05513/n;
      } else{
        Bz = 0.0;
      }
      double delta = .001*( (double)rand()/(double)RAND_MAX - .5 );

      for (j = 0; j < grid_N_p(theGrid,i); j++) {
        theCells[k][i][j].prim[RHO] = rho;
        theCells[k][i][j].prim[PPP] = Pp;
        theCells[k][i][j].prim[URR] = 0.0;
        theCells[k][i][j].prim[UPP] = omega*(1.+delta);
        theCells[k][i][j].prim[UZZ] = delta;
        theCells[k][i][j].prim[BRR] = 0.0;
        theCells[k][i][j].prim[BPP] = 0.0;
        theCells[k][i][j].prim[BZZ] = Bz;
        theCells[k][i][j].prim[PSI] = 0.0;
        theCells[k][i][j].wiph = omega*r;
        theCells[k][i][j].divB = 0.0;
        theCells[k][i][j].GradPsi[0] = 0.0;
        theCells[k][i][j].GradPsi[1] = 0.0;
        theCells[k][i][j].GradPsi[2] = 0.0;
      }
    }
  }

}
