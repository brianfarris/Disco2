#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"

void cell_add_src_grdisc( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses , double dt )
{
    int i,j,k;
  
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k )
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
    
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i )
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);

            for( j=0 ; j<sim_N_p(theSim,i) ; ++j )
            {
                struct Cell * c = &(theCells[k][i][j]);
                double phi = c->tiph-.5*c->dphi;
                double dphi = c->dphi;

                double rho = c->prim[RHO];
                double Pp  = c->prim[PPP];
                double r   = .5*(rp+rm);
                double vr  = c->prim[URR];
                double vz  = c->prim[UZZ];
                double vp  = c->prim[UPP]*r;
                double dz = zp-zm;
                double dV = dphi*.5*(rp*rp-rm*rm)*dz;

                double z  = .5*(zp+zm);
            }
        }
    }
}

