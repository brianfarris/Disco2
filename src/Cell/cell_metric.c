#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/Metric.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"


void cell_metric_horizon_jazz(struct Cell ***theCells, struct Sim *theSim)
{
    int NUM_Q = sim_NUM_Q(theSim);

    int nh,q;
    int i,j,k;

    nh = sim_N(theSim, R_DIR);

    for(i=0; i<sim_N(theSim, R_DIR); i++)
    {
        double rp = sim_FacePos(theSim, i, R_DIR);
        if(rp > metric_horizon(theSim))
        {
            nh = i-2;
            break;
        }
    }

    for(k=0; k<sim_N(theSim, Z_DIR); k++)
    {
        double zp = sim_FacePos(theSim, k, Z_DIR);
        double zm = sim_FacePos(theSim, k-1, Z_DIR);
        double z = 0.5*(zm+zp);

        for(i=0; i<nh; i++)
        {
            double rp = sim_FacePos(theSim, i, R_DIR);
            double rm = sim_FacePos(theSim, i-1, R_DIR);
            double r = 0.5*(rm+rp);
            
            int nans = 0;
            int attempt = 0;
            
            do
            {
                nans = 0;
                for(j=0; j<sim_N_p(theSim, i); j++)
                {
                    double dphi = theCells[k][i][j].dphi;
                    double phi = theCells[k][i][j].tiph - 0.5*dphi;

                    for(q=0; q<NUM_Q; q++)
                    {
                        //Check for NaNs
                        if(theCells[k][i][j].prim[q]
                                != theCells[k][i][j].prim[q])
                        {
                            nans++;
                            printf("Got a NaN here! (q=%d, r=%.12lg)\n", q, r);
                            (*cell_single_init_ptr(theSim))(&(theCells[k][i][j]), theSim, i,j,k);
                            /*
                            int jL = j>0 ? j-1 : sim_N_p(theSim,i)-1;
                            int jR = j<sim_N_p(theSim,i)-1 ? j+1 : 0;
                            double qL = theCells[k][i][jL].prim[q];
                            double qR = theCells[k][i][jR].prim[q];

                            if(qL == qL && qR == qR)
                            {
                                printf("   Fixed with phi-average.\n");
                                theCells[k][i][j].prim[q] = 0.5*(qL+qR);
                            }
                            else if(qL == qL)
                            {
                                printf("   Fixed with left.\n");
                                theCells[k][i][j].prim[q] = qL;
                            }
                            else if(qR == qR)
                            {
                                printf("   Fixed with right.\n");
                                theCells[k][i][j].prim[q] = qR;
                            }
                            else
                                printf("   No Fix, hopefully next time.\n");
                            */
                        }
                    }
                }
                attempt++;
            }
            while(nans > 0 && attempt < sim_N_p(theSim, i));
        }
    }
}
