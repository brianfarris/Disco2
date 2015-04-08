#define CELL_PRIVATE_DEFS
#define EOS_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

//Input Initial condition from file

void cell_init_file_radtxt(struct Cell *c, double r, double phi, double z,
                                struct Sim *theSim)
{
    //Read input from text file "input_r.txt".  Each row (after the first) of
    //the file is to have whitespace separated columns in the order:
    // r RHO PPP URR UPP (UZZ QQ1 QQ2...)
    // The first row contains two integers: the number of following rows and
    // the number of variables (must be at least 4).  It is assumed values of 
    // r are unique and increasing.

    int i, atleast2;
    double r1, *dat1, r2, *dat2, w1, w2;
    int nr, nq;
    double rho, P, vr, vp, vz, *q;

    FILE *f = fopen("input_r.txt", "r");
    fscanf(f, "%d %d\n", &nr, &nq);

    dat1 = (double *)malloc(nq * sizeof(double));
    dat2 = (double *)malloc(nq * sizeof(double));
    q = nq>5 ? (double*)malloc((nq-5)*sizeof(double)) : NULL;

    r2 = -1;
    atleast2 = 0;
    while(r2 < r || atleast2 < 2)
    {
        r1 = r2;
        for(i=0; i<nq; i++)
            dat1[i] = dat2[i];
        fscanf(f, "%lg ", &r2);
        for(i=0; i<nq; i++)
            fscanf(f, "%lg ", &dat2[i]);
        atleast2++;
    }
    fclose(f);

    w1 = (r2-r)/(r2-r1);
    w2 = (r-r1)/(r2-r1);

    rho = w1*dat1[0] + w2*dat2[0];
    P = w1*dat1[1] + w2*dat2[1];
    vr = w1*dat1[2] + w2*dat2[2];
    vp = w1*dat1[3] + w2*dat2[3];
    
    vz = nq>4 ? w1*dat1[4]+w2*dat2[4] : 0.0;

    if(nq > 5)
        for(i=0; i<nq-5; i++)
            q[i] = w1*dat1[i+5] + w2*dat2[i+5];


    free(dat1);
    free(dat2);
    if(q != NULL)
        free(q);

    c->prim[RHO] = rho;
    c->prim[PPP] = P;
    c->prim[URR] = vr;
    c->prim[UPP] = vp;
    c->prim[UZZ] = vz;

    for(i=sim_NUM_C(theSim); i<sim_NUM_Q(theSim); i++)
        c->prim[i] = q[i-5];
}

void cell_init_file_calc(struct Cell *c, double r, double phi, double z,
                            struct Sim *theSim)
{
    int opt = sim_InitPar0(theSim);

    if(opt == 0)
        cell_init_file_radtxt(c, r, phi, z, theSim);
    else
        printf("ERROR: cell_init_file given bad option.\n");
}

void cell_single_init_file(struct Cell *theCell, struct Sim *theSim,
                            int i, int j, int k)
{
    double rm = sim_FacePos(theSim,i-1,R_DIR);
    double rp = sim_FacePos(theSim,i,R_DIR);
    double r = 0.5*(rm+rp);
    double zm = sim_FacePos(theSim,k-1,Z_DIR);
    double zp = sim_FacePos(theSim,k,Z_DIR);
    double z = 0.5*(zm+zp);
    double t = theCell->tiph-.5*theCell->dphi;

    cell_init_file_calc(theCell, r, t, z, theSim);
}

void cell_init_file(struct Cell ***theCells,struct Sim *theSim,
                    struct MPIsetup *theMPIsetup)
{
    int i, j, k;
    for (k = 0; k < sim_N(theSim,Z_DIR); k++)
    {
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double z = 0.5*(zm+zp);

        for (i = 0; i < sim_N(theSim,R_DIR); i++)
        {
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double r = 0.5*(rm+rp);

            for (j = 0; j < sim_N_p(theSim,i); j++)
            {
                double t = theCells[k][i][j].tiph-.5*theCells[k][i][j].dphi;

                cell_init_file_calc(&(theCells[k][i][j]), r, t, z, theSim);

                if(PRINTTOOMUCH)
                {
                    printf("(%d,%d,%d) = (%.12lg, %.12lg, %.12lg): (%.12lg, %.12lg, %.12lg, %.12lg, %.12lg)\n",
                          i,j,k,r,t,z,theCells[k][i][j].prim[RHO], theCells[k][i][j].prim[URR],
                          theCells[k][i][j].prim[UPP], theCells[k][i][j].prim[UZZ],
                          theCells[k][i][j].prim[PPP]);
                }
            }
        }
    }
}
