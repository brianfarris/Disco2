#define RIEMANN_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Cell.h"
#include "../Headers/EOS.h"
#include "../Headers/Face.h"
#include "../Headers/Metric.h"
#include "../Headers/header.h"


// ********************************************************************************************
// WE REALLY SHOULD IMPROVE THE COMMENTING OF ALL OF THESE ROUTINES. 
// THERE IS SOME COMPLICATED STUFF HERE. 
// LETS CHOOSE A REFERENCE SUCH AS TORO AND IDENTIFY LINES OF CODE WITH EQUATIONS IN THE BOOK.
// ********************************************************************************************

// this routine is only called by riemann_set_vel.
// It is used to find various L/R quantities. 
// TODO: WRITE THIS.
void LR_speed_grdisc(double *prim, double r, int *n ,double GAMMALAW, double *p_vn, double *p_cf2, double *Fm, double *p_mn)
{

}

// Find velocities needed for the Riemann problem
// TODO: WRITE THIS.
void riemann_set_vel_grdisc(struct Riemann *theRiemann, struct Sim *theSim, double r, double GAMMALAW)
{

}

//TODO: WRITE THIS.
void riemann_set_flux_grdisc(struct Riemann *theRiemann, struct Sim *theSim, double GAMMALAW, int SetState)
{
}

void riemann_visc_flux_grdisc(struct Riemann *theRiemann, struct Sim *theSim)
{
}
