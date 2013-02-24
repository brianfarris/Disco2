#ifndef ONESTEP_H
#define ONESTEP_H
struct cell;
struct grid;
struct gravMass;

void onestep(struct Cell ***,struct Grid *,struct GravMass * ,double,double);
void onestep_Psi( struct Cell *** , struct Grid *, double );
#endif
