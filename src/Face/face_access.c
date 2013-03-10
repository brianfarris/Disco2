#define FACE_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/Cell.h"
#include "../Headers/header.h"

struct Cell *face_L_pointer( struct Face * theFaces,int n){
  return(theFaces[n].L) ;
}
struct Cell *face_R_pointer( struct Face * theFaces,int n){
  return(theFaces[n].R) ;
}
double face_deltaL( struct Face * theFaces,int n){
  return(theFaces[n].deltaL);
}
double face_deltaR( struct Face * theFaces,int n){
  return(theFaces[n].deltaR);
}
double face_cm( struct Face * theFaces,int n){
  return(theFaces[n].cm);
}
double face_dA( struct Face * theFaces,int n){
  return(theFaces[n].dA);
}
double face_r( struct Face * theFaces,int n){
  return(theFaces[n].r);
}



