#define FACE_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/Cell.h"
#include "../Headers/header.h"


//struct Face *face_pointer( struct Face * theFaces,int n){
//  return &( theFaces[n] );
//}
struct Cell *face_L_pointer( struct Face * theFaces,int n){
  return(theFaces[n].L) ;
}
struct Cell *face_R_pointer( struct Face * theFaces,int n){
  return(theFaces[n].R) ;
}
double face_deltaL( struct Face * thisface, int n){
  return(thisface[n].deltaL);
}
double face_deltaR( struct Face * thisface, int n){
  return(thisface[n].deltaR);
}
double face_cm( struct Face * thisface, int n){
  return(thisface[n].cm);
}
double face_dA( struct Face * thisface, int n){
  return(thisface[n].dA);
}
double face_r( struct Face * thisface, int n){
  return(thisface[n].r);
}



