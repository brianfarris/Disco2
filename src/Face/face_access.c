#define FACE_PRIVATE_DEFS
#include <stdlib.h>
#include <stdIO.h>
#include <math.h>
#include "../Headers/Grid.h"
#include "../Headers/Face.h"
#include "../Headers/Cell.h"
#include "../Headers/header.h"


struct Face *face_pointer( struct Face * theFaces,int n){
  return &( theFaces[n] );
}
struct Cell *face_L_pointer( struct Face * theFaces,int n){
  return(theFaces[n].L) ;
}
struct Cell *face_R_pointer( struct Face * theFaces,int n){
  return(theFaces[n].R) ;
}
double face_deltaL( struct Face * thisface){
  return(thisface->deltaL);
}
double face_deltaR( struct Face * thisface){
  return(thisface->deltaR);
}
double face_cm( struct Face * thisface){
  return(thisface->cm);
}
double face_dA( struct Face * thisface){
  return(thisface->dA);
}



