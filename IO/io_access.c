#define IO_PRIVATE_DEFS
#include <string.h>
#include <stdlib.h>
#include "../Headers/IO.h"
#include "../Headers/header.h"

double io_tcheck(struct IO * theIO){
  return(theIO->tcheck);
}
int io_nfile(struct IO * theIO){
  return(theIO->nfile);
}

