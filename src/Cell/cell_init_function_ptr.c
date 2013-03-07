#include <stdio.h>
#include <stdlib.h>
#include "../Headers/Cell.h"
#include "../Headers/Grid.h"
#include "../Headers/header.h"
void (*cell_init_ptr(struct Grid * theGrid))(struct Cell *** , struct Grid * ){
  if (grid_InitialDataType(theGrid)==FLOCK){
    return(&cell_init_flock);
  } else if (grid_InitialDataType(theGrid)==SHEAR){
    return(&cell_init_shear);
  } else{
    printf("ERROR\n");
    exit(0);
  }
}

void (*cell_single_init_ptr(struct Grid * theGrid))(struct Cell *** , struct Grid *,int,int,int ){
  if (grid_InitialDataType(theGrid)==FLOCK){
    return(&cell_single_init_flock);
  } else if (grid_InitialDataType(theGrid)==SHEAR){
    return(&cell_single_init_shear);
  } else{
    printf("ERROR\n");
    exit(0);
  }
}


