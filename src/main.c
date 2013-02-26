//                         ...WWn+,....
//                     ..JWyyZyZZXyyyZWL
//                   .XyZyZZZyyXyyZZyZyWWl
//                 .dWyZyZZyyZZyZZZyZyZyyWh.                       .J7!J
//                .,WyyZZyyZyZyZyyZyyZZyZyX\                    .JC..J=
//               .HZZZyZZyZyZyZZyZZyZyyZZZyW.                  .=?^^J&..
//               jZyyZyZyZyZZyyZyyyZyZyZyyZW`                ..:.JZWC.J\
//              .WyZyZyZyZyXkXXWWVWWVWWQkyZyF              .CJJVn.`3J=^J\
//              .HZyZyZyZyyW!`^^^^^^^^^?HZyX'            .CJ!    1,^?X!
//              KyyZyyZyZZX%.JgJJ.^^^.J+dkX%           .C`J:      j2!
//              ?WyZyXyZyyW+XHWkpHJJ+jHWHHF          .V76.r        l
//               ,kZyyZyWYY+WHHHWK`^^^?HkV!       .,Cv  ,=         j
//                .7HyZWIr^^^?77!^^^^^^J%       .J:J^              ,
//                  TWkXkJJ:^^^.JJ:`+1z6      .?.C7               .;
//                    .4WWD^^^J1.^^^^^^^?l .,!                ..?`
//                      .Z$^^^!^^??77Tf7!.?`               .?!
//                      ydt^^^^^^^^^J..J!              .,?`
//                     Jwdr^^^^^^^^.S?n.            .,^
//                  .,v!w1S.^^^^^^^`R?1Wl,       .,^
//               .,^  2 0?J2^^^^^^^^?k??4e..  .,^
//             .?7?i, ?J3??vo^^^^...^jAy&X.h.?`
//          .,^     JC?K??1udh.^J$j$^J#| ?tj
//        .?...      .qCuv= jvk`!!J,^dMr   ?
//      .Oz1?`        dP`    4de^`T:J#d$   j.
//    .?,:J:          .Y,    ,ZMN,^.M3dw   2l
//  .J.J:^J`         ,^ `,   .SWJWNdb,dw  .\.,
// .6CJ,J?        .J\    .,   0vhH1??dkd. ,  l
// v I?!?+4`   .?! .:     .+  w?JNm&z60d:.\  ?
// u.C^^^J\   ``   .`      .o j=?vh,^JId;J   .!
// j^`S:J'    .|   .      d.`4jz???ZGD?zyCJ$  |
// ,.^^J\      ?.  J      `+^^Jz???????zP^:r  1
// .cJ,?.C!     t  r       `+:JI???????zP^J`  ,
//  $?:^^J`     .,.:        1:JI???????zPJ:.  .,
//  j.^^J\       j,`   `  J7!^JI???????1$jJ4.  |
//  .xA+z         .`      .l^^JI???????1P^^j   l
//   6.^?i   ..JZ!2...     ,x^JI???????1$^.r  Jo
//    .JC??!`^^^I.$J`    .c I^JI????jx?1$^J, .CI
//    $^^^^^^^^.I2^?...,=?C,C^Jk???j=jxz$^^?7:^z
//    ?&.^.J`J.JJ:^^^^^^^^^^:^J;7wv'  ?7$^^^^JJ^
//       .777=! `1+J.:^:^^:^:.J!     .  7777?t
//                   `?P777?`       .^       r

#include <stdio.h>
#include <stdlib.h>
#include "Headers/MPIsetup.h"
#include "Headers/Cell.h"
#include "Headers/Grid.h"
#include "Headers/Face.h"
#include "Headers/GravMass.h"
#include "Headers/IO.h"
#include "Headers/TimeStep.h"
#include "Headers/header.h"

int main(int argc, char **argv) {

  // start MPI 
  struct MPIsetup * theMPIsetup = mpisetup_create(argc,argv);
  mpisetup_setprocs(theMPIsetup);
  mpisetup_cart_create(theMPIsetup);
  mpisetup_left_right(theMPIsetup);

  //grid
  struct Grid * theGrid = grid_create(theMPIsetup);
  grid_set_N_p(theGrid);
  grid_set_rz(theGrid,theMPIsetup);
  grid_set_Ncells_and_offset(theGrid,theMPIsetup);

  // gravMass
  struct GravMass *theGravMasses = gravMass_create(grid_NumGravMass(theGrid));
  gravMass_initialize(theGravMasses);
  gravMass_clean_pi(theGravMasses,theGrid);

  // allocate memory for data 
  struct Cell ***theCells = cell_create(theGrid,theMPIsetup);
  cell_clean_pi(theCells,theGrid);
  // set initial data 
  cell_init(theCells,theGrid,theMPIsetup);

  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid,theMPIsetup);
  cell_syncproc_z(theCells,theGrid,theMPIsetup);

  //set conserved quantities
  cell_calc_cons(theCells,theGrid);

  struct TimeStep * theTimeStep = timestep_create();

  double dtcheck = timestep_get_T_MAX(theTimeStep)/timestep_NUM_CHECKPOINTS(theTimeStep);
  double tcheck = dtcheck;
  printf("dtcheck: %f\n",dtcheck);
  int nfile=0;
  char filename[256];

  while( timestep_get_t(theTimeStep) < timestep_get_T_MAX(theTimeStep) ){
    timestep_set_dt(theTimeStep,theCells,theGrid);
    cell_copy(theCells,theGrid);
    gravMass_copy(theGravMasses,theGrid);
    timestep_set_RK(theTimeStep,0.0);
    timestep_substep(theTimeStep,theCells,theGrid,theGravMasses,theMPIsetup,1.0);
    timestep_set_RK(theTimeStep,0.5);
    timestep_substep(theTimeStep,theCells,theGrid,theGravMasses,theMPIsetup,0.5);
    timestep_update_Psi(theTimeStep,theCells,theGrid,theMPIsetup);
    timestep_update_t(theTimeStep); 
    if( timestep_get_t(theTimeStep)>tcheck){
      sprintf(filename,"checkpoint_%04d.h5",nfile);
      struct IO *theIO = io_create(theGrid);
      io_flattened_prim(theIO,theCells,theGrid);
      io_hdf5_out(theIO,theGrid,filename);
      io_destroy(theIO,theGrid);
      tcheck += dtcheck;
      ++nfile;
    }
  }
  /*
     cell_printscreen(theCells,theGrid);
     struct io *theIO = io_new(theGrid);
     io_flattened_prim(theIO,theCells,theGrid);

     io_hdf5_out(theIO,theGrid);
     io_delete(theIO,theGrid);

     struct io *input_prims = io_new(theGrid);
     io_hdf5_in(input_prims,theGrid);
     io_unflattened_prim(theIO,theCells,theGrid);
     io_delete(input_prims,theGrid);
     */
  //inter-processor syncs
  cell_syncproc_r(theCells,theGrid,theMPIsetup);
  cell_syncproc_z(theCells,theGrid,theMPIsetup);

  // clean up
  cell_destroy(theCells,theGrid);
  grid_destroy(theGrid);
  gravMass_destroy(theGravMasses);

  // exit MPI 
  mpisetup_destroy(theMPIsetup);
  return(0);
}
