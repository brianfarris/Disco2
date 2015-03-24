/*
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
*/

#include <stdio.h>
#include <stdlib.h>
#include "Headers/MPIsetup.h"
#include "Headers/Cell.h"
#include "Headers/Sim.h"
#include "Headers/EOS.h"
#include "Headers/Face.h"
#include "Headers/GravMass.h"
#include "Headers/IO.h"
#include "Headers/TimeStep.h"
#include "Headers/Diagnostics.h"
#include "Headers/Metric.h"
#include "Headers/header.h"


int main(int argc, char **argv) {
  char* inputfilename = argv[1];

  printf("\nWelcome to Disco!\n\nInitializing MPI.\n");

  // start MPI 
  struct MPIsetup * theMPIsetup = mpisetup_create(argc,argv);
  mpisetup_setprocs(theMPIsetup,inputfilename);
  mpisetup_cart_create(theMPIsetup);
  mpisetup_left_right(theMPIsetup);

  printf("Setting user parameters.\n");

  // Set up everything that will not change throughout the simulation
  struct Sim * theSim = sim_create(theMPIsetup);
  sim_read_par_file(theSim,theMPIsetup,inputfilename);
  sim_alloc_arr(theSim,theMPIsetup); 
  sim_set_rz(theSim,theMPIsetup);
  sim_set_N_p(theSim);
  sim_set_misc(theSim,theMPIsetup);

  printf("Building simulation structures.\n");
  //allocate memory for gravitating masses
  struct GravMass *theGravMasses = gravMass_create(sim_NumGravMass(theSim));

  // allocate memory for data 
  struct Cell ***theCells = cell_create(theSim,theMPIsetup);

  // make sure that phi is between 0 and 2 pi. I'm not sure why this matters.
  cell_clean_pi(theCells,theSim);

  printf("Performing initial sync.\n");
  // inter-processor syncs. cell_create currently sets up tiph_0 randomly, 
  // so we need to make sure ghost-zones are set appropriately
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  cell_syncproc_z(theCells,theSim,theMPIsetup);

  printf("Building more structures.\n");
  // create TimeStep struct
  struct TimeStep * theTimeStep = timestep_create(theSim);

  // create IO struct
  struct IO *theIO = io_create(theSim);

  //Point Cell & Riemann to Newtonian/GR versions
  metric_init_background(theSim);
  if(sim_Background(theSim) != NEWTON)
    metric_init_metric(theSim);
  eos_init(theSim);

  printf("Setting initial conditions.\n");
  // set initial data 
  if (sim_Restart(theSim)==1){ // getting initial data from checkpoint file
    io_allocbuf(theIO,theSim); // allocate memory for a buffer to store checkpoint data
    io_hdf5_in(theIO,theSim,theTimeStep); // read from hdf5 file into buffer
    io_readbuf(theIO,theCells,theSim,theGravMasses); // read from buffer into theCells and theGravMasses
    io_deallocbuf(theIO); // get rid of buffer
  }else{
    (*gravMass_init_ptr(theSim))(theGravMasses); // set up gravitating masses. 
    (*cell_init_ptr(theSim))(theCells,theSim,theMPIsetup); // setup initial data using routine specified in .par file
  }
  gravMass_clean_pi(theGravMasses,theSim); // make sure GravMasses have phi between 0 and 2 pi.

  printf("Setting up I/O.\n");
  // set tcheck, dtcheck, and nfile. Needs to be done after ID because it depends on current time.
  io_setup(theIO,theSim,theTimeStep);

  printf("Last sync.\n");
  // inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  cell_syncproc_z(theCells,theSim,theMPIsetup);

  printf("Calculating initial conserved quantities.\n");
  // inter-processor syncs
  // set conserved quantities
  cell_calc_cons(theCells,theSim);
  
  printf("Setting up diagnostics.\n");
  // set up diagnostics struct
  struct Diagnostics * theDiagnostics = diagnostics_create(theSim,theTimeStep,theMPIsetup);
 
  // Initial Diagnostic dump
  diagnostics_set(theDiagnostics,theCells,theSim,theTimeStep,theMPIsetup,theGravMasses);
  MPI_Barrier(sim_comm);    
  diagnostics_print(theDiagnostics,theTimeStep,theSim,theMPIsetup);

  // Initial Checkpoint
  if( timestep_get_t(theTimeStep)>=io_tcheck(theIO))
  {
    io_allocbuf(theIO,theSim);
    io_setbuf(theIO,theCells,theSim,theGravMasses);
    io_hdf5_out(theIO,theSim,theTimeStep);
    io_deallocbuf(theIO);
  }

  int mu,nu,la;
  struct Metric *g = metric_create(0,3,0,1,theSim);
  double gmn[16], dgmn[64], igmn[16], a, da[4], b[3], db[12];
  a = metric_lapse(g);
  for(mu=0; mu<3; mu++)
      b[mu] = metric_shift_u(g,mu);
  for(mu=0; mu<4; mu++)
      for(nu=0; nu<4; nu++)
      {
          gmn[4*mu+nu] = metric_g_dd(g,mu,nu);
          igmn[4*mu+nu] = metric_g_uu(g,mu,nu);
      }
  for(la=0; la<4; la++)
  {
      da[la] = metric_dlapse(g,la);
      for(mu=0; mu<3; mu++)
        db[3*la+mu] = metric_dshift_u(g,la,mu);
      for(mu=0; mu<4; mu++)
          for(nu=0; nu<4; nu++)
              dgmn[16*la+4*mu+nu] = metric_dg_dd(g,la,mu,nu);
  }
  int k[4];
  for(mu = 0; mu < 4; mu++)
      k[mu] = metric_killcoord(g, mu);
  metric_destroy(g);
  printf("***** Metric Test *****\n");
  printf("a = %.6f\n", a);
  printf("da = (%.6f, %.6f, %.6f, %.6f)\n", da[0], da[1], da[2], da[3]);
  printf("\n");
  printf("b = (%.6f, %.6f, %.6f)\n", b[0], b[1], b[2]);
  printf("dbdt = (%.6f, %.6f, %.6f)\n", db[0], db[1], db[2]);
  printf("dbdr = (%.6f, %.6f, %.6f)\n", db[3], db[4], db[5]);
  printf("dbdp = (%.6f, %.6f, %.6f)\n", db[6], db[7], db[8]);
  printf("dbdz = (%.6f, %.6f, %.6f)\n", db[9], db[10], db[11]);
  printf("\n");
  printf("     (%.6f, %.6f, %.6f, %.6f)\n", gmn[0],gmn[1],gmn[2],gmn[3]);
  printf(" g = (%.6f, %.6f, %.6f, %.6f)\n", gmn[4],gmn[5],gmn[6],gmn[7]);
  printf("     (%.6f, %.6f, %.6f, %.6f)\n", gmn[8],gmn[9],gmn[10],gmn[11]);
  printf("     (%.6f, %.6f, %.6f, %.6f)\n", gmn[12],gmn[13],gmn[14],gmn[15]);
  printf("\n");
  printf("     (%.6f, %.6f, %.6f, %.6f)\n", igmn[0],igmn[1],igmn[2],igmn[3]);
  printf("ig = (%.6f, %.6f, %.6f, %.6f)\n", igmn[4],igmn[5],igmn[6],igmn[7]);
  printf("     (%.6f, %.6f, %.6f, %.6f)\n", igmn[8],igmn[9],igmn[10],igmn[11]);
  printf("     (%.6f, %.6f, %.6f, %.6f)\n", igmn[12],igmn[13],igmn[14],igmn[15]);
  printf("\n");
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[0],dgmn[1],dgmn[2],dgmn[3]);
  printf("dgdt = (%.6f, %.6f, %.6f, %.6f)\n", dgmn[4],dgmn[5],dgmn[6],dgmn[7]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[8],dgmn[9],dgmn[10],dgmn[11]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[12],dgmn[13],dgmn[14],dgmn[15]);
  printf("\n");
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[16],dgmn[17],dgmn[18],dgmn[19]);
  printf("dgdr = (%.6f, %.6f, %.6f, %.6f)\n", dgmn[20],dgmn[21],dgmn[22],dgmn[23]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[24],dgmn[25],dgmn[26],dgmn[27]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[28],dgmn[29],dgmn[30],dgmn[31]);
  printf("\n");
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[32],dgmn[33],dgmn[34],dgmn[35]);
  printf("dgdp = (%.6f, %.6f, %.6f, %.6f)\n", dgmn[36],dgmn[37],dgmn[38],dgmn[39]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[40],dgmn[41],dgmn[42],dgmn[43]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[44],dgmn[45],dgmn[46],dgmn[47]);
  printf("\n");
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[48],dgmn[49],dgmn[50],dgmn[51]);
  printf("dgdz = (%.6f, %.6f, %.6f, %.6f)\n", dgmn[52],dgmn[53],dgmn[54],dgmn[55]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[56],dgmn[57],dgmn[58],dgmn[59]);
  printf("       (%.6f, %.6f, %.6f, %.6f)\n", dgmn[60],dgmn[61],dgmn[62],dgmn[63]);
  printf("\n");
  printf(" k = (%d, %d, %d, %d)\n", k[0], k[1], k[2], k[3]);
  printf("***** Metric Test End *****\n");

  //Begin the run
  printf("Let's go!\nStarting Loop.\n");
  while( timestep_get_t(theTimeStep) < sim_get_T_MAX(theSim) ){
    // here the actual timestep is taken
    //timestep_forward_euler(theTimeStep,theSim,theCells,theGravMasses,theMPIsetup);
    timestep_rk2(theTimeStep,theSim,theCells,theGravMasses,theMPIsetup);

    // calculate diagnostics
    diagnostics_set(theDiagnostics,theCells,theSim,theTimeStep,theMPIsetup,theGravMasses);
    //write diagnostics to file
    MPI_Barrier(sim_comm);    
    diagnostics_print(theDiagnostics,theTimeStep,theSim,theMPIsetup);

    // checkpointing
    if( timestep_get_t(theTimeStep)>=io_tcheck(theIO)){ // time to write checkpoint file
      io_allocbuf(theIO,theSim); // allocate memory for a buffer to store simulation data
      io_setbuf(theIO,theCells,theSim,theGravMasses); // fill the buffer
      io_hdf5_out(theIO,theSim,theTimeStep); // write contents to file
      io_deallocbuf(theIO); // get rid of buffer
    }
  }

  //Make sure final checkpoint exists.
  if(io_nfile(theIO) < sim_NUM_CHECKPOINTS(theSim))
  {
    io_allocbuf(theIO,theSim);
    io_setbuf(theIO,theCells,theSim,theGravMasses);
    io_hdf5_out(theIO,theSim,theTimeStep);
    io_deallocbuf(theIO);
  }

  //inter-processor syncs
  cell_syncproc_r(theCells,theSim,theMPIsetup);
  cell_syncproc_z(theCells,theSim,theMPIsetup);

  // clean up
  diagnostics_destroy(theDiagnostics,theSim);
  cell_destroy(theCells,theSim);
  sim_destroy(theSim);
  gravMass_destroy(theGravMasses);
  io_destroy(theIO);
  timestep_destroy(theTimeStep);
  mpisetup_destroy(theMPIsetup);
  return(0);
}
