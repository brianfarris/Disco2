#define SIM_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../Headers/Sim.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"

int readvar( char * filename , char * varname , int vartype , void * ptr ){
  FILE * inFile = fopen( filename , "r" );
  char s[512];
  char nm[512];
  char s1[512];
  int found = 0;
  while( (fgets(s,512,inFile) != NULL) && found==0 ){
    sscanf(s,"%s ",nm);
    if( strcmp(nm,varname)==0 ){
      strcpy(s1,s);
      found=1;
    }
  }

  fclose( inFile );
  if( found==0 ) {
    printf("cant find %s\n",varname);
    return(1);
  }
  char * s2 = s1+strlen(nm)+strspn(s1+strlen(nm),"\t :=>_");

  double temp;
  char stringval[256];

  sscanf(s2,"%lf",&temp);
  sscanf(s2,"%256s",stringval);

  if( vartype == VAR_INT ){
    *((int *)   ptr) = (int)temp;
  }else if( vartype == VAR_DOUB ){
    *((double *)ptr) = (double)temp;
  }else{
    strcpy( ptr , stringval );
  }

  return(0);
}

int sim_read_par_file(struct Sim * theSim, struct MPIsetup * theMPIsetup, char * inputfilename){

  char * pfile = inputfilename;
  int err=0;  
  int nrank;
  for( nrank=0 ; nrank<mpisetup_NumProcs(theMPIsetup) ; ++nrank ){
    if( mpisetup_MyProc(theMPIsetup)==nrank ){
      err += readvar( pfile , "Restart"              , VAR_INT  , &(theSim->Restart)  );
      err += readvar( pfile , "InitialDataType"      , VAR_INT  , &(theSim->InitialDataType)  );
      err += readvar( pfile , "GravMassType"      , VAR_INT  , &(theSim->GravMassType)  );
      err += readvar( pfile , "MassRatio"      , VAR_DOUB  , &(theSim->MassRatio)  );
      err += readvar( pfile , "OrbShrinkTscale"      , VAR_DOUB  , &(theSim->OrbShrinkTscale)  );
      err += readvar( pfile , "OrbShrinkT0"      , VAR_DOUB  , &(theSim->OrbShrinkT0)  );
      err += readvar( pfile , "BoundTypeR"         , VAR_INT  , &(theSim->BoundTypeR)  );
      err += readvar( pfile , "BoundTypeZ"         , VAR_INT  , &(theSim->BoundTypeZ)  );
      err += readvar( pfile , "ZeroPsiBndry"         , VAR_INT  , &(theSim->ZeroPsiBndry)  );
      err += readvar( pfile , "NoInnerBC"         , VAR_INT  , &(theSim->NoInnerBC)  );
      err += readvar( pfile , "NumR"              , VAR_INT  , &(theSim->N_global[R_DIR]) );
      err += readvar( pfile , "NumZ"              , VAR_INT  , &(theSim->N_global[Z_DIR]) );
      err += readvar( pfile , "ng"              , VAR_INT  , &(theSim->ng));
      err += readvar( pfile , "R_Min"             , VAR_DOUB , &(theSim->MIN[R_DIR])  );
      err += readvar( pfile , "R_Max"             , VAR_DOUB , &(theSim->MAX[R_DIR])  );
      err += readvar( pfile , "Z_Min"             , VAR_DOUB , &(theSim->MIN[Z_DIR])  );
      err += readvar( pfile , "Z_Max"             , VAR_DOUB , &(theSim->MAX[Z_DIR])  );
      err += readvar( pfile , "NP_CONST"             , VAR_INT , &(theSim->NP_CONST)  );
      err += readvar( pfile , "NPCAP"             , VAR_INT , &(theSim->NPCAP)  );
      err += readvar( pfile , "aspect"             , VAR_DOUB , &(theSim->aspect)  );
      err += readvar( pfile , "NUM_C"              , VAR_INT  , &(theSim->NUM_C) );
      err += readvar( pfile , "NUM_N"              , VAR_INT  , &(theSim->NUM_N) );
      err += readvar( pfile , "Time_Max"       , VAR_DOUB , &(theSim->T_MAX)  );
      err += readvar( pfile , "Num_Checkpoints"   , VAR_INT  , &(theSim->NUM_CHECKPOINTS)  );
      err += readvar( pfile , "Num_Diag_Dump"   , VAR_INT  , &(theSim->NUM_DIAG_DUMP)  );
      err += readvar( pfile , "Num_Diag_Measure"   , VAR_INT  , &(theSim->NUM_DIAG_MEASURE)  );
      err += readvar( pfile , "Move_Cells"        , VAR_INT  , &(theSim->MOVE_CELLS)  );
      err += readvar( pfile , "RiemannSolver"   , VAR_INT , &(theSim->Riemann)  );
      err += readvar( pfile , "Adiabatic_Index"   , VAR_DOUB , &(theSim->GAMMALAW)  );
      err += readvar( pfile , "Cooling"   , VAR_INT , &(theSim->COOLING)  );
      err += readvar( pfile , "Explicit_Viscosity" , VAR_DOUB  , &(theSim->EXPLICIT_VISCOSITY)  );
      err += readvar( pfile , "ViscConst"         , VAR_INT  , &(theSim->VISC_CONST)  );
      err += readvar( pfile , "DivB_Ch"           , VAR_DOUB , &(theSim->DIVB_CH)  );
      err += readvar( pfile , "DivB_l"            , VAR_DOUB , &(theSim->DIVB_L)  );
      err += readvar( pfile , "CFL"               , VAR_DOUB , &(theSim->CFL)  );
      err += readvar( pfile , "PLM"               , VAR_DOUB , &(theSim->PLM)  );
      err += readvar( pfile , "POWELL"               , VAR_INT , &(theSim->POWELL)  );
      err += readvar( pfile , "Grav_2D"            , VAR_INT  , &(theSim->GRAV2D)  );
      err += readvar( pfile , "G_EPS"             , VAR_DOUB , &(theSim->G_EPS)  );
      err += readvar( pfile , "RhoSinkOn"   , VAR_INT, &(theSim->RhoSinkOn)  );
      err += readvar( pfile , "PHI_ORDER"             , VAR_DOUB , &(theSim->PHI_ORDER)  );
      err += readvar( pfile , "Rho_Floor"         , VAR_DOUB , &(theSim->RHO_FLOOR)  );
      err += readvar( pfile , "Cs_Floor"          , VAR_DOUB , &(theSim->CS_FLOOR)  );
      err += readvar( pfile , "Cs_Cap"            , VAR_DOUB , &(theSim->CS_CAP)  );
      err += readvar( pfile , "Vel_Cap"           , VAR_DOUB , &(theSim->VEL_CAP)  );
      err += readvar( pfile , "SET_T"           , VAR_INT , &(theSim->SET_T)  );
      err += readvar( pfile , "HoR"           , VAR_DOUB , &(theSim->HoR)  );
      err += readvar( pfile , "runtype"           , VAR_INT , &(theSim->runtype)  );
      err += readvar( pfile , "DAMP_TIME"           , VAR_DOUB , &(theSim->DAMP_TIME)  );
      err += readvar( pfile , "RDAMP_INNER"           , VAR_DOUB , &(theSim->RDAMP_INNER)  );
      err += readvar( pfile , "RDAMP_OUTER"           , VAR_DOUB , &(theSim->RDAMP_OUTER)  );
      err += readvar( pfile , "RLogScale"           , VAR_DOUB , &(theSim->RLogScale)  );
      err += readvar( pfile , "ZLogScale"           , VAR_DOUB , &(theSim->ZLogScale)  );
      err += readvar( pfile , "HiResSigma"           , VAR_DOUB , &(theSim->HiResSigma)  );
      err += readvar( pfile , "HiResR0a"           , VAR_DOUB , &(theSim->HiResR0a)  );
      err += readvar( pfile , "HiResR0b"           , VAR_DOUB , &(theSim->HiResR0b)  );
      err += readvar( pfile , "HiResFac"           , VAR_DOUB , &(theSim->HiResFac)  );
      err += readvar( pfile , "w_a_type"            , VAR_INT , &(theSim->W_A_TYPE));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  int errtot;
  MPI_Allreduce( &err , &errtot , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );

  if( errtot > 0 ){
    printf("Read Failed\n");
    printf("there were %d errors\n",errtot);
    return(1);
  }

  return(0);

}



