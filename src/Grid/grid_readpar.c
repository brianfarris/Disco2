#define GRID_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../Headers/Grid.h"
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

int grid_read_par_file(struct Grid * theGrid, struct MPIsetup * theMPIsetup, char * inputfilename){

  char * pfile = inputfilename;
  int err=0;  
  int nrank;
  for( nrank=0 ; nrank<mpisetup_NumProcs(theMPIsetup) ; ++nrank ){
    if( mpisetup_MyProc(theMPIsetup)==nrank ){
      err += readvar( pfile , "Restart"              , VAR_INT  , &(theGrid->Restart)  );
      err += readvar( pfile , "InitialDataType"      , VAR_INT  , &(theGrid->InitialDataType)  );
      err += readvar( pfile , "GravMassType"      , VAR_INT  , &(theGrid->GravMassType)  );
      err += readvar( pfile , "BoundTypeR"         , VAR_INT  , &(theGrid->BoundTypeR)  );
      err += readvar( pfile , "BoundTypeZ"         , VAR_INT  , &(theGrid->BoundTypeZ)  );
      err += readvar( pfile , "NumR"              , VAR_INT  , &(theGrid->N_r_global) );
      err += readvar( pfile , "NumZ"              , VAR_INT  , &(theGrid->N_z_global  ));
      err += readvar( pfile , "ng"              , VAR_INT  , &(theGrid->ng));
      err += readvar( pfile , "R_Min"             , VAR_DOUB , &(theGrid->RMIN)  );
      err += readvar( pfile , "R_Max"             , VAR_DOUB , &(theGrid->RMAX)  );
      err += readvar( pfile , "Z_Min"             , VAR_DOUB , &(theGrid->ZMIN)  );
      err += readvar( pfile , "Z_Max"             , VAR_DOUB , &(theGrid->ZMAX)  );
      err += readvar( pfile , "NP_CONST"             , VAR_INT , &(theGrid->NP_CONST)  );
      err += readvar( pfile , "aspect"             , VAR_DOUB , &(theGrid->aspect)  );
      err += readvar( pfile , "NUM_Q"              , VAR_INT  , &(theGrid->NUM_Q) );
      err += readvar( pfile , "Time_Max"       , VAR_DOUB , &(theGrid->T_MAX)  );
      err += readvar( pfile , "Num_Checkpoints"   , VAR_INT  , &(theGrid->NUM_CHECKPOINTS)  );
      err += readvar( pfile , "Move_Cells"        , VAR_INT  , &(theGrid->MOVE_CELLS)  );
      //     err += readvar( pfile , "NumGravMass"        , VAR_INT  , &(theGrid->NumGravMass)  );
      err += readvar( pfile , "Adiabatic_Index"   , VAR_DOUB , &(theGrid->GAMMALAW)  );
      err += readvar( pfile , "Include_Viscosity" , VAR_INT  , &(theGrid->INCLUDE_VISCOSITY)  );
      err += readvar( pfile , "Explicit_Viscosity" , VAR_DOUB  , &(theGrid->EXPLICIT_VISCOSITY)  );
      err += readvar( pfile , "DivB_Ch"           , VAR_DOUB , &(theGrid->DIVB_CH)  );
      err += readvar( pfile , "DivB_l"            , VAR_DOUB , &(theGrid->DIVB_L)  );
      err += readvar( pfile , "CFL"               , VAR_DOUB , &(theGrid->CFL)  );
      err += readvar( pfile , "PLM"               , VAR_DOUB , &(theGrid->PLM)  );
      err += readvar( pfile , "POWELL"               , VAR_DOUB , &(theGrid->POWELL)  );
      err += readvar( pfile , "Grav_2D"            , VAR_INT  , &(theGrid->GRAV2D)  );
      err += readvar( pfile , "G_EPS"             , VAR_DOUB , &(theGrid->G_EPS)  );
      err += readvar( pfile , "PHI_ORDER"             , VAR_DOUB , &(theGrid->PHI_ORDER)  );
      err += readvar( pfile , "Rho_Floor"         , VAR_DOUB , &(theGrid->RHO_FLOOR)  );
      err += readvar( pfile , "Cs_Floor"          , VAR_DOUB , &(theGrid->CS_FLOOR)  );
      err += readvar( pfile , "Cs_Cap"            , VAR_DOUB , &(theGrid->CS_CAP)  );
      err += readvar( pfile , "Vel_Cap"           , VAR_DOUB , &(theGrid->VEL_CAP)  );
      err += readvar( pfile , "runtype"           , VAR_INT , &(theGrid->runtype)  );
      //err += readvar( pfile , "BC_Damp"           , VAR_INT  , &(theGrid->BC_DAMP));
      //err += readvar( pfile , "Num_Reports"       , VAR_INT  , &(theGrid->NUM_REPORTS)  );
      //err += readvar( pfile , "Num_Outputs"       , VAR_INT  , &(theGrid->NUM_OUTPUTS)  );
      //err += readvar( pfile , "Restart"           , VAR_INT  , &(theGrid->RESTART)  );
      //err += readvar( pfile , "Mp_Thresh"         , VAR_DOUB , &(theGrid->PLANET_MASS)  );
      //err += readvar( pfile , "Disk_Mach"         , VAR_DOUB , &(theGrid->DISK_MACH)  );
      //err += readvar( pfile , "Free_Planet"       , VAR_INT  , &(theGrid->FREE_PLANET)  );
      //err += readvar( pfile , "Damp_Timescale"    , VAR_DOUB , &(theGrid->DAMP_TIME)  );
      //err += readvar( pfile , "Np_Const"          , VAR_INT  , &(theGrid->NP_CONST)           );
      //err += readvar( pfile , "Cell_Aspect_Ratio" , VAR_DOUB , &(theGrid->CELL_ASPECT_RATIO)  );
      //err += readvar( pfile , "Logarithmic_Zoning", VAR_DOUB , &(theGrid->LOG_ZONING)  );
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



