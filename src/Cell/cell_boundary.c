#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/TimeStep.h"
#include "../Headers/header.h"


void cell_boundary_outflow_r_inner( struct Cell *** theCells , struct Face * theFaces ,struct Sim * theSim,struct MPIsetup * theMPIsetup, struct TimeStep * theTimeStep ){
  int NUM_Q = sim_NUM_Q(theSim);

  int n,q;
  int i,j,k;
  double r_face,r_face_m1,r_cell;

  if( mpisetup_check_rin_bndry(theMPIsetup) ){

    if (sim_NoInnerBC(theSim)!=1){ // if the global inner radius is set negative, we don't apply an inner BC
      for( i=0; i>=0 ; --i ){
      //for( i=sim_Nghost_min(theSim,R_DIR)-1 ; i>=0 ; --i ){   //TODO: Good idea/bad idea?
        r_face=sim_FacePos(theSim,i,R_DIR);
        r_face_m1=sim_FacePos(theSim,i-1,R_DIR);
        r_cell = 0.5*(r_face+r_face_m1);
        //printf("r_cell inner: %e\n",r_cell);

        for( n=timestep_n(theTimeStep,i,R_DIR) ; n<timestep_n(theTimeStep,i+1,R_DIR) ; ++n ){
          for( q=0 ; q<NUM_Q ; ++q ){
            face_L_pointer(theFaces,n)->prim[q] = 0.0;
          }
        } 
        for( n=timestep_n(theTimeStep,i,R_DIR) ; n<timestep_n(theTimeStep,i+1,R_DIR) ; ++n ){
          struct Cell * cL = face_L_pointer(theFaces,n);
          struct Cell * cR = face_R_pointer(theFaces,n);
          for( q=0 ; q<NUM_Q ; ++q ){
            cL->prim[q] += cR->prim[q]*face_dA(theFaces,n);
          }
        }

        for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
          double zp = sim_FacePos(theSim,k,Z_DIR);
          double zm = sim_FacePos(theSim,k-1,Z_DIR);
          double dz = zp-zm;
          for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
            double dA = dz*r_face*theCells[k][i][j].dphi;
            for( q=0 ; q<NUM_Q ; ++q ){
              theCells[k][i][j].prim[q] /= dA;
            }
            //printf("r: %e, theCells[k][i][j].prim[URR]: %e\n",r_face,theCells[k][i][j].prim[URR]);
            if( theCells[k][i][j].prim[URR] > 0.0 ) theCells[k][i][j].prim[URR] = 0.0;
            if (KEP_BNDRY==1){
              theCells[k][i][j].prim[UPP] = pow(r_cell,-1.5);          
            }
          }
        }
      }
    }
  }
}

void cell_boundary_outflow_r_outer( struct Cell *** theCells , struct Face * theFaces ,struct Sim * theSim,struct MPIsetup * theMPIsetup, struct TimeStep * theTimeStep ){
  int Nf = timestep_n(theTimeStep,sim_N(theSim,R_DIR)-1,R_DIR);
  int NUM_Q = sim_NUM_Q(theSim);
  int n1 = timestep_n(theTimeStep,sim_N(theSim,R_DIR)-2,R_DIR);

  int n,q;
  int j,k;
  double r_face,r_face_p1,r_cell;

  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( n=n1 ; n<Nf ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_R_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=n1 ; n<Nf ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(theFaces,n);
      }
    }
    r_face = sim_FacePos(theSim,sim_N(theSim,R_DIR)-2,R_DIR);
    r_face_p1 = sim_FacePos(theSim,sim_N(theSim,R_DIR)-1,R_DIR);
    r_cell = 0.5*(r_face+r_face_p1);
    //printf("r_cell outer: %e\n",r_cell);
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
      double zp = sim_FacePos(theSim,k,Z_DIR);
      double zm = sim_FacePos(theSim,k-1,Z_DIR);
      double dz = zp-zm;
      for( j=0 ; j<sim_N_p(theSim,sim_N(theSim,R_DIR)-1) ; ++j ){
        double dA = dz*r_face*theCells[k][sim_N(theSim,R_DIR)-1][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[k][sim_N(theSim,R_DIR)-1][j].prim[q] /= dA;
        }
        //theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] *= -1.;
        if( theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] < 0.0 ) theCells[k][sim_N(theSim,R_DIR)-1][j].prim[URR] = 0.0;
        if (KEP_BNDRY==1){
          theCells[k][sim_N(theSim,R_DIR)-1][j].prim[UPP] = pow(r_cell,-1.5);
        }
      }
    }
  }
}

void cell_boundary_outflow_z_bot( struct Cell *** theCells , struct Face * theFaces, struct Sim * theSim,struct MPIsetup * theMPIsetup,struct TimeStep * theTimeStep ){
  int NUM_Q = sim_NUM_Q(theSim);

  int j,i;
  int n0 = timestep_n(theTimeStep,1,Z_DIR);

  int n,q;

  if(mpisetup_check_zbot_bndry(theMPIsetup)){
    for( n=0 ; n<n0 ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_L_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=0 ; n<n0 ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cL->prim[q] += cR->prim[q]*face_dA(theFaces,n);
      }
    }

    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double dA = .5*(rp*rp-rm*rm)*theCells[0][i][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[0][i][j].prim[q] /= dA;
        }
        theCells[0][i][j].prim[UZZ] *= -1.0;
      }
    }
  }
}


void cell_boundary_outflow_z_top( struct Cell *** theCells , struct Face * theFaces, struct Sim * theSim,struct MPIsetup * theMPIsetup,struct TimeStep * theTimeStep ){
  int NUM_Q = sim_NUM_Q(theSim);

  int j,i;
  int Nf = timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-1,Z_DIR);
  int n1 = timestep_n(theTimeStep,sim_N(theSim,Z_DIR)-2,Z_DIR);

  int n,q;

  if(mpisetup_check_ztop_bndry(theMPIsetup)){
    for( n=n1 ; n<Nf ; ++n ){
      for( q=0 ; q<NUM_Q ; ++q ){
        face_R_pointer(theFaces,n)->prim[q] = 0.0;
      }
    }
    for( n=n1 ; n<Nf ; ++n ){
      struct Cell * cL = face_L_pointer(theFaces,n);
      struct Cell * cR = face_R_pointer(theFaces,n);
      for( q=0 ; q<NUM_Q ; ++q ){
        cR->prim[q] += cL->prim[q]*face_dA(theFaces,n);
      }
    }


    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      double rp = sim_FacePos(theSim,i,R_DIR);
      double rm = sim_FacePos(theSim,i-1,R_DIR);
      for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
        double dA = .5*(rp*rp-rm*rm)*theCells[sim_N(theSim,Z_DIR)-1][i][j].dphi;
        for( q=0 ; q<NUM_Q ; ++q ){
          theCells[sim_N(theSim,Z_DIR)-1][i][j].prim[q] /= dA;
        }
        theCells[sim_N(theSim,Z_DIR)-1][i][j].prim[UZZ] *= -1.0;
      }
    }
  }
}

void cell_boundary_fixed_r_inner( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup, 
    void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int) ){
  int i,j,k;
  if (sim_NoInnerBC(theSim)!=1){ // if the global inner radius is set negative, we don't apply an inner BC
    if(mpisetup_check_rin_bndry(theMPIsetup)){
      for( i=0 ; i<sim_Nghost_min(theSim,R_DIR) ; ++i ){
        for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
          for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
            (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
          }
        }
      }
    }
  }

  //TODO: Find a better way of doing this.
  if(sim_InitialDataType(theSim)==GBONDI && sim_NoInnerBC(theSim))
      if(mpisetup_check_rin_bndry(theMPIsetup))
          for(k=0; k<sim_N(theSim,Z_DIR); k++)
          {
              double zm = sim_FacePos(theSim,k-1,Z_DIR);
              double zp = sim_FacePos(theSim,k,Z_DIR);
              if(zm < 0 && zp > 0)
                for(j=0; j<sim_N_p(theSim,0); j++)
                    (*single_init_ptr)(&(theCells[k][0][j]),theSim,0,j,k);
          }
}

void cell_boundary_fixed_r_outer( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup, 
    void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int) ){
    int i,j,k;
  if( mpisetup_check_rout_bndry(theMPIsetup) ){
    for( i=sim_N(theSim,R_DIR)-1 ; i>sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR)-1 ; --i ){
      for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
        }
      }
    }
  }

}

void cell_boundary_fixed_z_bot( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup,
    void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int)  ){
  int i,j,k;
  if(mpisetup_check_zbot_bndry(theMPIsetup)){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( k=0 ; k<sim_Nghost_min(theSim,Z_DIR) ; ++k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
        }
      }
    }
  }
}

void cell_boundary_fixed_z_top( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup,
    void (*single_init_ptr)(struct Cell *,struct Sim *,int,int,int)  ){
  int i,j,k;
  if( mpisetup_check_ztop_bndry(theMPIsetup) ){
    for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
      for( k=sim_N(theSim,Z_DIR)-1 ; k>sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR)-1 ; --k ){
        for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
          (*single_init_ptr)(&(theCells[k][i][j]),theSim,i,j,k);
        }
      }
    }
  }

}

void cell_boundary_ssprofile_r_inner( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup)
{
    int i,j,k,q;
    if (sim_NoInnerBC(theSim)!=1) // if the global inner radius is set negative, we don't apply an inner BC
    {
        if(mpisetup_check_rin_bndry(theMPIsetup))
        {
            int iIn = sim_Nghost_min(theSim, R_DIR);
            double rIn = 0.5*(sim_FacePos(theSim,iIn,R_DIR) + sim_FacePos(theSim,iIn-1,R_DIR));
            double *primIn = theCells[0][iIn][0].prim;
            for( i=0 ; i<iIn ; ++i )
            {
                double rp = sim_FacePos(theSim,i,R_DIR);
                double rm = sim_FacePos(theSim,i-1,R_DIR);
                double r = 0.5*(rp+rm);
                for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k )
                    for( j=0 ; j<sim_N_p(theSim,i) ; ++j )
                    {
                        if(sim_InitialDataType(theSim) == SSDISC)
                        {
                            if(sim_InitPar0(theSim)==0)
                            {
                                theCells[k][i][j].prim[RHO] = primIn[RHO] * pow(r/rIn, -0.6);
                                theCells[k][i][j].prim[PPP] = primIn[PPP] * pow(r/rIn, -1.5);
                                theCells[k][i][j].prim[URR] = primIn[URR] * pow(r/rIn, -0.4);
                                theCells[k][i][j].prim[UPP] = primIn[UPP] * pow(r/rIn, -1.5);
                                theCells[k][i][j].prim[UZZ] = 0.0;
                                for(q=sim_NUM_C(theSim); q < sim_NUM_Q(theSim); q++)
                                    theCells[k][i][j].prim[q] = primIn[q] * pow(r/rIn, -0.6);
                            }
                            else if(sim_InitPar0(theSim)==1)
                            {
                                theCells[k][i][j].prim[RHO] = primIn[RHO] * pow(r/rIn, -1.5);
                                theCells[k][i][j].prim[PPP] = primIn[PPP] * pow(r/rIn, -1.5);
                                theCells[k][i][j].prim[URR] = primIn[URR] * pow(r/rIn, 0.5);
                                theCells[k][i][j].prim[UPP] = primIn[UPP] * pow(r/rIn, -1.5);
                                theCells[k][i][j].prim[UZZ] = 0.0;
                                for(q=sim_NUM_C(theSim); q < sim_NUM_Q(theSim); q++)
                                    theCells[k][i][j].prim[q] = primIn[q] * pow(r/rIn, -1.5);
                            }
                        }
                        else if(sim_InitialDataType(theSim) == NTDISC)
                        {
                            double M = sim_GravM(theSim);
                            double rs = 6*M;
                            double C = (1.0-3*M/r)/(1.0-3*M/rIn);
                            double D = (1.0-2*M/r)/(1.0-2*M/rIn);
                            double P = (1.0 - sqrt(rs/r) + sqrt(3*M/r)*(atanh(sqrt(3*M/r))-atanh(sqrt(3*M/rs)))) / (1.0 - sqrt(rs/rIn) + sqrt(3*M/rIn)*(atanh(sqrt(3*M/rIn))-atanh(sqrt(3*M/rs))));
                            if(sim_InitPar0(theSim)==0)
                            {
                                theCells[k][i][j].prim[RHO] = primIn[RHO] * pow(r/rIn,-0.6) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
                                theCells[k][i][j].prim[PPP] = primIn[PPP] * pow(r/rIn,-1.5) * sqrt(C) * P / (D*D);
                                theCells[k][i][j].prim[URR] = primIn[URR] * pow(r/rIn,-0.4) * pow(C,-0.1) * pow(D,1.6) * pow(P,-0.6);
                                theCells[k][i][j].prim[UPP] = primIn[UPP] * pow(r/rIn,-1.5);
                                theCells[k][i][j].prim[UZZ] = 0.0;
                                for(q=sim_NUM_C(theSim); q < sim_NUM_Q(theSim); q++)
                                    theCells[k][i][j].prim[q] = primIn[q] * pow(r/rIn,-0.6) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
                            }
                        }
                        else if (sim_InitialDataType(theSim) == ADAF)
                        {
                            theCells[k][i][j].prim[RHO] = primIn[RHO] * pow(r/rIn,-0.5);
                            theCells[k][i][j].prim[PPP] = primIn[PPP] * pow(r/rIn,-1.5);
                            theCells[k][i][j].prim[URR] = primIn[URR] * pow(r/rIn,-0.5);
                            theCells[k][i][j].prim[UPP] = primIn[UPP] * pow(r/rIn,-1.5);
                        }

                    }
            }
        }
    }
}

void cell_boundary_ssprofile_r_outer( struct Cell *** theCells, struct Sim *theSim,struct MPIsetup * theMPIsetup)
{
    int i,j,k,q;
   if(mpisetup_check_rout_bndry(theMPIsetup))
    {
        int iOut = sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR)-1;
        double rOut = 0.5*(sim_FacePos(theSim,iOut,R_DIR) + sim_FacePos(theSim,iOut-1,R_DIR));
        double *primOut = theCells[0][iOut][0].prim;
        for(i=iOut+1; i<sim_N(theSim,R_DIR); i++)
        {
            double rp = sim_FacePos(theSim,i,R_DIR);
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double r = 0.5*(rp+rm);
            for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k )
                for( j=0 ; j<sim_N_p(theSim,i) ; ++j )
                {
                    if(sim_InitialDataType(theSim) == SSDISC)
                    {
                        if(sim_InitPar0(theSim)==0)
                        {
                            theCells[k][i][j].prim[RHO] = primOut[RHO] * pow(r/rOut, -0.6);
                            theCells[k][i][j].prim[PPP] = primOut[PPP] * pow(r/rOut, -1.5);
                            theCells[k][i][j].prim[URR] = primOut[URR] * pow(r/rOut, -0.4);
                            theCells[k][i][j].prim[UPP] = primOut[UPP] * pow(r/rOut, -1.5);
                            theCells[k][i][j].prim[UZZ] = 0.0;
                            for(q=sim_NUM_C(theSim); q < sim_NUM_Q(theSim); q++)
                                theCells[k][i][j].prim[q] = primOut[q] * pow(r/rOut, -0.6);
                        }
                        else if(sim_InitPar0(theSim)==1)
                        {
                            theCells[k][i][j].prim[RHO] = primOut[RHO] * pow(r/rOut, -1.5);
                            theCells[k][i][j].prim[PPP] = primOut[PPP] * pow(r/rOut, -1.5);
                            theCells[k][i][j].prim[URR] = primOut[URR] * pow(r/rOut, 0.5);
                            theCells[k][i][j].prim[UPP] = primOut[UPP] * pow(r/rOut, -1.5);
                            theCells[k][i][j].prim[UZZ] = 0.0;
                            for(q=sim_NUM_C(theSim); q < sim_NUM_Q(theSim); q++)
                                theCells[k][i][j].prim[q] = primOut[q] * pow(r/rOut, -1.5);
                        }
                    }
                    
                    else if(sim_InitialDataType(theSim) == NTDISC)
                    {
                        double M = sim_GravM(theSim);
                        double rs = 6*M;
                        double C = (1.0-3*M/r)/(1.0-3*M/rOut);
                        double D = (1.0-2*M/r)/(1.0-2*M/rOut);
                        double P = (1.0 - sqrt(rs/r) + sqrt(3*M/r)*(atanh(sqrt(3*M/r))-atanh(sqrt(3*M/rs)))) / (1.0 - sqrt(rs/rOut) + sqrt(3*M/rOut)*(atanh(sqrt(3*M/rOut))-atanh(sqrt(3*M/rs))));
                        if(sim_InitPar0(theSim)==0)
                        {
                            theCells[k][i][j].prim[RHO] = primOut[RHO] * pow(r/rOut,-0.6) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
                            theCells[k][i][j].prim[PPP] = primOut[PPP] * pow(r/rOut,-1.5) * sqrt(C) * P / (D*D);
                            theCells[k][i][j].prim[URR] = primOut[URR] * pow(r/rOut,-0.4) * pow(C,-0.1) * pow(D,1.6) * pow(P,-0.6);
                            theCells[k][i][j].prim[UPP] = primOut[UPP] * pow(r/rOut,-1.5);
                            theCells[k][i][j].prim[UZZ] = 0.0;
                            for(q=sim_NUM_C(theSim); q < sim_NUM_Q(theSim); q++)
                                theCells[k][i][j].prim[q] = primOut[q] * pow(r/rOut,-0.6) * pow(C,0.6) * pow(D,-1.6) * pow(P,0.6);
                        }
                    }
                    else if (sim_InitialDataType(theSim) == ADAF)
                    {
                        theCells[k][i][j].prim[RHO] = primOut[RHO] * pow(r/rOut,-0.5);
                        theCells[k][i][j].prim[PPP] = primOut[PPP] * pow(r/rOut,-1.5);
                        theCells[k][i][j].prim[URR] = primOut[URR] * pow(r/rOut,-0.5);
                        theCells[k][i][j].prim[UPP] = primOut[UPP] * pow(r/rOut,-1.5);
                    }
                }
        }
    }
}
