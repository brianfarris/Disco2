#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

double maxvel(double * prim , double w , double r ,struct Sim * theSim){

    double Pp  = prim[PPP];
    double rho = prim[RHO];
    double vp;
    if (NO_W_IN_CFL==1){ 
        vp  = prim[UPP]*r;//-w;
    } else{
        //vp = prim[UPP]*r-sim_W_A(theSim,r);
        vp  = prim[UPP]*r-w;
    }
    double vr  = prim[URR];
    double vz  = prim[UZZ];
    double cs  = sqrt(sim_GAMMALAW(theSim)*Pp/rho);

    double Br  = prim[BRR];
    double Bp  = prim[BPP];
    double Bz  = prim[BZZ];
    double b2  = (Br*Br+Bp*Bp+Bz*Bz)/rho;

    double cf2 = cs*cs + b2;

    double maxv = sqrt(cf2) + sqrt( vr*vr + vp*vp + vz*vz );

    return(maxv);

}

double cell_mindt( struct Cell *** theCells, struct Sim * theSim, struct GravMass * theGravMasses ){
    int i_m,j_m,k_m;
    double dt_m = 1.e100;//HUGE_VAL;
    double a_m,r_m,dx_m;
    double mag_vel_m;
    int i,j,k;
    for( k=sim_Nghost_min(theSim,Z_DIR) ; k<sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR) ; ++k ){
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        double dz = zp-zm;
        for( i=sim_Nghost_min(theSim,R_DIR) ; i<sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR) ; ++i ){
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            double dr = rp-rm;
            double r = .5*(rp+rm);
            for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
                int jm = j-1;
                if( j==0 ) jm = sim_N_p(theSim,i)-1;
                double w = .5*(theCells[k][i][j].wiph+theCells[k][i][jm].wiph);
                double dx = dr;
                double rdphi = .5*(rp+rm)*theCells[k][i][j].dphi;
                if( rdphi<dr ) {
                    dx = rdphi;
                }
                if( dx>dz ) {
                    dx = dz;
                }
                double a  = maxvel( theCells[k][i][j].prim , w , r ,theSim);
                double rho = theCells[k][i][j].prim[RHO]; 
                double Pp  = theCells[k][i][j].prim[PPP]; 


                double dt = sim_CFL(theSim)*dx/a;
                if( sim_EXPLICIT_VISCOSITY(theSim)>0.0 ){
                    double nu;
                    if (sim_VISC_CONST(theSim)==1){
                        nu = sim_EXPLICIT_VISCOSITY(theSim);
                    } else{
                        double tiph = theCells[k][i][j].tiph - 0.5*theCells[k][i][j].dphi;
                        if (sim_InitialDataType(theSim)==SHEAR){
                            double HoR = 0.1;
                            //nu = sim_EXPLICIT_VISCOSITY(theSim)*HoR*HoR*pow(fabs((r*cos(tiph))),1.5);
                            nu = sim_EXPLICIT_VISCOSITY(theSim)*sim_GAMMALAW(theSim)*Pp/rho*pow(fabs(r*cos(tiph)),2.0);
                            if (r*cos(tiph)>20.) nu=0.000000001;
                        } else{
                            double M0 = gravMass_M(theGravMasses,0);
                            double M1 = gravMass_M(theGravMasses,1);
                            double dist_bh0 = gravMass_dist(theGravMasses,0,r,tiph,0.);
                            double dist_bh1 = gravMass_dist(theGravMasses,1,r,tiph,0.);
                            double alpha = sim_EXPLICIT_VISCOSITY(theSim);
                            double eps = sim_G_EPS(theSim);
                            //nu = sim_EXPLICIT_VISCOSITY(theSim)*sim_GAMMALAW(theSim)*Pp/rho*pow(r,1.5);
                            nu = alpha*Pp/rho/sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
                        }
                    }
                    double dt_visc = .25*dx*dx/nu;
                    dt = dt/( 1. + dt/dt_visc );
                }
                if( dt_m > dt ) {
                    dt_m = dt;
                }
            } 
        }
    }
    double dt2;
    MPI_Allreduce( &dt_m , &dt2 , 1 , MPI_DOUBLE , MPI_MIN , sim_comm );
    return( dt2 );

}



