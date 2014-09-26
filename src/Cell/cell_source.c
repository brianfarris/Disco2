#define CELL_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/header.h"

double fgrav( double M , double r , double eps, double n ){
    return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
}

void get_rho_sink( struct GravMass * theGravMasses, struct Sim * theSim, int p, double dt, double r0,double r1,double M0, double M1, double rho, double P, double * drho_dt_sink){

    double Mtotal = M0+M1;
    double alpha = sim_EXPLICIT_VISCOSITY(theSim);

    double one_o_nu = 1.0/(alpha*P/rho)*sqrt(pow(r0,-3)*M0/Mtotal + pow(r1,-3)*M1/Mtotal);
    double t_visc0 = TVISC_FAC * 2./3. * r0*r0*one_o_nu;
    double t_visc1 = TVISC_FAC * 2./3. * r1*r1*one_o_nu;    

    // Don't let the viscous timescale get too short, things could get unstable.
    if (t_visc0 < 10.* dt){
        t_visc0 = 10.*dt;
    }
    if (t_visc1 < 10.* dt){
        t_visc1 = 10.*dt;
    }

    *drho_dt_sink = 0.0;
    double sink_size = 0.5; //only apply the sink within a given radius.

    if (p==0){
        if (r0<sink_size){
            *drho_dt_sink = rho / t_visc0;
        }
    } else if (p==1){
        if (r1<sink_size){
            *drho_dt_sink = rho / t_visc1;
        }

    } else{
        exit(1);
    }
}

void gravMassForce( struct GravMass * theGravMasses ,struct Sim * theSim, int p , double r , double phi , double * fr , double * fp ){
    double G_EPS=sim_G_EPS(theSim);
    double PHI_ORDER=sim_PHI_ORDER(theSim);

    double rp = gravMass_r(theGravMasses,p);
    double pp = gravMass_phi(theGravMasses,p);
    double cosp = cos(phi);
    double sinp = sin(phi);
    double dx = r*cosp-rp*cos(pp);
    double dy = r*sinp-rp*sin(pp);
    double script_r = sqrt(dx*dx+dy*dy);

    double cosa = dx/script_r;
    double sina = dy/script_r;

    double cosap = cosa*cosp+sina*sinp;
    double sinap = sina*cosp-cosa*sinp;

    double f1;
    f1 = -fgrav( gravMass_M(theGravMasses,p) , script_r , G_EPS, PHI_ORDER );
    *fr = cosap*f1;
    *fp = sinap*f1;
}

void cell_add_src( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses , double dt ){
    int GRAV2D=sim_GRAV2D(theSim);
    int i,j,k;
    double Mdot_temp[2];
    Mdot_temp[0] = 0.0;
    Mdot_temp[1] = 0.0;

    int imin_noghost;
    if (sim_Nghost_min(theSim,R_DIR) == 1){
        imin_noghost = 0;
    } else{
        imin_noghost = sim_Nghost_min(theSim,R_DIR);
    }
    int kmin_noghost = sim_Nghost_min(theSim,Z_DIR);
    int imax_noghost = sim_N(theSim,R_DIR) - sim_Nghost_max(theSim,R_DIR);
    int kmax_noghost = sim_N(theSim,Z_DIR) - sim_Nghost_max(theSim,Z_DIR);

    double Mtotal = 1.0;
    double sep = 1.0;
    double M0 = gravMass_M(theGravMasses,0);
    double M1 = gravMass_M(theGravMasses,1);

    double rbh0 = gravMass_r(theGravMasses,0);
    double rbh1 = gravMass_r(theGravMasses,1);

    double a = rbh0+rbh1;

    double phi_bh0 = gravMass_phi(theGravMasses,0);
    double phi_bh1 = gravMass_phi(theGravMasses,1);

    double xbh0 = rbh0*cos(phi_bh0);
    double ybh0 = rbh0*sin(phi_bh0);

    double xbh1 = rbh1*cos(phi_bh1);
    double ybh1 = rbh1*sin(phi_bh1);

    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
                struct Cell * c = &(theCells[k][i][j]);
                double phi = c->tiph-.5*c->dphi;
                double dphi = c->dphi;

                double rho = c->prim[RHO];
                double Pp  = c->prim[PPP];
                double r   = .5*(rp+rm);
                double vr  = c->prim[URR];
                double vz  = c->prim[UZZ];
                double vp  = c->prim[UPP]*r;
                double dz = zp-zm;
                double dV = dphi*.5*(rp*rp-rm*rm)*dz;

                double xpos = r*cos(phi);
                double ypos = r*sin(phi);

                double r0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
                double r1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));

                double z  = .5*(zp+zm);
                double R  = sqrt(r*r+z*z);
                double sint;
                double cost;
                double gravdist;
                if (GRAV2D==1){
                    sint = 1.0;
                    cost = 0.0;
                    gravdist = r;
                }else{
                    sint = r/R;
                    cost = z/R;
                    gravdist = R;
                }
                double fr,fp;
                double Fr = 0.0;
                double Fp = 0.0;
                int p;
                for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
                    gravMassForce( theGravMasses , theSim , p , gravdist , phi , &fr , &fp );
                    Fr += fr;
                    Fp += fp;
                }
                double w_a = sim_rOm_a(theSim,r,a);
                double rdrOm_a = sim_rdrOm_a(theSim,r,a);
                double F_centrifugal_r = w_a*w_a/r;
                double F_coriolis_r =  2.*w_a*vp/r;
                double F_coriolis_phi = -2.*w_a*vr/r;
                double dtOm_a = sim_dtOm_a(theSim,r,a);
                double F_euler_phi =  -(r*dtOm_a + vr*rdrOm_a);

                c->cons[SRR] += dt*dV*( rho*vp*vp + Pp )/r;
                c->cons[SRR] += dt*dV*rho*(Fr*sint+F_centrifugal_r/*+F_coriolis_r*/);
                c->cons[LLL] += dt*dV*rho*(Fp/*+F_coriolis_phi+F_euler_phi*/ -r*dtOm_a)*r;
                c->cons[SZZ] += dt*dV*rho*Fr*cost;
                c->cons[TAU] += dt*dV*rho*( (Fr*sint+F_centrifugal_r)*vr+Fr*vz*cost + (Fp+F_euler_phi)*vp);

                if (sim_RhoSinkOn(theSim)==1){
                    double drho_dt_sink;
                    int num_sinks = sim_GravMassType(theSim);
                    for( p=0 ; p<num_sinks; ++p ){
                        get_rho_sink( theGravMasses, theSim,p,dt,r0,r1,M0,M1,rho,Pp, &drho_dt_sink);
                        c->cons[RHO] -= drho_dt_sink * dt * dV;
                        if ((i>=imin_noghost) && (i<imax_noghost) && (k>=kmin_noghost) && (k<kmax_noghost)){
                            Mdot_temp[p] += drho_dt_sink*dV;
                        }
                    }
                }

                //cooling
                if (sim_COOLING(theSim)==1){
                    if (sim_InitialDataType(theSim)==SHEAR){
                        double cooling_term = (Pp-.01)/(sim_GAMMALAW(theSim)-1.)*dV/1.e-2;
                        if (cooling_term<0.0) cooling_term = 0.0;
                        c->cons[TAU] -= cooling_term*dt;
                        c->Cool = cooling_term/dV;
                    }else{
                        double numfac = 1.0;
                        double alpha = sim_EXPLICIT_VISCOSITY(theSim);
                        double cooling_cap = c->cons[TAU]*numfac;
                        double cooling_term = 9./4.*alpha*pow(sim_PoRho_r1(theSim),-3)/rho*pow(Pp/rho,4.)*dt*dV;
                        if (cooling_term>cooling_cap) cooling_term = cooling_cap;
                        c->cons[TAU] -= cooling_term;
                        c->Cool = cooling_term/dt/dV;
                    }
                }
            }
        }
    }

    double Mdot_reduce[2];
    MPI_Allreduce( Mdot_temp,Mdot_reduce , 2, MPI_DOUBLE, MPI_SUM, sim_comm);
    int p;
    for( p=0 ; p<sim_NumGravMass(theSim); ++p ){
        gravMass_set_Mdot(theGravMasses,Mdot_reduce[p],p);
    }
}

void cell_add_visc_src( struct Cell *** theCells ,struct Sim * theSim, struct GravMass * theGravMasses, double dt ){
    double Mtotal = 1.0;
    double sep = 1.0;
    double M0 = gravMass_M(theGravMasses,0);
    double M1 = gravMass_M(theGravMasses,1);

    double r0 = gravMass_r(theGravMasses,0);
    double r1 = gravMass_r(theGravMasses,1);

    double a = r0 + r1;

    double phi_bh0 = gravMass_phi(theGravMasses,0);
    double phi_bh1 = gravMass_phi(theGravMasses,1);

    double xbh0 = r0*cos(phi_bh0);
    double ybh0 = r0*sin(phi_bh0);

    double xbh1 = r1*cos(phi_bh1);
    double ybh1 = r1*sin(phi_bh1);


    int i,j,k;
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
                struct Cell * c = &(theCells[k][i][j]);
                double dphi = c->dphi;
                double rho = c->prim[RHO];
                double P = c->prim[PPP];
                double r   = .5*(rp+rm);
                double vr  = c->prim[URR];
                double dz = zp-zm;
                double dV = dphi*.5*(rp*rp-rm*rm)*dz;

                double phi = c->tiph - 0.5*c->dphi;

                double xpos = r*cos(phi);
                double ypos = r*sin(phi);

                double dist_bh0 = sqrt((xpos-xbh0)*(xpos-xbh0)+(ypos-ybh0)*(ypos-ybh0));
                double dist_bh1 = sqrt((xpos-xbh1)*(xpos-xbh1)+(ypos-ybh1)*(ypos-ybh1));



                double alpha = sim_EXPLICIT_VISCOSITY(theSim);
                double rdrOm_a = sim_rdrOm_a(theSim,r,a);
                double eps = sim_G_EPS(theSim);
                double Sigma_nu;
                if (sim_VISC_CONST(theSim)==1){
                    Sigma_nu = rho*alpha;
                }else{
                    Sigma_nu = alpha*P/sqrt(pow(dist_bh0*dist_bh0+eps*eps,-1.5)*M0+pow(dist_bh1*dist_bh1+eps*eps,-1.5)*M1);
                }
                c->cons[TAU] += dt*dV*Sigma_nu*rdrOm_a*(rdrOm_a + r*c->grad[UPP] + c->gradp[URR] );
            }
        }
    }
}

void cell_add_visc_src_old( struct Cell *** theCells ,struct Sim * theSim, double dt ){
    int i,j,k;
    for( k=0 ; k<sim_N(theSim,Z_DIR) ; ++k ){
        double zm = sim_FacePos(theSim,k-1,Z_DIR);
        double zp = sim_FacePos(theSim,k,Z_DIR);
        for( i=0 ; i<sim_N(theSim,R_DIR) ; ++i ){
            double rm = sim_FacePos(theSim,i-1,R_DIR);
            double rp = sim_FacePos(theSim,i,R_DIR);
            for( j=0 ; j<sim_N_p(theSim,i) ; ++j ){
                struct Cell * c = &(theCells[k][i][j]);
                double dphi = c->dphi;
                double rho = c->prim[RHO];
                double r   = .5*(rp+rm);
                double vr  = c->prim[URR];
                double dz = zp-zm;
                double dV = dphi*.5*(rp*rp-rm*rm)*dz;
                double nu = sim_EXPLICIT_VISCOSITY(theSim);
                c->cons[SRR] += -dt*dV*nu*rho*vr/r/r;
            }
        }
    }
}


