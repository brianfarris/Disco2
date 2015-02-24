#include <math.h>
#include "../Headers/Metric.h"
#include "../Headers/Sim.h"

//Kerr metric in cylindrical Kerr-Schild coordinates.
//
// z_cyl = R_ks cos(theta_ks)
// r_cyl = sqrt(R_ks^2 - z_cyl^2)

#define SIG2 (R*R + A*A*z*z/(R*R))
#define DEL (R*R - 2*M*R + A*A)
#define AAA ((R*R+A*A)*(R*R+A*A) - A*A*(R*R-2*M*R+A*A)*r*r/(R*R))
#define BBB (2.0*M*R*R*R/(R*R*R*R + A*A*z*z))
#define DISIG2DR (-2*r*(R*R*R*R-A*A*z*z) / ((R*R*R*R+A*A*z*z)*(R*R*R*R+A*A*z*z)))
#define DISIG2DZ (-2*z*(R*R*R*R+A*A*r*r) / ((R*R*R*R+A*A*z*z)*(R*R*R*R+A*A*z*z)))
#define DDELDR (2.0*r - 2.0*M*r/R)
#define DDELDZ (2.0*z - 2.0*M*z/R)
#define DAAADR (2.0*r*((2.0*R*R+A*A)*R*R*R*R + A*A*(M*R*(R*R+z*z)-A*A*z*z))/(R*R*R*R))
#define DAAADZ (2.0*z*(2.0*(R*R+A*A)*R*R*R*R - A*A*r*r*(M*R-A*A))/(R*R*R*R))

double metric_g_dd_exact_kerr_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    double g;
    
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double R = sqrt(r*r+z*z);

    if(mu == 0 && nu == 0)
        g =  -1.0 + 2*M*R/SIG2;

    else if(mu == 1 && nu == 1)
        g = 1.0 + 2*M*r*r/(R*SIG2) + A*A*z*z*z*z/(R*R*R*R*R*R);

    else if(mu == 2 && nu == 2)
        g = AAA*r*r / (R*R*SIG2);

    else if(mu == 3 && nu == 3)
        g = 1.0 + 2*M*z*z/(R*SIG2) + A*A*r*r*z*z/(R*R*R*R*R*R);

    else if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
        g = 2*M*r/SIG2;

    else if((mu == 0 && nu == 2) || (mu == 2 && nu == 0))
        g = -2*A*M*r*r / (R*SIG2);

    else if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
        g = 2*M*z/SIG2;

    else if((mu == 1 && nu == 2) || (mu == 2 && nu == 1))
        g = -A * (1.0+2*M*R/SIG2) * r*r*r / (R*R*R);

    else if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        g = 2*M*r*z/(R*SIG2) - A*A*r*z*z*z/(R*R*R*R*R*R);

    else if((mu == 2 && nu == 3) || (mu == 3 && nu == 2))
        g = -A * (1.0+2*M*R/SIG2) * r*r*z / (R*R*R);

    else
        g = 0.0;

    return g;
}

double metric_g_uu_exact_kerr_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    double g;
    
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double R = sqrt(r*r+z*z);

    if(mu == 0 && nu == 0)
        g = -1.0 - 2*M*R/SIG2;

    else if(mu == 1 && nu == 1)
        g = (R*R*z*z + r*r*DEL)/(R*R*SIG2);

    else if(mu == 2 && nu == 2)
        g = R*R / (r*r*SIG2);

    else if(mu == 3 && nu == 3)
        g = 1.0 - 2*M*z*z/(R*SIG2);

    else if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
        g = 2*M*r/SIG2;

    else if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
        g = 2*M*z/SIG2;

    else if((mu == 1 && nu == 2) || (mu == 2 && nu == 1))
        g = A*r/(R*SIG2);

    else if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        g = r*z * (A*A-2*M*R) / (R*R*SIG2);

    else if((mu == 2 && nu == 3) || (mu == 3 && nu == 2))
        g = A*z/(R*SIG2);

    else
        g = 0.0;

    return g;
}

double metric_dg_dd_exact_kerr_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    double dg;
    
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double R = sqrt(r*r+z*z);

    if(k == 1)
    {
        if(mu == 0 && nu == 0)
            dg =  2*M*(r/(R*SIG2) + R*DISIG2DR);

        else if(mu == 1 && nu == 1)
            dg = 2*M*((2*r*R-r*r*r/R)/(R*R*SIG2) + r*r/R*DISIG2DR)
                - 6.0*A*A*pow(z,4)*r/pow(R,8);

        else if(mu == 2 && nu == 2)
            dg = (DAAADR*r*r+2*r*AAA) / (R*R*SIG2)
                    - AAA*r*r*(4*r*R*R) / (R*R*R*R*SIG2*SIG2);

        else if(mu == 3 && nu == 3)
            dg = 2*M*z*z*(-r/(R*R*R*SIG2) + DISIG2DR/R)
                 + A*A*z*z*(2*r/pow(R,6)-6*r*r*r/pow(R,8));

        else if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
            dg = 2*M*(1.0/SIG2 + r*DISIG2DR);

        else if((mu == 0 && nu == 2) || (mu == 2 && nu == 0))
            dg = -2*M*A*((2*r*R-r*r*r/R)/(R*R*SIG2) + r*r/R*DISIG2DR);

        else if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
            dg = 2*M*z*DISIG2DR;

        else if((mu == 1 && nu == 2) || (mu == 2 && nu == 1))
            dg = -A * (3*r*r/(R*R*R) - 3*pow(r/R,4)/R + 2*M*(3*r*r/(R*R*SIG2)
                    - 2*pow(r/R,4)/SIG2 + r*r*r/(R*R)*DISIG2DR));

        else if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
            dg = 2*M*z*(z*z/(R*R*R*SIG2) + r/R*DISIG2DR)
                     - A*A*z*z*z*(R*R - 6*r*r)/pow(R,8);

        else if((mu == 2 && nu == 3) || (mu == 3 && nu == 2))
            dg = -A*z * (2*r/(R*R*R) - 3*pow(r/R,3)/(R*R) + 2*M*(2*r/(R*R*SIG2)
                    - 2*pow(r/R,3)/(R*SIG2) + r*r/(R*R)*DISIG2DR));

        else
            dg = 0.0;
    }
    else if(k == 3)
    {
        if(mu == 0 && nu == 0)
            dg =  2*M*(z/(R*SIG2) + R*DISIG2DZ);

        else if(mu == 1 && nu == 1)
            dg = 2*M*r*r*(-z/(R*R*R*SIG2) + DISIG2DZ/R)
                + A*A*(4*z*z*z*R*R - 6*pow(z,5))/pow(R,8);

        else if(mu == 2 && nu == 2)
            dg = (DAAADZ*r*r) / (R*R*SIG2)
                - AAA*r*r*(4*z*R*R+2*A*A*z) / (R*R*R*R*SIG2*SIG2);

        else if(mu == 3 && nu == 3)
            dg = 2*M*((z*(2*R*R-z*z))/(R*R*R*SIG2) + z*z*DISIG2DZ/R)
                 + A*A*r*r*(2*z*R*R-6*z*z*z)/pow(R,8);

        else if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
            dg = 2*M*r*DISIG2DZ;

        else if((mu == 0 && nu == 2) || (mu == 2 && nu == 0))
            dg = -2*M*A*r*r*(-z/(R*R*R*SIG2) + DISIG2DZ/R);

        else if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
            dg = 2*M*(1.0/SIG2 + z*DISIG2DZ);

        else if((mu == 1 && nu == 2) || (mu == 2 && nu == 1))
            dg = -A*r*r*r * (-3*z/pow(R,5) + 2*M*(-2*z/(pow(R,4)*SIG2)
                     + DISIG2DZ/(R*R)));

        else if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
            dg = 2*M*r*(r*r/(R*R*R*SIG2) + z/R*DISIG2DZ)
                     - A*A*r*z*z*(3*R*R - 6*z*z)/pow(R,8);

        else if((mu == 2 && nu == 3) || (mu == 3 && nu == 2))
            dg = -A*r*r * ((R*R-3*z*z)/pow(R,5) + 2*M*(z*DISIG2DZ/(R*R)
                    + (R*R-2*z*z)/(R*R*R*R*SIG2)));

        else
            dg = 0.0;
    }
    else
        dg = 0.0;

    return dg;
}

double metric_dg_uu_exact_kerr_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    double dg;
    
    double M = sim_GravM(theSim);
    double A = M*sim_GravA(theSim);
    double R = sqrt(r*r+z*z);

    if(k == 1)
    {
        if(mu == 0 && nu == 0)
            dg =  -2*M*(r/(R*SIG2) + R*DISIG2DR);

        else if(mu == 1 && nu == 1)
            dg = ((2*r*DEL+r*r*DDELDR+2*r*z*z) * R*R*SIG2
                     - (r*r*DEL+R*R*z*z) * (4*R*R*r)) / (R*R*R*R*SIG2*SIG2);

        else if(mu == 2 && nu == 2)
            dg = -2*z*z/(r*r*r*SIG2) + R*R*DISIG2DR/(r*r);

        else if(mu == 3 && nu == 3)
            dg = -2*M*z*z*(-r/(R*R*R*SIG2) + DISIG2DR/R);

        else if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
            dg = 2*M*(1.0/SIG2 + r*DISIG2DR);

        else if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
            dg = 2*M*z*DISIG2DR;

        else if((mu == 1 && nu == 2) || (mu == 2 && nu == 1))
            dg = A * (z*z/(R*R*R*SIG2) + r*DISIG2DR/R);

        else if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
            dg = z*((A*A-2*M*(R+r*r/R))*R*R*SIG2
                     - (r*A*A-2*M*R*r)*(4*R*R*r)) / (R*R*R*R*SIG2*SIG2);

        else if((mu == 2 && nu == 3) || (mu == 3 && nu == 2))
            dg = A*z * (-r/(R*R*R*SIG2) + DISIG2DR/R);

        else
            dg = 0.0;
    }
    else if(k == 3)
    {
        if(mu == 0 && nu == 0)
            dg =  -2*M*(z/(R*SIG2) + R*DISIG2DZ);

        else if(mu == 1 && nu == 1)
            dg = ((r*r*DDELDZ+2*r*r*z+4*z*z*z) * R*R*SIG2
                     - (r*r*DEL+R*R*z*z) * (4*R*R*z+2*A*A*z))
                      / (R*R*R*R*SIG2*SIG2);

        else if(mu == 2 && nu == 2)
            dg = 2*z/(r*r*SIG2) + R*R*DISIG2DZ/(r*r);

        else if(mu == 3 && nu == 3)
            dg = -2*M*((z*(2*R*R-z*z))/(R*R*R*SIG2) + z*z*DISIG2DZ/R);

        else if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
            dg = 2*M*r*DISIG2DZ;

        else if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
            dg = 2*M*(1.0/SIG2 + z*DISIG2DZ);

        else if((mu == 1 && nu == 2) || (mu == 2 && nu == 1))
            dg = A*r * (-z/(R*R*R*SIG2) + DISIG2DZ/R);

        else if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
            dg = r*((A*A-2*M*(R+z*z/R))*R*R*SIG2
                 - (z*A*A-2*M*R*z)*(4*R*R*z+2*A*A*z)) / (R*R*R*R*SIG2*SIG2);

        else if((mu == 2 && nu == 3) || (mu == 3 && nu == 2))
            dg = A * (r*r/(R*R*R*SIG2) + z*DISIG2DZ/R);

        else
            dg = 0.0;
    }
    else
        dg = 0.0;

    return dg;
}

void metric_killing_exact_kerr_ks(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 0;
}

double metric_horizon_exact_kerr_ks(struct Sim *theSim)
{
    double a = sim_GravA(theSim);
    double M = sim_GravM(theSim);
    return M * (1.0 + sqrt(1.0 - a*a));
}
