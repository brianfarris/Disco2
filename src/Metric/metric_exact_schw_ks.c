#include "../Headers/Metric.h"
#include "../Headers/Sim.h"
#include <math.h>

//Schwarzschild metric in Kerr-Schild coordinates.

double metric_g_dd_exact_schw_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(mu == 0 && nu == 0)
    {
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return -1.0+2.0*M/R;
    }
    if((mu==0 && nu==1) || (mu==1 && nu==0))
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 2.0*M*r/(R*R);
    }
    if((mu==0 && nu==3) || (mu==3 && nu==0))
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 2.0*M*z/(R*R);
    }
    if(mu == 1 && nu == 1)
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 1.0 + 2.0*M*r*r/(R*R*R);
    }
    if(mu == 2 && nu == 2)
        return r*r;
    if(mu == 3 && nu == 3)
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 1.0 + 2.0*M*z*z/(R*R*R);
    }
    if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
    {
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 2.0*M*r*z/(R*R*R);
    }
    return 0.0;
}

double metric_g_uu_exact_schw_ks(int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(mu == 0 && nu == 0)
    {
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return -1.0-2.0*M/R;
    }
    if((mu==0 && nu==1) || (mu==1 && nu==0))
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 2.0*M*r/(R*R);
    }
    if((mu==0 && nu==3) || (mu==3 && nu==0))
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 2.0*M*z/(R*R);
    }
    if(mu == 1 && nu == 1)
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 1.0 - 2.0*M*r*r/(R*R*R);
    }
    if(mu == 2 && nu == 2)
        return 1.0/(r*r);
    if(mu == 3 && nu == 3)
    {   
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return 1.0 - 2.0*M*z*z/(R*R*R);
    }
    if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
    {
        double M = sim_GravM(theSim);
        double R = sqrt(r*r+z*z);
        return -2.0*M*r*z/(R*R*R);
    }
    return 0.0;
}

double metric_dg_dd_exact_schw_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(k==0 || k==2)
        return 0.0;

    if(k == 1)
    {
        if(mu == 0 && nu == 0)
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -2.0*M*r/(R*R*R);
        }
        if((mu==0 && nu==1) || (mu==1 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*(z*z-r*r)/(R*R*R*R);
        }
        if((mu==0 && nu==3) || (mu==3 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -4.0*M*r*z/(R*R*R*R);
        }
        if(mu == 1 && nu == 1)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*r*(2.0*z*z-r*r)/(R*R*R*R*R);
        }
        if(mu == 2 && nu == 2)
            return 2.0*r;
        if(mu == 3 && nu == 3)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -6.0*M*r*z*z/(R*R*R*R*R);
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*z*(z*z-2.0*r*r)/(R*R*R*R*R);
        }
        return 0.0;
    }
    if(k == 3)
    {
        if(mu == 0 && nu == 0)
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -2.0*M*z/(R*R*R);
        }
        if((mu==0 && nu==1) || (mu==1 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -4.0*M*r*z/(R*R*R*R);
        }
        if((mu==0 && nu==3) || (mu==3 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*(r*r-z*z)/(R*R*R*R);
        }
        if(mu == 1 && nu == 1)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -6.0*M*r*r*z/(R*R*R*R*R);
        }
        if(mu == 3 && nu == 3)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*z*(2.0*r*r-z*z)/(R*R*R*R*R);
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*r*(r*r-2.0*z*z)/(R*R*R*R*R);
        }
        return 0.0;
    }

    return 0.0;
}

double metric_dg_uu_exact_schw_ks(int k, int mu, int nu, double t, double r, double p, double z, struct Sim *theSim)
{
    if(k==0 || k==2)
        return 0.0;

    if(k == 1)
    {
        if(mu == 0 && nu == 0)
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*r/(R*R*R);
        }
        if((mu==0 && nu==1) || (mu==1 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*(z*z-r*r)/(R*R*R*R);
        }
        if((mu==0 && nu==3) || (mu==3 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -4.0*M*r*z/(R*R*R*R);
        }
        if(mu == 1 && nu == 1)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -2.0*M*r*(2.0*z*z-r*r)/(R*R*R*R*R);
        }
        if(mu == 2 && nu == 2)
            return -2.0/(r*r*r);
        if(mu == 3 && nu == 3)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 6.0*M*r*z*z/(R*R*R*R*R);
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -2.0*M*z*(z*z-2.0*r*r)/(R*R*R*R*R);
        }
        return 0.0;
    }
    if(k == 3)
    {
        if(mu == 0 && nu == 0)
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*z/(R*R*R);
        }
        if((mu==0 && nu==1) || (mu==1 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -4.0*M*r*z/(R*R*R*R);
        }
        if((mu==0 && nu==3) || (mu==3 && nu==0))
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 2.0*M*(r*r-z*z)/(R*R*R*R);
        }
        if(mu == 1 && nu == 1)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return 6.0*M*r*r*z/(R*R*R*R*R);
        }
        if(mu == 3 && nu == 3)
        {   
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -2.0*M*z*(2.0*r*r-z*z)/(R*R*R*R*R);
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double M = sim_GravM(theSim);
            double R = sqrt(r*r+z*z);
            return -2.0*M*r*(r*r-2.0*z*z)/(R*R*R*R*R);
        }
        return 0.0;
    }

    return 0.0;
}

void metric_killing_exact_schw_ks(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 0;
}
