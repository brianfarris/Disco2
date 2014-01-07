#include <math.h>

//Schwarzschild metric in Kerr-Schild coordinates.

double RG_KS = 1.0;

double metric_g_dd_exact_schw_ks(int mu, int nu, double t, double r, double p, double z)
{
    if(mu == 0 && nu == 0)
        return -1.0+RG_KS/sqrt(r*r+z*z);
    if(mu == 1 && nu == 1)
    {   
        double R = sqrt(r*r+z*z);
        return 1.0+RG_KS*r*r/(R*R*R);
    }
    if(mu == 2 && nu == 2)
        return r*r;
    if(mu == 3 && nu == 3)
    {   
        double R = sqrt(r*r+z*z);
        return 1.0+RG_KS*z*z/(R*R*R);
    }
    if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
        return RG_KS*r/(r*r+z*z);
    if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
        return RG_KS*z/(r*r+z*z);
    if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
    {
        double R = sqrt(r*r+z*z);
        return RG_KS*r*z/(R*R*R);
    }
    return 0.0;
}

double metric_g_uu_exact_schw_ks(int mu, int nu, double t, double r, double p, double z)
{
    if(mu == 0 && nu == 0)
        return -1.0+RG_KS/sqrt(r*r+z*z);
    if(mu == 1 && nu == 1)
    {   
        double R = sqrt(r*r+z*z);
        return 1.0+RG_KS*r*r/(R*R*R);
    }
    if(mu == 2 && nu == 2)
        return r*r;
    if(mu == 3 && nu == 3)
    {   
        double R = sqrt(r*r+z*z);
        return 1.0+RG_KS*z*z/(R*R*R);
    }
    if((mu == 0 && nu == 1) || (mu == 1 && nu == 0))
        return RG_KS*r/(r*r+z*z);
    if((mu == 0 && nu == 3) || (mu == 3 && nu == 0))
        return RG_KS*z/(r*r+z*z);
    if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
    {
        double R = sqrt(r*r+z*z);
        return RG_KS*r*z/(R*R*R);
    }
}

double metric_dg_dd_exact_schw_ks(int k, int mu, int nu, double t, double r, double p, double z)
{
    if(k!=1)
        return 0.0;
    if(mu==0 && nu==0)
        return -RG_KS/(r*r);
    if((mu==0 && nu==1) || (mu==1 && nu==0))
        return -2.0*RG_KS/(r*r);
    if(mu==1 && nu==1)
        return -RG_KS/(r*r);
    if(mu==2 && nu==2)
        return 2.0*r;
    return 0.0;
}

double metric_dg_uu_exact_schw_ks(int k, int mu, int nu, double t, double r, double p, double z)
{
    if(k==1 && mu==2 && nu==2)
        return -2.0/(r*r*r);
    return 0.0;
}

void metric_killing_exact_schw_ks(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 0;
}
