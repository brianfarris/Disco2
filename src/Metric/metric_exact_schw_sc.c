#include <math.h>

//Schwarzschild metric in cylindrical Scwarzschild coordinates.

double RG = 1.0;

double metric_g_dd_exact_schw_sc(int mu, int nu, double t, double r, double p, double z)
{
    if(mu == 0 && nu == 0)
        return -1.0+RG/sqrt(r*r+z*z);
    if(mu == 1 && nu == 1)
    {   
        double R = sqrt(r*r+z*z);
        return r*r/(R*(R-RG)) + z*z/(R*R);
    }
    if(mu == 2 && nu == 2)
        return r*r;
    if(mu == 3 && nu == 3)
    {   
        double R = sqrt(r*r+z*z);
        return z*z/(R*(R-RG)) + r*r/(R*R);
    }
    if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
    {
        double R = sqrt(r*r+z*z);
        return RG*r*z/(R*R*(R-RG));
    }
    return 0.0;
}

double metric_g_uu_exact_schw_sc(int mu, int nu, double t, double r, double p, double z)
{
    if(mu == 0 && nu == 0)
    {
        double R = sqrt(r*r+z*z);
        return -R/(R-RG);
    }
    if(mu == 1 && nu == 1)
    {   
        double R = sqrt(r*r+z*z);
        return (r*r*(R-RG) + z*z*R)/(R*R*R);
    }
    if(mu == 2 && nu == 2)
        return 1.0/(r*r);
    if(mu == 3 && nu == 3)
    {   
        double R = sqrt(r*r+z*z);
        return (z*z*(R-RG) + r*r*R)/(R*R*R);
    }
    if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
    {
        double R = sqrt(r*r+z*z);
        return -RG*r*z/(R*R*R);
    }

    return 0.0;
}

double metric_dg_dd_exact_schw_sc(int k, int mu, int nu, double t, double r, double p, double z)
{
    if(k == 1)
    {
        if(mu == 0 && nu == 0)
        {
            double R = sqrt(r*r+z*z);
            return -r*RG/(R*R*R);
        }
        if(mu == 1 && nu == 1)
        {   
            double R = sqrt(r*r+z*z);
            return -r*RG*(r*r*R - 2.0*z*z*(R-RG)) / (R*R*R*R*(R-RG)*(R-RG));
        }
        if(mu == 2 && nu == 2)
            return 2.0*r;
        if(mu == 3 && nu == 3)
        {   
            double R = sqrt(r*r+z*z);
            return -r*z*z*RG*(3.0*R-2*RG) / (R*R*R*R*(R-RG)*(R-RG));
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double R = sqrt(r*r+z*z);
            return z*RG*(R*R*(R-RG)-r*r*(3.0*R-2.0*RG)) / (R*R*R*R*(R-RG)*(R-RG));
        }
    }
    
    else if(k == 3)
    {
        if(mu == 0 && nu == 0)
        {
            double R = sqrt(r*r+z*z);
            return -z*RG/(R*R*R);
        }
        if(mu == 1 && nu == 1)
        {   
            double R = sqrt(r*r+z*z);
            return -z*r*r*RG*(3.0*R-2.0*RG) / (R*R*R*R*(R-RG)*(R-RG));
        }
        if(mu == 3 && nu == 3)
        {   
            double R = sqrt(r*r+z*z);
            return -z*RG*(R*z*z - 2.0*r*r*(R-RG)) / (R*R*R*R*(R-RG)*(R-RG));
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double R = sqrt(r*r+z*z);
            return r*RG*(R*R*(R-RG)-z*z*(3.0*R-2.0*RG)) / (R*R*R*R*(R-RG)*(R-RG));
        }
    }
    
    return 0.0;
}

double metric_dg_uu_exact_schw_sc(int k, int mu, int nu, double t, double r, double p, double z)
{
    if(k == 1)
    {
        if(mu == 0 && nu == 0)
        {
            double R = sqrt(r*r+z*z);
            return r*RG / (R*(R-RG)*(R-RG));
        }
        if(mu == 1 && nu == 1)
        {   
            double R = sqrt(r*r+z*z);
            return r*RG*(R*R-3.0*z*z) / (R*R*R*R*R);
        }
        if(mu == 2 && nu == 2)
            return -2.0/(r*r*r);
        if(mu == 3 && nu == 3)
        {   
            double R = sqrt(r*r+z*z);
            return 3.0*r*z*z*RG / (R*R*R*R*R);
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double R = sqrt(r*r+z*z);
            return -z*RG*(R*R-3.0*r*r) / (R*R*R*R*R);
        }
    }
    
    else if(k == 3)
    {
        if(mu == 0 && nu == 0)
        {
            double R = sqrt(r*r+z*z);
            return z*RG / (R*(R-RG)*(R-RG));
        }
        if(mu == 1 && nu == 1)
        {   
            double R = sqrt(r*r+z*z);
            return 3.0*r*r*z*RG / (R*R*R*R*R);
        }
        if(mu == 3 && nu == 3)
        {   
            double R = sqrt(r*r+z*z);
            return z*RG*(R*R-3.0*r*r) / (R*R*R*R*R);
        }
        if((mu == 1 && nu == 3) || (mu == 3 && nu == 1))
        {
            double R = sqrt(r*r+z*z);
            return -r*RG*(R*R-3.0*z*z) / (R*R*R*R*R);
        }
    }
    return 0.0;
}

void metric_killing_exact_schw_sc(int *k)
{
    k[0] = 1;
    k[1] = 0;
    k[2] = 1;
    k[3] = 0;
}
