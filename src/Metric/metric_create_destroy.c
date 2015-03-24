#define METRIC_PRIVATE_DEFS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../Headers/header.h"
#include "../Headers/Metric.h"

struct Metric* metric_create_exact(double t, double r, double p, double z, 
                                    struct Sim *theSim)
{
    int i,j;
    struct Metric *g = malloc(sizeof(struct Metric));
   
    g->x[0] = t;
    g->x[1] = r;
    g->x[2] = p;
    g->x[3] = z;
    for(i=0; i<4; i++)
        for(j=i; j<4; j++)
        {
            g->g_dd[j+4*i-i*(i+1)/2] = metric_g_dd_exact(i,j,t,r,p,z,theSim);
            g->g_uu[j+4*i-i*(i+1)/2] = metric_g_uu_exact(i,j,t,r,p,z,theSim);
        }
    metric_killing_exact(g->killing);
    for(i=0; i<4; i++)
        g->da[i] = 0.0;
    for(i=0; i<12; i++)
        g->db[i] = 0.0;
    g->dg_dd = NULL;
    g->dg_uu = NULL;
    g->length_dg = 0;
    g->theSim = theSim;

    return g;
}

void metric_create_der_exact(struct Metric *g)
{
    if(g->length_dg == 0)
    {
        int mu,nu,la,k;
        for(mu=0; mu<4; mu++)
            if(!g->killing[mu])
                (g->length_dg)++;
        if(g->length_dg == 0)
            return;
        g->dg_dd = malloc(10*(g->length_dg)*sizeof(double));
        g->dg_uu = malloc(10*(g->length_dg)*sizeof(double));
        
        k = 0;
        for(mu=0; mu<4; mu++)
            if(!g->killing[mu])
            {
                for(nu=0; nu<4; nu++)
                    for(la=nu; la<4; la++)
                    {
                        g->dg_dd[k*10 + la+4*nu-nu*(nu+1)/2] = metric_dg_dd_exact(mu,nu,la, g->x[0],g->x[1],g->x[2],g->x[3],g->theSim);
                        g->dg_uu[k*10 + la+4*nu-nu*(nu+1)/2] = metric_dg_uu_exact(mu,nu,la, g->x[0],g->x[1],g->x[2],g->x[3],g->theSim);
                    }
                
                double dg00uu = metric_dg_uu_exact(mu, 0, 0, g->x[0], g->x[1],
                                            g->x[2],g->x[3],g->theSim);
                double a = 1.0/sqrt(-g->g_uu[0]);
                double da = 0.5*a*a*a*dg00uu;
                double b[3], db[3];
                for(nu=0; nu<3; nu++)
                {
                    b[nu] = a*a * g->g_uu[nu+1];
                    double dg0iuu = metric_dg_uu_exact(mu,0,nu+1,g->x[0],
                                        g->x[1],g->x[2],g->x[3],g->theSim);
                    db[nu] = (a*a*a*a*dg0iuu + 2*a*b[nu]*da) / (a*a);
                }

                g->da[mu] = da;
                for(nu=0; nu<3; nu++)
                    g->db[3*mu+nu] = db[nu];

                k++;
            }
    }
}

struct Metric* metric_create_adm(double t, double r, double p, double z, 
                                    struct Sim *theSim)
{
    int i,j;
    struct Metric *g = malloc(sizeof(struct Metric));
   
    g->x[0] = t;
    g->x[1] = r;
    g->x[2] = p;
    g->x[3] = z;
    double a, b[3], gam[9], igam[9], b2, bd[3], ia2;

    a = metric_lapse_adm(t, r, p, z, theSim);
    for(i=0; i<3; i++)
    {
        b[i] = metric_shift_adm(i+1, t, r, p, z, theSim);
        gam[3*i+i] = metric_spatial_adm(i+1, i+1, t, r, p, z, theSim);
        igam[3*i+i] = metric_ispatial_adm(i+1, i+1, t, r, p, z, theSim);

        for(j=i+1; j<3; j++)
        {
            gam[3*i+j] = metric_spatial_adm(i+1, j+1, t, r, p, z, theSim);
            gam[3*j+i] = gam[3*i+j];
            igam[3*i+j] = metric_ispatial_adm(i+1, j+1, t, r, p, z, theSim);
            igam[3*j+i] = igam[3*i+j];
        }
    }

    ia2 = 1.0/(a*a);
    for(i=0; i<3; i++)
    {
        bd[i] = 0.0;
        for(j=0; j<3; j++)
            bd[i] += b[j]*gam[3*i+j];
    }
    b2 = b[0]*bd[0]+b[1]*bd[1]+b[2]*bd[2];

    g->g_dd[0] = -a*a+b2;
    g->g_dd[1] = bd[0];
    g->g_dd[2] = bd[1];
    g->g_dd[3] = bd[2];
    g->g_dd[4] = gam[0];
    g->g_dd[5] = gam[1];
    g->g_dd[6] = gam[2];
    g->g_dd[7] = gam[4];
    g->g_dd[8] = gam[5];
    g->g_dd[9] = gam[8];

    g->g_uu[0] = -ia2;
    g->g_uu[1] = b[0]*ia2;
    g->g_uu[2] = b[1]*ia2;
    g->g_uu[3] = b[2]*ia2;
    g->g_uu[4] = igam[0] - b[0]*b[0]*ia2;
    g->g_uu[5] = igam[1] - b[0]*b[1]*ia2;
    g->g_uu[6] = igam[2] - b[0]*b[2]*ia2;
    g->g_uu[7] = igam[4] - b[1]*b[1]*ia2;
    g->g_uu[8] = igam[5] - b[1]*b[2]*ia2;
    g->g_uu[9] = igam[8] - b[2]*b[2]*ia2;

    for(i=0; i<4; i++)
        g->da[i] = 0.0;
    for(i=0; i<12; i++)
        g->db[i] = 0.0;
    metric_killing_exact(g->killing);
    metric_killing_boost(g->killing);
    g->dg_dd = NULL;
    g->dg_uu = NULL;
    g->length_dg = 0;
    g->theSim = theSim;

    return g;
}

void metric_create_der_adm(struct Metric *g)
{
    if(g->length_dg == 0)
    {
        int mu,k,i,j;
        for(mu=0; mu<4; mu++)
            if(!g->killing[mu])
                (g->length_dg)++;
        if(g->length_dg == 0)
            return;
        g->dg_dd = malloc(10*(g->length_dg)*sizeof(double));

        double a, b[3], bd[3], gam[9], ia2, dg00;
        double da, db[3], dbd[3], dgam[9];
        double t, r, p, z;

        if(PRINTTOOMUCH)
        {
            for(i=0; i<10; i++)
                printf("%.6e ", g->g_dd[i]);
            printf("\n");
            for(i=0; i<10; i++)
                printf("%.6e ", g->g_uu[i]);
            printf("\n");
        }

        ia2 = -g->g_uu[0];
        a = 1.0/sqrt(ia2);
        b[0] = g->g_uu[1] * a*a;
        b[1] = g->g_uu[2] * a*a;
        b[2] = g->g_uu[3] * a*a;
        gam[0] = g->g_dd[4];
        gam[1] = g->g_dd[5];
        gam[2] = g->g_dd[6];
        gam[3] = g->g_dd[5];
        gam[4] = g->g_dd[7];
        gam[5] = g->g_dd[8];
        gam[6] = g->g_dd[6];
        gam[7] = g->g_dd[8];
        gam[8] = g->g_dd[9];
        t = g->x[0];
        r = g->x[1];
        p = g->x[2];
        z = g->x[3];
        
        k = 0;
        for(mu=0; mu<4; mu++)
            if(!g->killing[mu])
            {
                da = metric_dlapse_adm(mu, t, r, p, z, g->theSim);
                if(PRINTTOOMUCH)
                    printf("Making derivatives for x[%d], da = %.6e\n", mu, da);
                for(i=0; i<3; i++)
                {
                    db[i] = metric_dshift_adm(mu, i+1, t, r, p, z, g->theSim);
                    dgam[4*i] = metric_dspatial_adm(mu,i+1,i+1,t,r,p,z,
                                                    g->theSim);
                    for(j=i; j<3; j++)
                    {
                        dgam[3*i+j] = metric_dspatial_adm(mu,i+1,j+1,t,r,p,z,
                                                            g->theSim);
                        dgam[3*j+i] = dgam[3*i+j];
                    }
                }
                for(i=0; i<3; i++)
                {
                    dbd[i] = 0.0;
                    for(j=0; j<3; j++)
                        dbd[i] += gam[3*i+j]*db[j] + dgam[3*i+j]*b[j];
                }
                dg00 = -2*a*da;
                for(i=0; i<3; i++)
                    dg00 += bd[i]*db[i] + dbd[i]*b[i];

                g->da[mu] = da;
                for(i=0; i<3; i++)
                    g->db[3*mu+i] = db[i];
                g->dg_dd[10*k] = dg00;
                for(i=0; i<3; i++)
                    g->dg_dd[10*k+i+1] = dbd[i];
                g->dg_dd[10*k+4] = dgam[0];
                g->dg_dd[10*k+5] = dgam[1];
                g->dg_dd[10*k+6] = dgam[2];
                g->dg_dd[10*k+7] = dgam[4];
                g->dg_dd[10*k+8] = dgam[5];
                g->dg_dd[10*k+9] = dgam[8];
            
                if(PRINTTOOMUCH)
                {
                    printf("%.6e\n", a);
                    printf("%.6e\n", da);
                    printf("%.6e %.6e %.6e\n", b[0], b[1], b[2]);
                    printf("%.6e %.6e %.6e\n", db[0], db[1], db[2]);
                    printf("%.6e %.6e %.6e\n", dbd[0], dbd[1], dbd[2]);
                    for(i=0; i<9; i++)
                        printf("%.6e ", gam[i]);
                    printf("\n");
                    for(i=0; i<9; i++)
                        printf("%.6e ",dgam[i]);
                    printf("\n");
                    for(i=0; i<10; i++)
                        printf("%.6e ", g->dg_dd[10*k+i]);
                    printf("\n");
                }

                k++;
            }
    }
}

void metric_destroy(struct Metric *g)
{
    if (g->dg_dd != NULL)
        free(g->dg_dd);
    if (g->dg_uu != NULL)
        free(g->dg_uu);
    free(g);
}
