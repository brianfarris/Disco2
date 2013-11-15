#define METRIC_PRIVATE_DEFS
#include <stdio.h>
#include <stdlib.h>
#include "../Headers/Metric.h"

struct Metric* metric_create(double t, double r, double p, double z)
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
            g->g_dd[j+4*i-i*(i+1)/2] = metric_g_dd_exact(i,j,t,r,p,z);
            g->g_uu[j+4*i-i*(i+1)/2] = metric_g_uu_exact(i,j,t,r,p,z);
        }
    metric_killing_exact(g->killing);
    g->dg_dd = NULL;
    g->dg_uu = NULL;
    g->length_dg = 0;

    return g;
}

void metric_create_der(struct Metric *g)
{
    if(g->length_dg == 0)
    {
        int a,b,c,k;
        for(a=0; a<4; a++)
            if(!g->killing[a])
                (g->length_dg)++;
        if(g->length_dg == 0)
            return;
        g->dg_dd = malloc(10*(g->length_dg)*sizeof(double));
        g->dg_uu = malloc(10*(g->length_dg)*sizeof(double));
        
        k = 0;
        for(a=0; a<4; a++)
            if(!g->killing[a])
            {
                for(b=0; b<4; b++)
                    for(c=b; c<4; c++)
                    {
                        g->dg_dd[k*10 + c+4*b-b*(b+1)/2] = metric_dg_dd_exact(a,b,c, g->x[0],g->x[1],g->x[2],g->x[3]);
                        g->dg_uu[k*10 + c+4*b-b*(b+1)/2] = metric_dg_uu_exact(a,b,c, g->x[0],g->x[1],g->x[2],g->x[3]);
                    }
                k++;
            }
    }
}

void metric_destroy(struct Metric *g)
{
    if(g->length_dg > 0)
    {
        free(g->dg_dd);
        free(g->dg_uu);
    }
    free(g);
}
