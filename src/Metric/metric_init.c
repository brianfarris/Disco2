#define METRIC_PRIVATE_DEFS
#include <stdio.h>
#include <stdlib.h>
#include "../Headers/header.h"
#include "../Headers/Cell.h"
#include "../Headers/Riemann.h"
#include "../Headers/Sim.h"
#include "../Headers/Metric.h"

void metric_init_background(struct Sim *theSim)
{
    if(sim_Background(theSim) == NEWTON)
    {
        cell_prim2cons = &cell_prim2cons_newt;
        cell_cons2prim = &cell_cons2prim_newt;
        riemann_set_flux = &riemann_set_flux_newt;
        riemann_set_vel = &riemann_set_vel_newt;
    }
    else if(sim_Background(theSim) == GR)
    {
        cell_prim2cons = &cell_prim2cons_gr;
        cell_cons2prim = &cell_cons2prim_gr;
        riemann_set_flux = &riemann_set_flux_gr;
        riemann_set_vel = &riemann_set_vel_gr;
    }
}

void metric_init_metric(struct Sim *theSim)
{
    if(sim_Metric(theSim) == SR)
    {
        metric_g_dd_exact = &metric_g_dd_exact_sr;
        metric_g_uu_exact = &metric_g_uu_exact_sr;
        metric_dg_dd_exact = &metric_dg_dd_exact_sr;
        metric_dg_uu_exact = &metric_dg_uu_exact_sr;
        metric_killing_exact = &metric_killing_exact_sr;
    }
    else
    {
        printf("ERROR: Trying to initialize with an unknown metric: %d\n", sim_Metric(theSim));
        exit(0);
    }
}
