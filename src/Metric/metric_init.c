#define METRIC_PRIVATE_DEFS
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

    }

}
