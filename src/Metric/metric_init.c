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
        cell_add_src = &cell_add_src_newt;
        cell_mindt = &cell_mindt_newt;
        riemann_set_flux = &riemann_set_flux_newt;
        riemann_set_vel = &riemann_set_vel_newt;
    }
    else if(sim_Background(theSim) == GR)
    {
        cell_prim2cons = &cell_prim2cons_gr;
        cell_cons2prim = &cell_cons2prim_gr;
        cell_add_src = &cell_add_src_gr;
        cell_mindt = &cell_mindt_gr;
        riemann_set_flux = &riemann_set_flux_gr;
        riemann_set_vel = &riemann_set_vel_gr;
    }
    else if(sim_Background(theSim) == GRVISC1)
    {
        cell_prim2cons = &cell_prim2cons_gr;
        cell_cons2prim = &cell_cons2prim_gr;
        cell_add_src = &cell_add_src_gr;
        cell_mindt = &cell_mindt_gr;
        riemann_set_flux = &riemann_set_flux_gr;
        riemann_set_vel = &riemann_set_vel_gr;
    }
    else if(sim_Background(theSim) == GRDISC)
    {
        cell_prim2cons = &cell_prim2cons_grdisc;
        cell_cons2prim = &cell_cons2prim_grdisc;
        cell_add_src = &cell_add_src_grdisc;
        cell_mindt = &cell_mindt_grdisc;
        riemann_set_flux = &riemann_set_flux_grdisc;
        riemann_set_vel = &riemann_set_vel_grdisc;
    }

    if(sim_Frame(theSim) == FR_EULER)
    {
        printf("Euler!\n");
        metric_frame_U_u = &metric_frame_U_u_euler;
        metric_frame_dU_du = &metric_frame_dU_du_euler;
    }
    else if(sim_Frame(theSim) == FR_KEP)
    {
        printf("Kepler!\n");
        metric_frame_U_u = &metric_frame_U_u_kep;
        metric_frame_dU_du = &metric_frame_dU_du_kep;
    }
    else if(sim_Frame(theSim) == FR_ACC)
    {
        printf("Accreting!\n");
        metric_frame_U_u = &metric_frame_U_u_acc;
        metric_frame_dU_du = &metric_frame_dU_du_acc;
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
    else if(sim_Metric(theSim) == SCHWARZSCHILD_SC)
    {
        metric_g_dd_exact = &metric_g_dd_exact_schw_sc;
        metric_g_uu_exact = &metric_g_uu_exact_schw_sc;
        metric_dg_dd_exact = &metric_dg_dd_exact_schw_sc;
        metric_dg_uu_exact = &metric_dg_uu_exact_schw_sc;
        metric_killing_exact = &metric_killing_exact_schw_sc;
    }
    else if(sim_Metric(theSim) == SCHWARZSCHILD_KS)
    {
        metric_g_dd_exact = &metric_g_dd_exact_schw_ks;
        metric_g_uu_exact = &metric_g_uu_exact_schw_ks;
        metric_dg_dd_exact = &metric_dg_dd_exact_schw_ks;
        metric_dg_uu_exact = &metric_dg_uu_exact_schw_ks;
        metric_killing_exact = &metric_killing_exact_schw_ks;
    }
    else if(sim_Metric(theSim) == SR_CART)
    {
        metric_g_dd_exact = &metric_g_dd_exact_sr_cart;
        metric_g_uu_exact = &metric_g_uu_exact_sr_cart;
        metric_dg_dd_exact = &metric_dg_dd_exact_sr_cart;
        metric_dg_uu_exact = &metric_dg_uu_exact_sr_cart;
        metric_killing_exact = &metric_killing_exact_sr_cart;
    }
    else
    {
        printf("ERROR: Trying to initialize with an unknown metric: %d\n", sim_Metric(theSim));
        exit(0);
    }
}
