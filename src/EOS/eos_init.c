#define EOS_PRIVATE_DEFS
#include <stdio.h>
#include <stdlib.h>
#include "../Headers/header.h"
#include "../Headers/EOS.h"
#include "../Headers/Sim.h"
#include "../Headers/Metric.h"

void eos_init(struct Sim *theSim)
{
    int coolType = sim_CoolingType(theSim);
    int eosType = sim_EOSType(theSim);

    if(eosType == EOS_IDEAL_NEWT)
    {
        eos_cs2_prim  = &eos_cs2_prim_ideal_newt;
        eos_eps_prim  = &eos_eps_prim_ideal_newt;
        eos_temp_prim = &eos_temp_prim_ideal_newt;
    }
    else if(eosType == EOS_IDEAL_GR)
    {
        eos_cs2_prim  = &eos_cs2_prim_ideal_gr;
        eos_eps_prim  = &eos_eps_prim_ideal_gr;
        eos_temp_prim = &eos_temp_prim_ideal_gr;
    }

    if(coolType == COOL_NONE)
        eos_cool = &eos_cool_none;
    else if(coolType == COOL_ISOTHERM)
        eos_cool = &eos_cool_isotherm;
    else if(sim_Background(theSim) == COOL_BB_ES)
        eos_cool = &eos_cool_bb_es;
    else if(sim_Background(theSim) == COOL_BB_FF)
        eos_cool = &eos_cool_bb_ff;
    else if(sim_Background(theSim) == COOL_NU)
        eos_cool = &eos_cool_neutrino;
}

