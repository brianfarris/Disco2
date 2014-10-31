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

    if(coolType == COOL_NONE)
    {
        eos_cool = &eos_cool_none;
    }
    else if(coolType == COOL_ISOTHERM)
    {
        eos_cool = &eos_cool_isotherm;
    }
    else if(sim_Background(theSim) == COOL_BB_ES)
    {
        eos_cool = &eos_cool_bb_es;
    }
    else if(sim_Background(theSim) == COOL_VISC)
    {
        eos_cool = &eos_cool_visc;
    }
}

