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

    if(eosType == EOS_GAMMALAW)
    {
        eos_cs2  = &eos_cs2_gammalaw;
        eos_eps  = &eos_eps_gammalaw;
        eos_ppp  = &eos_ppp_gammalaw;
        eos_dpppdrho = &eos_dpppdrho_gammalaw;
        eos_dpppdttt = &eos_dpppdttt_gammalaw;
        eos_depsdrho = &eos_depsdrho_gammalaw;
        eos_depsdttt = &eos_depsdttt_gammalaw;
    }
    else if(eosType == EOS_GASRAD)
    {
        eos_cs2  = &eos_cs2_gasrad;
        eos_eps  = &eos_eps_gasrad;
        eos_ppp  = &eos_ppp_gasrad;
        eos_dpppdrho = &eos_dpppdrho_gasrad;
        eos_dpppdttt = &eos_dpppdttt_gasrad;
        eos_depsdrho = &eos_depsdrho_gasrad;
        eos_depsdttt = &eos_depsdttt_gasrad;
    }
    else if(eosType == EOS_GASRADDEG)
    {
        eos_cs2  = &eos_cs2_gasraddeg;
        eos_eps  = &eos_eps_gasraddeg;
        eos_ppp  = &eos_ppp_gasraddeg;
        eos_dpppdrho = &eos_dpppdrho_gasraddeg;
        eos_dpppdttt = &eos_dpppdttt_gasraddeg;
        eos_depsdrho = &eos_depsdrho_gasraddeg;
        eos_depsdttt = &eos_depsdttt_gasraddeg;
    }
    else if(eosType == EOS_PWF)
    {
        eos_cs2  = &eos_cs2_pwf;
        eos_eps  = &eos_eps_pwf;
        eos_ppp  = &eos_ppp_pwf;
        eos_dpppdrho = &eos_dpppdrho_pwf;
        eos_dpppdttt = &eos_dpppdttt_pwf;
        eos_depsdrho = &eos_depsdrho_pwf;
        eos_depsdttt = &eos_depsdttt_pwf;
    }

    else
    {
        printf("ERROR: Unrecognized EOS.\n");
        exit(0);
    }

    if(coolType == COOL_NONE)
        eos_cool = &eos_cool_none;
    else if(coolType == COOL_ISOTHERM)
        eos_cool = &eos_cool_isotherm;
    else if(coolType == COOL_BB_ES)
        eos_cool = &eos_cool_bb_es;
    else if(coolType == COOL_BB_FF)
        eos_cool = &eos_cool_bb_ff;
    else if(coolType == COOL_NU)
        eos_cool = &eos_cool_neutrino;
    else
    {
        printf("ERROR: Unrecognized Cooling.\n");
        exit(0);
    }
}

