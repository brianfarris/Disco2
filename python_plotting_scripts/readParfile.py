import sys

parnames = ['InitialDataType', 
            'GravMassType', 
            'Background', 
            'Metric',
            'Frame', 
            'Restart', 
            'NumR', 
            'NumZ', 
            'ng', 
            'R_Min', 
            'R_Max', 
            'Z_Min', 
            'Z_Max',
            'NP_CONST', 
            'aspect',
            'NUM_C', 
            'NUM_N', 
            'Time_Max', 
            'Num_Checkpoints', 
            'Num_Diag_Dump',
            'Num_Diag_Measure',
            'NoInnerBC', 
            'BoundTypeRIn', 
            'BoundTypeROut', 
            'BoundTypeZIn', 
            'BoundTypeZBot',
            'ZeroPsiBndry',
            'BoundTypeSource',
            'BoundPar1',
            'BoundPar2',
            'BoundPar3',
            'BoundPar4',
            'Move_Cells',
            'RiemannSolver',
            'Adiabatic_Index',
            'CFL',
            'PLM',
            'Grav_2D',
            'G_EPS',
            'GravM',
            'GravA',
            'AlphaVisc',
            'BoostType',
            'BinA',
            'BinW',
            'BinM',
            'EOSType',
            'EOSPar1',
            'EOSPar2',
            'EOSPar3',
            'EOSPar4',
            'CoolingType',
            'CoolPar1',
            'CoolPar2',
            'CoolPar3',
            'CoolPar4',
            'InitPar0',
            'InitPar1',
            'InitPar2',
            'InitPar3',
            'InitPar4',
            'InitPar5',
            'InitPar6',
            'PHI_ORDER',
            'Rho_Floor',
            'Cs_Floor',
            'Cs_Cap',
            'Vel_Cal',
            'runtype',
            'DAMP_TIME',
            'RDAMP_INNER',
            'RDAMP_OUTER',
            'RLogScale',
            'ZLogScale',
            'HiResSigma',
            'HiResR0',
            'HiResFac']

def readParfile(filename):
    # Read a parameter file and load the contents into a dict.
    f = open(filename, "r")

    pars = dict()

    for line in f:
        
        words = line.split()
        if len(words) < 2:
            continue

        if words[0] in parnames:
            key = words[0]
            sval = words[1]
            try:
                val = int(sval)
            except ValueError:
                try:
                    val = float(sval)
                except ValueError:
                    val = None
            pars[key] = val

    f.close()

    return pars


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("I need a parfile to nom.\n")

    pars = readParfile(sys.argv[1])
    print(pars)


