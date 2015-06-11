## discopy 

This is a python module for the analysis, transformation, and plotting of Disco
data.  It contains standalone functions for reading checkpoint files and
parfiles for basic data analysis.  More advanced features use the `Grid` object
to transform data between grids, make grid-accurate plots, and save archives.

### Basic Use

First import the package.
```python
import discopy as dp
```

Read in checkpoint files with `readCheckpoint()`.  It returns a tuple of 1D 
`numpy` arrays containing the checkpoint data.  Duplicate rows (e.g. output from a parallel run) are removed.  Given positions are cell-centered values.
```python
dat = dg.readCheckpoint("checkpoint_0100.h5")
t = dat[0]      # Time
r = dat[1]      # r of cell centres
phi = dat[2]    # phi of cell centres
z = dat[3]      # z of cell centres
rho = dat[4]    # prim[RHO]
P = dat[5]      # prim[PPP]
vr = dat[6]     # prim[URR]
vp = dat[7]     # prim[UPP]
vz = dat[8]     # prim[UZZ]
dV = dat[9]     # Volume of each cell
q = dat[10]     # First passive scalar if present, np.zeros(r.shape) otherwise.
piph = dat[11]  # phi value for face in positive phi direction of each cell.
```

Read in parameter files with `readParfile()`.  It returns a dictionary
containing all parameters.  Keys are defined in `dp.parnames`, using the 
same strings the parfile uses.

```python
pars = dg.readParfile("vortex.par")
pars['NumR']   #Number of radial non-ghost zones.
```

### Advanced Use

`Grid` objects contain all the necessary grid data: face locations, numbers of
zones and ghost zones in each dimension, and primitive values for each cell.

Make a grid, from a parfile, or a parfile and checkpoint:
```python
import discopy as dp

pars = dg.readParfile("vortex.par")
checkpointName = "checkpoint_0100.h5"

g1 = dg.Grid(pars)
g2 = dg.Grid(pars, checkpointName)
```

Grid objects can be saved as "archives".  These files are similar to, but they 
contain all the information needed to construct the grid, as well as a copy
of the parameter file.  That is, a single archive file contains all the
information necessary to start a `disco` run.

Grids can also output checkpoint files, with boundary zones properly copied
to enable running on a parallel machine.

```python
g2.saveArchive("archive.h5")
g2.saveCheckpoint("input.h5", numProc=64)
```




