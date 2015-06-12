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

#### Remapping Checkpoints

The main reason `discopy` exists is to map `disco` output onto different grids 
suitable for restarting a calculation.  This could be used, for instance, to 
bootstrap convergence at high resolution.  Data from a run that has already
converged to a quasi steady steady state can be remapped to a higher
resolution so the initial transients do not have to be recalculated.

To remap a checkpoint into `disco` input you need a checkpoint file (that you
wish to remap), the parfile used to create that checkpoint, and the parfile
for the new run.

```python
import discopy as dp

#Load the relevant parfiles.
par1 = dp.readParfile("lores.par")
par2 = dp.readParfile("hires.par")

# Make Grid objects for each.
g1 = dp.Grid(par1, "checkpoint_0100.h5")
g2 = dp.Grid(par2)

# Do the remap
dp.remapGrid(g1, g2)

# Save the new grid, and archives for good measure
g1.saveArchive("archive_lores.h5")
g2.saveArchive("archive_hires.h5")
g2.saveCheckpoint("input.h5", numProc=64)
```

Saving the archive files is useful, I highly recommend it.  Since they contain 
a copy of the parfile, they can be used to completely restart a run if you
lose the other checkpoint data.  You can also load the archive in `discopy` 
and immediately save a checkpoint for an arbitrary number of processes.  This
is very useful if you want to run on a different number of cores!


