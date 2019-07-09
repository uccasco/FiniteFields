# FiniteFields
Implementation of finite field methods in LAMMPS

This is an implementation of the constant D ensemble in LAMMPS, based
on the work of Stengel, Spaldin and Vanderbilt [Nature Phys. 5,
304–308 (2009)], which was later adopted for finite temperature MD
simulation by Zhang and Sprik [Phys. Rev. B 93, 144201 (2016)].

The code is directly based upon the fix_efield routines that come as
standard in LAMMPS. However, an important difference is that the
itinerant polarization must be supplied to the fix from the input
file. This can be achieved with the following commands:

```
# Store the initial atom positions
fix initconfig all store/state 0 xu yu zu

compute displ all displace/atom

variable OmegaPxi atom (c_displ[1]+f_initconfig[1])*q
variable OmegaPyi atom (c_displ[2]+f_initconfig[2])*q
variable OmegaPzi atom (c_displ[3]+f_initconfig[3])*q

compute OmegaPx all reduce sum v_OmegaPxi
compute OmegaPy all reduce sum v_OmegaPyi
compute OmegaPz all reduce sum v_OmegaPzi

fix Densemble all dfield 0.0 0.0 0.0 c_OmegaPx c_OmegaPy c_OmegaPz

```

In the above example, the D field has been set to zero in all three
Cartesian directions. If we wanted to use hybrid boundary conditions,
e.g. setting D<sub>x</sub> = 2.0 V/A but not setting D<sub>y</sub> or
D<sub>z</sub>, then we could use:

```
fix Densemble all dfield 2.0 NULL NULL c_OmegaPx c_OmegaPy c_OmegaPz
```

It has been tested with the Sept 2017 version of the LAMMPS
distribution, and with real units. Use of the code comes with no
warranty of any kind. As the code is a direct modification of existing
LAMMPS source code, it is distributed under the terms of the GNU
Public License (GPL). Please see the [LAMMPS
website](https://lammps.sandia.gov/doc/Intro_opensource.html) for
further details. To see where modifications have been made, comments
with the initials `// SJC` have been added.

Functionality to calculate the energy has been implemented, but it
will only work for real units. Initial tests seem to suggest it is
working, but this is still preliminary.

To use the constant D ensemble, it should be sufficient to simply copy
the relevant `.cpp` and `.h` files into the main LAMMPS src directory
and recompile.

