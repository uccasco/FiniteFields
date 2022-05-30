# FiniteFields

Implementation of finite field methods in LAMMPS. The tagging systems
is as follows:

tag = 03Mar2020; should be compatible with stable 03Mar2020 version of LAMMPS
tag = 29Oct2020; should be compatible with stable 29Oct2020 version of LAMMPS

Source files in the master branch should be compatible with the
29Sept2021 stable version of LAMMPS. Changes between 29Oct2020 and
29Sept2021 are minimal, and the "fix_dfield" implementation gives
identical trajectories between the two trajectories for short
simulations of pure water, an electrolyte solution, and a NaCl (111)
surface. <b>The "fix_dfield_tip4p" and "fix_efield_tip4p" have not been
tested between versions.</b>


This is an implementation of the constant D ensemble in LAMMPS, based
on the work of Stengel, Spaldin and Vanderbilt [Nature Phys. 5,
304–308 (2009)], which was later adopted for finite temperature MD
simulation by Zhang and Sprik [Phys. Rev. B 93, 144201 (2016)]. In
addition to the original citations, if you use this code, it would be
greatly appreciated if you cite:

> S. J. Cox and M. Sprik, <i>J. Chem. Phys.</i> <b>151</b>, 064506 (2019); doi: 10.1063/1.5099207

The full article can be found
[here](https://aip.scitation.org/doi/10.1063/1.5099207). If you use
the TIP4P routines (see below), please additionally cite:

> T. Sayer and S. J. Cox, <i> Phys. Chem. Chem. Phys. </i> <b>21</b>, 14546 (2019); doi: 10.1039/c9cp02193k

The full article can be found
[here](https://pubs.rsc.org/en/content/articlepdf/2019/cp/c9cp02193k). The
code is directly based upon the fix_efield routines that come as
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

## Use with TIP4P style potentials

To apply fields to TIP4P-style models requires appropriate
distribution of the forces from the virtual site to the O and H
atoms. The files "*_tip4p.cpp/.h" provide the capability for
this. Assuming O is atom type 1 and H is atom type 2, and that the
bond and angle IDs are both 1, an electric field would be applied with
the following command in the LAMMPS input file:

```
fix Eensemble fluid efield/tip4p 0.0 0.0 -0.031 1 2 1 1 0.1546 -1.112800
```

In the above, the virtual site is separated from the oxygen by 0.1546
Angstrom, and carries a charge -1.1128. The electric field along z is
-0.031 V/Ang.

To apply a displacement field is a bit more involved. Like above, the
command is altered to:

```
fix Densemble fluid dfield/tip4p NULL NULL -1.499 c_OmegaPx c_OmegaPy c_OmegaPz 1 2 1 1 0.1546 -1.112800
```

However, to calculate the polarization we now need to treat the water
separately from anything else in our system. For TIP4P/2005 this can
be done as follows:

```
# Store the initial atom positions
fix initconfig_wat       wat store/state 0 xu yu zu
fix initconfig_nonwat nonwat store/state 0 xu yu zu

compute displ_wat       wat displace/atom
compute displ_nonwat nonwat displace/atom

variable OmegaPxi_wat atom (c_displ_wat[1]+f_initconfig_wat[1])*q*0.73612446361690741440
variable OmegaPyi_wat atom (c_displ_wat[2]+f_initconfig_wat[2])*q*0.73612446361690741440
variable OmegaPzi_wat atom (c_displ_wat[3]+f_initconfig_wat[3])*q*0.73612446361690741440

variable OmegaPxi_nonwat atom (c_displ_nonwat[1]+f_initconfig_nonwat[1])*q
variable OmegaPyi_nonwat atom (c_displ_nonwat[2]+f_initconfig_nonwat[2])*q
variable OmegaPzi_nonwat atom (c_displ_nonwat[3]+f_initconfig_nonwat[3])*q

variable OmegaPxi atom v_OmegaPxi_wat+v_OmegaPxi_nonwat
variable OmegaPyi atom v_OmegaPyi_wat+v_OmegaPyi_nonwat
variable OmegaPzi atom v_OmegaPzi_wat+v_OmegaPzi_nonwat

compute OmegaPx all reduce sum v_OmegaPxi
compute OmegaPy all reduce sum v_OmegaPyi
compute OmegaPz all reduce sum v_OmegaPzi
```

where groups "wat" and "nonwat" have been defined previously in the
LAMMPS input file.
