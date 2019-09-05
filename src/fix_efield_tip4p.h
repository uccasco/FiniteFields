/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Modified for TIP4P capability (2017-2019): Stephen J. Cox
   (University of Cambridge)

   Copyright (2019) Stephen J. Cox
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(efield/tip4p,FixEfieldTip4p)

#else

#ifndef LMP_FIX_EFIELD_TIP4P_H
#define LMP_FIX_EFIELD_TIP4P_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEfieldTip4p : public Fix {
 public:
  FixEfieldTip4p(class LAMMPS *, int, char **);
  ~FixEfieldTip4p();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double memory_usage();
  double compute_scalar();
  double compute_vector(int);

 private:
  double ex,ey,ez;
  int varflag,iregion;
  char *xstr,*ystr,*zstr,*estr;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  int ilevel_respa;
  double qe2f;
  int qflag,muflag;

  int maxatom;
  double **efield_tip4p;

  int force_flag;
  double fsum[4],fsum_all[4];

  //SJC:
  int typeO, typeH, typeB, typeA;
  double qdist;
  double qM;
  double alpha;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix efield/tip4p does not exist

Self-explanatory.

E: Fix efield/tip4p requires atom attribute q or mu

The atom style defined does not have this attribute.

E: Variable name for fix efield/tip4p does not exist

Self-explanatory.

E: Variable for fix efield/tip4p is invalid style

The variable must be an equal- or atom-style variable.

E: Region ID for fix aveforce does not exist

Self-explanatory.

E: Fix efield/tip4p with dipoles cannot use atom-style variables

This option is not supported.

W: The minimizer does not re-orient dipoles when using fix efield/tip4p

This means that only the atom coordinates will be minimized,
not the orientation of the dipoles.

E: Cannot use variable energy with constant efield/tip4p in fix efield/tip4p

LAMMPS computes the energy itself when the E-field is constant.

E: Must use variable energy with fix efield/tip4p

You must define an energy when performing a minimization with a
variable E-field.

*/
