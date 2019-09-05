/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Christina Payne (Vanderbilt U)
                        Stan Moore (Sandia) for dipole terms

   Modified from constant E to constant D, adding in TIP4P
   capability (2017-2019): Stephen J. Cox (University of Cambridge)

   Copyright (2019) Stephen J. Cox
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_dfield_tip4p.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "region.h"
#include "memory.h"
#include "error.h"
#include "math_const.h" // SJC: access to 4*pi
#include "compute.h"

//SJC:
#include "angle.h"
#include "bond.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

#define INVOKED_SCALAR 1

/* ---------------------------------------------------------------------- */

FixDfieldTip4p::FixDfieldTip4p(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), xstr_d(NULL), ystr_d(NULL), zstr_d(NULL),
  xstr_p(NULL), ystr_p(NULL), zstr_p(NULL),
  estr(NULL), idregion(NULL), dfield_tip4p(NULL)
{
  // SJC: At the moment, I've only worked with real units. Possible it
  // might work with other unit systems, but will need testing. Note
  // that the energy conversion factors at hard-coded in at the
  // moment, so the energy computed by this fix won't give sensible
  // numbers for other unit systems yet.
  if (strcmp(update->unit_style,"real") != 0) error->warning(FLERR,"Energy computed with fix dfield/tip4p only compatible with 'real' energy units");
    
  // SJC: Note that unlike the LAMMPS fix_efield, we require the
  // intinerant polarization (Px, Py and Pz, multiplied by the volume)
  // to be passed to the fix too.
  if (narg < 15) error->all(FLERR,"Illegal fix dfield/tip4p command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  scalar_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  extscalar = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  dxflag = dyflag = dzflag = 1; // SJC: default to d-field being applied
                                // in all directions.

  qqr2e = force->qqr2e;
  // SJC: This converts q*p to a force The following discussion on the
  // mailing list may be useful:
  // http://lammps.sandia.gov/threads/msg10607.html.  Note that qqr2e
  // is initially for converting energies. However, it also works out
  // that this properly converts 4*pi*P to a force. As the D field is
  // in V/Ang, it has a different conversion factor.
  
  qe2f = force->qe2f; 
  xstr_d = ystr_d = zstr_d = NULL;

  // SJC: In order to calculate the energy, I need to define some of
  // my own conversion factors. If this works, then I could consider
  // adding these to update.cpp in time.
  Dconv_P = 4803204.25; // StatcoulPer_e/CmPerAng**2
  Dconv_D = 333564.09519815206; // StatVoltPerVolt/CmPerAng
  Dconv   = 1.4393260000000001e-11; // KcalMolPerErg*CmPerAng**3.
  
  // SJC: my understanding of this is that these three clauses
  // determine if a number or a LAMMPS variable for the field is being
  // passed. In the next three clauses, I should just be able to
  // change "ex->dx" etc. For the polarization, this must be a
  // "compute" defined in the LAMMPS input script and not a constant
  // number.
  if (strstr(arg[3],"v_") == arg[3]) {
    error->all(FLERR,"fix dfield/tip4p only constant displacement fields at the moment");
    int n = strlen(&arg[3][2]) + 1;
    xstr_d = new char[n];
    strcpy(xstr_d,&arg[3][2]);
  } else {
    if (strcmp(arg[3],"NULL") == 0){dxflag = 0; dx=0.0;} else {dx = force->numeric(FLERR,arg[3]);}
    xstyle = CONSTANT;
  }

  if (strstr(arg[4],"v_") == arg[4]) {
    error->all(FLERR,"fix dfield/tip4p only constant displacement fields at the moment");
    int n = strlen(&arg[4][2]) + 1;
    ystr_d = new char[n];
    strcpy(ystr_d,&arg[4][2]);
  } else {
    if (strcmp(arg[4],"NULL") == 0){dyflag = 0; dy=0.0;} else{dy = force->numeric(FLERR,arg[4]);}
    ystyle = CONSTANT;
  }

  if (strstr(arg[5],"v_") == arg[5]) {
    error->all(FLERR,"fix dfield/tip4p only constant displacement fields at the moment");
    int n = strlen(&arg[5][2]) + 1;
    zstr_d = new char[n];
    strcpy(zstr_d,&arg[5][2]);
  } else {
    //    dz = qe2f * force->numeric(FLERR,arg[5]);
    if (strcmp(arg[5],"NULL") == 0){dzflag = 0; dz=0.0;} else{dz = force->numeric(FLERR,arg[5]);}
    zstyle = CONSTANT;
  }

  // SJC: the below three clauses get the name of the compute for the
  // polarization components.
  if (strstr(arg[6],"c_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    xstr_p = new char[n];
    strcpy(xstr_p,&arg[6][2]);
  } else {
    error->all(FLERR,"Polarization in fix dfield/tip4p must be a compute.");
  }

  if (strstr(arg[7],"c_") == arg[7]) {
    int n = strlen(&arg[7][2]) + 1;
    ystr_p = new char[n];
    strcpy(ystr_p,&arg[7][2]);
  } else {
    error->all(FLERR,"Polarization in fix dfield/tip4p must be a compute.");
  }

  if (strstr(arg[8],"c_") == arg[8]) {
    int n = strlen(&arg[8][2]) + 1;
    zstr_p = new char[n];
    strcpy(zstr_p,&arg[8][2]);
  } else {
    error->all(FLERR,"Polarization in fix dfield/tip4p must be a compute.");
  }

  typeO = force->inumeric(FLERR,arg[9]);
  typeH = force->inumeric(FLERR,arg[10]);
  typeB = force->inumeric(FLERR,arg[11]);
  typeA = force->inumeric(FLERR,arg[12]);
  qdist = force->numeric(FLERR,arg[13]);
  qM    = force->numeric(FLERR,arg[14]);

  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (cos(0.5*theta) * blen);
  
  // optional args

  iregion = -1;
  idregion = NULL;
  estr = NULL;

  int iarg = 15;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      error->all(FLERR,"fix dfield/tip4p not yet implemented with 'region' keyword");
    } else if (strcmp(arg[iarg],"energy") == 0) {
      error->all(FLERR,"fix dfield/tip4p not yet implemented with 'energy' keyword");
    } else error->all(FLERR,"Illegal fix dfield/tip4p command");
  }

  force_flag = 0;
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;

  maxatom = atom->nmax;
  memory->create(dfield_tip4p,maxatom,4,"dfield_tip4p:dfield_tip4p");
}

/* ---------------------------------------------------------------------- */

FixDfieldTip4p::~FixDfieldTip4p()
{
  // SJC: added xstr_p etc.
  delete [] xstr_d;
  delete [] ystr_d;
  delete [] zstr_d;
  delete [] xstr_p;
  delete [] ystr_p;
  delete [] zstr_p;
  delete [] estr;
  delete [] idregion;
  memory->destroy(dfield_tip4p);
}

/* ---------------------------------------------------------------------- */

int FixDfieldTip4p::setmask()
{
  int mask = 0;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDfieldTip4p::init()
{
  qflag = muflag = 0;
  if (atom->q_flag) qflag = 1;
  if (atom->mu_flag && atom->torque_flag) muflag = 1;
  if (!qflag && !muflag)
    error->all(FLERR,"Fix efield requires atom attribute q or mu");

  // SJC: warn that energy hasn't been implemented for dipoles yet...
  if(muflag){error->warning(FLERR,"Energy calculation not implemented propely for dipoles yet");}
  
  // check variables

  if (xstr_d) {
    xvar = input->variable->find(xstr_d);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix dfield/tip4p does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix dfield/tip4p is invalid style");
  }
  if (ystr_d) {
    yvar = input->variable->find(ystr_d);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix dfield/tip4p does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix dfield/tip4p is invalid style");
  }
  if (zstr_d) {
    zvar = input->variable->find(zstr_d);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix dfield/tip4p does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix dfield/tip4p is invalid style");
  }
  if (estr) {
    evar = input->variable->find(estr);
    if (evar < 0)
      error->all(FLERR,"Variable name for fix efield does not exist");
    if (input->variable->atomstyle(evar)) estyle = ATOM;
    else error->all(FLERR,"Variable for fix efield is invalid style");
  } else estyle = NONE;

  // SJC: maybe need to do something similar for the polarization?
  int iOmegaPx = modify->find_compute(xstr_p);
  int iOmegaPy = modify->find_compute(ystr_p);
  int iOmegaPz = modify->find_compute(zstr_p);

  if (iOmegaPx < 0 || iOmegaPy < 0 || iOmegaPz < 0)
    error->all(FLERR,"Could not find compute for polaraztion in fix dfield/tip4p");
  c_OmegaPx = modify->compute[iOmegaPx];
  c_OmegaPy = modify->compute[iOmegaPy];
  c_OmegaPz = modify->compute[iOmegaPz];

  // set index and check validity of region

  // SJC: I think all the things below regarding region and energy
  // keywords should be OK, as I call an error during construction if
  // these keywords have been used.

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix aveforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (muflag && varflag == ATOM)
    error->all(FLERR,"Fix dfield/tip4p with dipoles cannot use atom-style variables");

  if (muflag && update->whichflag == 2 && comm->me == 0)
    error->warning(FLERR,
                   "The minimizer does not re-orient dipoles "
                   "when using fix dfield/tip4p");

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant efield in fix dfield/tip4p");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix dfield/tip4p");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixDfieldTip4p::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixDfieldTip4p::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   SJC: apply the force according the vanderbilt hamiltonian
------------------------------------------------------------------------- */

void FixDfieldTip4p::post_force(int vflag)
{
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  //SJC:
  int *type = atom->type;
  
  // reallocate efield array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(dfield_tip4p);
    memory->create(dfield_tip4p,maxatom,4,"dfield_tip4p:dfield_tip4p");
  }

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  double **x = atom->x;
  double fx,fy,fz;

  // constant efield

  if (varflag == CONSTANT) {
    double unwrap[3];

    double Omega = domain->xprd * domain->yprd * domain->zprd;
    double Omegainv = 1.0/Omega;

    // SJC: do I need to worry about invoking? c.f. compute_heat_flux.cpp line 109
    modify->clearstep_compute(); // SJC: this is required to clear the
				 // invoked flags, otherwise the
				 // polarization isn't properly updated

    // SJC: This is a bit of a hack. Not sure if it's entirely needed
    // but it seems to work.
    if (!(c_OmegaPx->invoked_flag & INVOKED_SCALAR)) {
      c_OmegaPx->compute_scalar();
      c_OmegaPx->invoked_flag |= INVOKED_SCALAR;
    }
    if (!(c_OmegaPy->invoked_flag & INVOKED_SCALAR)) {
      c_OmegaPy->compute_scalar();
      c_OmegaPy->invoked_flag |= INVOKED_SCALAR;
    }
    if (!(c_OmegaPz->invoked_flag & INVOKED_SCALAR)) {
      c_OmegaPz->compute_scalar();
      c_OmegaPz->invoked_flag |= INVOKED_SCALAR;
    }

    double Px = Omegainv*c_OmegaPx->scalar;
    double Py = Omegainv*c_OmegaPy->scalar;
    double Pz = Omegainv*c_OmegaPz->scalar;

    //SJC: This computes the energy.
    double DminusP_x = 0.0;
    double DminusP_y = 0.0;
    double DminusP_z = 0.0;
    if(dxflag){DminusP_x += Dconv_D*dx - Dconv_P*MY_4PI*Px;}
    if(dyflag){DminusP_y += Dconv_D*dy - Dconv_P*MY_4PI*Py;}
    if(dzflag){DminusP_z += Dconv_D*dz - Dconv_P*MY_4PI*Pz;}
    fsum[0] = 0.5*Dconv*Omega*(DminusP_x*DminusP_x + DminusP_y*DminusP_y + DminusP_z*DminusP_z)/MY_4PI;
    
    if (qflag) {
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

	  //SJC:
	  if(type[i]==typeO){
	    if(dxflag){fx = qM*(1-alpha)*(qe2f*dx - qqr2e*MY_4PI*Px);} else {fx = 0.0;}
	    if(dyflag){fy = qM*(1-alpha)*(qe2f*dy - qqr2e*MY_4PI*Py);} else {fy = 0.0;}
	    if(dzflag){fz = qM*(1-alpha)*(qe2f*dz - qqr2e*MY_4PI*Pz);} else {fz = 0.0;}
	    f[i][0] += fx;
	    f[i][1] += fy;
	    f[i][2] += fz;
	  }
	  else if(type[i]==typeH){
	    if(dxflag){fx = (0.5*alpha*qM + q[i])*(qe2f*dx - qqr2e*MY_4PI*Px);} else {fx = 0.0;}
	    if(dyflag){fy = (0.5*alpha*qM + q[i])*(qe2f*dy - qqr2e*MY_4PI*Py);} else {fy = 0.0;}
	    if(dzflag){fz = (0.5*alpha*qM + q[i])*(qe2f*dz - qqr2e*MY_4PI*Pz);} else {fz = 0.0;}
	    f[i][0] += fx;
	    f[i][1] += fy;
	    f[i][2] += fz;
	  }
	  else{
	    if(dxflag){fx = q[i]*(qe2f*dx - qqr2e*MY_4PI*Px);} else {fx = 0.0;}
	    if(dyflag){fy = q[i]*(qe2f*dy - qqr2e*MY_4PI*Py);} else {fy = 0.0;}
	    if(dzflag){fz = q[i]*(qe2f*dz - qqr2e*MY_4PI*Pz);} else {fz = 0.0;}	    
	    f[i][0] += fx;
	    f[i][1] += fy;
	    f[i][2] += fz;
	  }
	  
          domain->unmap(x[i],image[i],unwrap);

          fsum[1] += fx;
          fsum[2] += fy;
          fsum[3] += fz;
        }
    }

    // dipole interactions
    // no force, torque = mu cross E, potential energy = -mu dot E

    if (muflag) {
      double **mu = atom->mu;
      double **t = atom->torque;
      double tx,ty,tz;
      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit) {
	  if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;

	  double ex = qe2f*dx - qqr2e*MY_4PI*Px;
	  double ey = qe2f*dy - qqr2e*MY_4PI*Py;
	  double ez = qe2f*dz - qqr2e*MY_4PI*Pz;
	  
	  tx = ez*mu[i][1] - ey*mu[i][2];
	  ty = ex*mu[i][2] - ez*mu[i][0];
	  tz = ey*mu[i][0] - ex*mu[i][1];
	  t[i][0] += tx;
	  t[i][1] += ty;
	  t[i][2] += tz;
	  fsum[0] -= mu[i][0]*ex + mu[i][1]*ey + mu[i][2]*ez;
	}
    }

  // variable efield, wrap with clear/add
  // potential energy = evar if defined, else 0.0

  } else {

    error->all(FLERR,"fix dfield/tip4p only constant displacement fields at the moment");
    
  }
}

/* ---------------------------------------------------------------------- */

void FixDfieldTip4p::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDfieldTip4p::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixDfieldTip4p::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

double FixDfieldTip4p::compute_scalar(void)
{
  // SJC: Note the difference in structure from fix_efield
  return fsum[0];
}

/* ----------------------------------------------------------------------
   return total extra force due to fix
------------------------------------------------------------------------- */

double FixDfieldTip4p::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return fsum_all[n+1];
}
