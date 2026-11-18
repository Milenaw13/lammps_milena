// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "fix_langevin_spin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "modify.h"
#include "random_mars.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixLangevinSpin::FixLangevinSpin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), random(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix langevin/spin", error);

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  ecouple_flag = 1;

  temp = utils::numeric(FLERR,arg[3],false,lmp);
  alpha_t = utils::numeric(FLERR,arg[4],false,lmp);
  seed = utils::inumeric(FLERR,arg[5],false,lmp);

  // optional keywords

  tallyflag = 0;
  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "tally") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix langevin/spin tally", error);
      tallyflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else error->all(FLERR, iarg, "Unknown fix langevin/spin keyword: {}", arg[iarg]);
  }

  energy = 0.0;
  energy_avg[0] = energy_avg[1] = 0.0;
  ilast = -1;

  if (alpha_t < 0.0) {
    error->all(FLERR,"Illegal langevin/spin command");
  } else if (alpha_t == 0.0) {
    tdamp_flag = 0;
  } else {
    tdamp_flag = 1;
  }

  if (temp < 0.0) {
    error->all(FLERR,"Illegal langevin/spin command");
  } else if (temp == 0.0) {
    temp_flag = 0;
  } else {
    temp_flag = 1;
  }

  // initialize Marsaglia RNG with processor-unique seed

  // random = new RanPark(lmp,seed + comm->me);
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixLangevinSpin::~FixLangevinSpin()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixLangevinSpin::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::init()
{
  // fix_langevin_spin has to be the last defined fix

  int flag_force = 0;
  int flag_lang = 0;
  const auto &fixes = modify->get_fix_list();
  for (int i = 0; i < (int) fixes.size(); i++) {
     if (strcmp("precession/spin",fixes[i]->style)==0) flag_force = MAX(flag_force,i);
     if (strcmp("langevin/spin",fixes[i]->style)==0) flag_lang = i;
  }
  if (flag_force >= flag_lang) error->all(FLERR,"Fix langevin/spin has to come after all other spin fixes");

  gil_factor = 1.0/(1.0+(alpha_t)*(alpha_t));
  dts = 0.25 * update->dt;

  double hbar = force->hplanck/MY_2PI;  // eV/(rad.THz)
  double kb = force->boltz;             // eV/K

  D = (alpha_t*(1.0+(alpha_t)*(alpha_t))*kb*temp);
  D /= (hbar*dts);
  sigma = sqrt(2.0*D);
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^respa")) {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa-1);
  } else post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_tdamping(double spi[3], double fmi[3])
{
  double cpx = fmi[1]*spi[2] - fmi[2]*spi[1];
  double cpy = fmi[2]*spi[0] - fmi[0]*spi[2];
  double cpz = fmi[0]*spi[1] - fmi[1]*spi[0];

  // adding the transverse damping

  fmi[0] -= alpha_t*cpx;
  fmi[1] -= alpha_t*cpy;
  fmi[2] -= alpha_t*cpz;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::add_temperature(double fmi[3])
{
  double rx = sigma*random->gaussian();
  double ry = sigma*random->gaussian();
  double rz = sigma*random->gaussian();

  // adding the random field

  fmi[0] += rx;
  fmi[1] += ry;
  fmi[2] += rz;

  // adding gilbert's prefactor

  fmi[0] *= gil_factor;
  fmi[1] *= gil_factor;
  fmi[2] *= gil_factor;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpin::compute_single_langevin(int i, double spi[3], double fmi[3])
{
  int *mask = atom->mask;
  ilast = -1;
  if (mask[i] & groupbit) {
    if (tallyflag) {
      ilast = i;
      fmlast[0] = fmi[0];
      fmlast[1] = fmi[1];
      fmlast[2] = fmi[2];
      tally_average(spi,fmi);
    }
    if (tdamp_flag) add_tdamping(spi,fmi);
    if (temp_flag) add_temperature(fmi);
  }
}

/* ----------------------------------------------------------------------
   tally the energy removed from spin i by the bath in the spin update that
   follows compute_single_langevin(i), called by fix nve/spin with the spin
   before (sold) and after (snew) the update.
   The energy of spin i is -hbar fm.s, with fm its conservative precession
   vector (rad/ps) saved in compute_single_langevin(). The update rotates s
   around the total vector (conservative + damping + random field); the
   conservative part alone keeps fm.s constant, so the change of the spin
   energy, -hbar fm.(snew - sold), is the energy exchanged with the bath.
   This is exact for interactions linear in each spin (exchange, DMI, Zeeman)
   and holds to first order in the step otherwise.
------------------------------------------------------------------------- */

void FixLangevinSpin::tally_spin_update(int i, double *sold, double *snew)
{
  if (!tallyflag || (i != ilast)) return;

  const double hbar = force->hplanck/MY_2PI;    // eV/(rad.THz)
  energy += hbar*(fmlast[0]*(snew[0] - sold[0]) + fmlast[1]*(snew[1] - sold[1]) +
                  fmlast[2]*(snew[2] - sold[2]));
  ilast = -1;
}

/* ----------------------------------------------------------------------
   average power exchanged with the bath during one spin update of length
   dts, given s and the conservative precession vector fm (rad/ps):
     P = alpha/(1+alpha^2) [2 kB T s.fm - hbar |s x fm|^2]
   damping part [0] and random field part [1], tallied as energy REMOVED
   from the spins. Their sum agrees with the exact tally on average, but
   misses the work of the actual random field, so it diffuses away from it.
------------------------------------------------------------------------- */

void FixLangevinSpin::tally_average(double spi[3], double fmi[3])
{
  const double hbar = force->hplanck/MY_2PI;    // eV/(rad.THz)
  const double kb = force->boltz;               // eV/K
  const double pref = alpha_t*gil_factor*dts;

  const double cpx = fmi[1]*spi[2] - fmi[2]*spi[1];
  const double cpy = fmi[2]*spi[0] - fmi[0]*spi[2];
  const double cpz = fmi[0]*spi[1] - fmi[1]*spi[0];
  const double sdotf = spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2];

  energy_avg[0] += pref*hbar*(cpx*cpx + cpy*cpy + cpz*cpz);
  if (temp_flag) energy_avg[1] -= pref*2.0*kb*temp*sdotf;
}

/* ----------------------------------------------------------------------
   cumulative energy removed from the spins of the group by the bath
------------------------------------------------------------------------- */

double FixLangevinSpin::compute_scalar()
{
  if (!tallyflag) return 0.0;

  double energy_all = 0.0;
  MPI_Allreduce(&energy, &energy_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy_all;
}

/* ----------------------------------------------------------------------
   cumulative energy removed from the spins from the average power of the
   damping (n = 0) and of the random field (n = 1)
------------------------------------------------------------------------- */

double FixLangevinSpin::compute_vector(int n)
{
  if (!tallyflag) return 0.0;

  double energy_all = 0.0;
  MPI_Allreduce(&energy_avg[n], &energy_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return energy_all;
}
