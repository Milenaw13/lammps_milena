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
  if (mask[i] & groupbit) {
    if (tdamp_flag) add_tdamping(spi,fmi);
    if (temp_flag) add_temperature(fmi);
  }
}

/* ----------------------------------------------------------------------
   draw the random field of spin i for a spin update of length step,
   which is then kept fixed during the update (fix nve/spin midpoint).
   The field is constant during the update, so its integral (a rotation
   angle) has the variance sigma^2 step^2 = 2 D step^2 with D ~ 1/step:
   sigma was set for a step of dts, so it is scaled by sqrt(dts/step).
   Zero for atoms not in the group.
------------------------------------------------------------------------- */

void FixLangevinSpin::draw_noise(int i, double step, double *noise)
{
  noise[0] = noise[1] = noise[2] = 0.0;
  if (!temp_flag || !(atom->mask[i] & groupbit)) return;
  const double amp = sigma * sqrt(dts / step);
  noise[0] = amp * random->gaussian();
  noise[1] = amp * random->gaussian();
  noise[2] = amp * random->gaussian();
}

/* ----------------------------------------------------------------------
   same as compute_single_langevin(), with a given random field
------------------------------------------------------------------------- */

void FixLangevinSpin::compute_single_langevin_noise(int i, double spi[3], double fmi[3],
                                                    const double *noise)
{
  if (atom->mask[i] & groupbit) {
    if (tdamp_flag) add_tdamping(spi,fmi);
    if (temp_flag) {
      fmi[0] += noise[0];
      fmi[1] += noise[1];
      fmi[2] += noise[2];
      fmi[0] *= gil_factor;
      fmi[1] *= gil_factor;
      fmi[2] *= gil_factor;
    }
  }
}

/* ----------------------------------------------------------------------
   tally the energy removed from spin i by the bath in one spin update of
   length step, called by fix nve/spin after the update with the spin
   before (sold) and after (snew) it, and with fm the conservative
   precession vector (rad/ps, without damping and random field) used for
   the update: at sold (sequential update of single spins) or at the
   midpoint of sold and snew (midpoint integrator).
   The energy of the spins changes by -hbar fm.(snew - sold): exactly for
   interactions linear in each spin (sequential) or quadratic in the spins
   (midpoint: the gradient at the midpoint gives the exact difference),
   to first order in the step otherwise. The update rotates s around the
   total vector (conservative + damping + random field); the conservative
   part alone conserves this energy, so the change is the energy exchanged
   with the bath.
   In addition, the expected exchange from the average power, given s and
   fm (damping part [0] and random field part [1]):
     P = alpha/(1+alpha^2) [2 kB T s.fm - hbar |s x fm|^2]
   Its sum agrees with the exact tally on average, but misses the work of
   the actual random field, so it diffuses away from it.
   All tallied as energy REMOVED from the spins (fix langevin convention).
------------------------------------------------------------------------- */

void FixLangevinSpin::tally_spin_update(int i, const double *fm, const double *sold,
                                        const double *snew, double step)
{
  if (!tallyflag || !(atom->mask[i] & groupbit)) return;

  const double hbar = force->hplanck/MY_2PI;    // eV/(rad.THz)
  const double kb = force->boltz;               // eV/K

  energy += hbar*(fm[0]*(snew[0] - sold[0]) + fm[1]*(snew[1] - sold[1]) +
                  fm[2]*(snew[2] - sold[2]));

  const double pref = alpha_t*gil_factor*step;
  const double cpx = fm[1]*sold[2] - fm[2]*sold[1];
  const double cpy = fm[2]*sold[0] - fm[0]*sold[2];
  const double cpz = fm[0]*sold[1] - fm[1]*sold[0];
  const double sdotf = sold[0]*fm[0] + sold[1]*fm[1] + sold[2]*fm[2];

  if (tdamp_flag) energy_avg[0] += pref*hbar*(cpx*cpx + cpy*cpy + cpz*cpz);
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
