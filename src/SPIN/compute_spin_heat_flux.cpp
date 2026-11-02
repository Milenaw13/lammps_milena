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

/* ----------------------------------------------------------------------
   Contributing author: M. A. Wittermans (Utrecht University)
   Based on fix spin/heat/flux by Yanguang Zhou (UCLA), extended by
   M. A. Wittermans to all spin components.

   Spin heat flux across the interface between two groups, as in:
   Y. Zhou, J. Tranchida, Y. Ge, J. Murthy, T. S. Fisher,
   Phys. Rev. B 101, 224303 (2020), eq. (13).
------------------------------------------------------------------------- */

#include "compute_spin_heat_flux.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "math_const.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_spin_exchange.h"
#include "pair_spin_exchange_spline.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

ComputeSpinHeatFlux::ComputeSpinHeatFlux(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), id_a(nullptr), id_b(nullptr), list(nullptr), channel(nullptr),
    qall(nullptr), last_step(-1)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "compute spin/heat/flux", error);
  if (narg > 5) error->all(FLERR, 5, "Too many arguments for compute spin/heat/flux");
  if (igroup) error->all(FLERR, 1, "Compute spin/heat/flux must use group all");
  if (!atom->sp_flag)
    error->all(FLERR, "Compute spin/heat/flux requires atom style spin");

  delete[] id_a;
  id_a = utils::strdup(arg[3]);
  delete[] id_b;
  id_b = utils::strdup(arg[4]);
  if (group->find(id_a) < 0)
    error->all(FLERR, 3, "Could not find compute spin/heat/flux group ID {}", id_a);
  if (group->find(id_b) < 0)
    error->all(FLERR, 4, "Could not find compute spin/heat/flux group ID {}", id_b);

  // one channel per unordered pair of atom types

  const int ntypes = atom->ntypes;
  nchannel = ntypes * (ntypes + 1) / 2;
  memory->create(channel, ntypes + 1, ntypes + 1, "spin/heat/flux:channel");
  int k = 0;
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) channel[i][j] = channel[j][i] = k++;

  memory->create(vector, nchannel, "spin/heat/flux:vector");
  memory->create(qall, nchannel + 1, "spin/heat/flux:qall");

  scalar_flag = 1;
  extscalar = 1;
  vector_flag = 1;
  size_vector = nchannel;
  extvector = 1;
}

/* ---------------------------------------------------------------------- */

ComputeSpinHeatFlux::~ComputeSpinHeatFlux()
{
  delete[] id_a;
  delete[] id_b;
  memory->destroy(channel);
  memory->destroy(vector);
  memory->destroy(qall);
}

/* ---------------------------------------------------------------------- */

void ComputeSpinHeatFlux::init()
{
  const int ia = group->find(id_a);
  const int ib = group->find(id_b);
  if (ia < 0) error->all(FLERR, Error::NOLASTLINE, "Compute spin/heat/flux group {} does not exist", id_a);
  if (ib < 0) error->all(FLERR, Error::NOLASTLINE, "Compute spin/heat/flux group {} does not exist", id_b);
  bit_a = group->bitmask[ia];
  bit_b = group->bitmask[ib];

  // values of a previous run may no longer be valid

  last_step = -1;

  hbar = force->hplanck / MY_2PI;

  // collect all exchange pair styles, also as hybrid/overlay sub-styles

  if (!force->pair)
    error->all(FLERR, Error::NOLASTLINE, "Compute spin/heat/flux requires a pair style");
  if (force->pair_match("^hybrid/scaled", 0))
    error->all(FLERR, Error::NOLASTLINE,
               "Compute spin/heat/flux does not support pair style hybrid/scaled");

  styles.clear();
  find_exchange_styles("spin/exchange");
  find_exchange_styles("spin/exchange/spline");
  if (styles.empty())
    error->all(FLERR, Error::NOLASTLINE,
               "Compute spin/heat/flux requires pair style spin/exchange or "
               "spin/exchange/spline");

  // other magnetic pair styles change the spin precession, which is included,
  // but their own coupling between the two groups is not part of the flux

  if (comm->me == 0) {
    // pair_match() with nsub = 0 returns a style only if exactly one matches

    int nspin = 0;
    if (force->pair_match("^spin", 0, 0))
      nspin = 1;
    else
      while (force->pair_match("^spin", 0, nspin + 1)) nspin++;
    if (nspin > (int) styles.size())
      error->warning(FLERR,
                     "Compute spin/heat/flux only includes the exchange interaction "
                     "in the heat flux, other magnetic pair styles are ignored");
  }

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ----------------------------------------------------------------------
   add every instance of an exchange pair style with this exact name
------------------------------------------------------------------------- */

void ComputeSpinHeatFlux::find_exchange_styles(const std::string &name)
{
  std::vector<Pair *> found;
  Pair *p = force->pair_match(name, 1, 0);
  if (p) {
    found.push_back(p);
  } else {
    for (int n = 1; (p = force->pair_match(name, 1, n)); n++) found.push_back(p);
  }

  for (auto *pair : found) {
    ExchangeStyle s;
    s.pair = pair;
    s.bs = dynamic_cast<PairSpinExchange *>(pair);
    s.spline = dynamic_cast<PairSpinExchangeSpline *>(pair);
    int dim = 0;
    s.cut = (double **) pair->extract("cut_spin_exchange", dim);
    if ((!s.bs && !s.spline) || !s.cut || (dim != 2))
      error->all(FLERR, Error::NOLASTLINE,
                 "Compute spin/heat/flux could not access pair style {}", name);
    styles.push_back(s);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpinHeatFlux::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeSpinHeatFlux::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  compute_all();
  scalar = qall[nchannel];
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeSpinHeatFlux::compute_vector()
{
  invoked_vector = update->ntimestep;
  compute_all();
  for (int k = 0; k < nchannel; k++) vector[k] = qall[k];
}

/* ----------------------------------------------------------------------
   Q_AB = sum_{i in A} sum_{j in B} J(r_ij) (w_i x s_i) . s_j
   w_i = atom->fm[i], the precession vector (rad/ps), ds_i/dt = w_i x s_i
   sum over all periodic images of j within the cutoff of each exchange style
   qall[k] = contribution of the pairs in channel k, qall[nchannel] = total
   the pair loop runs once per time step for the scalar and the vector
------------------------------------------------------------------------- */

void ComputeSpinHeatFlux::compute_all()
{
  if (last_step == update->ntimestep) return;
  last_step = update->ntimestep;

  neighbor->build_one(list);

  double **x = atom->x;
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *type = atom->type;
  int *mask = atom->mask;

  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  std::vector<double> qlocal(nchannel + 1, 0.0);

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (!(mask[i] & bit_a)) continue;
    const int itype = type[i];

    // ds_i/dt = w_i x s_i

    double dsi[3];
    dsi[0] = fm[i][1] * sp[i][2] - fm[i][2] * sp[i][1];
    dsi[1] = fm[i][2] * sp[i][0] - fm[i][0] * sp[i][2];
    dsi[2] = fm[i][0] * sp[i][1] - fm[i][1] * sp[i][0];

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & bit_b)) continue;

      const double delx = x[i][0] - x[j][0];
      const double dely = x[i][1] - x[j][1];
      const double delz = x[i][2] - x[j][2];
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];
      const int tlo = MIN(itype, jtype);
      const int thi = MAX(itype, jtype);

      double spj[3] = {sp[j][0], sp[j][1], sp[j][2]};
      const double dsdots = dsi[0] * spj[0] + dsi[1] * spj[1] + dsi[2] * spj[2];
      double qij = 0.0;

      for (const auto &s : styles) {
        if (!s.pair->setflag[tlo][thi]) continue;
        const double cut = s.cut[itype][jtype];
        if (rsq > cut * cut) continue;

        if (s.bs) {
          // compute_exchange() adds J(r)/hbar * s_j

          double fij[3] = {0.0, 0.0, 0.0};
          s.bs->compute_exchange(i, j, rsq, fij, spj);
          qij += hbar * (dsi[0] * fij[0] + dsi[1] * fij[1] + dsi[2] * fij[2]);
        } else {
          double J, dJ;
          s.spline->exchange_J(itype, jtype, sqrt(rsq), J, dJ);
          qij += J * dsdots;
        }
      }
      qlocal[channel[itype][jtype]] += qij;
      qlocal[nchannel] += qij;
    }
  }

  MPI_Allreduce(qlocal.data(), qall, nchannel + 1, MPI_DOUBLE, MPI_SUM, world);
}
