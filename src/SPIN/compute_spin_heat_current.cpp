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

   Volume-integrated heat current carried by the exchange interaction
   between magnetic spins (Hardy form, energy-weighted position derivative):
     J_prec = 1/2 sum_i sum_j r_ij J(r_ij) s_j . ds_i/dt
   optionally the current of the work done by the mechanical exchange forces:
     J_mech = 1/2 sum_i sum_j r_ij 1/2 F_ij . (v_i + v_j)
   The precession term is the volume form of the spin heat flux of
   Y. Zhou, J. Tranchida, Y. Ge, J. Murthy, T. S. Fisher,
   Phys. Rev. B 101, 224303 (2020), eq. (13).
------------------------------------------------------------------------- */

#include "compute_spin_heat_current.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_nve_spin.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_spin_exchange.h"
#include "pair_spin_exchange_spline.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

ComputeSpinHeatCurrent::ComputeSpinHeatCurrent(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), mech_flag(0), lattice_flag(0), list(nullptr), channel(nullptr)
{
  if (!atom->sp_flag) error->all(FLERR, "Compute spin/heat/current requires atom style spin");

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mech") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "compute spin/heat/current mech", error);
      mech_flag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown compute spin/heat/current keyword: {}", arg[iarg]);
    }
  }

  if (igroup && (comm->me == 0))
    error->warning(FLERR,
                   "Compute spin/heat/current with group {} includes only the pairs with "
                   "the first atom in the group, this is not the heat current of the system",
                   arg[1]);

  // one channel per unordered pair of atom types

  const int ntypes = atom->ntypes;
  nchannel = ntypes * (ntypes + 1) / 2;
  memory->create(channel, ntypes + 1, ntypes + 1, "spin/heat/current:channel");
  int k = 0;
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) channel[i][j] = channel[j][i] = k++;

  // vector = J_prec (3), [J_mech (3)], J_prec per channel (3 each)

  offset_ch = mech_flag ? 6 : 3;
  vector_flag = 1;
  size_vector = offset_ch + 3 * nchannel;
  extvector = 1;
  memory->create(vector, size_vector, "spin/heat/current:vector");
}

/* ---------------------------------------------------------------------- */

ComputeSpinHeatCurrent::~ComputeSpinHeatCurrent()
{
  memory->destroy(channel);
  memory->destroy(vector);
}

/* ---------------------------------------------------------------------- */

void ComputeSpinHeatCurrent::init()
{
  hbar = force->hplanck / MY_2PI;

  // collect all exchange pair styles, also as hybrid/overlay sub-styles

  if (!force->pair)
    error->all(FLERR, Error::NOLASTLINE, "Compute spin/heat/current requires a pair style");
  if (force->pair_match("^hybrid/scaled", 0))
    error->all(FLERR, Error::NOLASTLINE,
               "Compute spin/heat/current does not support pair style hybrid/scaled");

  styles.clear();
  find_exchange_styles("spin/exchange");
  find_exchange_styles("spin/exchange/spline");
  if (styles.empty())
    error->all(FLERR, Error::NOLASTLINE,
               "Compute spin/heat/current requires pair style spin/exchange or "
               "spin/exchange/spline");

  // other magnetic pair styles change the spin precession,
  // but their own energy transfer is not part of the current

  if (comm->me == 0) {
    // pair_match() with nsub = 0 returns a style only if exactly one matches

    int nspin = 0;
    if (force->pair_match("^spin", 0, 0))
      nspin = 1;
    else
      while (force->pair_match("^spin", 0, nspin + 1)) nspin++;
    if (nspin > (int) styles.size())
      error->warning(FLERR,
                     "Compute spin/heat/current only includes the exchange interaction, "
                     "other magnetic pair styles are ignored");
  }

  // the exchange forces act on the atoms only if fix nve/spin moves the lattice,
  // same test as in PairSpin::init_style()

  lattice_flag = 0;
  auto fixes = modify->get_fix_by_style("^nve/spin");
  if (fixes.size() == 1) lattice_flag = (dynamic_cast<FixNVESpin *>(fixes.front()))->lattice_flag;
  if (mech_flag && !lattice_flag && (comm->me == 0))
    error->warning(FLERR,
                   "Compute spin/heat/current mech yes: the exchange forces do not act on "
                   "the atoms (no fix nve/spin with lattice moving), the mech current is zero");

  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
}

/* ----------------------------------------------------------------------
   add every instance of an exchange pair style with this exact name
------------------------------------------------------------------------- */

void ComputeSpinHeatCurrent::find_exchange_styles(const std::string &name)
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
      error->all(FLERR, Error::NOLASTLINE, "Compute spin/heat/current could not access pair style {}",
                 name);
    styles.push_back(s);
  }
}

/* ---------------------------------------------------------------------- */

void ComputeSpinHeatCurrent::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   J_prec = 1/2 sum_i sum_j r_ij q_ij,  q_ij = J(r_ij) s_j . ds_i/dt
   J_mech = 1/2 sum_i sum_j r_ij F_ij . v_i
   r_ij = x_i - x_j, ds_i/dt = w_i x s_i with w_i = atom->fm[i] (rad/ps),
   F_ij = force on i due to j from the exchange (2x the half force of one
   neighbor list entry). The sum over the full neighbor list visits each
   pair from both atoms, so the v_j part of 1/2 F_ij.(v_i+v_j) is counted
   from the owner of j and ghost velocities are not needed.
------------------------------------------------------------------------- */

void ComputeSpinHeatCurrent::compute_vector()
{
  invoked_vector = update->ntimestep;

  neighbor->build_one(list);

  double **x = atom->x;
  double **v = atom->v;
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *type = atom->type;
  int *mask = atom->mask;

  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  const int domech = mech_flag && lattice_flag;
  std::vector<double> jlocal(size_vector, 0.0);

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    const int itype = type[i];

    double spi[3] = {sp[i][0], sp[i][1], sp[i][2]};

    // ds_i/dt = w_i x s_i

    double dsi[3];
    dsi[0] = fm[i][1] * spi[2] - fm[i][2] * spi[1];
    dsi[1] = fm[i][2] * spi[0] - fm[i][0] * spi[2];
    dsi[2] = fm[i][0] * spi[1] - fm[i][1] * spi[0];

    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

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
      double fi[3] = {0.0, 0.0, 0.0};
      double eij[3] = {0.0, 0.0, 0.0};
      if (domech) {
        const double inorm = 1.0 / sqrt(rsq);
        eij[0] = -inorm * delx;
        eij[1] = -inorm * dely;
        eij[2] = -inorm * delz;
      }

      for (const auto &s : styles) {
        if (!s.pair->setflag[tlo][thi]) continue;
        const double cut = s.cut[itype][jtype];
        if (rsq > cut * cut) continue;

        if (s.bs) {
          // compute_exchange() adds J(r)/hbar * s_j

          double fmij[3] = {0.0, 0.0, 0.0};
          s.bs->compute_exchange(i, j, rsq, fmij, spj);
          qij += hbar * (dsi[0] * fmij[0] + dsi[1] * fmij[1] + dsi[2] * fmij[2]);
          if (domech) s.bs->compute_exchange_mech(i, j, rsq, eij, fi, spi, spj);
        } else {
          double J, dJ;
          s.spline->exchange_J(itype, jtype, sqrt(rsq), J, dJ);
          qij += J * dsdots;
          if (domech) s.spline->compute_exchange_mech(i, j, rsq, eij, fi, spi, spj);
        }
      }

      // 1/2 r_ij q_ij, total and per channel

      const double hq = 0.5 * qij;
      const int k = offset_ch + 3 * channel[itype][jtype];
      jlocal[0] += hq * delx;
      jlocal[1] += hq * dely;
      jlocal[2] += hq * delz;
      jlocal[k] += hq * delx;
      jlocal[k + 1] += hq * dely;
      jlocal[k + 2] += hq * delz;

      // fi = 1/2 F_ij -> 1/2 r_ij F_ij.v_i = r_ij fi.v_i

      if (domech) {
        const double fv = fi[0] * v[i][0] + fi[1] * v[i][1] + fi[2] * v[i][2];
        jlocal[3] += fv * delx;
        jlocal[4] += fv * dely;
        jlocal[5] += fv * delz;
      }
    }
  }

  MPI_Allreduce(jlocal.data(), vector, size_vector, MPI_DOUBLE, MPI_SUM, world);
}
