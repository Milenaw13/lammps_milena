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
   Contributing author: Milena Wittermans (Utrecht University)
   Based on pair style spin/exchange by Julien Tranchida (SNL)
   and Aidan Thompson (SNL).
------------------------------------------------------------------------- */

#include "pair_spin_exchange_spline.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "info.h"
#include "memory.h"
#include "neigh_list.h"
#include "table_file_reader.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

// tolerances for the requirement J(rc) = 0 and dJ/dr(rc) = 0,
// relative to the largest |J| and |dJ/dr| in the table

static constexpr double TOL_JCUT = 1.0e-4;
static constexpr double TOL_DJCUT = 1.0e-3;

// relative tolerances for equally spaced points and for r_last = rc

static constexpr double TOL_UNIFORM = 1.0e-3;
static constexpr double TOL_RCUT = 1.0e-6;

/* ----------------------------------------------------------------------
   derivative at x[0] of the parabola through (x[0..2], y[0..2])
------------------------------------------------------------------------- */

static double end_slope(double x0, double x1, double x2, double y0, double y1, double y2)
{
  return y0 * (2.0 * x0 - x1 - x2) / ((x0 - x1) * (x0 - x2)) +
      y1 * (x0 - x2) / ((x1 - x0) * (x1 - x2)) + y2 * (x0 - x1) / ((x2 - x0) * (x2 - x1));
}

/* ---------------------------------------------------------------------- */

PairSpinExchangeSpline::PairSpinExchangeSpline(LAMMPS *lmp) :
    PairSpin(lmp), tabindex(nullptr), cut_spin_exchange(nullptr)
{
  e_offset = 0;
  nextrapolate = 0;
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairSpinExchangeSpline::~PairSpinExchangeSpline()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_exchange);
    memory->destroy(tabindex);
    memory->destroy(cutsq);
    memory->destroy(emag);
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::settings(int narg, char **arg)
{
  PairSpin::settings(narg, arg);

  if (narg != 1) error->all(FLERR, "Illegal pair_style spin/exchange/spline command");

  cut_spin_exchange_global = utils::numeric(FLERR, arg[0], false, lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs
   pair_coeff I J exchange rc file keyword [offset yes/no]
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if ((narg != 6) && (narg != 8))
    error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (strcmp(arg[2], "exchange") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  const double rc = utils::numeric(FLERR, arg[3], false, lmp);
  if (rc <= 0.0) error->all(FLERR, "Exchange cutoff must be > 0.0");
  if (rc > cut_spin_exchange_global)
    error->all(FLERR, "Exchange cutoff {} is larger than the pair style cutoff {}", rc,
               cut_spin_exchange_global);

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "offset") == 0) {
      e_offset = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  }

  // read and spline the table on proc 0, then broadcast

  tables.emplace_back();
  Table &tb = tables.back();
  if (comm->me == 0) read_table(tb, arg[4], arg[5], rc);
  bcast_table(tb);
  const int itable = (int) tables.size() - 1;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut_spin_exchange[i][j] = rc;
      tabindex[i][j] = itable;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   read one table section, check it, and compute the spline
   only called by proc 0
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::read_table(Table &tb, const std::string &file,
                                        const std::string &keyword, double rc)
{
  TableFileReader reader(lmp, file, "spin/exchange/spline");

  char *line = reader.find_section_start(keyword);
  if (!line) error->one(FLERR, "Did not find keyword {} in table file {}", keyword, file);

  // parameter line: N n

  int n = 0;
  line = reader.next_line();
  try {
    ValueTokenizer values(line);
    while (values.has_next()) {
      std::string word = values.next_string();
      if (word == "N")
        n = values.next_int();
      else
        error->one(FLERR, "Invalid keyword {} in table {} parameters", word, keyword);
    }
  } catch (TokenizerException &e) {
    error->one(FLERR, "Error parsing parameters of table {}: {}", keyword, e.what());
  }
  if (n < 3) error->one(FLERR, "Table {} must set N >= 3", keyword);

  tb.r.resize(n);
  tb.J.resize(n);
  for (int i = 0; i < n; i++) {
    line = reader.next_line();
    if (!line)
      error->one(FLERR, "Data missing when reading table {} line {} of {}", keyword, i + 1, n);
    try {
      ValueTokenizer values(line);
      values.next_int();
      tb.r[i] = values.next_double();
      tb.J[i] = values.next_double();
    } catch (TokenizerException &e) {
      error->one(FLERR, "Error parsing table {} line {} of {}: {}\nLine was: {}", keyword, i + 1,
                 n, e.what(), line);
    }
  }

  // distances must increase and the table must end at the cutoff

  if (tb.r[0] <= 0.0) error->one(FLERR, "Distances in table {} must be > 0.0", keyword);
  for (int i = 1; i < n; i++)
    if (tb.r[i] <= tb.r[i - 1])
      error->one(FLERR, "Distances in table {} must be strictly increasing (line {})", keyword,
                 i + 1);
  if (fabs(tb.r[n - 1] - rc) > TOL_RCUT * rc)
    error->one(FLERR, "Table {} ends at r = {} but the exchange cutoff is {}: the last point "
               "of the table must be at the cutoff", keyword, tb.r[n - 1], rc);

  // J and dJ/dr must vanish at the cutoff, otherwise energy and forces jump there

  double jmax = 0.0, djmax = 0.0;
  for (int i = 0; i < n; i++) jmax = MAX(jmax, fabs(tb.J[i]));
  for (int i = 1; i < n; i++)
    djmax = MAX(djmax, fabs((tb.J[i] - tb.J[i - 1]) / (tb.r[i] - tb.r[i - 1])));

  const double djhi = end_slope(tb.r[n - 1], tb.r[n - 2], tb.r[n - 3], tb.J[n - 1], tb.J[n - 2],
                                tb.J[n - 3]);
  if (fabs(tb.J[n - 1]) > TOL_JCUT * jmax)
    error->one(FLERR, "Table {}: J(rc) = {} eV is not zero; the exchange curve must go "
               "smoothly to zero at the cutoff", keyword, tb.J[n - 1]);
  if (fabs(djhi) > TOL_DJCUT * djmax)
    error->one(FLERR, "Table {}: dJ/dr(rc) = {} eV/Angstrom is not zero; the exchange curve "
               "must go smoothly to zero at the cutoff", keyword, djhi);

  // check for equally spaced points, which allows a direct bin lookup

  const double dr = (tb.r[n - 1] - tb.r[0]) / (n - 1);
  tb.uniform = 1;
  for (int i = 1; i < n; i++)
    if (fabs(tb.r[i] - tb.r[i - 1] - dr) > TOL_UNIFORM * dr) tb.uniform = 0;
  tb.invdr = 1.0 / dr;

  // spline with the slope at r[0] from the data and dJ/dr = 0 imposed at rc

  tb.dJlo = end_slope(tb.r[0], tb.r[1], tb.r[2], tb.J[0], tb.J[1], tb.J[2]);
  spline_table(tb, tb.dJlo, 0.0);
}

/* ----------------------------------------------------------------------
   cubic spline with given end slopes: sets the second derivatives J2
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::spline_table(Table &tb, double yp1, double ypn)
{
  const int n = tb.r.size();
  const std::vector<double> &x = tb.r;
  const std::vector<double> &y = tb.J;
  std::vector<double> &y2 = tb.J2;
  std::vector<double> u(n);
  y2.resize(n);

  y2[0] = -0.5;
  u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  for (int i = 1; i < n - 1; i++) {
    const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    const double p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }
  const double qn = 0.5;
  const double un =
      (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (int k = n - 2; k >= 0; k--) y2[k] = y2[k] * y2[k + 1] + u[k];
}

/* ----------------------------------------------------------------------
   broadcast a table from proc 0 to all procs
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::bcast_table(Table &tb)
{
  int n = tb.r.size();
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  tb.r.resize(n);
  tb.J.resize(n);
  tb.J2.resize(n);
  MPI_Bcast(tb.r.data(), n, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb.J.data(), n, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb.J2.data(), n, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb.invdr, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb.dJlo, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb.uniform, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinExchangeSpline::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, Error::NOLASTLINE,
               "All pair coeffs are not set. Status\n" + Info::get_pair_coeff_status(lmp));

  tabindex[j][i] = tabindex[i][j];
  cut_spin_exchange[j][i] = cut_spin_exchange[i][j];

  return cut_spin_exchange_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinExchangeSpline::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "cut") == 0) return (void *) &cut_spin_exchange_global;
  return nullptr;
}

/* ----------------------------------------------------------------------
   report pairs that came closer than the start of their table
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::finish()
{
  bigint all = 0;
  MPI_Allreduce(&nextrapolate, &all, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if ((comm->me == 0) && (all > 0))
    error->warning(FLERR,
                   "Pair spin/exchange/spline: {} pair distances were below the start of "
                   "the exchange table; J(r) was extrapolated linearly", all);
  nextrapolate = 0;
}

/* ----------------------------------------------------------------------
   exchange coupling J(r) (eV) and dJ/dr (eV/Angstrom) for a type pair
   J and dJ/dr come from the same cubic polynomial of the spline bin,
   below the first table point J(r) is extrapolated linearly
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::exchange_J(int itype, int jtype, double r, double &J, double &dJ)
{
  const Table &tb = tables[tabindex[itype][jtype]];
  const int n = tb.r.size();
  const double *x = tb.r.data();

  if (r < x[0]) {
    nextrapolate++;
    dJ = tb.dJlo;
    J = tb.J[0] + dJ * (r - x[0]);
    return;
  }

  // find bin k with x[k] <= r <= x[k+1]

  int k;
  if (tb.uniform) {
    k = static_cast<int>((r - x[0]) * tb.invdr);
    k = MIN(k, n - 2);
    while ((k > 0) && (r < x[k])) k--;
    while ((k < n - 2) && (r > x[k + 1])) k++;
  } else {
    int klo = 0, khi = n - 1;
    while (khi - klo > 1) {
      const int kmid = (khi + klo) >> 1;
      if (x[kmid] > r)
        khi = kmid;
      else
        klo = kmid;
    }
    k = klo;
  }

  const double h = x[k + 1] - x[k];
  const double a = (x[k + 1] - r) / h;
  const double b = (r - x[k]) / h;
  J = a * tb.J[k] + b * tb.J[k + 1] +
      ((a * a * a - a) * tb.J2[k] + (b * b * b - b) * tb.J2[k + 1]) * (h * h) / 6.0;
  dJ = (tb.J[k + 1] - tb.J[k]) / h +
      (-(3.0 * a * a - 1.0) * tb.J2[k] + (3.0 * b * b - 1.0) * tb.J2[k + 1]) * h / 6.0;
}

/* ---------------------------------------------------------------------- */

void PairSpinExchangeSpline::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double evdwl, ecoul;
  double xi[3], eij[3];
  double delx, dely, delz;
  double spi[3], spj[3];
  double fi[3];
  double local_cut2;
  double rsq, r, inorm, Jex, dJex, sdots;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **fm = atom->fm;
  double **sp = atom->sp;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // checking size of emag

  if (nlocal_max < nlocal) {    // grow emag lists if necessary
    nlocal_max = nlocal;
    memory->grow(emag, nlocal_max, "pair/spin:emag");
  }

  // computation of the exchange interaction
  // loop over atoms and their neighbors

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];
    emag[i] = 0.0;

    // loop on neighbors

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      local_cut2 = cut_spin_exchange[itype][jtype] * cut_spin_exchange[itype][jtype];
      if (rsq > local_cut2) continue;

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      r = sqrt(rsq);
      inorm = 1.0 / r;
      eij[0] = -inorm * delx;
      eij[1] = -inorm * dely;
      eij[2] = -inorm * delz;

      exchange_J(itype, jtype, r, Jex, dJex);

      // precession: fm_i += J/hbar s_j

      fm[i][0] += Jex / hbar * spj[0];
      fm[i][1] += Jex / hbar * spj[1];
      fm[i][2] += Jex / hbar * spj[2];

      sdots = spi[0] * spj[0] + spi[1] * spj[1] + spi[2] * spj[2];
      if (e_offset) sdots -= 1.0;

      // mechanical force from dJ/dr, half per neighbor list entry

      fi[0] = fi[1] = fi[2] = 0.0;
      if (lattice_flag) {
        fi[0] = -0.5 * dJex * sdots * eij[0];
        fi[1] = -0.5 * dJex * sdots * eij[1];
        fi[2] = -0.5 * dJex * sdots * eij[2];
      }

      if (eflag) {
        evdwl = -0.5 * Jex * sdots;
        emag[i] += evdwl;
      } else
        evdwl = 0.0;

      f[i][0] += fi[0];
      f[i][1] += fi[1];
      f[i][2] += fi[2];
      if (newton_pair || j < nlocal) {
        f[j][0] -= fi[0];
        f[j][1] -= fi[1];
        f[j][2] -= fi[2];
      }

      if (evflag)
        ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, ecoul, fi[0], fi[1], fi[2], delx, dely,
                     delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   update the pair interactions fmi acting on the spin ii
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2, delx, dely, delz, rsq, Jex, dJex;
  int j, jnum, jtype;
  int *jlist;

  // check if interaction applies to type of ii

  const int itype = type[ii];
  int locflag = 0;
  for (int k = 1; k <= atom->ntypes; k++) {
    if (((k <= itype) && (setflag[k][itype] == 1)) || ((k > itype) && (setflag[itype][k] == 1))) {
      locflag = 1;
      break;
    }
  }
  if (locflag == 0) return;

  jlist = list->firstneigh[ii];
  jnum = list->numneigh[ii];

  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];
    j &= NEIGHMASK;
    jtype = type[j];
    local_cut2 = cut_spin_exchange[itype][jtype] * cut_spin_exchange[itype][jtype];

    delx = x[ii][0] - x[j][0];
    dely = x[ii][1] - x[j][1];
    delz = x[ii][2] - x[j][2];
    rsq = delx * delx + dely * dely + delz * delz;

    if (rsq <= local_cut2) {
      exchange_J(itype, jtype, sqrt(rsq), Jex, dJex);
      fmi[0] += Jex / hbar * sp[j][0];
      fmi[1] += Jex / hbar * sp[j][1];
      fmi[2] += Jex / hbar * sp[j][2];
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinExchangeSpline::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cut_spin_exchange, n + 1, n + 1, "pair/spin/exchange/spline:cut_spin_exchange");
  memory->create(tabindex, n + 1, n + 1, "pair/spin/exchange/spline:tabindex");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}
