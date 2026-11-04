/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(spin/exchange/spline,PairSpinExchangeSpline);
// clang-format on
#else

#ifndef LMP_PAIR_SPIN_EXCHANGE_SPLINE_H
#define LMP_PAIR_SPIN_EXCHANGE_SPLINE_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinExchangeSpline : public PairSpin {
 public:
  PairSpinExchangeSpline(class LAMMPS *);
  ~PairSpinExchangeSpline() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;
  void finish() override;

  void compute(int, int) override;
  void compute_single_pair(int, double *) override;

  // exchange coupling J(r) in eV and its derivative dJ/dr in eV/Angstrom

  void exchange_J(int, int, double, double &, double &);

  double cut_spin_exchange_global;    // global exchange cutoff distance

 protected:
  // cubic spline through the tabulated (r, J) points of one table section

  struct Table {
    std::vector<double> r;     // distances (Angstrom), strictly increasing
    std::vector<double> J;     // J(r) values (eV)
    std::vector<double> J2;    // spline second derivatives
    double invdr;              // 1/spacing, used for equally spaced tables
    double dJlo;               // dJ/dr at r[0], used to extrapolate below r[0]
    int uniform;               // 1 if points are equally spaced
    std::string keyword;       // section keyword in the table file
    bigint nbelow;             // pair distances below r[0] seen by compute()
  };

  int e_offset;                    // apply energy offset
  std::vector<Table> tables;       // all tables read by pair_coeff
  int **tabindex;                  // table index for each type pair
  double **cut_spin_exchange;      // cutoff distance exchange

  void allocate() override;
  void read_table(Table &, const std::string &, const std::string &, double);
  void bcast_table(Table &);
  void spline_table(Table &, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
