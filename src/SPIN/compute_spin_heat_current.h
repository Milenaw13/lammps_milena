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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(spin/heat/current,ComputeSpinHeatCurrent);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPIN_HEAT_CURRENT_H
#define LMP_COMPUTE_SPIN_HEAT_CURRENT_H

#include "compute.h"

#include <vector>

namespace LAMMPS_NS {

class ComputeSpinHeatCurrent : public Compute {
 public:
  ComputeSpinHeatCurrent(class LAMMPS *, int, char **);
  ~ComputeSpinHeatCurrent() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_vector() override;

 private:
  // one exchange pair style (or hybrid sub-style) that contributes J(r)

  struct ExchangeStyle {
    class Pair *pair;                        // the style, for setflag
    class PairSpinExchange *bs;              // non-null for spin/exchange
    class PairSpinExchangeSpline *spline;    // non-null for spin/exchange/spline
    double **cut;                            // per type pair cutoffs
  };

  int mech_flag;       // 1 if the work of the mechanical exchange forces is output
  int lattice_flag;    // 1 if the exchange forces act on the atoms (fix nve/spin)
  double hbar;         // Planck constant / 2pi (eV.ps.rad-1)
  std::vector<ExchangeStyle> styles;
  class NeighList *list;

  // channels = unordered type pairs (I <= J), in the order
  // (1,1) (1,2) ... (1,N) (2,2) ... (N,N)

  int nchannel;     // number of channels = N*(N+1)/2
  int **channel;    // channel index of each type pair
  int offset_ch;    // index of the first channel entry in vector

  void find_exchange_styles(const std::string &);
};

}    // namespace LAMMPS_NS

#endif
#endif
