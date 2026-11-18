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

#ifdef FIX_CLASS
// clang-format off
FixStyle(langevin/spin,FixLangevinSpin);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_SPIN_H
#define LMP_FIX_LANGEVIN_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinSpin : public Fix {
 public:
  int tdamp_flag, temp_flag;    // damping and temperature flags

  FixLangevinSpin(class LAMMPS *, int, char **);
  ~FixLangevinSpin() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void add_tdamping(double *, double *);    // add transverse damping
  void add_temperature(double[3]);
  void compute_single_langevin(int, double *, double *);
  void tally_spin_update(int, double *, double *);
  double compute_scalar() override;
  double compute_vector(int) override;

 protected:
  double alpha_t;       // transverse mag. damping
  double dts;           // magnetic timestep
  double temp;          // spin bath temperature
  double D, sigma;      // bath intensity var.
  double gil_factor;    // gilbert's prefactor

  int tallyflag;             // 1 if energy exchange with the bath is tallied
  double energy;             // cumulative energy removed from the local spins
  double energy_avg[2];      // same, from the average power of the damping [0]
                             // and of the random field [1]
  int ilast;                 // spin of the last call of compute_single_langevin
  double fmlast[3];          // its conservative precession vector (rad/ps)

  void tally_average(double *, double *);

  int nlevels_respa;
  class RanMars *random;
  int seed;
};

}    // namespace LAMMPS_NS

#endif
#endif
