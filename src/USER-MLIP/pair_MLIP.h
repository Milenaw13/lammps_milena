/* ----------------------------------------------------------------------
 *   This is the MLIP-LAMMPS interface
 *   MLIP is a software for Machine Learning Interatomic Potentials
 *   MLIP is released under the "New BSD License", see the LICENSE file.
 *   Contributors: Evgeny Podryabinkin

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   LAMMPS is distributed under a GNU General Public License
   and is not a part of MLIP.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Evgeny Podryabinkin (Skoltech)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(mlip,PairMLIP)

#else

#ifndef LMP_PAIR_MLIP_H
#define LMP_PAIR_MLIP_H

#include <stdio.h>
#include <string>
#include <map>
#include "pair.h"

struct CallbackWrapper
{
    void* pair_ptr = nullptr;
    void  callback_call(double*);
};

extern void* MLIP_init(std::map<std::string, std::string>, double&, MPI_Comm&, void(*CallbackComm)(void*, double*));
extern void  MLIP_CalcCfgForLammps(void*, void*, int, int, int*, int*, int**, double*, double**, int*, bool, double**, double&, double*, double*, double**);
extern void  MLIP_finalize(void*);


namespace LAMMPS_NS {

class PairMLIP : public Pair {
 public:
  void* mlip_wrapper_ptr = nullptr;
  double cutoff;
  double* comm_arr;

  PairMLIP(class LAMMPS *);
  virtual ~PairMLIP();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);

 protected:
  bool inited;
  double cutoffsq;
  double **site_virial9;       // 9-comp per-atom virial buffer when only vatom is requested
  int nmax_site_virial9;
  std::map<std::string, std::string> mlip_settings;

  void allocate();
  
  int pack_reverse_comm(int, int, double*);
  void unpack_reverse_comm(int, int*, double*);

 public:
  static void CallbackComm(void*, double*);
  void reverse_comm_callback(double*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
