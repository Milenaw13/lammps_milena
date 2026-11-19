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

   Modified (branch ferh-production-mlip, 2026-11): per-atom centroid
   virial (9 components, many-body form) for compute centroid/stress/atom
   and compute heat/flux; CENTROID_AVAIL. Needs the mlip-3 patch in
   lib/mlip/mlip-3_site_virial.patch. See lib/mlip/README_ferh.md.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "pair_MLIP.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"


using namespace LAMMPS_NS;


#define MAXLINE 1024


void PairMLIP::CallbackComm(void* p_pair, double* p_arr)
{
  PairMLIP* p_pair_mlip = static_cast<PairMLIP*>(p_pair);
  p_pair_mlip->reverse_comm_callback(p_arr);
}


/* ---------------------------------------------------------------------- */

PairMLIP::PairMLIP(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;

  single_enable = 0;

  inited = false;
  allocated = 0;

  comm_reverse = 3;

  // per-atom virial: many-body (Fan) form, 9 components -> cvatom
  centroidstressflag = CENTROID_AVAIL;
  site_virial9 = nullptr;
  nmax_site_virial9 = 0;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairMLIP::~PairMLIP()
{
  if (copymode) return;

  memory->destroy(site_virial9);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
  
  if (inited && mlip_wrapper_ptr != nullptr) 
    MLIP_finalize(mlip_wrapper_ptr);

  mlip_wrapper_ptr = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairMLIP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  {
    double lattice[9];
    lattice[0] = domain->xprd;
    lattice[1] = 0.0;
    lattice[2] = 0.0;
    lattice[3] = domain->xy;
    lattice[4] = domain->yprd;
    lattice[5] = 0.0;
    lattice[6] = domain->xz;
    lattice[7] = domain->yz;
    lattice[8] = domain->zprd;
                          
    double en = 0.0;
    double virstr[9];
    double *p_site_en = NULL;
    double **p_site_virial = NULL;

    if (eflag_atom) p_site_en = &eatom[0];

    // the library writes the 9-component centroid virial (ghosts included);
    // use cvatom if requested, otherwise a private buffer that is
    // symmetrized into vatom below
    const int nall = atom->nlocal + atom->nghost;
    if (cvflag_atom) {
      p_site_virial = cvatom;
    } else if (vflag_atom) {
      if (atom->nmax > nmax_site_virial9) {
        memory->destroy(site_virial9);
        nmax_site_virial9 = atom->nmax;
        memory->create(site_virial9,nmax_site_virial9,9,"pair:site_virial9");
      }
      for (int i = 0; i < nall; i++)
        for (int k = 0; k < 9; k++) site_virial9[i][k] = 0.0;
      p_site_virial = site_virial9;
    }
                        
//                                printf("PairMLIP::compute begin\n");
//                                //printf("\tPairMLIP:: CutOff = %f\n", cutoff);
//                                printf("\tPairMLIP:: inum = %d\n", list->inum);
//                                printf("\tPairMLIP:: x[0] = %e %e %e\n", x[0][0], x[0][1], x[0][2]);
//                                printf("\tPairMLIP:: lat = %e %e %e   %e %e %e   %e %e %e\n", lattice[0], lattice[1], lattice[2], lattice[3], lattice[4], lattice[5], lattice[6], lattice[7], lattice[8]);
    
    MLIP_CalcCfgForLammps(mlip_wrapper_ptr,
                          this, 
                          list->inum,              // input parameter: number of neighborhoods
                          atom->nghost,
                          list->ilist,                             // input parameter: 
                          list->numneigh,                     // input parameter: number of neighbors in each neighborhood (inum integer numbers)
                          list->firstneigh,                  // input parameter: pointer to the first neighbor
                          lattice,                  // input parameter: lattice (9 double numbers)
                          atom->x,                         // input parameter: array of coordinates of atoms
                          atom->type,                     // input parameter: array of atom types (inum of integer numbers)
                          (comm->nprocs > 1),
                          atom->f,                    // output parameter: forces on atoms (cartesian, n x 3 double numbers)
                          en,                    // output parameter: summ of site energies 
                          virstr,            // output parameter: stresses in eV (9 double numbers)
                          p_site_en,       // output parameter: array of site energies (inum double numbers). if =nullptr while call no site energy calculation is done
                          p_site_virial);                             
    

//                                  printf("\tPairMLIP:: en = %f\n", en);
//                                  printf("\tPairMLIP:: inum = %d\n", list->inum);
//                                  printf("\tPairMLIP:: x[0] = %f %f %f\n", x[0][0], x[0][1], x[0][2]);
//                                  printf("\tPairMLIP:: f[0] = %f %f %f\n", f[0][0], f[0][1], f[0][2]);
//                                  printf("\tPairMLIP:: lat = %e %e %e   %e %e %e   %e %e %e\n", lattice[0], lattice[1], lattice[2], lattice[3], lattice[4], lattice[5], lattice[6], lattice[7], lattice[8]);
//                                  printf("PairMLIP::compute end\n");

    if (vflag_atom && p_site_virial != nullptr) {
      double **w = p_site_virial;
      for (int i = 0; i < nall; i++) {
        vatom[i][0] += w[i][0];
        vatom[i][1] += w[i][1];
        vatom[i][2] += w[i][2];
        vatom[i][3] += 0.5*(w[i][3] + w[i][6]);
        vatom[i][4] += 0.5*(w[i][4] + w[i][7]);
        vatom[i][5] += 0.5*(w[i][5] + w[i][8]);
      }
    }

    if (eflag_global) eng_vdwl += en;
    if (vflag_fdotr) 
      virial_fdotr_compute();
    else if (vflag)
    {
      virial[0] += virstr[0];
      virial[1] += virstr[4];
      virial[2] += virstr[8];
      virial[3] += (virstr[1] + virstr[3]) / 2;
      virial[4] += (virstr[2] + virstr[6]) / 2;
      virial[5] += (virstr[5] + virstr[7]) / 2;
    }
  }

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMLIP::allocate()
{
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  allocated = 1;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMLIP::settings(int narg, char **arg)
{
  if (narg < 1) 
    error->all(FLERR, "Illegal pair_style command");

  mlip_settings.clear();
  mlip_settings["mlip"] = "true";
  mlip_settings["calculate_efs"] = "true";
  mlip_settings["extrapolation_control:save_extrapolative_to"] = "preselected.cfg";

  for (int i=0; i<narg; i++)
  {
    if (strlen(arg[i]) > 999)
      error->all(FLERR, "MLIP pair_style argument is too long");
    std::string s(arg[i]);
    int count = 0;
    for (int j=0; j<s.size(); j++)
       if (s[j] == '=') count++;
    if (count != 1)
      error->all(FLERR, "MLIP pair_style argument must contain only one \"=\" simbol");

    if (s.substr(0, s.find('=')) == "load_from")
      mlip_settings["mlip:load_from"] = s.substr(s.find('=')+1);
    else
    {
      mlip_settings[s.substr(0, s.find('='))] = s.substr(s.find('=')+1);
      mlip_settings["extrapolation_control:" + s.substr(0, s.find('='))] = s.substr(s.find('=')+1);
    }
  }
}

/* ----------------------------------------------------------------------
   set flags for type pairs
------------------------------------------------------------------------- */

void PairMLIP::coeff(int narg, char **arg)
{
  if (strcmp(arg[0],"*") || strcmp(arg[1],"*") )
    error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMLIP::init_style()
{
  if (force->newton_pair != 1)
      error->all(FLERR, "Pair style MLIP requires Newton pair on");

  if (inited)
    MLIP_finalize(mlip_wrapper_ptr);
  
  mlip_wrapper_ptr = MLIP_init(mlip_settings, cutoff, world, CallbackComm);

  cutoffsq = cutoff*cutoff;
  int n = atom->ntypes;
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      cutsq[i][j] = cutoffsq;

  inited = true;

  //int irequest = neighbor->request(this,instance_me);
  //neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMLIP::init_one(int i, int j)
{
  return cutoff;
}


void PairMLIP::reverse_comm_callback(double* p_arr)
{
  comm_arr = p_arr;
  comm->reverse_comm(this);
}


int PairMLIP::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) 
  {
    buf[m++] = comm_arr[3*i+0];
    buf[m++] = comm_arr[3*i+1];
    buf[m++] = comm_arr[3*i+2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairMLIP::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) 
  {
    j = list[i];
    comm_arr[3*j+0] += buf[m++];
    comm_arr[3*j+1] += buf[m++];
    comm_arr[3*j+2] += buf[m++];
  }
}
