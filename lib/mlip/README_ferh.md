# MLIP-3 interface with per-atom centroid virial (branch ferh-production-mlip)

This branch = ferh-production (spline exchange, spin heat computes, midpoint
fix nve/spin) + the MLIP-3 LAMMPS interface (`src/USER-MLIP`, from
https://gitlab.com/ivannovikov/interface-lammps-mlip-3, commit 68f02dd) with a fix of
the per-atom virial. MLIP is not a LAMMPS package: keep this branch out of pull requests.

## What was changed and why

Stock mlip-3 (7fe598d) never fills the per-atom virial, so `compute stress/atom`,
`compute centroid/stress/atom` and `compute heat/flux` miss the whole MTP virial term.

- `mlip-3_site_virial.patch` (for mlip-3 7fe598d, `src/external/interface.cpp`):
  fills the 9-component per-atom virial in the many-body (Fan) form,
  W_j += (x_i - x_j) (x) dU_i/dx_j, for local and ghost neighbours.
- `src/USER-MLIP/pair_MLIP.cpp/.h`: `centroidstressflag = CENTROID_AVAIL`; passes
  `cvatom` (9 columns) to the library for centroid/stress/atom; for stress/atom a private
  9-column buffer is symmetrized into `vatom` (6 columns). `ev_setup` -> `ev_init`.

Use `compute heat/flux ke pe <centroid/stress/atom compute>`. The symmetric
stress/atom is not exact for a many-body potential (about 7 % error in the virial
term of the test below).

Validation (local, 2026-11-13, `check/`): sum of the per-atom virials = global virial
(1e-14, 1 and 4 ranks); free FeRh cluster with MTP + spin/exchange/spline + midpoint
fix nve/spin: dG/dt = J(heat/flux, centroid) + J_prec(spin/heat/current) to O(dt^2)
(1.3e-4 at dt 0.1 fs, 3.2e-5 at 0.05 fs); 1 vs 4 ranks 1e-14.

## Build on Snellius (traditional make, as for the earlier lmp_mpi builds)

The two parts (mlip-3 library and LAMMPS) MUST be built with the same modules, in a
shell where the modules are loaded. Do every step in ONE shell:

```bash
module purge
module load 2025 GCC/14.2.0 OpenMPI/5.0.7-GCC-14.2.0 OpenBLAS/0.3.29-GCC-14.2.0
which mpicc mpicxx                      # must be the OpenMPI/5.0.7 module, not /usr/bin
```

### 1. LAMMPS source: this branch

In the existing clone (e.g. `~/lammps-ferh`, branch ferh-production). An earlier
`install.sh` copied the unpatched interface into `src/`; remove that first, otherwise
`git switch` refuses (untracked files in the way) or the old pair_MLIP stays in use:

```bash
cd ~/lammps-ferh/src
make no-user-mlip                       # removes src/pair_MLIP.* installed earlier
git status --short                      # untracked src/USER-MLIP/ from install.sh? then:
mv USER-MLIP ~/USER-MLIP.old            #   (keep a copy, do not delete blindly)
cd ..
git fetch origin
git switch ferh-production-mlip         # tracks origin/ferh-production-mlip
git log --oneline -2
grep -c CENTROID_AVAIL src/USER-MLIP/pair_MLIP.cpp      # must be 1
```

### 2. mlip-3 library with the patch (fresh clone, so no old patch is mixed in)

```bash
git clone https://gitlab.com/ashapeev/mlip-3.git ~/mlip-3-ferh
cd ~/mlip-3-ferh
git checkout 7fe598d
git apply ~/lammps-ferh/lib/mlip/mlip-3_site_virial.patch
git status --short                      # " M src/external/interface.cpp"
./configure --compiler=gnu              # MPI is autodetected from mpicc in PATH (step 0!)
grep -E "^USE_MPI|^CXX_MPI" make/config.mk   # must say USE_MPI = 1, CXX_MPI = mpicxx
make -j 16 libinterface CPPFLAGS_MPI="-DMLIP_MPI -DOMPI_SKIP_MPICXX"
```

Check the library BEFORE building LAMMPS:

```bash
nm -C lib/lib_mlip_interface.a | grep -c "MPI::"                       # must be 0
nm -C lib/lib_mlip_interface.a | grep " T MLIP_init" | grep -c ompi_communicator_t   # must be 1
```

If the second check gives 0, the library was built serially (configure did not see
mpicc): load the modules, `make clean`, `./configure --compiler=gnu` again.

### 3. Build LAMMPS

```bash
cp ~/mlip-3-ferh/lib/lib_mlip_interface.a ~/lammps-ferh/lib/mlip/
cp ~/lammps-ferh/lib/mlip/Makefile.lammps.template ~/lammps-ferh/lib/mlip/Makefile.lammps
cd ~/lammps-ferh/src
make yes-user-mlip yes-spin yes-manybody yes-molecule   # + any other package you use
make ps | grep YES
grep -c CENTROID_AVAIL pair_MLIP.cpp                    # installed copy patched: 1
make -j 16 mpi
./lmp_mpi -h | grep -E " mlip |spin/exchange/spline|spin/heat/current"
```

Do NOT run `install.sh` / `preinstall.sh` of interface-lammps-mlip-3 on this tree: they
copy the unpatched pair_MLIP from the interface repository over the patched one.

### 4. Quick check on Snellius (a few minutes, small allocation or interactive)

```bash
cd check                                # copy lib/mlip/check somewhere to run it
ln -s <path to pot.mtp> pot.mtp
srun -n 1 ~/lammps-ferh/src/lmp_mpi -in in.sumrule -log log.np1 | grep -E "^(GLOBAL|CENTROID)"
srun -n 4 ~/lammps-ferh/src/lmp_mpi -in in.sumrule -log log.np4 | grep -E "^(GLOBAL|CENTROID)"
```

The GLOBAL numbers and the first 6 CENTROID numbers must agree to ~12 digits, on 1 and
4 ranks (the local values: GLOBAL -1.38946810273932e10 -2.20720103098402e10 ...).
Optional exactness test: `in.cluster_spin` with `-var start read -var integ midpoint
-var dt 0.0001 -var nrun 400 -var tag mid_dt1_np1` and with `-var dt 0.00005 -var nrun
800 -var tag mid_dt05_np1`, then `python3 analyse_cluster_spin.py` (edit the list of tags
at the bottom to the runs you made); "centroid + J_prec" must drop ~4x with dt.

## Link errors and their cause

| Error | Cause | Fix |
|---|---|---|
| undefined reference to `MPI::Comm::Comm()`, `MPI::Datatype::Free()`, `ompi_mpi_cxx_op_intercept` | library compiled against an MPI with C++ bindings (other/older OpenMPI headers) | same modules for both, rebuild library with `-DOMPI_SKIP_MPICXX` |
| undefined reference to `MLIP_init(..., int&, ...)` | serial library (USE_MPI = 0) linked into lmp_mpi | modules loaded before `./configure`, check `make/config.mk` |
| undefined reference to `MLIP_init(..., ompi_communicator_t*&, ...)` | MPI library linked into a serial LAMMPS (`make serial`) | build `make mpi` |
| multiple definition of `MPI_Init(int*, char***)` | serial library + serial LAMMPS: two sets of MPI stubs | do not use; build both with MPI |
| The Open MPI wrapper compiler was unable to find the specified compiler | wrapper from another MPI installation than the loaded compiler | `which mpicxx`; module purge + load |

## Local build (WSL, CMake)

`~/SANDBOXES/lammps_sandbox_claude_ferh_production`: library in `mlip/mlip-3` (patched
clone of 7fe598d), built with the miniforge MPI:

```bash
cd mlip/mlip-3
OMPI_CXX=g++ OMPI_CC=gcc OMPI_FC=gfortran make -j8 libinterface USE_MPI=1 \
  CC_MPI=$HOME/miniforge3/bin/mpicc CXX_MPI=$HOME/miniforge3/bin/mpicxx \
  FC_MPI=$HOME/miniforge3/bin/mpif90 CPPFLAGS_MPI="-DMLIP_MPI -DOMPI_SKIP_MPICXX"
```

LAMMPS (worktree lammps-mlip): symlinks `src/pair_MLIP.{cpp,h} -> USER-MLIP/` (in
.git/info/exclude), then the usual CMake configure with
`-D CMAKE_CXX_STANDARD_LIBRARIES="<abs path>/mlip/mlip-3/lib/lib_mlip_interface.a -lgfortran"`.
