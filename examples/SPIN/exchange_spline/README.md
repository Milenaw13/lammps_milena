# pair_style spin/exchange/spline example

Exchange interaction J(r) given as a table and interpolated with a cubic spline.

- `write_exchange_table.py`: writes table files from any function J(r); run as a
  script it writes `Fe_bethe_slater.table` (Bethe-Slater J(r) of bcc iron from
  `../iron`, shifted so that J and dJ/dr vanish at the 5 Angstrom cutoff).
  Use it as a module to tabulate your own fitted J(r) (see the docstring for
  units and the factor 2 for codes that count each pair twice, e.g. TB2J).
- `in.spin.iron_spline`: spin-lattice dynamics of bcc iron (EAM + spline exchange)
  without thermostat; the total energy is conserved.
- `in.dimer_scan` + `check_exchange_table.py`: two atoms with parallel spins at
  distances r; compares the LAMMPS energy (-J) and force (dJ/dr) with the table
  spline and with a finite difference of the energy:

```
lmp -in in.dimer_scan -var table Fe_bethe_slater.table -var key FE_BS \
    -var rc 5.0 -var r0 1.5 -var dr 0.0013 -var n 2690
python3 check_exchange_table.py Fe_bethe_slater.table FE_BS
```
