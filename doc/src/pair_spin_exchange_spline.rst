.. index:: pair_style spin/exchange/spline

pair_style spin/exchange/spline command
=======================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style spin/exchange/spline cutoff

* cutoff = global cutoff for the exchange interaction (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style spin/exchange/spline 7.0
   pair_coeff * * exchange 7.0 J_FeRh.table FE_FE
   pair_coeff 1 2 exchange 7.0 J_FeRh.table FE_RH offset yes

   pair_style hybrid/overlay eam/alloy spin/exchange/spline 7.0
   pair_coeff * * eam/alloy FeRh.eam.alloy Fe Rh
   pair_coeff 1 1 spin/exchange/spline exchange 7.0 J_FeRh.table FE_FE offset yes
   pair_coeff 1 2 spin/exchange/spline exchange 7.0 J_FeRh.table FE_RH offset yes

Description
"""""""""""

.. versionadded:: TBD

Style *spin/exchange/spline* computes the exchange interaction between
pairs of magnetic spins with the same energy, precession vectors, and
forces as the :doc:`spin/exchange <pair_spin_exchange>` pair style:

.. math::

   H_{ex} = -\sum_{i<j}^N J(r_{ij}) \,\vec{s}_i \cdot \vec{s}_j

where :math:`\vec{s}_i` and :math:`\vec{s}_j` are unit vectors along the
magnetic spins of two atoms and :math:`r_{ij}` is the distance between
them.  The sum counts each pair of atoms *once*.  The only difference to
*spin/exchange* is the distance dependence of :math:`J(r)`: instead of
the Bethe-Slater function, :math:`J(r)` is a cubic spline through
:math:`(r, J)` points read from a tabulated file.  This allows exchange
interactions that change sign several times with distance, as is common
for interactions computed from first principles.

As for *spin/exchange*, each spin :math:`i` precesses around

.. math::

   \vec{\omega}_{i} = \frac{1}{\hbar} \sum_{j} J(r_{ij})\,\vec{s}_{j}

and, for spin-lattice calculations (see :doc:`fix nve/spin
<fix_nve_spin>` with the *lattice moving* option), atom :math:`i`
experiences the force

.. math::

   \vec{F}_{i} = \sum_{j} \frac{d J(r_{ij})}{d r_{ij}}
   \left( \vec{s}_{i}\cdot \vec{s}_{j} \right) \vec{e}_{ij}

which couples the magnetic and the lattice degrees of freedom.  The
energy :math:`J(r)` and its derivative :math:`dJ/dr` are evaluated
from the same cubic polynomial of the spline, so forces are exactly
consistent with the energy.

The following coefficients must be defined for each pair of atom types
via the :doc:`pair_coeff <pair_coeff>` command as in the examples above:

* the word *exchange*
* :math:`R_c` = cutoff of this pair (distance units), must not be larger
  than the global cutoff
* file name of the table
* keyword of the table section in the file

followed by the optional keyword *offset* with the value *yes* or *no*
(see below).

----------

**Table file format**\ :

The file can contain several sections, for example one for each pair of
atom types.  Blank lines and text following a "#" character are
ignored.  Each section has the format:

.. parsed-literal::

   # J(r) for Fe-Fe pairs, r in Angstrom, J in eV   (comment)
   FE_FE                                            (keyword)
   N 601                                            (number of points)

   1 1.000 -1.25000e-02                             (index, r, J)
   2 1.010 -1.24100e-02
   ...
   601 7.000 0.0

The line after the keyword must be "N" followed by the number of
points, at least 3.  Each of the following *N* lines contains an index
(ignored), a distance :math:`r` (distance units), and the value
:math:`J(r)` (energy units).  The distances must be strictly increasing.
They do not need to be equally spaced, but equally spaced points allow
a faster lookup of the spline interval during a run.

The table must describe a *smooth* curve, typically obtained by fitting
the raw exchange data (e.g. from first-principles calculations) and
sampling the fitted curve densely.  The spline passes exactly through
all points, so a table of raw data with scatter between neighbor shells
leads to strongly oscillating :math:`J(r)` and :math:`dJ/dr`.
The *examples/SPIN/exchange_spline* folder contains a Python script that
writes a table from any fitted function.

The last point of the table must be at the cutoff :math:`R_c` of the
pair, and :math:`J(r)` must go smoothly to zero there:
:math:`J(R_c) = 0` and :math:`dJ/dr(R_c) = 0`.  Otherwise the energy
and the forces would jump when atoms cross the cutoff, and the total
energy would not be conserved.  LAMMPS stops with an error if
:math:`|J(R_c)|` exceeds :math:`10^{-4}` times the largest
:math:`|J|` in the table, or if the slope at :math:`R_c` estimated from
the last three points exceeds :math:`10^{-3}` times the largest slope in
the table.  When building the spline, :math:`dJ/dr(R_c) = 0` is used as
boundary condition.

For distances shorter than the first point of the table,
:math:`J(r)` is extrapolated linearly with the slope at the first
point.  At the end of each run, LAMMPS prints a warning with the number
of pair distances for which this happened.  If this occurs often, the
table should start at a shorter distance.

.. note::

   The sign and prefactor convention of :math:`J` must match the
   Hamiltonian above, where each pair is counted once and the spins are
   unit vectors.  Several codes that compute exchange interactions from
   first principles (e.g. TB2J) write the Hamiltonian as a sum over
   :math:`i \neq j`, which counts each pair twice.  Their values must
   then be multiplied by 2 for this pair style.

----------

**Offsetting magnetic forces and energies**\ :

The *offset* keyword works exactly as for the :doc:`spin/exchange
<pair_spin_exchange>` pair style: with *offset yes*,
:math:`\vec{s}_i \cdot \vec{s}_j` is replaced by
:math:`\vec{s}_i \cdot \vec{s}_j - 1` in the energy and the mechanical
forces, so that both vanish for parallel spins.  The precession vectors
are not affected.  As for *spin/exchange*, the value of the last
pair_coeff command that sets *offset* applies to all pairs.

----------

Restrictions
""""""""""""

This pair style is part of the SPIN package.  It is only enabled if
LAMMPS was built with that package, and requires the atom style "spin"
and metal units.  See the :doc:`Build package <Build_package>` page for
more info.

The tables are not stored in restart files; the pair_style and
pair_coeff commands must be given again after a
:doc:`read_restart <read_restart>` command.

Related commands
""""""""""""""""

:doc:`pair_style spin/exchange <pair_spin_exchange>`,
:doc:`atom_style spin <atom_style>`, :doc:`pair_coeff <pair_coeff>`,
:doc:`fix nve/spin <fix_nve_spin>`, :doc:`pair_style table <pair_table>`

Default
"""""""

The default *offset* keyword value is *no*.
