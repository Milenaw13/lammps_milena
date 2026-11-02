.. index:: compute spin/heat/flux

compute spin/heat/flux command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID all spin/heat/flux group-A group-B

* ID is documented in :doc:`compute <compute>` command; the group must be *all*
* spin/heat/flux = style name of this compute command
* group-A = ID of the group whose spin motion is used (see below)
* group-B = ID of the group on the other side of the interface

Examples
""""""""

.. code-block:: LAMMPS

   region   left  block INF 10.0 INF INF INF INF
   region   right block 10.0 INF INF INF INF INF
   group    left  region left
   group    right region right

   compute  qspin all spin/heat/flux right left
   variable qx equal c_qspin/(ly*lz)
   fix      gk all ave/correlate 1 2000 2000 v_qx type auto file Qspin.acf ave running

   # 3 atom types: 6 channels, all auto- and cross-correlations
   fix      gkc all ave/correlate 1 2000 2000 c_qspin[*] type full file Qchan.acf ave running

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the heat flow carried by the
magnetic spins across the interface between two groups of atoms, A and
B, due to the exchange interaction.  Following :ref:`(Zhou)
<Zhou2020>`, it is

.. math::

   Q_{AB} = \sum_{i \in A} \sum_{j \in B} J(r_{ij})
            \left( \vec{\omega}_i \times \vec{s}_i \right) \cdot \vec{s}_j
          = \sum_{i \in A} \sum_{j \in B} J(r_{ij})\,
            \frac{d\vec{s}_i}{dt} \cdot \vec{s}_j

where :math:`\vec{s}_i` is the unit vector along the spin of atom
:math:`i`, :math:`\vec{\omega}_i` is its precession vector (the
magnetic force of all magnetic pair styles and fixes, divided by
:math:`\hbar`), so that :math:`d\vec{s}_i/dt = \vec{\omega}_i \times
\vec{s}_i`, and :math:`J(r_{ij})` is the exchange interaction between
atoms :math:`i` and :math:`j` at distance :math:`r_{ij}`.  The sum
includes all periodic images of atom :math:`j` within the cutoff.

:math:`J(r)` and its cutoff are taken from the exchange pair styles that
are defined in the input, :doc:`spin/exchange <pair_spin_exchange>` and
:doc:`spin/exchange/spline <pair_spin_exchange_spline>`, also when they
are used as sub-styles of :doc:`pair_style hybrid/overlay
<pair_hybrid>`.  For each pair of atom types, the exchange style (or
styles) assigned to that pair with :doc:`pair_coeff <pair_coeff>` and
its cutoff for that pair are used.  :math:`J(r)` follows the convention
of these pair styles, in which every pair of atoms is counted once in
the energy.  No exchange parameters are given to this compute.

The interface between the groups is usually a plane normal to one of
the box directions, e.g. :math:`x`, with group A on one side and group B
on the other.  With periodic boundaries, the two groups share two such
planes and :math:`Q_{AB}` is the total heat flow through both.  The
heat flux (heat flow per area) is :math:`Q_{AB}/A_{\perp}`, where
:math:`A_{\perp}` is the area of the interface; the compute does not
divide by the area, so that this can be done with an :doc:`equal-style
variable <variable>`, as in the example above.  This is the same as
for :doc:`compute heat/flux <compute_heat_flux>`, which does not
divide by the volume.

In :ref:`(Zhou) <Zhou2020>`, eq. (13), group A is the right side
:math:`R` of the interface and group B is the left side :math:`L`, and
the magnon thermal conductivity is obtained from the Green-Kubo
relation

.. math::

   \kappa_{magnon} = \frac{V}{k_B T^2} \int_0^\infty
     \left\langle q(t)\, q(0) \right\rangle dt,
   \qquad q = Q_{AB}/A_{\perp}

where :math:`V` is the volume of the system.  The time correlation
function can be computed with :doc:`fix ave/correlate
<fix_ave_correlate>`.  The order of the two groups changes the sign of
:math:`Q_{AB}`, which matters for correlations with the lattice heat
flux, but not for the autocorrelation.  Because only the spin motion of
group A enters, :math:`Q_{AB}` and :math:`-Q_{BA}` are equal only on
average, not at every time step.

The compute also splits :math:`Q_{AB}` into exchange channels, one for
each pair of atom types :math:`\{I, J\}` (the order of the two types
does not matter):

.. math::

   Q_{AB}^{IJ} = \sum_{i \in A} \sum_{j \in B,\ \{t_i, t_j\} = \{I, J\}}
                 J(r_{ij})\,\frac{d\vec{s}_i}{dt} \cdot \vec{s}_j,
   \qquad \sum_{I \le J} Q_{AB}^{IJ} = Q_{AB}

where :math:`t_i` is the type of atom :math:`i`.  For :math:`N` atom
types there are :math:`N(N+1)/2` channels, in the order (1,1), (1,2),
..., (1,N), (2,2), ..., (N,N); pairs of types without an exchange
interaction give zero.  This shows which exchange coupling carries the
heat.  The thermal conductivity splits into

.. math::

   \kappa_{magnon} = \sum_{k,l} \kappa_{kl}, \qquad
   \kappa_{kl} = \frac{V}{k_B T^2} \int_0^\infty
     \left\langle q_k(t)\, q_l(0) \right\rangle dt

where :math:`k` and :math:`l` run over the channels and :math:`q_k` is
the heat flux of channel :math:`k`.  The cross terms :math:`k \neq l`
are in general not small, because the channels share spins: the motion
of one spin enters every channel it takes part in.  Therefore the sum
of the autocorrelations of the channels is *not* the thermal
conductivity; use the correlation of the total heat flux, or all
correlations between channels (:doc:`fix ave/correlate
<fix_ave_correlate>` with *type full*, see the example above).  A
useful measure of the contribution of channel :math:`k` is
:math:`\sum_l (\kappa_{kl} + \kappa_{lk})/2`: these contributions add
up to :math:`\kappa_{magnon}` and can be negative.  The directory
``examples/SPIN/spin_heat_flux`` contains an input and a Python script
that computes them.  Contributions of single channels need longer runs
to converge than the total thermal conductivity.

The lattice (phonon) heat flux is not included; it can be computed
with :doc:`compute heat/flux <compute_heat_flux>`.  Other magnetic pair
styles (e.g. *spin/dmi*, *spin/neel*) contribute to
:math:`\vec{\omega}_i`, but their own coupling between the two groups is
not included in :math:`Q_{AB}`; LAMMPS prints a warning if such pair
styles are defined.

The value is computed from the precession vectors of the current time
step, which are available after the forces have been computed.  Use the
compute in output commands that are invoked at the end of a time step
(:doc:`thermo <thermo_style>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`fix ave/correlate <fix_ave_correlate>`, :doc:`fix print
<fix_print>`, :doc:`variables <variable>` used by these).

----------

Output info
"""""""""""

This compute calculates a global scalar, the heat flow :math:`Q_{AB}`,
and a global vector of length :math:`N(N+1)/2` (:math:`N` = number of
atom types) with the heat flow of each channel, in the order given
above.  The elements of the vector add up to the scalar.  These values
can be accessed by any command that uses global scalar or vector values
from a compute as input.  See the :doc:`Howto output <Howto_output>`
page for an overview of LAMMPS output options.

The scalar and vector values are "extensive".  They are in units of
energy per time, i.e. eV/ps for :doc:`metal units <units>`.

Restrictions
""""""""""""

This compute is part of the SPIN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This compute requires atom style *spin*, the group *all*, and pair style
*spin/exchange* or *spin/exchange/spline* (also as sub-styles of
*hybrid* or *hybrid/overlay*).  Pair style *hybrid/scaled* is not
supported.

Related commands
""""""""""""""""

:doc:`compute spin <compute_spin>`,
:doc:`compute heat/flux <compute_heat_flux>`,
:doc:`pair_style spin/exchange <pair_spin_exchange>`,
:doc:`pair_style spin/exchange/spline <pair_spin_exchange_spline>`,
:doc:`fix ave/correlate <fix_ave_correlate>`

Default
"""""""

none

----------

.. _Zhou2020:

**(Zhou)** Zhou, Tranchida, Ge, Murthy, Fisher, Phys Rev B, 101, 224303 (2020).
