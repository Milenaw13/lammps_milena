.. index:: compute spin/heat/current

compute spin/heat/current command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID spin/heat/current keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* spin/heat/current = style name of this compute command
* zero or more keyword/value pairs may be appended
* keyword = *mech*

  .. parsed-literal::

       *mech* value = *yes* or *no*
         *yes* = also output the heat current of the mechanical exchange forces
         *no* = do not output it

Examples
""""""""

.. code-block:: LAMMPS

   compute  jspin all spin/heat/current
   variable jx equal c_jspin[1]
   variable jy equal c_jspin[2]
   variable jz equal c_jspin[3]
   fix      gk all ave/correlate 1 2000 2000 v_jx v_jy v_jz type auto file Jspin.acf ave running

   compute  jspin2 all spin/heat/current mech yes

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the heat current carried by the
magnetic spins through the exchange interaction, integrated over the
volume of the system.  It is the spin counterpart of :doc:`compute
heat/flux <compute_heat_flux>` and can be used in the same way to
compute the magnon thermal conductivity with the Green-Kubo relation.

The magnetic energy of each atom is half of the exchange energy of all
its pairs,

.. math::

   e_i = -\frac{1}{2} \sum_{j} J(r_{ij})\, \vec{s}_i \cdot \vec{s}_j

(with :math:`\vec{s}_i \cdot \vec{s}_j - 1` if the *offset* option of the
pair style is used).  The heat current is the time derivative of the
energy-weighted position :math:`\sum_i \vec{r}_i e_i` due to the
precession of the spins :ref:`(Hardy) <Hardy1963>`:

.. math::

   \vec{J}_{prec} = \frac{1}{2} \sum_i \sum_{j} \vec{r}_{ij}\, J(r_{ij})\,
     \frac{d\vec{s}_i}{dt} \cdot \vec{s}_j,
   \qquad \frac{d\vec{s}_i}{dt} = \vec{\omega}_i \times \vec{s}_i

where :math:`\vec{r}_{ij} = \vec{r}_i - \vec{r}_j`, :math:`\vec{s}_i` is
the unit vector along the spin of atom :math:`i`, :math:`\vec{\omega}_i`
is its precession vector (the magnetic force of all magnetic pair styles
and fixes, divided by :math:`\hbar`), and :math:`J(r_{ij})` is the
exchange interaction between atoms :math:`i` and :math:`j`.  The sum
over :math:`j` includes all periodic images within the cutoff.  The
pair term :math:`J(r_{ij})\,(d\vec{s}_i/dt) \cdot \vec{s}_j` is the
energy per time that the spin of atom :math:`i` takes from the spin of
atom :math:`j`; it is the same quantity that :doc:`compute
spin/heat/flux <compute_spin_heat_flux>` sums over the pairs across an
interface between two groups :ref:`(Zhou) <Zhou2020b>`.  Here every
pair is weighted with its distance vector instead, which is equivalent
to averaging the heat flow through a plane over all positions of the
plane.  The result does not depend on the choice of interfaces and is
the same kind of quantity as the lattice heat current of :doc:`compute
heat/flux <compute_heat_flux>`.

:math:`J(r)` and its cutoff are taken from the exchange pair styles that
are defined in the input, :doc:`spin/exchange <pair_spin_exchange>` and
:doc:`spin/exchange/spline <pair_spin_exchange_spline>`, also when they
are used as sub-styles of :doc:`pair_style hybrid/overlay
<pair_hybrid>`.  For each pair of atom types, the exchange style (or
styles) assigned to that pair with :doc:`pair_coeff <pair_coeff>` and
its cutoff for that pair are used.  No exchange parameters are given to
this compute.  All exchange interactions that act on the spins must be
included, otherwise the result is not a heat current; LAMMPS prints a
warning if other magnetic pair styles (e.g. *spin/dmi*, *spin/neel*)
are defined, whose energy transfer is not included.

As for :doc:`compute heat/flux <compute_heat_flux>`, the result is not
divided by the volume.  The thermal conductivity from the Green-Kubo
relation is

.. math::

   \kappa_{xx} = \frac{1}{V k_B T^2} \int_0^\infty
     \left\langle J_x(t)\, J_x(0) \right\rangle dt

and similarly for the *y* and *z* directions, where :math:`V` is the
volume of the system (the box volume for group *all*), including atoms
without a spin interaction.  In :doc:`metal units <units>`, a
conductivity in eV/(Angstrom ps K) is multiplied by 1602.18 to obtain
W/(m K).  The time correlation functions can be computed with
:doc:`fix ave/correlate <fix_ave_correlate>`.

The compute also splits :math:`\vec{J}_{prec}` into exchange channels,
one for each pair of atom types :math:`\{I, J\}` (the order of the two
types does not matter), by restricting the double sum to the pairs of
atoms with these types.  For :math:`N` atom types there are
:math:`N(N+1)/2` channels, in the order (1,1), (1,2), ..., (1,N), (2,2),
..., (N,N); pairs of types without an exchange interaction give zero.
The channels add up to :math:`\vec{J}_{prec}`.  As explained on the
:doc:`compute spin/heat/flux <compute_spin_heat_flux>` page, the cross
correlations between the channels are in general not small, so the
contribution of a channel to the thermal conductivity must be computed
from all correlations between the channels; all channels share the
same volume :math:`V`.

When the atoms move, :math:`J(r_{ij})` changes with time and the
exchange forces (the derivative :math:`dJ/dr`) exchange energy between
the spins and the lattice.  This energy is not carried by
:math:`\vec{J}_{prec}`.  It is carried by the work of the mechanical
exchange forces, a heat current of the same form as the virial part of
:doc:`compute heat/flux <compute_heat_flux>`:

.. math::

   \vec{J}_{mech} = \frac{1}{2} \sum_i \sum_{j} \vec{r}_{ij}\,
     \frac{1}{2} \vec{F}_{ij} \cdot (\vec{v}_i + \vec{v}_j),
   \qquad \vec{F}_{ij} = \frac{dJ}{dr}(r_{ij})\,
     (\vec{s}_i \cdot \vec{s}_j)\, \frac{\vec{r}_{ij}}{r_{ij}}

(again with :math:`\vec{s}_i \cdot \vec{s}_j - 1` for the *offset*
option), where :math:`\vec{F}_{ij}` is the exchange force on atom
:math:`i` due to atom :math:`j` and :math:`\vec{v}_i` is the velocity of atom
:math:`i`.  With *mech* = *yes*, it is output in addition to
:math:`\vec{J}_{prec}`.  It is zero if the exchange forces do not act on
the atoms, i.e. if no :doc:`fix nve/spin <fix_nve_spin>` with *lattice
moving* is defined.

.. note::

   :doc:`compute heat/flux <compute_heat_flux>`, used with
   :doc:`compute pe/atom <compute_pe_atom>` and :doc:`compute stress/atom
   <compute_stress_atom>` for all pair styles, already includes the
   magnetic energy :math:`e_i` in its convective part and
   :math:`\vec{J}_{mech}` in its virial part.  The total heat current of
   a spin-lattice simulation is therefore the sum of the heat current
   of :doc:`compute heat/flux <compute_heat_flux>` and
   :math:`\vec{J}_{prec}`.  Do *not* add :math:`\vec{J}_{mech}` to it a
   second time.  :math:`\vec{J}_{mech}` is meant to measure how much
   heat the exchange forces carry, e.g. through its correlations with
   :math:`\vec{J}_{prec}`, which cannot be separated from the output of
   :doc:`compute heat/flux <compute_heat_flux>` when other pair styles
   are used with :doc:`pair_style hybrid/overlay <pair_hybrid>`.

If a group other than *all* is used, only the pairs with atom :math:`i`
in the group are included (atom :math:`j` can be in any group).  This is
in general not the heat current of the atoms in the group, and LAMMPS
prints a warning.

The values are computed from the precession vectors of the current time
step, which are available after the forces have been computed.  Use the
compute in output commands that are invoked at the end of a time step
(:doc:`thermo <thermo_style>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`fix ave/correlate <fix_ave_correlate>`, :doc:`fix print
<fix_print>`, :doc:`variables <variable>` used by these).

----------

Output info
"""""""""""

This compute calculates a global vector.  Its first 3 values are the
:math:`x`, :math:`y`, :math:`z` components of :math:`\vec{J}_{prec}`.
With *mech* = *yes*, the next 3 values are the components of
:math:`\vec{J}_{mech}`.  They are followed by the 3 components of each
channel of :math:`\vec{J}_{prec}`, in the order of the channels given
above.  The length of the vector is :math:`3 + 3 N(N+1)/2` (*mech* =
*no*) or :math:`6 + 3 N(N+1)/2` (*mech* = *yes*), where :math:`N` is
the number of atom types.  These values can be accessed by any command
that uses global vector values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The vector values are "extensive".  They are in units of energy times
velocity, i.e. eV Angstrom/ps for :doc:`metal units <units>`.

Restrictions
""""""""""""

This compute is part of the SPIN package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This compute requires atom style *spin* and pair style *spin/exchange*
or *spin/exchange/spline* (also as sub-styles of *hybrid* or
*hybrid/overlay*).  Pair style *hybrid/scaled* is not supported.

Related commands
""""""""""""""""

:doc:`compute spin/heat/flux <compute_spin_heat_flux>`,
:doc:`compute heat/flux <compute_heat_flux>`,
:doc:`compute spin <compute_spin>`,
:doc:`pair_style spin/exchange <pair_spin_exchange>`,
:doc:`pair_style spin/exchange/spline <pair_spin_exchange_spline>`,
:doc:`fix ave/correlate <fix_ave_correlate>`

Default
"""""""

The keyword default is mech = no.

----------

.. _Hardy1963:

**(Hardy)** Hardy, Phys Rev, 132, 168 (1963).

.. _Zhou2020b:

**(Zhou)** Zhou, Tranchida, Ge, Murthy, Fisher, Phys Rev B, 101, 224303 (2020).
