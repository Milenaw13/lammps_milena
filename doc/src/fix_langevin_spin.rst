.. index:: fix langevin/spin

fix langevin/spin command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID langevin/spin T Tdamp seed keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin/spin = style name of this fix command
* T = desired temperature of the bath (temperature units, K in metal units)
* Tdamp = transverse magnetic damping parameter (adim)
* seed = random number seed to use for white noise (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *tally*

  .. parsed-literal::

       *tally* value = *no* or *yes*
         *no* = do not tally the energy exchanged with the bath
         *yes* = tally the energy exchanged with the bath

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all langevin/spin 300.0 0.01 21
   fix 2 hot langevin/spin 600.0 0.01 21 tally yes

Description
"""""""""""

Apply a Langevin thermostat as described in :ref:`(Mayergoyz) <Mayergoyz1>` to the
magnetic spins associated to the atoms.
Used with :doc:`fix nve/spin <fix_nve_spin>`, this command performs
Brownian dynamics (BD).
A random torque and a transverse dissipation are applied to each spin i according to
the following stochastic differential equation:

.. math::

   \frac{d \vec{s}_{i}}{dt} = \frac{1}{\left(1+\lambda^2 \right)} \left( \left(
   \vec{\omega}_{i} +\vec{\eta} \right) \times \vec{s}_{i} + \lambda\, \vec{s}_{i}
  \times\left( \vec{\omega}_{i} \times\vec{s}_{i} \right) \right)

with :math:`\lambda` the transverse damping, and :math:`\eta` a random vector.
This equation is referred to as the stochastic Landau-Lifshitz (sLL)
equation.

The components of :math:`\eta` are drawn from a Gaussian probability
law. Their amplitude is defined as a proportion of the temperature of
the external thermostat T (in K in metal units).

More details about this implementation are reported in :ref:`(Tranchida) <Tranchida2>`.

Note: due to the form of the sLL equation, this fix has to be defined just
before the nve/spin fix (and after all other magnetic fixes).
As an example:

.. code-block:: LAMMPS

   fix 1 all precession/spin zeeman 0.01 0.0 0.0 1.0
   fix 2 all langevin/spin 300.0 0.01 21
   fix 3 all nve/spin lattice moving

is correct, but defining a force/spin command after the langevin/spin command
would give an error message.

----------

.. versionadded:: TBD

   The *tally* keyword

The keyword *tally* enables the calculation of the cumulative energy
exchanged between the spins of the group and the bath.  This is useful,
for example, to measure the magnon heat flux in a non-equilibrium
simulation with a hot and a cold spin reservoir, in the same way as the
*tally* keyword of :doc:`fix langevin <fix_langevin>` is used for the
lattice.  The energy of spin :math:`i` is
:math:`-\hbar\, \vec{\omega}_{i} \cdot \vec{s}_{i}`, with
:math:`\vec{\omega}_{i}` its precession vector from all other
interactions (in rad/ps).  At each update of spin :math:`i` by
:doc:`fix nve/spin <fix_nve_spin>` (four per timestep), the spin is
rotated around the total precession vector, including damping and random
field.  Without damping and random field the rotation keeps
:math:`\vec{\omega}_{i} \cdot \vec{s}_{i}` constant, so the energy
exchanged with the bath during the update is

.. math::

   \Delta E_i = -\hbar\, \vec{\omega}_{i} \cdot \left( \vec{s}_{i}^{\,new}
   - \vec{s}_{i}^{\,old} \right)

This is exact for interactions that are linear in each spin (e.g. pair
styles *spin/exchange*, *spin/dmi* and the Zeeman term of
:doc:`fix precession/spin <fix_precession_spin>`), and holds to first
order in the timestep otherwise.  With the *midpoint* integrator of
:doc:`fix nve/spin <fix_nve_spin>`, all spins are updated at once and
:math:`\vec{\omega}_{i}` is evaluated at the midpoints of the spins
before and after the update; the tally is then exact for spin energies
that are quadratic in the spins, when the iterations are converged.  It is the spin equivalent of the
*tally* keyword of :doc:`fix langevin <fix_langevin>`: the sum of the
total energy and of the tallies of all thermostats (thermo keyword
*econserve*) is conserved up to the integration error.

In addition, the expected energy exchange is accumulated from the
average power of damping and random field, given the current spin
direction and precession vector:

.. math::

   P_i = \frac{\lambda}{1+\lambda^2} \left( 2 k_B T \, \vec{s}_{i} \cdot
   \vec{\omega}_{i} - \hbar \left| \vec{s}_{i} \times \vec{\omega}_{i}
   \right|^2 \right)

The first term is the energy given by the random field, the second the
energy taken by the damping.  :math:`P_i` summed over the group
vanishes when the spin temperature of the group, as defined by
:ref:`(Nurdin) <Nurdin1>` and used by :doc:`compute spin <compute_spin>`,
equals the bath temperature.  On average it agrees with the exact tally,
but it does not contain the work of the actual random field, so the
difference between the two grows like the square root of time.  For
heat fluxes, the exact tally (the scalar of this fix) should be used.

The groups of different *langevin/spin* fixes should not overlap when
*tally* is used.

Note: The random # *seed* must be a positive integer.  A Marsaglia random
number generator is used.  Each processor uses the input seed to
generate its own unique seed and its own stream of random numbers.
Thus the dynamics of the system will not be identical on two runs on
different numbers of processors.

With the *integrator midpoint* option of :doc:`fix nve/spin
<fix_nve_spin>`, the random vector of each spin is drawn once per spin
update and kept fixed during the iterations of the implicit midpoint
rule; see the :doc:`fix nve/spin <fix_nve_spin>` page.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  Because the state of the random number generator
is not saved in restart files, this means you cannot do "exact"
restarts with this fix, where the simulation continues on the same as
if no restart had taken place.  However, in a statistical sense, a
restarted simulation should produce the same behavior.

The :doc:`fix_modify <fix_modify>` *energy* option is not supported by
this fix, but this fix does add the cumulative energy removed from the
spins to the *ecouple* thermodynamic keyword (see the
:doc:`thermo_style <thermo_style>` command), so that *econserve* can be
used to monitor energy conservation, when the *tally* keyword is set to
*yes*.

This fix computes a global scalar and a global vector of length 2, which
can be accessed by various :doc:`output commands <Howto_output>`.  The
scalar is the cumulative energy removed from the spins of the group by
the bath (exact tally, see above), i.e. the negative of the energy added
to them.  The vector holds the cumulative energy removed from the spins
according to the average power, split into the part removed by the
damping (first element, positive) and the part removed by the random
field (second element, negative on average).  All are zero unless the
*tally* keyword is set to *yes*.  The scalar and the vector values are
"extensive" and in energy units.  They are not reset at the start of a
run and are not saved in restart files.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The *langevin/spin* fix is part of the SPIN package.  This style is
only enabled if LAMMPS was built with this package.  See the :doc:`Build package <Build_package>` page for more info.

The numerical integration has to be performed with *fix nve/spin*
when *fix langevin/spin* is enabled.

This fix has to be the last defined magnetic fix before the time
integration fix (e.g. *fix nve/spin*\ ).

Related commands
""""""""""""""""

:doc:`fix nve/spin <fix_nve_spin>`, :doc:`fix precession/spin <fix_precession_spin>`

Default
"""""""

The option default is tally = no.

----------

.. _Mayergoyz1:

**(Mayergoyz)** I.D. Mayergoyz, G. Bertotti, C. Serpico (2009). Elsevier (2009)

.. _Tranchida2:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).

.. _Nurdin1:

**(Nurdin)** Nurdin and Schotte, Phys. Rev. E, 61, 3579 (2000).
