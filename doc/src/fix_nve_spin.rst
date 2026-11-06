.. index:: fix nve/spin

fix nve/spin command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/spin keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/spin = style name of this fix command
* zero or more keyword/value pairs may be appended
* keyword = *lattice* or *integrator* or *iterations* or *tol*

  .. parsed-literal::

       *lattice* value = *moving* or *frozen*
         moving = integrate both spin and atomic degrees of freedom
         frozen = integrate spins on a fixed lattice
       *integrator* value = *sequential* or *midpoint*
         sequential = update the spins one after the other
         midpoint = update all spins at once with the implicit midpoint rule
       *iterations* value = N
         N = number of field evaluations per spin update (*midpoint* only)
       *tol* value = tolerance
         tolerance = stop the iterations when no spin component changes by more
                     than this value (*midpoint* only; 0.0 = always do N iterations)

Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all nve/spin lattice moving
   fix 1 all nve/spin lattice frozen
   fix 1 all nve/spin lattice moving integrator midpoint
   fix 1 all nve/spin integrator midpoint iterations 20 tol 1.0e-12

Description
"""""""""""

Perform a symplectic integration for the spin or spin-lattice system.

The *lattice* keyword defines if the spins are integrated on a lattice
of fixed atoms (lattice = frozen), or if atoms are moving
(lattice = moving).
The first case corresponds to a spin dynamics calculation, and
the second to a spin-lattice calculation.
By default a spin-lattice integration is performed (lattice = moving).

The *nve/spin* fix applies a Suzuki-Trotter decomposition to
the equations of motion of the spin lattice system, following the scheme:

.. image:: JPG/fix_integration_spin_stdecomposition.jpg
   :align: center

according to the implementation reported in :ref:`(Omelyan) <Omelyan1>`.

A sectoring method enables this scheme for parallel calculations.
The implementation of this sectoring algorithm is reported
in :ref:`(Tranchida) <Tranchida1>`; it is based on the synchronous
sublattice algorithm of :ref:`(Shim) <Shim1>`.
Each processor sub-domain is cut in two halves (sectors) along every
direction in which the simulation box is divided between several
processors (see the :doc:`processors <processors>` command), giving 2,
4, or 8 sectors.  Spins in sectors with the same index on different
processors do not interact and are updated at the same time.  This
requires the sub-domain of each processor to be at least twice as wide
as the largest magnetic pair cutoff plus the :doc:`neighbor skin
distance <neighbor>` in each of these directions.  The skin is needed
because atoms can move by up to half the skin distance between two
neighbor list builds.  There is no size requirement in directions that
are not divided between processors, so that, e.g., thin films can be run
in parallel with a processor grid like "processors \* \* 1".  With
:doc:`comm_style tiled <comm_style>`, all directions are cut in two.

.. versionchanged:: TBD

Sectors are only created in directions divided between processors
(before, 2x2x2 sectors were always required), the size requirement
includes the neighbor skin, and it is checked again whenever the
sub-domains change.

The sequential update of the spins described above is the default
(*integrator sequential*).  Since the spins of one sector have to be
updated before those of the next one, the spin update requires about 30
exchanges of data between processors per timestep, which limits the
parallel efficiency when there are few atoms per processor.

.. versionadded:: TBD

   *integrator midpoint*, *iterations*, and *tol* keywords

With *integrator midpoint*, each of the two spin updates per timestep
(over half a timestep, :math:`h = \Delta t/2`) is done for all spins at
once with the implicit midpoint rule

.. math::

   \vec{s}_i(t+h) = \vec{s}_i(t) + h \, \vec{\omega}_i(\vec{m}) \times \vec{m}_i,
   \qquad \vec{m} = \frac{\vec{s}(t) + \vec{s}(t+h)}{2}

where :math:`\vec{\omega}_i(\vec{m})` is the magnetic field (precession
vector) acting on spin *i*, computed from the midpoint values of all
spins.  The new spin is obtained as a rotation of the old one, so that
the length of each spin is preserved exactly.  The equation is solved by
iteration: each iteration computes the fields of all spins from the
current estimate of the midpoints and updates all spins.  The *iterations*
keyword sets the number of iterations (field evaluations), and the *tol*
keyword optionally stops the iterations earlier once no spin component
changes by more than *tol* (this requires an extra global sum per
iteration).  If the iterations do not converge to *tol*, a warning is
printed.  For timesteps of about 0.1 fs in bcc iron, each iteration
changes the spins about two orders of magnitude less than the previous
one, so that the default of 4 iterations solves the equation to nearly
round-off.  The change of the spins in the last iteration can be
monitored (see below) to check this for other systems and timesteps.
The remaining error of the iterations accumulates over long runs: with
the default settings and a timestep of 0.4 fs, the energy of 1024 spins
of bcc iron drifted by about 7e-4 eV over 100 ps (1.5e-5 eV at 0.2 fs).
For long runs at large timesteps, more iterations or a *tol* value
remove this drift.

When the equation is solved exactly, the implicit midpoint rule conserves
the magnetic energy for interactions that are quadratic in the spins
(e.g., exchange, DMI, anisotropy, Zeeman) and the total magnetization for
isotropic exchange :ref:`(Hairer) <Hairer1>`.  Like the sequential update,
it is second order in the timestep, and time-reversible when solved
exactly.  Unlike the spherical midpoint rule of :ref:`(McLachlan)
<McLachlan1>`, which evaluates the field at normalized midpoints, it is
not symplectic for spin systems in general, so that for interactions
that are not quadratic in the spins (e.g., cubic anisotropy,
biquadratic exchange, or *spin/neel*), a bounded energy error over long
runs is not guaranteed.  Its error at a given timestep is
larger than that of the sequential update: for bcc iron at 300 K with a
moving lattice, the spin trajectories deviated about 4 times more from a
converged reference, while energy conservation was the same.  It needs
more field evaluations per timestep than the sequential update (8
instead of 4 with the default settings), but only
one exchange of data between processors per iteration, and no sectors:
there is no restriction on the size of the processor sub-domains, any
processor grid or
:doc:`comm_style tiled <comm_style>` can be used, and the result does
not depend on the number of processors (except for round-off).  This
can make it faster than the sequential update when running with few atoms
per processor.

With :doc:`fix langevin/spin <fix_langevin_spin>` at a non-zero
temperature, the random field of each spin is drawn once per spin update
(of length :math:`h`) and kept fixed during the iterations, with an
amplitude for a time interval :math:`h`.  This corresponds to the
Stratonovich interpretation of the stochastic equation of motion
:ref:`(Garcia-Palacios) <GarciaPalacios1>`, as in the midpoint schemes of
:ref:`(Mentink) <Mentink1>`.  The random numbers then differ from those
of the sequential update, and with several processors also between
different numbers of processors.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.

With *integrator midpoint*, this fix computes a global vector of length
2, which can be accessed by various :doc:`output commands
<Howto_output>`.  The first element is the number of field evaluations
per spin during the last timestep (for both spin updates), the second
element is the largest change of a spin component in the last iteration.
These values are "intensive".  No global or per-atom quantities are
stored by this fix with *integrator sequential*.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

----------

Restrictions
""""""""""""

This fix style can only be used if LAMMPS was built with the SPIN
package.  See the :doc:`Build package <Build_package>` page for more
info.

To use the spin algorithm, it is necessary to define a map with
the atom_modify command. Typically, by adding the command:

.. code-block:: LAMMPS

   atom_modify map array

before you create the simulation box. Note that the keyword "hash"
instead of "array" is also valid.

When running on multiple processors with a moving lattice, the sectoring
algorithm relies on atoms moving less than half the neighbor skin
distance between two neighbor list builds.  This is only guaranteed
with the settings

.. code-block:: LAMMPS

   neigh_modify every 1 delay 0 check yes once no

With other settings (e.g. *every* > 1, *delay* > 0, or *once* yes), a
warning is printed, and another one if atoms did move too far during the
run (i.e. if there are "dangerous builds"), since then spins of atoms on
different processors may have been updated in the wrong order.  This
does not apply to *integrator midpoint*.

Related commands
""""""""""""""""

:doc:`atom_style spin <atom_style>`, :doc:`fix nve <fix_nve>`

Default
"""""""

The option defaults are lattice = moving, integrator = sequential,
iterations = 4, and tol = 0.0.

----------

.. _Omelyan1:

**(Omelyan)** Omelyan, Mryglod, and Folk, Phys. Rev. Lett.
86, 898 (2001).

.. _Tranchida1:

**(Tranchida)** Tranchida, Plimpton, Thibaudeau and Thompson,
Journal of Computational Physics, 372, 406-425, (2018).

.. _Hairer1:

**(Hairer)** Hairer, Lubich, and Wanner, Geometric Numerical
Integration, 2nd edition, Springer (2006), Chapter IV.

.. _Mentink1:

**(Mentink)** Mentink, Tretyakov, Fasolino, Katsnelson, and Rasing,
J Phys Condens Matter, 22, 176001 (2010).

.. _Shim1:

**(Shim)** Shim and Amar, Phys. Rev. B, 71, 125432 (2005).

.. _McLachlan1:

**(McLachlan)** McLachlan, Modin, and Verdier, Phys. Rev. E, 89,
061301 (2014).

.. _GarciaPalacios1:

**(Garcia-Palacios)** Garcia-Palacios and Lazaro, Phys. Rev. B, 58,
14937 (1998).
