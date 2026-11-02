# compute spin/heat/flux example

`in.spin.heat_flux`: bcc iron (EAM + spin/exchange), spin-lattice dynamics.
The two sites of the bcc cell have different atom types (1 = Fe1 at the cube
corners, 2 = Fe2 at the body centers) with identical interactions, so the heat
flow can be split into the exchange channels Fe1-Fe1, Fe1-Fe2 and Fe2-Fe2.

1. The spins are heated to 300 K with fix langevin/spin while the lattice
   starts at 600 K and equilibrates to about 300 K.
2. The thermostat is removed and an NVE run computes the heat flow carried by
   the spins from the right half of the box to the left half
   (`compute qs all spin/heat/flux right left`), as in
   Zhou et al., Phys. Rev. B 101, 224303 (2020), eq. (13).
   `c_qs` is the total heat flow, `c_qs[1]`, `c_qs[2]`, `c_qs[3]` are the
   channels 1-1, 1-2 and 2-2; they add up to `c_qs`.
3. The total flux q = Q/(ly*lz) is correlated with fix ave/correlate (file
   `Qspin.acf`) and integrated to a Green-Kubo estimate of the magnon thermal
   conductivity along x, kappa = V/(kB T^2) * integral <q(t) q(0)> dt (W/mK).
4. All auto- and cross-correlations of the channels are written to
   `Qspin_channels.acf` (fix ave/correlate, type full).  The last line of the
   output prints the command for `gk_spin_channels.py`, which splits kappa into
   channel contributions:

```
python3 gk_spin_channels.py Qspin_channels.acf --timestep 0.0001 \
    --volume 12059.40274 --area 525.876624 --temperature 300 --names Fe1,Fe2
```

The script prints kappa (equal to the value from the total flux, to the 6
digits that fix ave/correlate writes), the share of each channel (including its
cross terms with the other channels; the shares add up to kappa and can be
negative), the channel on its own, and the full matrix kappa_kl.

Why the cross terms matter: the Green-Kubo integral needs the correlation of
the TOTAL flux.  Adding the autocorrelations of the channels leaves out the
cross-correlations <q_k(t) q_l(0)>, which are large because the channels share
spins (the motion of one spin enters all channels it takes part in).

The runs are far too short for converged values: production calculations need
correlation times of several ps and many independent runs, and the channel
shares need more statistics than the total kappa.  In this example the shares
change completely between the 1 and 4 processor runs (Fe1-Fe1 and Fe2-Fe2 are
equivalent by symmetry here, which is a useful check of convergence).
With spin/exchange/spline instead of spin/exchange, the compute takes J(r)
from the table; nothing else changes in the input.

Log files for 1 and 4 processors differ because fix langevin/spin uses
different random numbers on each processor.
