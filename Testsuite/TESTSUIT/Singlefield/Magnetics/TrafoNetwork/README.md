TrafoNetwork
============

* model with 4 elements
* the outer ones are *Core*, i.e. the magnetic circuit
* the inner ones are *Coil* and *Resistor*
* all cross sections and lengths are 1
* conductivity = 0.1 to get a frequency spectrum from 1-10 Hz

Linear harmonic
---------------
* as comparison with analytic network model

Multiharmonic + Transient
-------------------------

* 3 frequency points (low=1, medium=2.5, high=10)
* `run.sh` creates all runs, choose `MAT` and `NHARMONICS`
* excitation amplitudes are adapted for material `analytic-poly-2` to give maximum amplitudes of B=1
* for convenience in postprocessing we do not change them for different materials

Tested for the following cases:
* `analytic-poly-2` and 5 harmonics - perfect correspondence with transient
* `analytic-bilinear` and 3,5,7 harmonics -> nice convergence
