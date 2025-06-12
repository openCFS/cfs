CurrentSheet3D multiharmonic
============================

* The H-distribution is linearly increasing with z from 0 to J
* conductivity = 0 gives a quasi-static solution
* several analytical material descriptions are provided in the material file
* the `run.sh` runs a transient simulation and some multiharmonic ones

Tested things
-------------

* `analytic-linear` works perfeclty
* `analytic-bilinear` fits well for the linear part, but the higher harmonics are much too small
* `analyic-poly-3` converging to something too high
