Eigensolver
=================

[MKL](https://software.intel.com/en-us/mkl-developer-reference-c-extended-eigensolver-naming-conventions)
implements [FEAST](http://www.feast-solver.org/).

Implementation from fwein into `Fabian_Wein_TopOpt_0415`, added in r14332, removed in r14409

(before remove)[https://cfs.mdmt.tuwien.ac.at/trac/browser/CFS%2B%2B/branches/Fabian_Wein_TopOpt_0415?rev=14395]

there is a method to estimate the search-interval for feast:
https://arxiv.org/abs/1308.4275

[brief-feastdescription](https://software.intel.com/en-us/mkl-developer-reference-c-the-feast-algorithm)

[feast-hcsrev](https://software.intel.com/en-us/mkl-developer-reference-c-feast-scsrev/-feast-hcsrev)

the indexes of ia and ja must be correct (number of elements PLUS 1)

EigenFrequencyDriver does the reading in of parameters

ToDo
==================

Implement differnt matrix types
-------------------------------

* just use the differnt versions of the FEAST function
* distingish between complex and real valued

Input Structure for the EigenFrequencyDriver
--------------------------------------------

# shift & num
 * shift point
 * number of modes

# interval
 * min
 * max

Parameters that should follow from the problem definition
 * isBloch -> determined by the presence of wave vectors
 * isQuadratic -> determined by the presence of damping

Make FEAST work with num-freq
-----------------------------
One needs to [estimate the search interval](https://arxiv.org/abs/1308.4275).

Use the community based FEAST
-----------------------------
* it offers routines for non symmetric/hermitian problems.
