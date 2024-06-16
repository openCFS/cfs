FEAST Eigensolver
=================

Intel oneAPI Math Kernel Library (oneMKL)
implements [FEAST](http://www.feast-solver.org/) in version 2.
We use the community version of FEAST with considerably more features.

There is a method to estimate the search interval for FEAST:
https://arxiv.org/abs/1308.4275

ToDo
==================

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