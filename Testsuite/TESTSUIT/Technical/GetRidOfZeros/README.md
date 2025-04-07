# Testing the `GetRidOfZeros` function

## Introduction

The function is implemented to get rid of zeros in non-zero entities of a system matrix.
The known reason for zeros ending up in a system matrix is a non-conforming interface (NCI).

The zeros may significantly affect the sparse solvers' performance if an NCI is large enough in comparison to the whole model. Both RAM and CPU time can increase.

## Behaviour of NCI issue

* For moving interfaces, CFS defines the sparsity pattern of the system matrix such that all DoFs of a non-conforming interface side can potentially be coupled to all DoFs of the other side
* For non-moving interfaces, there is still a large number of unnecessary zeros introduced in the system matrix
* Only in the assembly, the intersection grid is accounted for by creating non-zero entries for intersecting DoFs (i.e. for DoFs of elements that have an overlap)

## Implemented solution

We have implemented `GetRidOfZeros` function, which loops over all the elements of the system matrix and creates a copy containing only elements larger than the tolerance ($1e-20$ by default). Then the original system matrix is overwritten with the copy.

The implementation is tested for a limited amount of models. For example, the current implementation (24.07.2024) works only for Pardiso solver.

The function is applied automatically to all models with NCIs (if it's implemented for a specific case).
In the XML scheme, one can choose one of three options `auto` (default), `yes` or `no`. Additionally one can change the tolerance of zero removing procedure.

```
    <linearSystems getRidOfZeros="yes" getRidOfZerosTolerance="1e-14">
        <system>
            <solutionStrategy>
                <standard/>
            </solutionStrategy>
        </system>
    </linearSystems>
```

## Test

We test all the three options `auto` (default), `yes` or `no`. For `yes`, we select a tolerance.
The test checks only the functionality of `GetRidOfZeros`. **Namely, we compare the system matrices modyfied by `GetRidOfZeros` with the reference solutions validated in Python.** The correctness of the results is NOT validated. Additionally, the test checks neither the correctness of the implementation for all the PDEs, coupled systems or solvers, nor the performance increase (time decrease) for the solver and the whole computation!

To generate new reference matrices, run the simulation XML file and rename the exported files.
**After the generation of new reference matrices, verify them with `CheckRefMat.ipynb`.**
