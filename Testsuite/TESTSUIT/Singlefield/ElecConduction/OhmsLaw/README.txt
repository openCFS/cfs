To validate the results computed with cfs the power is computed 
analytically, based on the test OhmsLaw (Devel/CFS_SCR/cfs/Testsuite/TESTSUIT/Singlefield/ElecConduction/OhmsLaw).

The electric power is calculated with the formula 

	\int_Ω γ E \cdot E dΩ

with dΩ = dxdydz.

In the test OhmsLaw a 3d model of a cuboid wire is given. The wire consists of two parts with the following dimensions:

wire 1: x from [0,1], y from [0,1] and z from [0, 20]
wire 2: x from [0,1], y from [0,1] and z from [20, 30]

The value for the magnitude of the electric field intensity E in z-direction is obtained from the computation with openCFS. The value of γ
is defined in the mat.xml file.

When solving the integrals for both wires with the parameters
Ez = -0.0333333333333
γ = 30

the following results for the electric power are obtained:

elecPower on wire 1 = 0.6666
elecPower on wire 2 = 0.3333

They match the values obtained via the computation with cfs: 

elecPower on wire 1 = 0.6666
elecPower on wire 2 = 0.3333

For the corrected factor for the elecPowerDensity the formula is:

	e = γ E \cdot E

Calculating this with the upper values for γ and Ez the following electric power density e is calculated for all elements:

elecPowerDensity = 0.03333

which is also calculated by openCFS, seen in the elemResult in the history directory of the testcase.
