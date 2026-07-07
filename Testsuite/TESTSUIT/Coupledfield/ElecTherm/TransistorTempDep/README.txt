This testcase features a temperature dependant transistor in a short-circuit application in order to test the tripole model from Sebastian Eiser. This application represents the simplified model of a vertical DMOS transistor, which is relevant for the industry (e.g. KAI).
To convert specific resistance data provided by KAI in a .table format to the sliced data of the electric conductivity used by openCFS, use the provided file convert_table.m.


Comment dmayrhof 20240412: 

I disabled some results since they lead to quite large numerical problems.
The elecPowerDensity gets computed correctly for all regions which is enough for the test itself.
The reason why the power fails to produce reliable/testable results for vWire2 is the relatively high voltage of around 14V.
When computing the postProcResult via the ResultFunctor, vWire1 does not suffer a problem since the voltage is basically 0V - hence the numerical differences during the computation of the power do not play a significant role.
For vWire2 the computation looks different since the voltage level used for the computation is high while the entries for elemMatR used in the ResultFunctor are also relatively large (roughly 1e5 - see sim_log.txt).
Hence, getting accurate results which should lie somewhere way below 1e-10 is not possible due to numerical inaccuracies, leading to fluctuating results in the range of 1e-9 to 1e-8.
For vWire1 we get down to below 1e-20 due to the small voltage levels.
Hence, testing the result for vWire in a relative sense can be quite annyoing/non-stable.
The elecPowerDensity also suffers this problem but a bit differently since the values in general are way higher since the kernel used for the calculation in the CoefFunctionFormBased gives values in the range of 1e15.
One way to combat this problem would be to change the geometry.

Therefore, the following tests have been disabled:
-elecCurrentDensity on vGate 
-elecCurrentDensity on vWire2
-elecPowerDensity on vWire2
-elecPower on vWire2
Switching to absL2diff unfortunately does not help, hence, they are completely disabled.
