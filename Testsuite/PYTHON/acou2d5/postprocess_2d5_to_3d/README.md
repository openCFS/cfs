A Wavenumber spectrum for 2.5D analysis for 100 Hz excitation frequency is stored in three .cfs files, and we postprocess the wavenumber results to spatial sound pressure for z = "0.1,0.2,0.3,0.4,0.5"

Test with h5diff

h5diff usage:
### Use h5diff
-DTEST_H5DIFF:BOOL=ON
### Wet set the limit of relative difference to be 1e-4
-DEPSILON=1e-4
