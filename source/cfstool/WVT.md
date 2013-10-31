% Weight Vector Theory (WVT) Implementation in cfstool
% Simon Triebenbacher
% 29-10-2013

## Description

- Daten DVD mit
    + CFS++ binaries für Win/Lin
	+ WVT Beispieldaten
	+ sims Verezeichnis

- WVT arguments file.

- Surface WVT according to Hemp2 implemented but deactivated.

- Analytic $\mathbf{u}_1$ and $\mathbf{u}_2$

## Assumptions

- All input files containing meshes or data must be in HDF5 format.

- There must be at  least two input files with the same  grid and the harmonic
  fluid results  for the primary or  drive mode and the  secondary or Coriolis
  mode.

- There  may  be  an  optional  third input  file  containing  the  mean  flow
  data. Though  not strictly required  by the theory,  it is assumed  that the
  input files for  the primary and secondary modes and  the optional mean flow
  file, contain the same grid.

- There  may be  an  optional  fourth output  file.  The intermediate  element
  results  (e.g. weight  vector, force  densities, etc.)  get written  to this
  output file for postprocessing.

- The  actual WVT  results,  first  and foremost  the  phase difference,  gets
  written to a CSV file.

- If the fluid results (prim. &  sec.) have been calculated through a directly
  coupled FSI  simulation the harmonic  mechanic displacement for  the primary
  mode at  the sensor node  is read from the  input file of  the primary/drive
  mode (cf. (16) and (19) in Hemp1  $(u_y)_{P} = (u_y)_{P'}$. A named node for
  the sensor position  $P$ must be provided.  The  actual complex displacement
  value is either read from the history  or the mesh results branch inside the
  HDF5 file.

- If the fluid results (prim. & sec.) have been calculated through a decoupled
  FSI simulation, the  harmonic mech displacement for the primary  mode at the
  sensor node has to be provided through the WVT arguments file.

- Mean flow data can be provided as  results in respect to the grid inside the
  mean flow input  file, as analytic expression in respect  to the grid inside
  the mean flow input  file or as data with respect to  a point cloud provided
  as CSV file.

## Examples

- Straight Pipe with Stokes mean flow from HDF5

- Straight Pipe with Gersten-Herwig mean flow from profile.

~~~~ {.json}
{"wvt":
 {
  "integOrder":"3",
  "u_p":
  {
    "real":"0.00209",
    "imag":"0.0"
  },
  "V":
  {
    "meanVel":"1.0000000000e+01",
    "profile":"sample1D('profile.txt',sqrt(x^2+y^2),1)"
//    "profile":"1-(sqrt(x^2+y^2)/0.01315)^2"
//    "profile":"1"
  }
 }
}
~~~~

- Promass-F DN250 with mean flow from CSV.
