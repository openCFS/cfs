# Testcase for rayleigh damping in waterWave PDE

The test checks if the rayleigh damping works as expected. Currently only direct input of alpha and beta is supported. 
measuredFreq has to be defined in the mat file due to xml scheme, but is ignored for this PDE.

## Files:
- RectangularContainer.jou: Mesh with a rectangular water region enclosed in a container (mechanic region).
- Container2DrayleighDamping.xml applies an acceleration on the container in horizontal direction, with defined damping.
- Container2DrayleighDamping_undamped.xml applies the same acceleration with no damping defined 
- postprocess.ipynb: Pressure is measured on top right corner of the container and compared. Also the system matrices were exported to 
  verify the damping matrix.
