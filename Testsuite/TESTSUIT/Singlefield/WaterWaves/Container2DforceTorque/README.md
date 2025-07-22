# Testcase for Force/Torque (and respective densities) in waterWave PDE

## Files:
- RectangularContainer.jou: Mesh with a rectangular water region enclosed in a container (mechanic region).
- Container2DforceTorque.xml: The test checks if the postprocessing results work as expected, by applying and acceleration boundary condition on the water-container boundary in the waterWavePDE.
- Container2DforceTorque_waterWaveMech.xml: The results are verified via WaterWave-mechanic coupling in the file Container2Dacceleration_waterWaveMech.xml, by comparing the applied Force/Torque on a container filled with water to the postprocessing results in Postprocessing.ipynb
