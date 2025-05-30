# Testcase for Acceleration BC in waterWave PDE

The test checks if the acceleration boundary condition in the waterWave PDE works as expected, by comparing to the acceleration BC applied in the mechanic region in waterWaveMech.

## Files:
- RectangularContainer.jou: Mesh with a rectangular water region enclosed in a container (mechanic region).
- Container2Dacceleration.xml applies the acceleration in the waterWave PDE
- Container2Dacceleration_waterWaveMech.xml applies the acceleration on the container, coupled to the water region (WaterWaveMech).
- postprocess.ipynb: Pressure is measured on top right corner of container and compared
