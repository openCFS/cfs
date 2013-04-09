#!/usr/bin/env python
import os.path
import glob

# search all binary gid results
files = glob.glob('*.post.bin')
  
# create a gid batch file to be opened in gid bv Utilities->Tools->Read Batch Window  
batch = open("gid_screenshots.bch", "w")  
  
for file in files:  
  problem = file[0:file.find(".post.bin")]
  batch.write("Mescape Files Read " + os.getcwd() + "/" + problem + ".post.bin\n")
  batch.write("Mescape zoom frame\n")
  # batch.write("Mescape select volumes V_membrane V_simple_support V_line V_clamped V_hot V_gnd V_side\n")
  # batch.write("Mescape Results geometry original geometry Deformation mechDisplacement_(m)\n")
  # batch.write("Mescape Results Deformation mechDisplacement_(m)\n")
  batch.write("Mescape Results ContourFill penalizedStress_default_(plain)\n")
  batch.write("Mescape View HardCopy png " + os.getcwd() + "/" + problem + "-pen_stress.png\n")

  batch.write("Mescape Results ContourFill physicalPseudoDensity\n")
  batch.write("Mescape View HardCopy png " + os.getcwd() + "/" + problem + "-physical.png\n")

  batch.write("Mescape Results ContourFill mechPseudoDensity\n")
  batch.write("Mescape View HardCopy png " + os.getcwd() + "/" + problem + "-density.png\n")
