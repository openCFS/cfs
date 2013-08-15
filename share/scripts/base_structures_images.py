#!/usr/bin/env python

# this simple tool creates images as basis for meshes via create_mesh.py
# It would be easy to integrate mess_tool.py directly but it is nice to have the images anyway!

import argparse, Image, sys, os, copy

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="resolution of image", type=int, required = True )
parser.add_argument('--x_thickness', help="relative thickness for horizontal bar <= 1", type = float, required = True)
parser.add_argument('--y_thickness', help="relative thickness for vertical bar <= 1", type = float, required = True)
parser.add_argument('--file', help="optional give output file name")
parser.add_argument('--noshow', help="suppress showing the image", action='store_true')
arg = parser.parse_args()

img = Image.new("L", ((arg.res, arg.res)), "white")
pix = img.load() 

## draw cross
width = arg.res * arg.x_thickness
for x in range(int(0.5 * (arg.res - width)), int(0.5 * (arg.res + width))):
  for y in range(0, arg.res):
    pix[x,y] = 1

height = arg.res * arg.y_thickness
for y in range(int(0.5 * (arg.res - height)), int(0.5 * (arg.res + height))):
  for x in range(0, arg.res):
    pix[x,y] = 1
    
if not arg.noshow:    
  img.show()
  
if arg.file:
  img.save(arg.file)
  print "saved file '" + arg.file + "'" 