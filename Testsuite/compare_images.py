#!/usr/bin/env python
from PIL import Image
import numpy
import numpy.linalg
import os.path
import argparse 
import sys
# print sys.argv

parser = argparse.ArgumentParser(description = "Compares matviz generated images by the averaged Frobenuis norms for all channes normalized by size.")
parser.add_argument("reference", help="the reference image")
parser.add_argument("test", help="the image to compare against the reference")
parser.add_argument("--eps", help="normalized weighted Frobenius norm tolerance", type=float, default=1e-6)
args = parser.parse_args()


print("compare_images: ref='" + args.reference + "' test='" + args.test + "'")

if not os.path.exists(args.reference):
  print("error: reference image '" + args.reference + "' does not exist")
  os.sys.exit(1)
if not os.path.exists(args.test):
  print("error: test image '" + args.test + "' does not exist")
  os.sys.exit(1)

ref = Image.open(args.reference)
tst = Image.open(args.test)

if ref.size != tst.size: # if image sizes mismatch, resize to largest dimensions
   ref_width,ref_height = ref.size
   tst_width,tst_height = tst.size
   
   width = max(ref_width,tst_width)
   height = max(ref_height,tst_height)
   
   # used to be Image.ANTIALIAS but from PIL.Image 10.0 it needs to be LANCZOS
   ref = ref.resize((width,height),Image.LANCZOS)  
   tst = tst.resize((width,height),Image.LANCZOS)
   
   print("warning: image size mismatch! resized ref and test to width=",width, " and height=",height)
   
#    os.sys.exit(1)
  
ref.load()
tst.load()

a = numpy.asarray(ref).astype(float) * 1./255 
b = numpy.asarray(tst).astype(float) * 1./255

diff = 0.0
colors = 1 if len(a.shape) == 2 else min(3, a.shape[2]) # ignore potential alpha 
# rgb images
if colors > 1: # (1000, 1500, 3)
  for d in range(colors):
     diff += 1./colors * numpy.linalg.norm((a[:,:,d]-b[:,:,d]), ord='fro') 
     #print 'frobenius norm for d=' + str(d) + ' -> ' + str(diff)
else:
  diff = numpy.linalg.norm(a-b, ord='fro')        
  #print 'frobenius norm -> ' + str(diff)
diff /= float(a.shape[0]) * float(a.shape[1])

if diff <= args.eps:
  print('++ difference ' + str(diff) + ' is below tol ' + str(args.eps) + ' ++')        
  os.sys.exit(0)
else:
  print('** error: difference ' + str(diff) + ' above tol ' + str(args.eps) + ' **')        
  os.sys.exit(1) 
    