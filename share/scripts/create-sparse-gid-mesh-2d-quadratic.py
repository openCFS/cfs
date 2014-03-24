#!/usr/bin/python 

import Image, sys, os, bz2

# mimic conditional operator
def cond(test, trueval, falseval):
  if test:
    return trueval
  else:
    return falseval

if len(sys.argv) < 7 or len(sys.argv) > 9:
  print 'usage: ' + os.path.basename(sys.argv[0]) + ' <periodic|p|non-periodic|np> <mode=normal|n|quadratic|q> <infile.png> <threshold>  <out.mesh> <out.density.xml.bz2> [<scale>] [<x_pml>]'
  print 'with threshold from 255^-1..1'
  print 'with optional scale with 1.0 corresponding to a width of 1 meter'
  print 'with optional x_pml the size in fraction of an additional dense region added to the right'
  sys.exit(0)

periodic = False
if sys.argv[1] == "periodic" or sys.argv[1] == "p":
  periodic = True
else:
  if sys.argv[1] != "non-periodic" and sys.argv[1] != "np":     
    print "Invalid mode given, please choose either periodic|p|non-periodic|np"
    sys.exit(0)

regionstring = " 6 4 mech"
quad = False
if sys.argv[2] == "quadratic" or sys.argv[2] == "q":
  regionstring = " 7 8 mech"
  quad = True
else:
  if sys.argv[2] != "normal" and sys.argv[2] != "n":
    print "Invalid mode given, please choose either normal|n|quadratic|q!"
    sys.exit(0)

infile = sys.argv[3]
threshold = float(sys.argv[4])
gridfile = sys.argv[5]
densout  = sys.argv[6]

scale = 1.0
if len(sys.argv) >= 8:
  scale = float(sys.argv[7])

pml_fract = 0.0
if len(sys.argv) == 9:
  pml_fract = float(sys.argv[8])
  
# read the png into a list
im2 = Image.open(infile)
print "image mode: " + im2.mode
if im2.mode == 'I':
  print "Warning: mode is stupid, may give unusable results!"
im = im2.convert("L").transpose(Image.FLIP_TOP_BOTTOM)

nx, ny = im.size
imlist = list(im.getdata())

# in the PML case add negative values for pml elements which are not to be written to the density file and get another region
pml_x = int(nx * pml_fract) # we need if for the inner nodes
if pml_fract > 0:
  img_idx = 0
  data = list()
  for y in range(ny):
    for x in range(nx + pml_x):
      if x >= nx:
        data.append(512) # becomes negative when scaled later
      else:
        val = imlist[img_idx]
        data.append(val)
        img_idx += 1
  nx += pml_x
  imlist = data              
  

nnx = nx+1
nny = ny+1
offset = 0
if quad:
  offset = nnx*nny

numofallnodes = nnx*nny + 2*offset + 1

allelems = set()
pmlelems = set()
allnodes = set()


dout = bz2.BZ2File(densout, 'wb')
dout.write('<?xml version="1.0"?>\n <cfsErsatzMaterial infoWriteCounter="0" \
infoRejectCounter="0">\n <header>\n <mesh x="'+str(nx)+'" y="'+str(ny)+'" \
z="1"/>\n <design adapt_lower="true" \
constant="false" fixed="false" initial="0.5" lower="0" name="density" \
region="mech" scale="false" upper="1"/>\n <transferFunction \
application="mech" design="density" param="6" type="simp"/>\n \
</header>\n<set id="frompng">\n')
c = 0
el = 0
for y in range(ny):
  for x in range(nx):
# check for the density, if larger zero, we need the element and the
# corresponding nodes!
# CAUTION: we also have to take care of the periodic neighbors!
    val = 1.0 - float(imlist[el])/255.0
    el += 1
    needed = False
    # when non periodic needed stays False
    if val <= threshold and periodic: 
      # bottom elements
      if y == 0:
        elneighbor = x + (ny-1)*nx
        valneighbor = 1.0 - float(imlist[elneighbor])/255.0
        if valneighbor > threshold:
          needed = True
      # top elements
      if y == ny-1:
        elneighbor = x
        valneighbor = 1.0 - float(imlist[elneighbor])/255.0
        if valneighbor > threshold:
          needed = True
      # left elements
      if x == 0:
        elneighbor = y*nx+nx-1
        valneighbor = 1.0 - float(imlist[elneighbor])/255.0
        if valneighbor > threshold:
          needed = True
      # right elements
      if x == nx-1:
        elneighbor = y*nx
        valneighbor = 1.0 - float(imlist[elneighbor])/255.0
        if valneighbor > threshold:
          needed = True

    #print "el " + str(el).rjust(2) + ": val = " + str(val)
    if val < 0.0 or val > threshold or needed:
      c += 1
      if needed:
        val = 1e-9

      # we don't write pml elements to the density file
      if val >= 0.0: 
        dout.write('  <element nr="' + str(c) + '" type="density" design="' + str(val) + '"/>\n')
      else:
        pmlelems.add(x + nx*y + 1) # but mark them

      allelems.add(x + nx*y + 1)

      n1 = y*nnx + x + 1
      n2 = n1 + 1
      n3 = n1 + nnx + 1
      n4 = n1 + nnx
      allnodes.add(n1)
      allnodes.add(n2)
      allnodes.add(n3)
      allnodes.add(n4)
      if quad:
        allnodes.add(n1 + offset)       # n5
        allnodes.add(n1 + 2*offset + 1) # n6
        allnodes.add(n1 + nnx + offset) # n7
        allnodes.add(n1 + 2*offset)     # n8
dout.write('</set>\n</cfsErsatzMaterial>\n')
dout.close()


print "found " + str(len(allelems)) + " elements (incl. " + str(len(pmlelems)) + " pml) and " + str(len(allnodes)) + " nodes"
print "now writing mesh file " + gridfile + "..."


newnums = [0] * numofallnodes

c = 1
i = 1
for x in range(numofallnodes):
  if c in allnodes:
    newnums[c] = i
    i += 1
  c += 1


numnodebc = 0

# bottom (y=0)
for x in range(nnx):
  num = x+1
  if num in allnodes:
    numnodebc += 1
  if quad and (num + offset in allnodes):
    numnodebc += 1

# top (y=1)
for x in range(nnx):
  num = ny*nnx+x+1
  if num in allnodes:
    numnodebc += 1
  if quad and (num + offset in allnodes):
    numnodebc += 1

# left (x=0)
for y in range(nny):
  num = y*nnx+1
  if num in allnodes:
    numnodebc += 1
  if quad and (num + 2*offset in allnodes):
    numnodebc += 1

# right (x=1)
for y in range(nny):
  num = (y+1)*nnx
  if num in allnodes:
    numnodebc += 1
  if quad and (num + 2*offset in allnodes):
    numnodebc += 1

# inner nodes left from pml
if pml_x > 0:
  for y in range(nny):
    num = y * nnx + (nnx-pml_x)
    if num in allnodes:
      numnodebc += 1
    if quad and (num + 2*offset in allnodes):
      numnodebc += 1

gridf = open(gridfile, "w")
gridf.write("[Info]\n")
gridf.write("Version 1\n")
gridf.write("Dimension 2\n")
gridf.write("NumNodes " + str(len(allnodes)) + "\n")
gridf.write("Num3DElements 0\n")
gridf.write("Num2DElements " + str(len(allelems)) + "\n")
gridf.write("Num1DElements 0\n")
gridf.write("NumNodeBC " + str(numnodebc) + "\n")
gridf.write("NumSaveNodes 0\n")
gridf.write("NumSaveElements 0\n")
gridf.write("Num 2d-line      : 0\n")
gridf.write("Num 2d-line,quad : 0\n")
gridf.write("Num 3d-line      : 0\n")
gridf.write("Num 3d-line,quad : 0\n")
gridf.write("Num triangle     : 0\n")
gridf.write("Num triangle,quad: 0\n")
gridf.write("Num quadr        : " + str(len(allelems)) + "\n")
gridf.write("Num quadr,quad   : 0\n")
gridf.write("Num tetra        : 0\n")
gridf.write("Num tetra,quad   : 0\n")
gridf.write("Num brick        : 0\n")
gridf.write("Num brick,quad   : 0\n")
gridf.write("Num pyramid      : 0\n")
gridf.write("Num pyramid,quad : 0\n")
gridf.write("Num wedge        : 0\n")
gridf.write("Num wedge,quad   : 0\n\n")

gridf.write("[Nodes]\n")
gridf.write("#NodeNr x-coord y-coord z-coord\n")
c = 1
# we need to distort the node coordinates to handle not quadratic meshes
# we keep the maximal x-coordinate at 1 and only scale the y-coordinates
yfactor = float(ny)/float(nx-pml_x)
for y in range(nny):
  for x in range(nnx):
    if c in allnodes:
      #"%8d %20.13e %20.13e %20.13e\n", ++$c, 1.0*$x/$n, 1.0*$y/$n, 1.0*$z/$n
      gridf.write(str(newnums[c]).rjust(8) + " " + str(scale * float(x)/(nx-pml_x)) + " " + str(scale * yfactor*float(y)/ny) + " 0.000 \n")
    c += 1
if quad:
  # now comes the grid shifted to the right
  for y in range(nny):
    for x in range(nnx):
      if c in allnodes:
        #"%8d %20.13e %20.13e %20.13e\n", ++$c, 1.0*$x/$n, 1.0*$y/$n, 1.0*$z/$n
        gridf.write(str(newnums[c]).rjust(8) + " " + str(scale * float(x)/nx + 0.5/(nx-pml_x)) + " " + str(scale * yfactor*float(y)/ny) + " 0.000 \n")
      c += 1
  # now comes the grid shifted up
  for y in range(nny):
    for x in range(nnx):
      if c in allnodes:
        #"%8d %20.13e %20.13e %20.13e\n", ++$c, 1.0*$x/$n, 1.0*$y/$n, 1.0*$z/$n
        gridf.write(str(newnums[c]).rjust(8) + " " + str(scale * float(x)/(nx-pml_x)) + " " + str(scale * yfactor*(float(y)/ny + 0.5/ny)) + " 0.000 \n")
      c += 1

gridf.write("\n")
gridf.write("[1D Elements]\n")
gridf.write("#ElemNr  ElemType  NrOfNodes  Level\n")
gridf.write("#Node1 Node2 ... NodeNrOfNodes\n")
gridf.write("\n")
gridf.write("[2D Elements]\n")
gridf.write("#ElemNr  ElemType  NrOfNodes  Level\n")
gridf.write("#Node1 Node2 ... NodeNrOfNodes\n")
c = 0
for y in range(ny):
  for x in range(nx):
    elemnr   = x + nx*y + 1
    if elemnr in allelems:
      c += 1
      f = 4
      if quad:
        f = 8
      nodes    = [0.0]*f
      nodes[0] = newnums[y*nnx+x+1]
      nodes[1] = newnums[y*nnx+x+2]
      nodes[2] = newnums[(y+1)*nnx+x+2]
      nodes[3] = newnums[(y+1)*nnx+x+1]
      if quad:
        nodes[4] = newnums[y*nnx+x+1 + offset]
        nodes[5] = newnums[y*nnx+x+2 + 2*offset]
        nodes[6] = newnums[(y+1)*nnx+x+1 + offset]
        nodes[7] = newnums[y*nnx+x+1 + 2*offset]

      #gridf.write(str(elemnr) + regionstring + "\n")
      if elemnr in pmlelems:
        gridf.write(str(c) + regionstring + "1\n") # region mech1
      else:
        gridf.write(str(c) + regionstring + "\n")
      gridf.write(" ".join(map(str, nodes)) + "\n")
gridf.write("\n")
gridf.write("[3D Elements]\n")
gridf.write("#ElemNr  ElemType  NrOfNodes  Level\n")
gridf.write("#Node1 Node2 ... NodeNrOfNodes\n")
gridf.write("\n")


gridf.write("[Node BC]\n")
gridf.write("#NodeNr Level\n")

countdown = 0
# bottom (y=0)
for x in range(nnx):
  num = x+1
  if num in allnodes:
    gridf.write(str(newnums[num]) + " bottom\n")
    countdown += 1
  num += offset
  if quad and num in allnodes:
    gridf.write(str(newnums[num]) + " bottom\n")
    countdown += 1

counttop = 0
# top (y=1)
for x in range(nnx):
  num = ny*nnx+x+1
  if num in allnodes:
    gridf.write(str(newnums[num]) + " top\n")
    counttop += 1
  num += offset
  if quad and num in allnodes:
    gridf.write(str(newnums[num]) + " top\n")
    counttop += 1

countleft = 0
# left (x=0)
for y in range(nny):
  num = y*nnx+1
  if num in allnodes:
    gridf.write(str(newnums[num]) + " left\n")
    countleft += 1
  num += 2*offset
  if quad and num in allnodes:
    gridf.write(str(newnums[num]) + " left\n")
    countleft += 1

countright = 0
# right (x=1)
for y in range(nny):
  num = (y+1)*nnx
  if num in allnodes:
    gridf.write(str(newnums[num]) + " right\n")
    countright += 1
  num += 2*offset
  if quad and num in allnodes:
    gridf.write(str(newnums[num]) + " right\n")
    countright += 1

# inner bound between structure and pml - is 'right' w/o pml
countinner = 0
if pml_x > 0:
  for y in range(nny):
    num = y * nnx + (nnx-pml_x)
    if num in allnodes:
      gridf.write(str(newnums[num]) + " inner\n")
      countinner += 1
    num += 2*offset
    if quad and num in allnodes:
      gridf.write(str(newnums[num]) + " inner\n")
      countinner += 1


gridf.write("\n")


gridf.write("[Save Nodes]\n")
gridf.write("#NodeNr Level\n")
gridf.write("\n")

gridf.write("[Save Elements]\n")
gridf.write("#NodeNr Level\n")
gridf.write("\n")
gridf.close()

print 'found boundary nodes:'
print ' top:    ' + str(counttop)
print ' bottom: ' + str(countdown)
print ' left:   ' + str(countleft)
print ' right:  ' + str(countright)
print ' inner:  ' + str(countinner)
