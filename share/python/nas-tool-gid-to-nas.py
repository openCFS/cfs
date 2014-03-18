#!/usr/bin/python

import sys

# define a mapping from name in gid file to type of condition (boundary, force) in nastran file
# FORCE: needs an index and four double values: index is for multiload!
# SPC:   needs two values: a string defining type of bc (1: x, 2: y, 12: xy, 123: xyz, ...) and a value
table = {'load':  ['FORCE', 1, 1.0, 0.0, -1.0, 0.0],\
         'lfixed': ['SPC', '  12', 0.0],\
         'rfixed': ['SPC', '  2', 0.0]}










#####################################################
#####################################################
##          DO NOT EDIT BELOW THIS LINE            ##
#####################################################
#####################################################

# TODO:
# 1. check for nodes which have a bc and a force assigned at the same time 
#    and remove them in one of the settings




if len(sys.argv) < 3:
  print "usage: " + sys.argv[0] + " <grid-file> <outfile>"
  print "\033[93mInfo:\033[0m the beginning of the script has to be edited!"
  print "      see comment inside the script"
  sys.exit(1)

infile = sys.argv[1]
outfile = sys.argv[2]




#### helper function
def printWarning(message):
  print "\033[91m WARNING:\033[0m " +  message


print "\033[93mInfo:\033[0m"
print " found the following definitions:"
for i in table:
  print i.rjust(16) + ': ' + str(table[i])



# modify table so that it contains the right string
for i in table:
# must handle FORCE and SPC differently!
  out = ''
  if table[i][0] == 'FORCE':
    valuelist = table[i][2:]
    out = str(0).rjust(8)
    for l in valuelist:
      out += (str(l).ljust(8))
    table[i] = [table[i][0], table[i][1], out]

  if table[i][0] == 'SPC':
    valuelist = table[i][1:]
    for l in valuelist:
      out += (str(l).ljust(8))
    table[i] = [table[i][0], out]





# list containing the names of named entities in grid-file as indices, and the number as values
# this is just to check if all named entities in the mesh get a bc or if we missed something
names = {}


out = open(outfile, 'w')
sectinfo    = False
sectnodes   = False
sectelems   = False
sectnodesbc = False

griddimension = 0

# need for parsing the elements section consisting of two lines...
elemstmpline = ''

nodectr = 0
elemctr = 0

# write header into outfile
out.write('ENDCONTROL\n')

# create multiload header
multiloads = {}
for i in table:
  if table[i][0] == 'FORCE':
    k = str(table[i][1])
    if k in multiloads.keys():
      multiloads[k] += 1
    else:
      multiloads[k]  = 1

for i in multiloads.keys():
  out.write('SUBCASE       %s\n' % i)
  out.write('  LABEL= SUBCASE %s \n' % i)
  out.write('  LOAD =       %s\n' % i)



out.write('BEGIN BULK \n')


# we do everything in one go (i. e. parse the grid-file only once
for line in open(infile, 'r'):
  # parse through the grid-file
  if line[0] == '#': # skip comments
    continue
  if len(line.strip()) == 0:
    continue

  if line[0] == '[':
    sectinfo    = False
    sectnodes   = False
    sectelems   = False
    sectnodesbc = False

    if line.strip() == '[Info]':
      sectinfo = True
      print "\033[92mparsing header...\033[0m"
      continue

    if line.strip() == '[Nodes]':
      sectnodes = True
      print "\033[92mparsing nodes...\033[0m"
      continue

    if griddimension == 2:
      if line.strip() == '[2D Elements]':
        sectelems = True
        print " found " + str(nodectr) + " nodes"
        print "\033[92mparsing 2D elements...\033[0m"
        continue

    if griddimension == 3:
      if line.strip() == '[3D Elements]':
        sectelems = True
        print " found " + str(nodectr) + " nodes"
        print "\033[92mparsing 3D elements...\033[0m"
        continue

    if line.strip() == '[Node BC]':
      sectnodesbc = True
      print " found " + str(elemctr) + " elements"
      print "\033[92mparsing nodesbc...\033[0m"
      continue



######################
# file section: Info
######################
# mainly get grid dimension
  if sectinfo:
    if line.strip()[0:9] == 'Dimension':
      griddimension = int(line.split()[1])
      print " Dimension of grid: " + str(griddimension)


######################
# file section: Nodes
######################
  if sectnodes:
    ls = line.split()
    if not len(ls) == 4:
      print "error parsing nodes! offending line: " + line

    outstring = 'GRID'.ljust(8) + str(ls[0]).rjust(8) + '0'.rjust(8)
    for i in [1, 2, 3]:
      # convert to float
      d = float(ls[i])
      # format correctly
      o = '%3.4f' % d
      outstring += o.rjust(8)

    out.write(outstring + '\n')
    nodectr += 1


######################
# file section: Elements
######################
  if sectelems:
    # since the elements definitions' are spread to two lines,
    # we just remember the current line and skip to the next one
    if elemstmpline == '':
      elemstmpline = line 
      continue

    # at this point, tmpline must not be emtpy because it contains 
    # the data from the previous line which is needed here!
    if elemstmpline == '':
      print "error! tmpline is empty!"

    tls = elemstmpline.split() # contains info about element (number, type, numofnodes, region)
    ls = line.split() # contains node numbers forming the element

    # this check is mainly redundant because it should be correctly done by gid!
    if not (len(ls) == int(tls[2]) and len(tls) == 4):
      print "error parsing nodes! offending lines: " + elemstmpline + ' ' + line

    #############
    # 2D elements
    if griddimension == 2:
      # can only handle quad elements at the moment!
      outstring = 'CQUAD4'.ljust(8) + tls[0].rjust(8) + '1'.rjust(8)
      for i in range(4):
        outstring += str(ls[i]).rjust(8)
      out.write(outstring + '\n')


    #############
    # 3D elements
    if griddimension == 3:
      outstring = 'CHEXA'.ljust(8) + tls[0].rjust(8) + '1'.rjust(8)
      for i in range(6):
        outstring += str(ls[i]).rjust(8)

      outstring += '+'.ljust(8)
      out.write(outstring + '\n')

      outstring  = '+'.ljust(8)
      outstring += str(ls[6]).rjust(8)
      outstring += str(ls[7]).rjust(8)
      out.write(outstring + '\n')

    # make sure to delete contents of elemstmpline for next parse
    elemstmpline = ''
    elemctr += 1


######################
# file section: NodeBC
######################
  if sectnodesbc:
    num = line.split()[0]  # node number
    name = line.split()[1] # name given in grid-file
    # count the number of times a key appears in the grid-file
    if name in names.keys():
      names[name] += 1
    else:
      names[name] = 1
      print " new named entity in grid-file: " + name

    # check if named entity in the grid-file has been assigned
    if not name in table.keys():
      continue

    if table[name][0] == 'FORCE':
      out.write('%8s%8s%8s%40s\n' % ('FORCE'.ljust(8), str(table[name][1]).rjust(8), str(num).rjust(8), table[name][2]))
    if table[name][0] == 'SPC':
      out.write('%8s%8s%8s%16s\n' % ('SPC'.ljust(8), '1'.rjust(8), str(num).rjust(8), table[name][1]))

# here end the for loop over the lines of the input file


# check if all named entities in the grid-file have been assigned
for name in names.keys():
  if not name in table.keys():
    printWarning("name \033[93m" + name + "\033[0m in grid-file has not been assigned a bc/force!")
    continue

if griddimension == 2:
  out.write('PSHELL         1       1     1.0  \n')
  out.write('PSHELL         2       1     1.0  \n')
if griddimension == 3:
  out.write('PSOLID         1       1          \n')
  out.write('PSOLID         2       1          \n')
out.write('MAT1    1       2.87E-6 0.35    0.0     0.785E-5  12.E-6                +M1     \n')
out.write('+M1     100.    -100.   100.    \n')
out.write('ENDDATA\n')

out.close()

