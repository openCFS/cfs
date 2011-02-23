#!/usr/bin/env python

import sys
from optimization_tools import *



if len(sys.argv) < 5:
  print "usage: " + sys.argv[0] + " <nas-file> <density.xml> <outfile> <threshold>"
  print "\033[93mInfo:\033[0m"
  print " - <nas-file>:    the nas input file (from plato-n)"
  print " - <density.xml>: the density file from the cfs-simulation"
  print " - <outfile>:     the name of the output file for the new nas-file"
  print " - <threshold>:   a double value specifying where the density is cut off"
  sys.exit(1)

nasfile = sys.argv[1]
densfile = sys.argv[2]
outfile = sys.argv[3]
threshold = float(sys.argv[4])



#### helper function
def printWarning(message):
  print "\033[91m WARNING:\033[0m " +  message


def createSetString(set1, index):
  s1 = 'SET %d =' % index
  s1 = s1.ljust(8)
  l1 = len(set1)
  for i in range(l1):
    s1 += str(set1[i]).rjust(7)
    if i < l1-1:
      s1 += ','
    if i%5 == 4:
      s1 += '\n'
      s1 += ' '.ljust(8)
  return s1

# read the density values in a vector
dens = read_density_as_vector(densfile)
# read the element numbers in a vector
nums = read_density_as_vector(densfile, True)


set1 = [] # will contain nums where dens < threshold
set2 = [] # will contain rest

for i in range(len(dens)):
  if(dens[i] < threshold):
    set1.append(int(nums[i]))
  else:
    set2.append(int(nums[i]))

s1 = createSetString(set1, 1)
s2 = createSetString(set2, 2)


# now we have the correct strings for the sets
# we now need to adapt the nas-file

out = open(outfile, "w")
afterbulk = False
for line in open(nasfile, "r"):
  # we can simply copy everything up to BEGIN BULK
  if line.strip() == 'BEGIN BULK':
    # now insert the two set strings
    out.write('$HMSET        1        2 "xs_3D"\n')
    out.write(s1)
    out.write('\n')
    out.write('$HMSET        1        2 "xs_3D"\n')
    out.write(s2)
    out.write('\n')
    out.write('BEGIN BULK\n')
    afterbulk = True
    continue

  #if not afterbulk:
    #out.write(line)
    #continue

  # now we are after bulk and need to adapt the element regions
  if line[0:4] == 'GRID':
    out.write(line)
    continue


  # only 2D for now!!!
  if line[0:6] == 'CQUAD4':
    ls = line.split()
    newline = line 

    # check if number is in set1 
    # set1 contains nums where dens < threshold
    if int(ls[1]) in set1:
      newline = line[0:23]
      newline += "2"
      newline += line[24:len(line)]

    out.write(newline)
    continue

  # if we arrive here, we have a line that does not need to 
  # be modified, so we can just output it
  out.write(line)




out.close()
