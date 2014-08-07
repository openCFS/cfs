#!/usr/bin/python

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


################################################
## TODO
################################################


#### helper function
def printWarning(message):
  print "\033[91m WARNING:\033[0m " +  message


# takes a list of element numbers and creates the string needed
# for the nas-file
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
#dens = read_density_as_vector(densfile, "physical")
dens = read_density_as_vector(densfile)
# read the element numbers in a vector
nums = read_density_as_vector(densfile, "nr")




set1 = [] # will contain nums where dens > threshold
set2 = [] # will contain rest aka. void-region

# use density information to sort the element numbers into the 
# two sets
for i in range(len(dens)):
  if(dens[i] > threshold):
    set1.append(int(nums[i]))
  else:
    set2.append(int(nums[i]))

s1 = createSetString(set1, 1)
s2 = createSetString(set2, 2)


# now we have the correct strings for the sets
# we now need to adapt the nas-file

out = open(outfile, "w")
afterbulk = False
# numpoints = 0 # count number of points in the nas file
numelems = 0 # count number of elements in the nas file
for line in open(nasfile, "r"):
  # we can simply copy everything up to BEGIN BULK
  if line.strip() == 'BEGIN BULK':
    # now insert the two set strings
    out.write(s1)
    out.write('\n')
    out.write(s2)
    out.write('\n')
    out.write('BEGIN BULK\n')
    afterbulk = True
    continue


  # for check only: count number of points in the grid
  # if line[0:4] == 'GRID':
    # numpoints += 1
    #continue


  # 2D and 3D is basically the same
  # the definition of the elements in 3D spans over two lines.
  # however, the second line is not interesting and can simply be copied
  if line[0:6] == 'CQUAD4' or line[0:5] == 'CHEXA':
    # for check only: count number of elements in the grid
    numelems += 1

    ls = line.split()
    newline = line 

    # check if number is in set2 
    # set2 contains nums where dens <= threshold
    if int(ls[1]) in set2:
      newline = line[0:23]
      newline += "2"
      newline += line[24:len(line)]


    # newline is either the original line from the nas file if 
    # we do not have the void region 
    # else, we modify the line and put the element in the second (void) region
    out.write(newline)
    continue

  # if we arrive here, we have a line that does not need to 
  # be modified, so we can just output it
  out.write(line)


# now for the sanity check (can unfortunately only be done at the very end)
if not numelems == len(dens):
  printWarning("number of elements in nas file is not equal to number of elements in density file!")
  printWarning("nas: " + str(numelems) + " -- dens: " + str(len(dens)))
  printWarning("the newly created file \033[93m" + outfile + "\033[0m is probably incorrect!")
else:
  print "\033[93m Status:\033[0m process finished successfully (put " \
        + str(len(set2)) + " of " + str(numelems) + " elements in void region)"


out.close()
