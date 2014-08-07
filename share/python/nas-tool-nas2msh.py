#!/usr/bin/python

import sys


# materials:
# 1: where the structure is
# 2: void
# 3: where the force is applied
materials  = 'm        1      ISO    FIXED       1.000000000000E+00       3.428700000000E-01\n'
materials += 'm        2      ISO    FIXED       2.870000000000E-06       3.500000000000E-01\n'
materials += 'm        3      ISO    FIXED       5.000000000000E+00       3.500000000000E-01\n'

# elems where the force is applied
forceelems = [360,390,420,450,480,510,540,570]




if len(sys.argv) < 3:
  print "usage: " + sys.argv[0] + " <nas-file> <outfile>"
  print "\033[93mInfo:\033[0m"
  print " - <nas-file>:    the nas input file (from plato-n)"
  print " - <outfile>:     the name of the output file for the new nas-file"
  sys.exit(1)

nasfile = sys.argv[1]
outfile = sys.argv[2]


################################################
## TODO
################################################


#### helper function
def printWarning(message):
  print "\033[91m WARNING:\033[0m " +  message


def readSet(s, l, nf):
  ls = l[7:].split(',')
  ls.remove('\n')
  for it in ls:
    s.append(int(it))
  line = nf.readline()
  while(line[0:4] == '    '):
    ls = line.split(',')
    if '\n' in ls:
      ls.remove('\n')
    for it in ls:
      s.append(int(it))
    line = nf.readline()
  return line


set1 = [] # contains set1 from nasfile
set2 = [] # contains set2 from nasfile

# read nasfile
nf = open(nasfile)
line = "hh"
while line: # search for SET1
  line = nf.readline()
  if(line[0:5] == 'SET 1'):
    ls = line[7:].split(',')
    ls.remove('\n')
    for it in ls:
      set1.append(int(it))
    line = nf.readline()
    while(line[0:4] == '    '):
      ls = line.split(',')
      if '\n' in ls:
        ls.remove('\n')
      for it in ls:
        set1.append(int(it))
      line = nf.readline()
    break
  else: continue

if(line[0:5] == 'SET 2'):
  ls = line[7:].split(',')
  ls.remove('\n')
  for it in ls:
    set2.append(int(it))
  line = nf.readline()
  while(line[0:4] == '    '):
    ls = line.split(',')
    if '\n' in ls:
      ls.remove('\n')
    for it in ls:
      set2.append(int(it))
    line = nf.readline()
else:
  printWarning("wrong file format?")

out = open(outfile, "w")
sss = '$ PLATO-N MD Library 0.9\n$\n$\n$ Mesh file \n$\n$ c      1681 \n$ q      1600 \n$ m         2    \n$ grq       2    \n$ grc       0    \n$\n$ a         1         1\n$\n'
out.write(sss)

# lines look like this:
# GRID           1       0  1.0000  1.0000  0.0000
line = nf.readline()
ls = line.split()
while ls[0].strip() == 'GRID':
  lineout = 'c ' + ls[1].rjust(8) + '    '
  lineout += "%1.15E" % float(ls[3]) + '    '
  lineout += "%1.15E" % float(ls[4]) + '    '
  lineout += "%1.15E" % float(ls[5]) + '\n'
  out.write(lineout)


  line = nf.readline()
  ls = line.split()


out.write('$\n')


ls = line.split()
while ls[0].strip() == 'CQUAD4':
  lineout = 'q ' + ls[1].rjust(8) + '   CQUAD4 MEMBRANE'
  lineout += "%9d" % int(ls[3])
  lineout += "%9d" % int(ls[4])
  lineout += "%9d" % int(ls[5])
  lineout += "%9d" % int(ls[6])

  lineout += '  1.0000E+00        '

  numm = int(ls[1])
  if numm in forceelems:
    lineout += '3\n'
  else:
    if numm in set1:
      lineout += '1\n'
    else:
      lineout += '2\n'


  out.write(lineout)
  line = nf.readline()
  ls = line.split()

nf.close()
  
out.write('$\n')
out.write(materials)
out.write('$\n')

for i in set1:
  if not i in forceelems:
    out.write('grq        1%9d\n' % i)
for i in set2:
  out.write('grq        2%9d\n' % i)
for i in forceelems:
  out.write('grq        2%9d\n' % i)

out.write('$\n')

out.close()

