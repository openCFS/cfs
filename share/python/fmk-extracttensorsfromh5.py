#!/usr/bin/python

import os, sys, h5py, numpy, numpy.linalg
from numpy import sqrt, dot, sin, cos


def print_usage():
  print "usage: " + sys.argv[0] + " <mode> <fmo-h5-file> [<nr of rotangles> (mode tocfs)][<parametrization> (optional in mode tocfs)]"
  print "<mode> can be:"
  print " - tocfs    extract stiffness directions and write cfs density.xml-file"
  print "            additional parameters: <nr of rotangles>"
  print " - fmk:     produce txt-file suitable as warmstart file for fmk"
  print " - invhom:  produce txt-file for inverse homogenization"
  print " - fromcfs: read tensors from cfs h5-file and produce txt-file for fmk warmstart"
  print "<parametrization> is optional and can be:"
  print " - default:    parametrization [s1, 0, 0; 0, s2, 0; 0, 0, g]"
  print " - trans-iso:  parametrization [s1/(1-s12), sqrt(s12*s1*s2)/(1-s12), 0; sqrt(s12*s1*s2)/(1-s12), s2/(1-s12), 0; 0, 0, g]"
  print " - ortho:      parametrization [s1^2+s12^2, s1*s12+s12*s2, 0; s1*s12+s12*s2, s12^2+s2^2, 0; 0, 0, g]"
  





def read_cfs_h5(infile, outfile):
# read tensor entries from the h5-file
  f = h5py.File(infile, "r")
  data = f['Results/Mesh/MultiStep_1']
  laststep = len(data) - 2
  data = f['Results/Mesh/MultiStep_1/Step_' + str(laststep) + '/mechTensor/mech/Elements/Real']

  trace = 0.0
  out = open(outfile, 'w')
  for t in range(len(data)):
    tensor = data[t]
    # "e11", "e22", "e33", "e23", "e13", "e12";
    #    0      1      2      3      4      5
    tensor[2] *= 2.0
    tensor[3] *= sqrt(2.0)
    tensor[4] *= sqrt(2.0)
    trace += tensor[0] + tensor[1] + tensor[2]
    ts = map(str, tensor)
    out.write(ts[0] + '\n')
    out.write(ts[5] + '\n')
    out.write(ts[1] + '\n')
    out.write(ts[4] + '\n')
    out.write(ts[3] + '\n')
    out.write(ts[2] + '\n')
  print "trace = " + str(trace/len(data))
  f.close()
  return

  



def read_fmk_h5(infile):
# read tensor entries from the h5-file
  f = h5py.File(infile, "r")
  lastindex = f.keys()[-1]
  dim     = int(f[lastindex].get("ElementDimension/ElementDimension")[0])
  print "dimension of dataset = " + str(dim)
  tensors = list(f[lastindex].get("MaterialTensor/MaterialTensor"))
  complianceinit = list(f["0/Compliance/Compliance"])
  compliance = list(f[str(lastindex) + "/Compliance/Compliance"])
  f.close()

  numtens = len(tensors)/21.0

  if not numtens == round(numtens):
    print "\033[91m WARNING:\033[0m"
    print " the number of tensor entries could be wrong"
    print " found 21 entries, corresponding to " + str(numtens) + " tensors!"

  if not (dim == 2 or dim == 3):
    print "\033[91m ERROR:\033[0m"
    print " dimension of dataset must be 2 or 3!"
    sys.exit(1)


  lsstrlen = len(str(lastindex))
  print " compliance[" + "0".rjust(lsstrlen) + "] = " + str(complianceinit[0])
  print " compliance[" + str(lastindex) + "] = " + str(compliance[0])

  return tensors, dim, int(numtens)




def find_max_stiffness(tens, t, sint, cost, numsteps):
  # stiffness is the tensors entry 1,1
  # e11 e12 e13 e22 e23 e33
  # 0	1   2   3   4   5  ?????????
  stiffness  = tens[0]*cost**4.0 + tens[3]*sint**4.0
  stiffness += (2.0*tens[5] + 2.0*tens[1]) * cost**2.0 * sint**2.0
  stiffness += 2.0 * sqrt(2.0) * tens[2] * cost**3.0 * sint
  stiffness += 2.0 * sqrt(2.0) * tens[4] * cost * sint**3.0

  s1 = 0.0
  ind1 = 0
  for i in range(numsteps):
    if stiffness[i] > s1:
      s1 = stiffness[i]
      ind1 = i

  a1 = t[ind1]

  s2  = tens[0]*sin(a1)**4.0 + tens[3]*cos(a1)**4.0
  s2 += 2.0*(tens[1]+tens[5])*cos(a1)**2*sin(a1)**2
  s2 -= 2.0*sqrt(2)*tens[2]*cos(a1)*sin(a1)**3
  s2 -= 2.0*sqrt(2)*tens[4]*cos(a1)**3*sin(a1)

  g = tens[5]+.5*tens[0]*sin(2*a1)**2-tens[1]*sin(2*a1)**2
  g += .5*tens[3]*sin(2*a1)**2 - tens[5]*sin(2*a1)**2
  g += -.5*sqrt(2)*tens[2]*sin(4*a1) + .5*sqrt(2)*tens[4]*sin(4*a1)
  g *= .5
  
 
  alpha = s1
  beta = tens[1] + (tens[0]*sin(2*a1)**2)/4 - (tens[1]*sin(2*a1)**2)/2 + (tens[3]*sin(2*a1)**2)/4 - (tens[5]*sin(2*a1)**2)/2 - (2**.5*tens[2]*sin(4*a1))/4 + (2**.5*tens[4]*sin(4*a1))/4
  gamma = s2
  s12 = 0
  if parametrization == "ortho":
    lowereigenbound = 0.01
    g = max(0, -lowereigenbound + g)
    s12 = 0
    if abs(beta) > 1e-8:
      s1 = -lowereigenbound + (((2*beta**2*(- beta**2 + alpha*gamma)**.5 - alpha*beta**2 - 2*alpha*gamma**2 + alpha**2*gamma + 3*beta**2*gamma + gamma**3)/(alpha**2 - 2*alpha*gamma + 4*beta**2 + gamma**2))**.5*(beta**2 - alpha*(- beta**2 + alpha*gamma)**.5 + gamma*(- beta**2 + alpha*gamma)**.5))/(beta**2 + gamma**2 - alpha*gamma)
      s2 = -lowereigenbound + (gamma - (alpha*beta**2 - 2*beta**2*(- beta**2 + alpha*gamma)**.5 + beta**2*gamma)/(alpha**2 - 2*alpha*gamma + 4*beta**2 + gamma**2))**.5
      s12 = (beta*(gamma - (- beta**2 + alpha*gamma)**.5)*((2*beta**2*(- beta**2 + alpha*gamma)**.5 - alpha*beta**2 - 2*alpha*gamma**2 + alpha**2*gamma + 3*beta**2*gamma + gamma**3)/(alpha**2 - 2*alpha*gamma + 4*beta**2 + gamma**2))**.5)/(beta**2 + gamma**2 - alpha*gamma)
      if s1 < 0 or s2 < 0:
        s1 = -lowereigenbound + ((-(2*beta**2*(- beta**2 + alpha*gamma)**.5 + alpha*beta**2 + 2*alpha*gamma**2 - alpha**2*gamma - 3*beta**2*gamma - gamma**3)/(alpha**2 - 2*alpha*gamma + 4*beta**2 + gamma**2))**.5*(beta**2 + alpha*(- beta**2 + alpha*gamma)**.5 - gamma*(- beta**2 + alpha*gamma)**.5))/(beta**2 + gamma**2 - alpha*gamma)
        s2 = -lowereigenbound + (gamma - (2*beta**2*(- beta**2 + alpha*gamma)**.5 + alpha*beta**2 + beta**2*gamma)/(alpha**2 - 2*alpha*gamma + 4*beta**2 + gamma**2))**.5
        s12 = (beta*(-(2*beta**2*(- beta**2 + alpha*gamma)**.5 + alpha*beta**2 + 2*alpha*gamma**2 - alpha**2*gamma - 3*beta**2*gamma - gamma**3)/(alpha**2 - 2*alpha*gamma + 4*beta**2 + gamma**2))**.5*(gamma + (- beta**2 + alpha*gamma)**.5))/(beta**2 + gamma**2 - alpha*gamma)
      s1 = max(0, s1)
      s2 = max(0, s2)
    else:
      s1 = s1**.5
      s2 = s2**.5
  elif parametrization == "trans-iso" and beta > 0:# and beta < 0.45:
   # if gamma > 1e-2:
    s1 = alpha - beta**2/gamma
    #if alpha > 1e-2:
    s2 = gamma - beta**2/alpha
    #if alpha > 1e-2 or gamma > 1e-2:
    s12 = beta**2/(alpha*gamma)

    
  a1 = numpy.pi - t[ind1]
  a2 = a1 + 0.5*numpy.pi
  if a2 > numpy.pi:
    a2 -= numpy.pi

  a1 = ((a1+0.94)% numpy.pi)-0.94 # one load at lower right corner
 # a1 = ((a1+numpy.pi-1.7105)% numpy.pi)+1.7105-numpy.pi
  a1 = ((a1+.5*numpy.pi-1.625)% numpy.pi)+1.625-.5*numpy.pi # two loads at top and bottom
    
  return s1, a1, s2, g, s12



def to_cfs_xml_file(res, outfile, numtens):
  out = open(outfile, 'w')
  out.write('<?xml version="1.0"?>\n<cfsErsatzMaterial>\n')
  out.write(' <header>\n  <design constant="false" initial="0.5" name="density" region="all" scale="false"/>\n')
  out.write('  <transferFunction application="mech" design="density" param="1.0" type="simp"/>\n')
  out.write(' </header>\n <set id="1">\n')
  if parametrization == "trans-iso":
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="emodul-iso" design="%0.10e"/>\n' % res[i, 0])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="rotAngle" design="%1.4f"/>\n' % res[i, 1])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="emodul" design="%0.10e"/>\n' % res[i, 2])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="gmodul" design="%1.4f"/>\n' % res[i, 3])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="poisson" design="%1.4f"/>\n' % res[i, 4])
  elif parametrization == "ortho":
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="tensor11" design="%0.10e"/>\n' % res[i, 0])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="rotAngle" design="%1.4f"/>\n' % res[i, 1])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="tensor22" design="%0.10e"/>\n' % res[i, 2])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="tensor33" design="%1.4f"/>\n' % res[i, 3])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="tensor12" design="%1.4f"/>\n' % res[i, 4])
  else:
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="emodul-iso" design="%0.10e"/>\n' % res[i, 0])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="rotAngle" design="%1.4f"/>\n' % res[i, 1])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="emodul" design="%0.10e"/>\n' % res[i, 2])
    for i in range(numtens):
      out.write('  <element nr="' + str(i+1) + '" type="gmodul" design="%1.4f"/>\n' % res[i, 3])
      
      
  out.write(' </set>\n</cfsErsatzMaterial>\n')
  out.close()



def handle_mode_tocfs(outfile, tensors, dim, numtens, steps):
  # define a range of angles, make sure that pi is included
  t = numpy.arange(0.0, numpy.pi+1e-4, numpy.pi/steps)
  numsteps = len(t)
  cost = cos(t)
  sint = sin(t)
  res = numpy.zeros((numtens, 4))
  if parametrization == "ortho" or parametrization == "trans-iso":
    res = numpy.zeros((numtens, 5))
    
  for n in range(numtens):
    tens = tensors[n*21:(n+1)*21]
    s1, a1, s2, g, s12 = find_max_stiffness(tens, t, sint, cost, numsteps)
    res[n, 0] = s1
    res[n, 1] = a1
    res[n, 2] = s2
    res[n, 3] = g
    if parametrization == "ortho" or parametrization == "trans-iso":
      res[n, 4] = s12

  to_cfs_xml_file(res, outfile, numtens)









def handle_mode_fmk(outfile, tensors, dim, numtens):
# we must sort the entries obtained from the h5-file in a different manner
# to reread the correct values, we must transpose the matrix of entries
# in 2D, only entries 3 and 4 are swapped
# in 3D, things are a bit more complicated... ;)
  out = open(outfile, "w")
  for n in range(numtens):
    tmplist = map(str, tensors[n*21:(n+1)*21])
    if dim == 2:
      out.write(tmplist[0] + "\n")
      out.write(tmplist[1] + "\n")
      out.write(tmplist[3] + "\n")
      out.write(tmplist[2] + "\n")
      out.write(tmplist[4] + "\n")
      out.write(tmplist[5] + "\n")

    if dim == 3:
      out.write(tmplist[0] + "\n")
      out.write(tmplist[1] + "\n")
      out.write(tmplist[6] + "\n")
      out.write(tmplist[2] + "\n")
      out.write(tmplist[7] + "\n")

      out.write(tmplist[11] + "\n")
      out.write(tmplist[3] + "\n")
      out.write(tmplist[8] + "\n")
      out.write(tmplist[12] + "\n")
      out.write(tmplist[15] + "\n")

      out.write(tmplist[4] + "\n")
      out.write(tmplist[9] + "\n")
      out.write(tmplist[13] + "\n")
      out.write(tmplist[16] + "\n")
      out.write(tmplist[18] + "\n")

      out.write(tmplist[5] + "\n")
      out.write(tmplist[10] + "\n")
      out.write(tmplist[14] + "\n")
      out.write(tmplist[17] + "\n")
      out.write(tmplist[19] + "\n")

      out.write(tmplist[20] + "\n")

  out.close()


def handle_mode_invhom(outfile, tensors, dim, numtens):
  rescale = True
  out = open(outfile, "w")
  for n in range(numtens):
    tmplist = tensors[n*21:(n+1)*21]
    if rescale and dim == 2:
      tmplist[2] /= sqrt(2.0)
      tmplist[4] /= sqrt(2.0)
      tmplist[5] /= 2.0

    if rescale and dim == 3:
      print "ERROR! rescaling not implemented"
      sys.exit(1)

    out.write(' '.join( map(str, tmplist)) + "\n")

  out.close()







if len(sys.argv) < 3:
  print_usage()
  sys.exit(1)

mode = sys.argv[1]

# define infile and check for existence
infile = sys.argv[2]
if not os.path.isfile(infile):
  print "infile " + infile + " is not a real file!"
  sys.exit(1)
else:
  print "      reading infile -- " + infile

outfile = ""
cfssteps = 0

if mode == 'tocfs':
  print 'mode tocfs'
  outfile = infile + '.dens.xml'
  if len(sys.argv) < 4:
    print '\033[91m ERROR:\033[0m need to know number of steps for discrete rotation angles, exiting.'
    print 'call: ' + sys.argv[0] + ' <mode> <fmo-h5-file> <nr of rotangles>'
    sys.exit(1)
  cfssteps = int(sys.argv[3])
  if len(sys.argv) < 5:
    parametrization = "default"
  else:
    parametrization = sys.argv[4]
  
elif mode == 'fmk':
  print 'mode fmk'
  outfile = infile + '.fmk.txt'
elif mode == 'invhom':
  print 'mode invhom'
  outfile = infile + '.invhom.txt'
elif mode == 'fromcfs':
  outfile = infile + '.fromcfs.txt'
else:
  print 'mode ' + mode + ' not supported, exiting.'
  print_usage()
  sys.exit(1)


if outfile == '':
  print 'no outfile name set, exiting.'
  sys.exit(1)
else:
  print '  writing to outfile -- ' + outfile








if mode == 'tocfs':
  tensors, dim, numtens = read_fmk_h5(infile)
  handle_mode_tocfs(outfile, tensors, dim, numtens, cfssteps)
elif mode == 'fmk':
  tensors, dim, numtens = read_fmk_h5(infile)
  handle_mode_fmk(outfile, tensors, dim, numtens)
elif mode == 'invhom':
  tensors, dim, numtens = read_fmk_h5(infile)
  handle_mode_invhom(outfile, tensors, dim, numtens)
elif mode == 'fromcfs':
  # read cfs h5 file and convert to fmk txt file
  read_cfs_h5(infile, outfile)
