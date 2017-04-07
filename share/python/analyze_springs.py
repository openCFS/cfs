#!/usr/bin/env python

import os.path
from hdf5_tools import *
from cfs_utils import *
import argparse
from optimization_tools import *
import collections
import os.path

def read_spring_data(file,xml, spring, dic,ns,fix_force,detailed):
  # reads and calculates spring data
  sq_out = []
  max_f_tension = 0.
  max_f_shear = 0.
  max_u_spring_tension = 0.
  max_u_spring_shear = 0.
  max_u_spring = 0.
  disp = numpy.zeros((ns,3))
  scale = 0.5 if fix_force else 1.
  for i in range(ns):
    dx = scale * float(xpath(xml,'//result[@location="sup' + str(i+1) + '"]/item[last()]/@x'))
    dy = scale * float(xpath(xml,'//result[@location="sup' + str(i+1) + '"]/item[last()]/@y'))
    dz = scale * float(xpath(xml,'//result[@location="sup' + str(i+1) + '"]/item[last()]/@z'))
    disp[i][:] = [dx,dy,dz]
    sq_out.append(dx**2 + dy**2 + dz**2)
    max_f_tension = max(abs(max_f_tension),abs(spring*dy))
    max_f_shear = max(abs(max_f_shear),abs(numpy.sqrt((spring/3.57143*dx)**2 + (spring/3.57143*dz)**2)))
    max_u_spring = max(abs(max_u_spring), numpy.sqrt(dx**2+ dy**2 +dz**2))
    max_u_spring_tension = max(abs(max_u_spring_tension),abs(dy))
    max_u_spring_shear = max(abs(max_u_spring_shear),numpy.sqrt(dx**2+dz**2))
  a_all = max(sq_out)/min(sq_out)
  a_calc = sq_out[ns-1]/sq_out[0]
  f = numpy.sqrt((spring/3.57143 * disp[:,0])**2 + (spring * disp[:,1])**2 + (spring/3.57143 * disp[:,2])**2)
  dic['f_spring_sum'] = sum(f)
  if detailed:
    dic['max_u_spring'] =  max_u_spring
    dic['max_u_spring_tension'] = max_u_spring_tension
    dic['max_u_spring_shear'] = max_u_spring_shear
    dic['a_calc'] = a_calc
    dic['a_all'] = a_all
    dic['q_calc'] = f[ns-1]/f[0] 
  
  dic['max_f_tension'] = max_f_tension
  dic['max_f_shear'] =  max_f_shear
  dic['q_all'] = max(f)/min(f)
   
  

def process(file, sort,spring,fix_force,detailed):
  xml = open_xml(file)
  problem = file[:-9]
  if has(xml,'//iteration[last()]/@compliance'):
    dic = collections.OrderedDict()#{'total_wall' : float(xpath(xml,'/cfsInfo/summary/timer/@wall')), 'total_cpu' : float(xpath(xml,'/cfsInfo/summary/timer/@cpu'))})
    dic['obj'] = float(xpath(xml,'//iteration[last()]/@compliance'))
    dic['iter'] = int(xpath(xml,'//iteration[last()]/@number'))
    if has(xml,'//iteration[last()]/@volume'):
      dic['vol'] = float(xpath(xml,'//iteration[last()]/@volume'))
    else:
      dic['vol'] = float(xpath(xml,'//iteration[last()]/@globalTwoScaleVolume_allDesigns'))
    print(problem + ".h5")
    if os.path.exists(problem+".h5"):
      f = h5py.File(problem+".h5", 'r')
    else:
      print(problem+".h5 not found!")
      return None
    u_design = read_displacement(f,'design')
    u_ndesign = read_displacement(f,'non_design')
    scale = 0.5 if fix_force else 1.
    dic['max_u'] = scale * max(max(numpy.sqrt(u_design[:,1]**2+u_design[:,0]**2 + u_design[:,2]**2)),max(numpy.sqrt(u_ndesign[:,1]**2+u_ndesign[:,0]**2 + u_ndesign[:,2]**2)))
    read_spring_data(file,xml, spring, dic,5,fix_force,detailed)
    if has(xml,'//iteration[last()]/@sf_a_by_s'):
      dic['a'] = float(xpath(xml,'//iteration[last()]/@sf_a_by_s'))
    else:
      dic['a'] = -1.
    dic['unknowns'] = int(xpath(xml,'//system/@totalNumRows'))
    dic['problem'] = problem
    #dic[sort] = float(xpath(xml,'//ssor/@omega')) 
    return dic
  else:
    print("No optimization found for " + str(file))
    return None

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="info xml files (wildcards ok) to be ckecked")
parser.add_argument('--sort', help="how to sort the results (header label)", required=False)
parser.add_argument('--ref',help="reference info.xml without optimization",action="store_true")
parser.add_argument('--fix_force',help="fix force",action="store_true")
parser.add_argument('--d',help="detailed output")


args = parser.parse_args()

if len(args.input) == 0:
  print('give some .info.xml files as input')
  os.sys.exit(1)

res = []
spring = 100000.
if args.ref:
  xml = open_xml(args.input[0])
  dic = collections.OrderedDict()
  
  read_spring_data(args.input[0],xml, spring, dic,10,False,args.d)
  res.append(dic)
else:
  for file in args.input:
    #print("file = " + str(file))
    data =  process(file, args.sort,spring,args.fix_force,args.d)
    if data != None:
      res.append(data)

if args.sort:
  res = sorted(res, key=lambda k: k[args.sort])

header = "#"
for idx, key in enumerate(res[0].keys()):
  header += '(' + str(idx+1) + ') ' + key + ' \t'
print(header)
for item in res:
  line = ""
  for idx, key in enumerate(res[0].keys()):
    line += str(item[key])  + ' \t'
  print(line)



#args = parser.parse_args()


#if len(args.input) == 0:
#  print('usage: give .info.xml files as input')
#  os.sys.exit(1)
#for file in args.input:
#  output = process_file(file, args, [], recursive=False)
#  print("u_max \t max_spring_tension \t max_spring_disp")
#  problem = file.split('.info.xml')[0] 
#  f = h5py.File(problem+".h5", 'r')
#  print "problem = " + str(problem)
#  u_design = read_displacement(f,'design')
#  u_ndesign = read_displacement(f,'non_design')
#  u_max = max(max(numpy.sqrt(u_design[:,1]**2+u_design[:,0]**2 + u_design[:,2]**2)),max(numpy.sqrt(u_ndesign[:,1]**2+u_ndesign[:,0]**2 + u_ndesign[:,2]**2)))
#  xml = open_xml(problem+".info.xml")
#  file = open(problem+".spring_forces.csv","w")
#  file.write("x,y,z,dx,dy,dz,sx,sy,sz\n")
#  max_f_tension = 0.
#  max_u_tension = 0.
#  max_u_spring = 0.
#  spring = 100000.
#  for i in range(5): 
#          dx = float(xpath(xml,'//result[@location="sup1"]/item[last()]/@x'))
#          dy = float(xpath(xml,'//result[@location="sup1"]/item[last()]/@y'))
#          dz = float(xpath(xml,'//result[@location="sup1"]/item[last()]/@z'))
#          coord = xpath(xml,'//nodes[@name="sup' + str(i+1) + '"]/@coord')
#          vec = coord[1:len(coord)-1].split(',')
#          file.write(vec[0] + ',' + vec[1] + ',' + vec[2] + ','+ str(dx) + ',' + str(dy) + ',' + str(dz) + ',' + str(spring/3.57143*dx) + ',' + str(spring*dy) + ',' + str(spring/3.57143 * dz) +  '\n')
#          max_f_tension = max(abs(max_f_tension),abs(spring*dy))
#          max_u_tension = max(abs(max_u_tension),abs(dy))
#          max_u_spring = max(abs(max_u_spring),numpy.sqrt(dx**2+dy**2+dz**2))
#    
#  file.close()
#  print(str(u_max)+ " \t"+str(max_f_tension) + " \t"+ str(max_u_tension) + " \t" + str(max_u_spring))
