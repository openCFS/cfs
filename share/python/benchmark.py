#!/usr/bin/env python
import argparse
from cfs_utils import *
from mesh_tool import *
import platform

# the purposee of this script is to assist in running cfs benchmarks. It contains the probems.
# The resulting .info.xml shall be analysed with performance,py

matstr = """\
<?xml version="1.0" encoding="utf-8"?>  
<cfsMaterialDataBase xmlns="http://www.cfs++.org/material">
<material name="99lines">
  <mechanical>
    <density>
      <linear>
        <real>1e-8</real>
     </linear>
   </density>
     <elasticity>
       <linear>
         <isotropic>
           <elasticityModulus>
             <real>1</real>
           </elasticityModulus>
           <poissonNumber>
             <real>0.3</real>
           </poissonNumber>
         </isotropic>
       </linear>
     </elasticity>
  </mechanical>
</material>
</cfsMaterialDataBase>
"""

mech3dstr = """\
<?xml version="1.0" encoding="UTF-8"?>
<cfsSimulation xmlns="http://www.cfs++.org/simulation">
  <fileFormats>
    <output>
      <hdf5 directory="."/>
      <info/>
    </output>
    <materialData file="bench_mat.xml" format="xml" />
  </fileFormats>
  <domain geometryType="3d">
    <regionList>
      <region material="99lines" name="mech" />
    </regionList>
  </domain>
  <sequenceStep index="1">
    <analysis>
      <static/>
    </analysis>
    <pdeList>
      <mechanic subType="3d">
        <regionList>
          <region name="mech" />
        </regionList>
        <bcsAndLoads>
           <fix name="left"> 
              <comp dof="x"/> 
              <comp dof="y"/> 
              <comp dof="z"/>
           </fix>
           <force name="bottom_right">
             <comp dof="y" value="-1"/>
           </force>
        </bcsAndLoads>
        <storeResults>
          <elemResult type="physicalPseudoDensity">
            <allRegions/>
          </elemResult>
          <nodeResult type="mechDisplacement">
            <allRegions/>
          </nodeResult>
        </storeResults>
      </mechanic>
    </pdeList>
   <linearSystems>
      <system>
        <solverList>
          <cholmod/>
        </solverList>
      </system>
    </linearSystems> 
  </sequenceStep> 
  <optimization>
    <costFunction type="compliance" task="minimize" >
      <stopping queue="10" value="0.001" type="relativeCostChange"/>
    </costFunction>
    <constraint type="volume" value=".3" bound="upperBound" linear="true" mode="constraint" />
    <optimizer type="optimalityCondition" maxIterations="5"/>
    <ersatzMaterial region="mech" material="mechanic" method="simp">
      <filters>
        <filter neighborhood="maxEdge" value="1.3" />
      </filters>
      <design name="density" initial=".3" physical_lower="1e-9" upper="1.0" />
      <transferFunction type="simp" application="mech" param="3"/>
      <export save="last" write="iteration"/>
    </ersatzMaterial>
    <commit mode="forward" stride="999"/>
  </optimization>
</cfsSimulation>
"""
shapemapstr = """\
<?xml version="1.0" encoding="UTF-8"?>
<cfsSimulation xmlns="http://www.cfs++.org/simulation">
  <fileFormats>
    <output>
      <hdf5 directory="."/>
      <info/>
    </output>
    <materialData file="bench_mat.xml" format="xml"/>
  </fileFormats>
  <domain geometryType="3d">
    <regionList>
      <region material="99lines" name="mech"/>
    </regionList>
    <nodeList>
      <nodes name="load">
        <coord x="1" y="0.5" z="0.5"/>
      </nodes>
    </nodeList>
  </domain>
  <sequenceStep index="1">
    <analysis>
      <static/>
    </analysis>
    <pdeList>
      <mechanic subType="3d">
        <regionList>
          <region name="mech"/>
        </regionList>
        <bcsAndLoads>
           <fix name="left"> 
              <comp dof="x"/> 
              <comp dof="y"/> 
              <comp dof="z"/>
           </fix>
           <force name="load" >
             <comp dof="y" value="-1"/>
           </force>
        </bcsAndLoads>
        <storeResults>
          <nodeResult type="mechDisplacement">
            <allRegions/>
          </nodeResult>
          <elemResult type="mechPseudoDensity">
            <allRegions/>
          </elemResult>
        </storeResults>
      </mechanic>
    </pdeList>
    <linearSystems>
      <system>
        <solverList>
         <cg>
         <!-- intentionally this cannot solve the system -->
          <maxIter>10</maxIter>
         </cg>
        </solverList>
      </system>
    </linearSystems>
  </sequenceStep>
  <optimization>
    <costFunction type="compliance" sequence="1" linear="false">
      <stopping queue="55" value="0.0001" type="relativeCostChange" maxHours="48"/>
    </costFunction>
    <constraint type="volume" bound="upperBound" value="0.4" design="density" linear="false"/>
    <constraint type="slope" bound="upperBound" value="1.1/nx" design="node" linear="true"/>
    <constraint type="slope" bound="upperBound" value="1.1/nx" design="profile" linear="true"/>
    <constraint type="curvature" bound="upperBound" value="0.3/nx" design="node" linear="true"/>
    <constraint type="curvature" bound="upperBound" value="0.3/nx" design="profile" linear="true"/>
    <optimizer type="evaluate" maxIterations="1"/>
    <ersatzMaterial region="mech" material="mechanic" method="shapeMap">
      <shapeMap beta="20" overlap="tanh_sum" enforce_bounds="false" integration_order="11" shape="tanh" integration_strategy="tailored">
       <center bottom_up_sym="mirror"  front_back_sym="mirror" left_right_sym="mirror"  diagonal_sym="mirror" >
         <node lower="0" upper="1" initial=".25" dof="y" />
         <node slave="true" dof="z"/>
       </center> 
       <profile lower=".05" upper=".3" initial=".1"/>
      </shapeMap>
      <design name="density" initial=".5" physical_lower="1e-3" upper="1.0"/>
      <transferFunction type="identity" application="mech"/>
      <export save="last" write="iteration" density="true"/>
    </ersatzMaterial>
    <commit mode="forward" stride="1"/>
  </optimization>
</cfsSimulation>
"""

parser = argparse.ArgumentParser(description="Run cfs benchmarks with permutations of CFS/OMP/MKL_THREADS=VECLIB_MAXIMUM_THREADS for different problems.\n" 
                                           + "Not all cfs executables have mkl, but Setting MKL_THREADS and equally VECLIB_MAXIMUM_THREADS (Apple) works for all.\n"
                                           + "Use --dry to see what would happen")
parser.add_argument("cfs", nargs='+', help="cfs binary and optional second value the comment of the binary (e.g. cfs gcc_mkl)")
parser.add_argument("--res", help="edge resolution of bulk3d_<res>.mesh to scale problem size", type=int, default='30')
parser.add_argument("--threads", nargs='+', help="space separated list of CFS/OMP/MKL_THREADS=VECLIB_MAXIMUM_THREADS to be tested, see --permutate", default=['1','4'])
parser.add_argument("--permutate", help="perform  full permuation of threads = more calcuations", action='store_true')
parser.add_argument("--dry", help="don't execute but print what shall be run", action='store_true')
parser.add_argument("--log", help="if given, write a summary plot file with the given name")
parser.add_argument("--hostname", help="overwrite hostname detection")
parser.add_argument("--skip_cholmod", help="skip the BLAS benchmark which requies USE_SUITESPARSE", action='store_true')
parser.add_argument("--skip_pardiso", help="skip the BLAS benchmark which requires Intel's MKL", action='store_true')
parser.add_argument("--skip_shapemap", help="skip the shape map non-blas benchmark", action='store_true')

args = parser.parse_args()
WIN = platform.system() == 'Windows'

if len(args.cfs) > 2:
  print('Error: usage is benchmark.py <executable> [<comment>] [--<...>]')
  sys.exit(1)

if args.permutate and WIN:
  print('Error: --permutate does not work on Windows')
  sys.exit(1)

comment = 'default' if len(args.cfs) == 1 else args.cfs[1]  

mat = open("bench_mat.xml","w")
mat.write(matstr)
mat.close()

# prepare files. Use the python tools from the same location as we have this script from
mesh = 'bulk3d_' + str(args.res) + '.mesh'
if not os.path.exists(mesh):
  mo = create_3d_mesh(args.res)
  write_ansys_mesh(mo, mesh)
   
host = args.hostname if args.hostname else platform.uname().node.split('.')[0]
# host can be mojo.mi.uni-erlangen.de or ca8c49d1-e336-417a-ba4e-bb348a961a55.fritz.box
if host.count("-") >= 2:
  host = 'local' 
tests = []    
if not args.skip_cholmod:
  f = open("bench_cholmod.xml","w")
  f.write(mech3dstr)
  f.close()
  tests.append('cholmod')
if not args.skip_pardiso:
  f = open("bench_pardiso.xml","w")
  f.write(mech3dstr.replace('<cholmod/>','<pardiso/>'))
  f.close()
  tests.append('pardiso')
if not args.skip_shapemap:
  f = open("bench_shapemap.xml","w")
  f.write(shapemapstr)
  f.close()
  tests.append('shapemap')
  
# the threads we run with on linux/mac:. 0=CFS_NUM_THREADS, 1=OMP_NUM_THREADS, 2=MKL_NUM_THREADS=VECLIB_MAXIMUM_THREADS
threads = [] 
t0 = args.threads[0] # there is at least one entry due to nargs='+'
for i, t in enumerate(args.threads):
  if i > 0 and args.permutate:
    threads.append([t0,t0,t]) # only MKL high
    threads.append([t0,t,t0]) # only OMP high
    threads.append([t0,t,t])  # only OMP and MKL high
    threads.append([t,t0,t0]) # CFS and MKL high
  threads.append([t,t,t]) # the only possibility on Windows

# here we optinally write the log file to
log = open(args.log, "w") if args.log else False

# here we store the results
result = {'test' : [], 'wall' : [], 'cpu' : [], 'CFS_NUM_THREADS' : [], 'OMP_NUM_THREADS' : [], 'MKL/VECLIB_THREADS' : [], 'res' : []}
if comment != 'dummy':
  result['comment'] = []

for tst in tests:      
  for t in threads:
    if tst == 'shapemap' and t[2] != t0 and not (t[0] == t[1] == t[2]): # skip playing with mkl, we don't need it
      continue 

    problem = 'bench_' + host + '-' + comment + '-' + tst + '-res_' + str(args.res) + '-cfs_' + t[0] + '-omp_' + t[1] + '-mkl_' + t[2]
    env = 'CFS_NUM_THREADS=' + str(t[0]) + ' OMP_NUM_THREADS=' + str(t[1]) + ' MKL_NUM_THREADS=' + str(t[2]) + ' VECLIB_MAXIMUM_THREADS=' + str(t[2]) + ' '
    flag = ' -t ' + str(t[0]) if WIN else '' # on Windows no env and all threads same as we have no permutation
    cmd = '' if WIN else env
    cmd += args.cfs[0] + flag + ' -q -m ' + mesh + ' -p bench_' + tst + '.xml ' + problem
 
    if args.dry:
      print(cmd)
    else:
      execute(cmd,output=True)
      xml = open_xml(problem + '.info.xml')
      wall = xpath(xml, '//cfsInfo/summary/timer/@wall')
      cpu = xpath(xml, '//cfsInfo/summary/timer/@cpu')
      print('wall', wall, 'cpu', cpu, problem)
      print()
      result['test'].append(tst)
      result['wall'].append(wall)
      result['cpu'].append(cpu)
      result['CFS_NUM_THREADS'].append(str(t[0]))
      result['OMP_NUM_THREADS'].append(str(t[1]))
      result['MKL/VECLIB_THREADS'].append(str(t[2]))
      result['res'].append(str(args.res))
      if comment != 'dummy':
        result['comment'].append(comment)
 
# present the whole data 
hdr = '# host:' + host + ' mesh:' + str(args.res) + ' comment:' + comment + '\n'
hdr += '# '
for i, k in enumerate(result.keys()):
  hdr += str(i+1) + ':' + k + ' ' 
if log:
  print("write log file '" + args.log + "'")
  log.write(hdr + '\n')    
  
print(hdr)
  
for i in range(len(result['test'])):
  line = ''
  for v in result.values():
    line += v[i] + ' \t'   
  print(line)
  if log:
    log.write(line + '\n')    
  
    
