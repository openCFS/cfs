#!/usr/bin/env python
import os
import argparse
import platform
import cfs_utils as cu
import mesh_tool as mt


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
    <costFunction type="compliance" task="minimize" />
    <constraint type="volume" value=".3" bound="upperBound" linear="true" mode="constraint" />
    <optimizer type="ocm" maxIterations="5"/>
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
    <costFunction type="compliance" />
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
# generate with -c a bench_config,cmake file, parse it, and return the content as dict with bool values
def get_config(cfs): 
  cu.execute(cfs + ' -c bench_config.cmake', output = False)
  res = {}
  with open('bench_config.cmake') as f:
    for line in f:  
      # either a comment line with # or stuff like set(USE_PARDISO OFF)
      if line.startswith('set('):
        assert ' ' in line and ')' in line
        si = line.index(' ', 4)
        key = line[4:si]
        value = line[si+1:line.index(')',si)]
        assert value == 'ON' or value == 'OFF'
        res[key] = True if value == 'ON' else False
  return res 
        
# return the available cases based on the config dict
def get_cases(config):
  available = ['shapemap'] # even in minimal build
  if config['USE_SUITESPARSE']:
    available.append('cholmod')
  if config['USE_PARDISO']:
    available.append('pardiso')
  if config['USE_GINKGO']:  
    available.append('ginkgo')
  if config['USE_GINKGO'] and config['USE_CUDA']:        
    available.append('ginkgo-cuda')
  return available  

def prepare_tests(tests):
  for t in tests:
    f = open("bench_" + t + ".xml","w")
    if t == 'shapemap':
      f.write(shapemapstr)
    else:
      if t == 'ginkgo-cuda':
        f.write(mech3dstr.replace('<cholmod/>','<ginkgo cuda="on"/>'))
      elif t == 'ginkgo':
        f.write(mech3dstr.replace('<cholmod/>','<ginkgo cuda="off"/>'))
      else:            
        f.write(mech3dstr.replace('<cholmod/>','<' + t + '/>'))
    f.close()
    
parser = argparse.ArgumentParser(description="Run cfs benchmarks with permutations of CFS/OMP/MKL_THREADS=VECLIB_MAXIMUM_THREADS for different problems.\n" 
                                           + "Not all cfs executables have mkl, but Setting MKL_THREADS and equally VECLIB_MAXIMUM_THREADS (Apple) works for all.\n"
                                           + "Use --dry to see what would happen")
parser.add_argument("cfs", help="cfs binary")
parser.add_argument("--comment", help="optional comment, e.g. 'clang_openblas'")
parser.add_argument("--hostname", help="overwrite hostname detection")
parser.add_argument("--info", help="only give information about possible tests, then exit", action='store_true')
parser.add_argument("--res", help="edge resolution of bulk3d_<res>.mesh to scale problem size", type=int, default='30')
parser.add_argument("--threads", nargs='+', help="space separated list of CFS/OMP/MKL_THREADS=VECLIB_MAXIMUM_THREADS to be tested, see --permutate", default=['1','4'])
parser.add_argument("--permutate", help="perform  full permuation of threads = more calcuations", action='store_true')
parser.add_argument("--dry", help="don't execute but print what shall be run", action='store_true')
parser.add_argument("--log", help="if given, write a summary plot file with the given name")
parser.add_argument("--skip", nargs='*', help="do not run the given tests, check --info first, e.g. --skip shapemap pardiso")
parser.add_argument("--test", nargs='*', help="test only the given tests, check --info first, e.g. --test cholmod ginkgo-cuda")

args = parser.parse_args()
WIN = platform.system() == 'Windows'

if args.permutate and WIN:
  print('Error: --permutate does not work on Windows')
  sys.exit(1)  

if args.skip and args.test:
  print('Error: --skip and --test concurrently makes no sense')
  sys.exit(1)

# will print the cfs info and state that config.cmake is written
config = get_config(args.cfs) # dictionary with USE_PARDISO = True/False, ...
avail = get_cases(config)

print('\nThe following tests are available',avail)
if args.info:
  sys.exit(0)
else:
  print() # new line  

tests = []
if args.test:
  for t in args.test:
    if t not in avail:
      print("Error: invalid case in --test: '" + s + "'")
      sys.exit(1)
    tests.append(t)  
else:
  tests = [t for t in avail]

if args.skip: 
  for s in args.skip:
    if s not in avail:
      print("Error: invalid case in --skip: '" + s + "'")
      sys.exit(1)
    tests.remove(s)  

mat = open("bench_mat.xml","w")
mat.write(matstr)
mat.close()

# prepare files. Use the python tools from the same location as we have this script from
mesh = 'bulk3d_' + str(args.res) + '.mesh'
if not os.path.exists(mesh):
  mo = mt.create_3d_mesh(args.res)
  mt.write_ansys_mesh(mo, mesh)
   
host = args.hostname if args.hostname else platform.uname().node.split('.')[0]
# host can be mojo.mi.uni-erlangen.de or ca8c49d1-e336-417a-ba4e-bb348a961a55.fritz.box
if host.count("-") >= 2:
  host = 'local' 

# prepare all bech_*.xml files
prepare_tests(tests)
  
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
result = {'test' : [], 'wall' : [], 'cpu' : [], 'CFS_NUM_THREADS' : [], 'OMP_NUM_THREADS' : [], 'MKL/VECLIB_THREADS' : [], 'res' : [], 'mem in GB' : []}
if args.comment:
  result['comment'] = []

for tst in tests:      
  for t in threads:
    if (tst == 'shapemap' or tst.startswith('ginkgo')) and t[2] != t0 and not (t[0] == t[1] == t[2]): # skip playing with mkl, we don't need it
      continue 

    cmnt = '-' + args.comment if args.comment else ''
    problem = 'bench_' + host + cmnt + '-' + tst + '-res_' + str(args.res) + '-cfs_' + t[0] + '-omp_' + t[1] + '-mkl_' + t[2]
    env = 'CFS_NUM_THREADS=' + str(t[0]) + ' OMP_NUM_THREADS=' + str(t[1]) + ' MKL_NUM_THREADS=' + str(t[2]) + ' VECLIB_MAXIMUM_THREADS=' + str(t[2]) + ' '
    flag = ' -t ' + str(t[0]) if WIN else '' # on Windows no env and all threads same as we have no permutation
    cmd = '' if WIN else env
    cmd += args.cfs + flag + ' -q -m ' + mesh + ' -p bench_' + tst + '.xml ' + problem
 
    if args.dry:
      print(cmd)
    else:
      cu.execute(cmd,output=True)
      xml = cu.open_xml(problem + '.info.xml')
      wall = float(cu.xpath(xml, '//cfsInfo/summary/timer/@wall'))
      cpu = float(cu.xpath(xml, '//cfsInfo/summary/timer/@cpu'))
      if cu.has(xml, '//cfsInfo/summary/memory/@peak'):
        mem = cu.xpath(xml, '//cfsInfo/summary/memory/@peak') # Linux and Windows
      else:
        mem = cu.xpath(xml, '//cfsInfo/summary/memory/@final') # macOS
      m = float(mem)/1024 

      print('wall', wall, 'cpu', cpu, problem)
      print()
      max_name = max([len(t) for t in tests])
      print(max_name)
      result['test'].append(tst.ljust(max_name, ' ')) # full up to 'ginkgo-cuda' or shortet
      result['wall'].append(f"{wall:6.2f}")
      result['cpu'].append(f"{cpu:6.2f}")
      result['CFS_NUM_THREADS'].append(str(t[0]))
      result['OMP_NUM_THREADS'].append(str(t[1]))
      result['MKL/VECLIB_THREADS'].append(str(t[2]))
      result['res'].append(str(args.res))
      result['mem in GB'].append(str(round(m,4)))
 
# present the whole data 
cmnt = arg.comment if args.comment else '-'
hdr = '# host:' + host + ' mesh:' + str(args.res) + ' comment:' + cmnt + '\n'
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
  
    
