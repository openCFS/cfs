#!/usr/bin/env python
import argparse
import sys
import os.path
from cfs_utils import *
from performance import read_info

## execute cfs and return the timer with gap
# @param threads optional number of cores, given to cfs as -t
def run_and_read(binary, mesh, xml, problem, threads=None):
  cmd = binary + (' -t ' + str(threads) if threads else '') + (" -m " + mesh if mesh else '') + " -p " + xml + " " + problem
  execute(cmd, output=True)
  info = open_xml(problem + '.info.xml')
  timer = read_info(info, gap=True)
  return timer


## 'run': perform a series of cfs runs over the meshes, the core counts and the repetitions
def cmd_run(args):
  problem = args.input[:-4] if args.input.endswith('.xml') else args.input
  meshes = args.mesh if args.mesh else [None] # None -> don't give -m to cfs
  threads = args.threads if args.threads else [None] # None -> don't give -t to cfs

  infos = []
  for mesh in meshes:
    for t in threads:
      for i in range(args.repeat):
        # the job name becomes the name of the info.xml, hence it needs to be unique
        job = problem
        job += '-' + os.path.splitext(os.path.basename(mesh))[0] if len(meshes) > 1 else ''
        job += '-t' + str(t) if args.threads else ''
        job += '-' + str(i+1) if args.repeat > 1 else ''
        run_and_read(args.executable, mesh, problem + '.xml', job, threads=t)
        infos.append(job + '.info.xml')

  print('\nwrote ' + str(len(infos)) + ' info.xml file(s), analyze them with')
  print('  performance.py analyze ' + ' '.join(infos))
  return 0

# parsing structure to run performance test simulations
main = argparse.ArgumentParser(prog='measure.py', formatter_class=argparse.RawDescriptionHelpFormatter,
  description="Call cfs for the full product of the meshes, the core counts and the repetitions.\n"
              "Every job writes its own info.xml, named after the problem plus mesh, cores and run index.",
  epilog="""examples:
  measure.py bracket.xml
  measure.py bracket.xml -m coarse.mesh fine.mesh -t 1 2 4 -r 3
  measure.py bracket.xml --executable ~/code/cfs/bin/cfs

Analyse the resulting info.xml files with 'performance.py'""")
main.add_argument('input', help="the simulation input xml file")
main.add_argument('-t', '--threads', help="space separated list of core counts, given to cfs as -t", nargs='+', type=int)
main.add_argument('-m', '--mesh', help="space separated list of mesh files - default is <input>.mesh", nargs='+')
main.add_argument('--executable', help="what to call for cfs", default='cfs')
main.add_argument('-r', '--repeat', help="how often shall each simulation be repeated - default is 1", type=int, default=1)

args = main.parse_args()
sys.exit(cmd_run(args))

