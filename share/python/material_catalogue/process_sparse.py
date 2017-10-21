#!/usr/bin/env python
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction", type=int)
parser.add_argument("output", help="output job file")

args = parser.parse_args()
steps = args.stp
for i in range(steps+1):
  for j in range(steps+1):
    if i == 0 and j == 0:
      continue
    print('process_image.py --sparsemesh ' + str(i)+ '-'+str(j)+'.mesh' + ' --threshold 1e-8 --rhomin 1e-8 --densfile ' + str(i)+ '-'+str(j)+'_msfem.dens.xml ' + str(i)+ '-'+str(j)+'_msfem0_x.h5')
    #os.system('process_image.py --sparsemesh ' + str(i)+ '-'+str(j)+'.mesh' + ' --threshold 1e-8 --rhomin 1e-8 --densfile '+ str(i)+ '-'+str(j)+'_msfem.dens.xml '+ str(i)+ '-'+str(j)+'_msfem0_x.h5')
# write jobfile
jobfile = open(args.output, "w")
for i in range(steps+1):
  for j in range(steps+1):
    if i == 0 and j == 0:
      continue
    index = 0
    for k in range(8):
	    if k % 2 == 0:
	      jobfile.write('cfs.rel -m '+ str(i)+ '-'+str(j)+'.mesh' +' -x '+ str(i)+ '-'+str(j)+'_msfem.dens.xml ' +str(i)+ '-'+str(j)+'_msfem'+str(index)+'_x \n')
	    else:
	      jobfile.write('cfs.rel -m '+ str(i)+ '-'+str(j)+'.mesh' +' -x '+ str(i)+ '-'+str(j)+'_msfem.dens.xml ' +str(i)+ '-'+str(j)+'_msfem'+str(index)+'_y \n')
              index += 1
jobfile.close()

