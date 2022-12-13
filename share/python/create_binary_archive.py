#!/usr/bin/env python
import os
import shutil
import subprocess
import argparse
from datetime import datetime
from zipfile import ZIP_FILECOUNT_LIMIT

## the purpose of this file is to create a cfs binary arcive. 
# Similar to create_binary_archive.cmake but including share/python and less noise
parser = argparse.ArgumentParser(description = 'create a distributable cfs zip archive')
parser.add_argument('--archive', help="explicitly give archive name")
parser.add_argument('--password', help="optional encryption of arcive by password",)

args = parser.parse_args()

# is tested to be a build direcotry within cfs
pwd = os.getcwd()

# the cfs main directory
cfs = ''
with open("CMakeCache.txt") as config_file:
  for line in config_file:
    if 'CFS_SOURCE_DIR' in line:
      cfs = line.split('=')[1].rstrip()
      break
if not cfs:
  cfs = input('Could not find cfs source directory. Please enter:')

#cfs = os.path.abspath(os.path.join(pwd, os.pardir))

if not os.path.exists(cfs + '/share'): 
  print('The cfs directory seems to be wrong.')
  os.sys.exit()

if not os.path.exists(pwd + '/bin'):
  print('The current directory is no build directory.')
  os.sys.exit()

# here we copy the stuff, delete afterwards
tmp = os.path.abspath(os.path.join(pwd,'tmp_create_binary_archive'))
shutil.rmtree(tmp, ignore_errors=True)
os.mkdir(tmp)
# share


print('pwd=', pwd)
print('cfs=', cfs)
print('tmp=', tmp)

shutil.copytree(cfs + '/share', tmp + '/cfs/share')
# one might consider to remove the cfs script in /share/scripts/cfs   
os.remove(tmp + '/cfs/share/scripts/cfs')
shutil.rmtree(tmp + '/cfs/share/doc')
shutil.copytree('bin/', tmp + '/cfs/bin/')
shutil.copytree('lib/', tmp + '/cfs/lib/')

# construct zipfile name
name = os.path.basename(pwd)
zip_file = 'cfs_' + name + '_' + str(datetime.now().year) + '_' + str(datetime.now().month)  + '_' + str(datetime.now().day) + '.zip' if not args.archive else args.archive

call = ['zip', '-r']
if args.password:
  call.append('--password')
  call.append(args.password)
call.append(pwd + '/' + zip_file)
call.append('.')  

subprocess.call(call, cwd=tmp)
print('created cfs archive ', zip_file)
if args.password:
  print("password to decrypt: '" + args.password + "'")
