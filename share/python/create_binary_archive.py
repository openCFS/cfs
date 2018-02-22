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
cfs = os.path.abspath(os.path.join(pwd, os.pardir))

if not os.path.exists(cfs + '/share') and not os.path.exists(cfs + '/bin'):
  print('The current directory is no build directory within a CFS++ directory.')
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
shutil.copytree('bin/LINUX_X86_64', tmp + '/cfs/bin/LINUX_X64_64')
shutil.copytree('lib64/LINUX_X86_64', tmp + '/cfs/lib64/LINUX_X64_64')

# construct zipfile name
name = os.path.basename(pwd)
zip_file = 'cfs_' + name + '_' + str(datetime.now().year) + '_' + str(datetime.now().month)  + '_' + str(datetime.now().day) + '.zip' if not args.archive else args.archive

call = ['zip', '-r']
if args.password:
  call.append('--password')
  call.appen(args.password)
call.append(pwd + '/' + zip_file)
call.append('.')  

subprocess.call(call, cwd=tmp)
print('created cfs archive ', zip_file)
if args.password:
  print("password to decrypt '" + args.password + "'")
