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
parser.add_argument('--filename', help="optionally give the full archive name")
parser.add_argument('--arch', help="optionally give architecture info for the the archive name")
parser.add_argument('--target', help="root directory name within archive (default 'openCFS')", default='openCFS')
parser.add_argument('--password', help="optional encryption of arcive by password")

args = parser.parse_args()

if args.filename and args.arch:
  print("warning: --filename overwrites --arch")

date = str(datetime.now().year) + '_' + str(datetime.now().month)  + '_' + str(datetime.now().day)
arch = ('_' + args.arch + '_') if args.arch else '_' 
zip_file = args.target + arch + date if not args.filename else args.filename

# will be tested to be a build direcotry within cfs
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

if not os.path.exists(cfs + '/share'): 
  print('The cfs directory seems to be wrong.')
  os.sys.exit()

if not os.path.exists(pwd + '/bin'):
  print('The current directory is no build directory.')
  os.sys.exit()

# here we copy the stuff, delete afterwards
tmp = os.path.abspath(os.path.join(pwd,'tmp_' + zip_file))
shutil.rmtree(tmp, ignore_errors=True)
os.mkdir(tmp)

print('pwd=', pwd)
print('cfs=', cfs)
print('tmp=', tmp)

base = tmp + '/' + args.target

shutil.copytree(cfs + '/share', base + '/share')
if os.path.exists(base + '/share/python/__pycache__'):
  shutil.rmtree(base + '/share/python/__pycache__')
shutil.rmtree(base + '/share/doc')
shutil.copytree('bin/', base + '/bin/')
# do not copy lib - add copying shared libs when ipopt is enabled for windows
shutil.copytree('license/', base + '/license/')

# we use python when we need not encryption
if args.password:
  call = ['zip', '-r'] # requires the command zip!
  call.append('--password')
  call.append(args.password)
  call.append(pwd + '/' + zip_file + '.zip')
  call.append('.')  
  subprocess.call(call, cwd=tmp)
  print("password to decrypt: '" + args.password + "'")
else:
  shutil.make_archive(zip_file, 'zip', base)

print('created cfs archive', zip_file + '.zip')
