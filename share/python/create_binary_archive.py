#!/usr/bin/env python

# the purpose of this file is to create a cfs binary archive. 
# we include the original share, not the copy from the build directory
# includes some manual cleanup

import os
import shutil
import subprocess
import argparse
from datetime import datetime
import platform
import zipfile

# modify the root name of (cpacked) cfs archive
def modify_zip_root(org_zip, new_zip, root):
  source = zipfile.ZipFile(org_zip, 'r')
  target = zipfile.ZipFile(new_zip, 'w', zipfile.ZIP_DEFLATED)
  
  # find the org_root which is to be replaced by root
  org_root = source.filelist[0].filename
  sep = '/' if org_root.count('/') > 0 else '\\'
  org_root = org_root[:org_root.index(sep)]
  assert len(org_root) > 0
  
  for file in source.filelist:
    target.writestr(file.filename.replace(org_root, root), source.read(file.filename))
  target.close()
  source.close()   

# copy redistribuates, assumes Windows and oneAPI mkl/Fortran
def copy_redist(target):
  # check one of the directories
  def dir_cand(list):
    for l in list:
      if os.path.isdir(l):
        return l
    print('error: none of the directories exits:',list)
    os.sys.exit()    
  
  # extend list with one of the candidates
  def file_cand(files, base, cands):
    for f in cands:
      if os.path.exists(os.path.join(base,f)):
        files.append(f)
        return
    print('warning: none of the files',cands,'found in',base)  
  
  if platform.system() == 'Windows':
    cpl_base = dir_cand(["C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\windows\\redist\\intel64_win\\compiler", # e.g. oneAPI 2023.0.2
                         "C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\bin"]) # e.g. oneAPI 2024.2.2
    cpl_files = ['libiomp5md.dll','libifcoremd.dll','libmmd.dll','svml_dispmd.dll']
    copy_redist_helper(cpl_base, cpl_files, target)
  # on Linux we assume gcc Fortran and therefore gnu openmp (no libiomp.so needed)

def copy_redist_helper(dir, files, target):
  if not os.path.isdir(dir):
    print("error: directory to copy redistributables from does not exit:",dir)
    os.sys.exit()
  for f in files:
    # we assume purely windows!
    file = os.path.join(dir,f)
    if os.path.exists(file):
      shutil.copy(file, target)
    else:  
      print('warning: skip not existing redistributable',file)

# extract from the CMakeCache.txt in the build directory the cfs root directory
def cfs_root(build):
  if not os.path.isdir(build):
    print('error: neither directory or .zip:',build)
    os.sys.exit()
  
  if not os.path.exists(build + '/bin'):
    print('error: the current directory is no build directory.')
    os.sys.exit()
  
  cfs = ''
  with open(os.path.join(build,"CMakeCache.txt")) as config_file:
    for line in config_file:
      if 'CFS_SOURCE_DIR' in line:
        cfs = line.split('=')[1].rstrip()
        break
  if not cfs:
    print('error: could not find a valid CMakeCache.txt in',build)
    os.sys.exit()

  if not os.path.exists(cfs + '/share'): 
    print('error: the cfs directory seems to be wrong.')
    os.sys.exit()
     
  return cfs     

parser = argparse.ArgumentParser(description = 'create a or rename a distributable cfs zip archive')
parser.add_argument('input', nargs='?', help='build directory (empty = pwd) or cpack zip file (to rename root)')
parser.add_argument('--root', help="root directory name within archive (default 'openCFS')", default='openCFS')
parser.add_argument('--archive', help="optionally give the archive name w/o extension")
parser.add_argument('--keep', help="do not delete temp directory when creating a archive", action='store_true')
args = parser.parse_args()

date = str(datetime.now().year) + '_' + str(datetime.now().month)  + '_' + str(datetime.now().day)

# platform.system() : Darwin, Linux, Windows
# plarform.processor() : arm (Mac), 'Intel64 Family 6 Model 94 Stepping 3, GenuineIntel', 
# platform.machine(): 'AMD64' (Xenon), 'arm64' (Mac)
arch = 'macOS' if platform.system() == 'Darwin' else platform.system() # Linux/Windows
if platform.system() == 'Darwin': # extend when we support arm for Linux and Windows
  arch += '_' + platform.machine() # 'AMD64' (Xenon), 'arm64' (Mac) 
zip_file = args.root + '_' + arch + '_' + date if not args.archive else args.archive

if args.input and args.input.endswith('.zip'):
  print('rewrite',args.input,'to',zip_file + '.zip','setting root to',args.root)
  modify_zip_root(args.input, zip_file + '.zip', args.root)
else:  
  binary = os.getcwd() if not args.input else args.input
  cfs = cfs_root(binary)   
  
  # here we copy the stuff, delete afterwards when do not keep
  #tmp_dir = os.path.abspath(os.path.join(os.getcwd(),'tmp_' + zip_file))
  tmp_dir = 'tmp_' + zip_file
  shutil.rmtree(tmp_dir, ignore_errors=True)
  os.mkdir(tmp_dir)

  print('create',zip_file + '.zip','cfs',cfs,'build',binary,'root',args.root)

  base = tmp_dir + '/' + args.root


  shutil.copytree(cfs + '/share', base + '/share')
  if os.path.exists(base + '/share/python/__pycache__'):
    shutil.rmtree(base + '/share/python/__pycache__')
  shutil.rmtree(base + '/share/doc')
  shutil.copytree('bin/', base + '/bin/')
  shutil.copytree('license/', base + '/license/')
  # Window has the muparser.dll in bin, for UNIX we might have for Ipopt libhsl.dylib (likely) or libhsl.so (unlikely)
  suffix = 'dylib' if platform.system() == 'Darwin' else 'so' # no harm for Windows
  if os.path.exists('lib/libiomp5.so') or os.path.exists('lib/libhsl.' + suffix):
    os.makedirs(base + '/lib/')
    if os.path.exists('lib/libiomp5.so'): # we assume Apple Silicon macs, so no Intel oneAPI
      shutil.copy('lib/libiomp5.so', base + '/lib/')
    if os.path.exists('lib/libhsl.' + suffix):
      shutil.copy('lib/libhsl.' + suffix, base + '/lib/')

  if platform.system() in ['Windows', 'Linux']:
    copy_redist(base + '/bin/')

  shutil.make_archive(zip_file, 'zip', tmp_dir)

  if not args.keep:
    shutil.rmtree(base) # recursively removes content but not directory itself
    os.rmdir('tmp_' + zip_file)
