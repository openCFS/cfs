# -*- coding: utf-8 -*-
#
# This file is part of the NACS interface for Trelis.
#
# Copyright (c) 2011-2015 SIMetris GmbH, Erlangen
#             www.simetris.de
# All rights reserved.
#
# Permission is hereby granted, without written agreement to use and  
# copy this software as specified within license agreement the customer
# agreed to.
# The customer is NOT allowed to modify the content of this file nor to
# redistribute this file or part of it without written agreement of
# SIMetris GmbH, Erlangen.
#
# Revision information:
#  Revision: $Revision: 3120 $
#  Last Modified: $Date: 2015-10-09 18:00:47 +0200 (Fr, 09. Okt 2015) $
#  Last Author:   $Author: mmeiler $

def grepMethods(name, inspectItem=None,delim='\n'):
  res=[]
  if inspectItem is None:
    inspectItem=cubit
  for item in dir(inspectItem):
    if name in str(item).lower():
      res.append(item)
  print(delim.join(res))

def initNacsInterface():
  import os
  import sys
  nacsDir=r'C:\Program Files\SIMetris\NACS-2.2'
  nacsBinDir=nacsDir+r'\bin'
  nacsPyDir=nacsDir+r'\lib\site-packages'
  if not nacsBinDir in os.environ['PATH']:
    os.environ['PATH'] += ';'+nacsBinDir
  if not nacsPyDir in sys.path:
    sys.path.append(nacsPyDir)
  import h5py

if __name__=='__main__':
  pass