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
#  Revision: $Revision: 3128 $
#  Last Modified: $Date: 2015-10-16 19:55:19 +0200 (Fr, 16. Okt 2015) $
#  Last Author:   $Author: mmeiler $

from datetime import datetime
import h5py
import numpy as np
import getpass
import os
import socket

class NacsFileWriter():
  
  _vlenStr = h5py.special_dtype(vlen=bytes)
  
  def __init__(self, fileName=None, overwrite=True, 
               useCompression=False, compressionLevel=1):
    
    # flag for using compression or not 
    self._useCompression = useCompression
    self._compressionLevel = compressionLevel
    
    # dict with element type map
    self._elemAttrs = dict()
    self._elemAttrs['Num1DElems'] = range(2,4)
    self._elemAttrs['Num2DElems'] = range(4,9)
    self._elemAttrs['Num3DElems'] = range(9,18)
    self._elemAttrs['NumElems'] = range(2,18)
    self._elemAttrs['QuadraticElems'] = [5,7,8,10,12,13,15,17]
    self._elemAttrs[0] = 'UNDEF'
    self._elemAttrs[1] = 'POINT'
    self._elemAttrs[2] = 'LINE2'
    self._elemAttrs[3] = 'LINE3'
    self._elemAttrs[4] = 'TRIA3'
    self._elemAttrs[5] = 'TRIA6'
    self._elemAttrs[6] = 'QUAD4'
    self._elemAttrs[7] = 'QUAD8'
    self._elemAttrs[8] = 'QUAD9'
    self._elemAttrs[9] = 'TET4'
    self._elemAttrs[10] = 'TET10'
    self._elemAttrs[11] = 'HEXA8'
    self._elemAttrs[12] = 'HEXA20'
    self._elemAttrs[13] = 'HEXA27'
    self._elemAttrs[14] = 'PYRA5'
    self._elemAttrs[15] = 'PYRA13'
    self._elemAttrs[16] = 'WEDGE6'
    self._elemAttrs[17] = 'WEDGE15'
    
    #set fileName
    self.setFileName(fileName, overwrite)
    
  
  def __del__(self):
    self.close()
  
  
  def close(self):
    """
    This method will close the NRF file
    """
    if self._h5 is not None:    
      self._h5.flush()
      self._h5.close()
      self._h5 = None
      self._fileName = None
  
  
  def setFileName(self, fileName, overwrite=True):
    """
    This method will create a new NRF file or open and clear an existing
    NRF file
    """

    if os.path.exists(fileName) and not overwrite:
      self._fileName = None
      self._h5 = None
      raise IOError('File %s already exists at that location!' % os.path.basename(fileName))
    
    self._fileName = fileName
    self._h5 = h5py.File(self._fileName, 'w')
    self._initFileContents()
    
  
  def _initFileContents(self):
    """
    This method will create the basic groups for a NRF File
    """
    self._fileInfo = self._h5.create_group('FileInfo')
    if self._useCompression:
      self._content = self._fileInfo.create_dataset('Content',
                                                     chunks=True, 
                                                     maxshape=(5,),
                                                     data=np.array([]),
                                                     dtype='int32',
                                                     compression='gzip',
                                                     compression_opts=self._compressionLevel)
    else:
      self._content = self._fileInfo.create_dataset('Content',
                                                     chunks=True, 
                                                     maxshape=(5,),
                                                     data=np.array([]),
                                                     dtype='int32')
    
    creatorStr =  'NACS Interface for Trelis\n'
    creatorStr += 'Username: %s\n' % getpass.getuser()
    creatorStr += 'Hostname: %s' % socket.getfqdn()
    self._fileInfo.create_dataset('Creator', shape=(1,), dtype=self._vlenStr,
                                  data=creatorStr)
    self._fileInfo.create_dataset('Date', shape=(1,), dtype=self._vlenStr,
                                  data=str(datetime.now()))
    self._fileInfo.create_dataset('Version', shape=(1,), dtype=self._vlenStr,
                                  data='0.9')
    self._updateFileCapability(1)
    #self._updateFileCapability(2)
    
    self._mesh = self._h5.create_group('Mesh')
    
    self._nodes = self._mesh.create_group('Nodes')
    self._elements = self._mesh.create_group('Elements')
    self._groups = self._mesh.create_group('Groups')
    self._regions = self._mesh.create_group('Regions')
    
    self._results = self._h5.create_group('Results')
    
  def _updateFileCapability(self, cap):
    caps = self._h5['/FileInfo/Content']
    if cap not in caps:
      caps.resize(caps.shape[0]+1, 0)
      caps[-1] = cap    
    
    
  def setNodeCoordinates(self, nodeCoordinates):
    """
    This method will set and write the node coordinates
    """
    if self._useCompression:
      self._nodes.create_dataset('Coordinates', 
                                 data=nodeCoordinates,
                                 dtype='float64',
                                 compression='gzip',
                                 compression_opts=self._compressionLevel)
    else:
      self._nodes.create_dataset('Coordinates', 
                                 data=nodeCoordinates,
                                 dtype='float64')
    
    self._nodes.attrs.create('NumNodes', nodeCoordinates.shape[0], 
                             dtype='uint32')

  def setElementConnectivity(self, elemConnectivity):
    """
    This method will set and write the element connectivity
    """
    if self._useCompression:
      self._elements.create_dataset('Connectivity', 
                                    data=elemConnectivity,
                                    dtype='uint32',
                                    compression='gzip',
                                    compression_opts=self._compressionLevel)
    else:
      self._elements.create_dataset('Connectivity', 
                                    data=elemConnectivity,
                                    dtype='uint32')
                                  
  
  def setElementTypes(self, elemTypes):
    """
    This method will set and write the element types
    """
    if self._useCompression:
      self._elements.create_dataset('Types', 
                                    data=elemTypes,
                                    dtype='uint32',
                                    compression='gzip',
                                    compression_opts=self._compressionLevel)
    else:
      self._elements.create_dataset('Types', 
                                    data=elemTypes,
                                    dtype='uint32')
                                  
    if np.max(elemTypes) > 8:
      self._mesh.attrs.create('Dimension', 3 , dtype='uint32')
    else:
      self._mesh.attrs.create('Dimension', 2, dtype='uint32')
      
    for k, v in self._elemAttrs.iteritems():
      # summary entries      
      if isinstance(k, str):
        nElems = 0
        for elemType in v:
          nElems += len((elemTypes==elemType).nonzero()[0])
        if k == 'QuadraticElems':
          if nElems > 0:
            nElems=1
          self._elements.attrs.create(k, nElems, dtype='int32')
        else:
          self._elements.attrs.create(k, nElems, dtype='uint32')
      
      # single entries
      elif isinstance(k, int):
        nElems = len((elemTypes==k).nonzero()[0])
        self._elements.attrs.create('NUM_'+v, nElems, dtype='uint32')
                                  
  def addRegion(self, name, nodes, elements, dim=3):
    """
    This method will create a new region data set
    """
    newRegion = self._regions.create_group(name)
    
    if self._useCompression:
      newRegion.create_dataset('Nodes',
                               data=nodes,
                               dtype='int32',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
                               
      newRegion.create_dataset('Elements',
                               data=elements,
                               dtype='int32',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    else:
      newRegion.create_dataset('Nodes',
                               data=nodes,
                               dtype='int32')
                               
      newRegion.create_dataset('Elements',
                               data=elements,
                               dtype='int32')
    
    newRegion.attrs.create('Dimension', dim, dtype='uint32')
  
  def addGroup(self, name, nodes, elements=None, dim=2):
    """
    This method will create a new group data set
    """
    newGroup = self._groups.create_group(name)

    if self._useCompression:
      newGroup.create_dataset('Nodes',
                              data=nodes,
                              dtype='int32',
                              compression='gzip',
                              compression_opts=self._compressionLevel)
                               
      if elements is not None:
        newGroup.create_dataset('Elements',
                                data=elements,
                                dtype='int32',
                                compression='gzip',
                                compression_opts=self._compressionLevel)
    else:
      newGroup.create_dataset('Nodes',
                            data=nodes,
                            dtype='int32')
                             
      if elements is not None:
        newGroup.create_dataset('Elements',
                                data=elements,
                                dtype='int32')
    
    newGroup.attrs.create('Dimension', dim, dtype='uint32')
  
      
  def _getMultiStepStr(self, msStep):
    """
    This method returns a multisequence step string for the given msStep
    """
    return 'MultiStep_%.0f' % msStep
    
  def _getStepStr(self, step):
    """
    This method returns a step string for the given analysis step
    """
    return 'Step_%.0f' % step
    
  
  def addMeshResult(self, msStep, step, result, region, 
                    realData, imagData=None,
                    isNodeResult=True):
    """
    This method will create a new mesh result data set
    """

    meshResultGroup = self._results.require_group('Mesh')    
    if not 'ExternalFiles' in meshResultGroup.attrs.keys():
      meshResultGroup.attrs.create('ExternalFiles', 
                                   data=0, 
                                   shape=(1,), 
                                   dtype='uint32')    
    
    msStepGroup = meshResultGroup.require_group(self._getMultiStepStr(msStep))
    stepGroup = msStepGroup.require_group(self._getStepStr(step))
    resultGroup = stepGroup.require_group(result)
    regionGroup = resultGroup.require_group(region)
    if isNodeResult:
      entityGroup = regionGroup.require_group('Nodes')
    else:
      entityGroup = regionGroup.require_group('Elements')
    if self._useCompression:
      entityGroup.create_dataset('Real',
                                 data=realData,
                                 dtype='float64',
                                 compression='gzip',
                                 compression_opts=self._compressionLevel)
      if imagData is not None and len(imagData) != 0:
        entityGroup.create_dataset('Imag',
                                 data=imagData,
                                 dtype='float64',
                                 compression='gzip',
                                 compression_opts=self._compressionLevel)
    else:
      entityGroup.create_dataset('Real',
                                 data=realData,
                                 dtype='float64')
      if imagData is not None and len(imagData) != 0:
        entityGroup.create_dataset('Imag',
                                 data=imagData,
                                 dtype='float64')
  def addMeshResultDescription(self, msStep, result, dofNames,
                               entityNames, stepNumbers, stepValues, unit,
                               isNodeResult=True,
                               isHarmonic=True):
    """
    This method will add the result decription for a mesh result
    """
    meshResultGroup = self._results.require_group('Mesh')    
    msStepGroup = meshResultGroup.require_group(self._getMultiStepStr(msStep))
    if 'LastStepNum' in msStepGroup.attrs.keys():
      if stepNumbers[-1] > msStepGroup.attrs['LastStepNum']:
        msStepGroup.attrs['LastStepNum'] = stepNumbers[-1]
    else:
      msStepGroup.attrs.create('LastStepNum', 
                               data=stepNumbers[-1], 
                               dtype='uint32')
    
    if 'LastStepValue' in msStepGroup.attrs.keys():
      if stepValues[-1] > msStepGroup.attrs['LastStepValue']:
        msStepGroup.attrs['LastStepValue'] = stepValues[-1]
    else:
      msStepGroup.attrs.create('LastStepValue', 
                               data=stepValues[-1], 
                               dtype='float64')

                               
    resDescGroup = msStepGroup.require_group('ResultDescription')
    resultGroup = resDescGroup.require_group(result)
    
    if isNodeResult:
      data=np.array([1])
    else:
      data=np.array([4])      
    resultGroup.create_dataset('DefinedOn',
                               data=data,
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    nDofs = len(dofNames)
    if nDofs < 1:
      nDofs = 1
    resultGroup.create_dataset('NumDOFs',
                               data=np.array([nDofs]),
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    resultGroup.create_dataset('EntryType',
                               data=np.array([3]),
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    resultGroup.create_dataset('StepNumbers',
                               data=stepNumbers,
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    resultGroup.create_dataset('StepValues',
                               data=stepValues,
                               dtype='float64',
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    dt = h5py.special_dtype(vlen=unicode)
    resultGroup.create_dataset('EntityNames',
                               data=entityNames,
                               dtype=dt,
                               chunks=True,
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    resultGroup.create_dataset('DOFNames',
                               data=np.array(dofNames),
                               dtype=dt,
                               compression='gzip',
                               compression_opts=self._compressionLevel)
    resultGroup.create_dataset('Unit',
                               data=np.array(unit),
                               dtype=dt,
                               compression='gzip',
                               compression_opts=self._compressionLevel) 
    dtbytes = h5py.special_dtype(vlen=bytes)    
    
    if isHarmonic:
      aTypeStr = 'harmonic'
    else:
      aTypeStr = 'transient'
      
    if 'AnalysisType' in msStepGroup.attrs.keys():
      msStepGroup.attrs['AnalysisType'] = aTypeStr
    else:
      msStepGroup.attrs.create('AnalysisType', 
                               data=str(aTypeStr),
                               dtype=dtbytes)
                               
    for stepNr, stepValue in zip(stepNumbers, stepValues):
      stepStr = self._getStepStr(stepNr)
      stepGroup = msStepGroup.require_group(stepStr)
      if not 'StepValue' in stepGroup.attrs.keys():
        stepGroup.attrs.create('StepValue',
                               data=stepValue,
                               dtype='float64')
                               




if __name__=='__main__':
  pass
