# -*- coding: utf-8 -*-
#
# This file is part of the NACS-Python interface.
#
# Copyright (c) 2011-2014 SIMetris GmbH, Erlangen
#                         www.simetris.de
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
#    Revision: $Revision$
#    Last Modified: $Date$
#    Last Author:   $Author$
#
"""
Created on Fri Dec 11 09:35:12 2015
Last edited on Mar 08 2017

@author: martin, ivan lazarov, florian toth

Use Python 2.*, because Python 3.* is currently buggier!
"""

from __future__ import division, print_function

import os
from collections import OrderedDict
from datetime import datetime
import h5py
import numpy as np
from nrfcommon import NrfCommon

class NrfReader(NrfCommon):
  """
  This is a class that can handle NACS HDF5 result files, read, extract
  and provide information.
  """

  _vlenStr = h5py.special_dtype(vlen=bytes)
  
  _unknownStr2Int = {'Nodes': 1, 'Edges': 2, 'Faces': 3, 'Elements': 4,
                     'Regions': 7, 'ElementGroup': 8, 'NodeGroup': 9,
                     'Coils': 10, 'Unknowns': 11}
  _unknownInt2Str = {i: s for (s, i) in _unknownStr2Int.items()}
  # "5" = surface elements, but ParaView PlugIn handles only "4" at the moment 
  _unknownInt2Str[5] = 'Elements'
  
  _entryStr2Int = {'Unknown': 0, 'Scalar': 1, 'Vector': 3, 'Tensor': 6, 'String': 32}
  
  _entryInt2Str = {i: s for (s, i) in _entryStr2Int.items()}
  
  def __init__(self, fileName="", writeable = False, template=None,
               externalFiles=False):
    super(NrfReader, self).__init__()
    
    self._h5 = None
    
    if not os.path.exists(fileName):
      if not template:
        msg="Filename %s does not exist!" % fileName
        raise IOError(msg)
      else:
        self._createFile(fileName, template)
    
    self.fileName = fileName
    self.useExtFiles = externalFiles
    
    self.isWriteable = writeable
    if not self._h5:
      if writeable:
        self._h5 = h5py.File(self.fileName,'r+')
      else:
        self._h5 = h5py.File(self.fileName,'r')
    

    self._readInformation()
    
  def __del__(self):
    '''
    This method should clean up ...
    '''
    # self._h5.close() # leads to weird ImportError: sys.meta_path is None, Python is likely shutting down
  
  def _createFile(self, fileName, template):
    """
    The method creates a HDF5 file by copying the '/Mesh'-group from
    a given template file
    """
    # create new file, raise an error if file already exists
    self._h5 = h5py.File(fileName, 'w-')
    
    # copy mesh from template file
    self._h5.copy(template._h5['/Mesh'], '/Mesh')
    
    # create FileInfo
    path = '/FileInfo'
    self._h5.create_group(path)
    content = [1]
    self._h5[path].create_dataset('Content', shape=(1,), dtype=np.uint32,
                                 data=content, chunks=True, maxshape=(5,))
    self._h5[path].create_dataset('Creator', shape=(1,), dtype=self._vlenStr,
                                 data=os.path.basename(__file__)+' $Revision$')
    self._h5[path].create_dataset('Date', shape=(1,), dtype=self._vlenStr,
                                 data=str(datetime.now()))
    self._h5[path].create_dataset('Version', shape=(1,), dtype=self._vlenStr,
                                 data='0.9')
    
    # create mandatory groups
    self._h5.create_group('/Results')
    self._h5.create_group('/UserData')
  
  def _readInformation(self):
    """
    The method reads the result file and sets some variables as shown
    in the examples
    
    Parameters
    ----------
    nElements:  int
        the number of elements
    nNodes: int
        the number of nodes
    regions: list
        the name of a result regions
    groups: list
        the name of a result groups
    
    Examples
    --------
    >>> f._readInformation()
    >>> f.nElements
    4
    >>> f.nNodes
    9
    >>> [key for key in f.groups]
    []
    >>> [key for key in f.regions]
    ['domain']
    """
    
    self.hasHistory = 'History' in self._h5['Results'].keys()
    self.hasMesh = 'Mesh' in self._h5['Results'].keys()
    try:
      self.useExtFiles = (self._h5['/Results/Mesh'].attrs['ExternalFiles'] > 0)
    except:
      pass

    self.groups = self._h5['Mesh/Groups'].keys()
    self.regions = self._h5['Mesh/Regions'].keys()
    
    self.nElements = len(self._h5['Mesh/Elements/Types'])
    self.nNodes = len(self._h5['Mesh/Nodes/Coordinates'])
    
  def _updateLastStepInfo(self, msStep, stepNum, stepVal):
    """
    This method will create or update the last step attributes
    (info) for a given multisequence step in an existing HDF5 File

    Parameters
    ----------
    msStep: int
        the number of the multisequence step that has to be edited or
        created        
    stepNum: int
        the step number that has to be appended or created
    stepVal: float
        the step value that has to be appended or created

    Examples
    ----------
    f._updateLastStepInfo(1, 3, 2.0)
    """
    try:
      msGroup = self._h5['/Results/Mesh/'+self._getMultiStepStr(msStep)]
    except:
      raise RuntimeError('Multisequence step %d does not exist.' % msStep)
    
    if 'LastStepNum' in msGroup.attrs.keys():
      if stepNum > msGroup.attrs['LastStepNum']:
        msGroup.attrs['LastStepNum'] = stepNum
    else:
      msGroup.attrs.create('LastStepNum', data=stepNum, shape=(1,), dtype=np.uint32)
    
    if 'LastStepValue' in msGroup.attrs.keys():
      if stepVal > msGroup.attrs['LastStepValue']:
        msGroup.attrs['LastStepValue'] = stepVal
    else:
      msGroup.attrs.create('LastStepValue', data=stepVal, shape=(1,), dtype=np.float64)
    
  def createMultisequenceStep(self, step, analysisType, mesh=True, history=False):
    """
    This method will create a multisequence step in an existing
    HDF5 File

    Parameters
    ----------
    step: int
        the number of the multisequence step that has to be created        
    analysisType: str
        the name of the analysis type
    mesh: bool, optional
        defines if an mesh group should be created
    history: bool, optional
        defines if an history group should be created
        
    Examples
    ----------
    f.createMultisequenceStep(2, 'transient', mesh=True, history=False)
    """
    if not self.isWriteable:
      raise IOError('H5 File is not writeable! Open file using writable flag!')

    msStr = self._getMultiStepStr(step)

    resultGroup = self._h5.require_group('Results')
      
    if mesh:
      meshGroup = resultGroup.require_group('Mesh')
      self.hasMesh = True

      meshGroup.attrs.create('ExternalFiles', int(self.useExtFiles), dtype=np.uint32)
      
      msGroup = meshGroup.require_group(msStr)

      if 'AnalysisType' not in msGroup.attrs.keys():
        msGroup.attrs.create('AnalysisType', analysisType, dtype=self._vlenStr)
      elif msGroup.attrs['AnalysisType'] != analysisType:
        raise ValueError('Multisequence step %d already exists with different analysis type!' % step)
      
      msGroup.require_group('ResultDescription')
      
      self._updateFileCapability(2)
    
    if history:
      histGroup = resultGroup.require_group('History')
      self.hasHistory = True

      msGroup = histGroup.require_group(msStr)

      if 'AnalysisType' not in msGroup.attrs.keys():
        msGroup.attrs.create('AnalysisType', analysisType, dtype=self._vlenStr)
      elif msGroup.attrs['AnalysisType'] != analysisType:
        raise ValueError('Multisequence step %d already exists with different analysis type!' % step)
      
      msGroup.require_group('ResultDescription')
      
      self._updateFileCapability(3)
      
    
  def createMeshResult(self, msStep, name, dofNames, definedOn, entryType, unit):
    """
    This method will create a "ResultDescription" group for a given
    multisequence step in an existing HDF5 File

    Parameters
    ----------
    msStep: int
        the number of the multisequence step that has to be created        
    name: str
        the result name
    dofNames: np.array
        the name of the DOF
    definedOn: str
        the type of the entry on which the result is defined on
    entryType: str
        the type of the result entry
    unit: str
        the name of the units
        
    Examples
    ----------
    f.createMeshResult(2, 'mechDisplacement', np.array(['x','y']), 'Nodes', 'Vector', 'm')
    """
    if not self.isWriteable:
      raise IOError('H5 File is not writeable! Open file using writable flag!')
      
    path = '/Results/Mesh'
    msStr = self._getMultiStepStr(msStep)

    if not self.hasMesh or msStr not in self._h5[path].keys():
      raise RuntimeError('MultiSequence step %d does not exist!' % msStep)
    
    path += '/' + msStr + '/ResultDescription'
    
    if name not in self._h5[path].keys():
      self._h5[path].create_group(name)
    else:
      raise RuntimeError('Result "%s" already exists in multisequence step %d' % (name, msStep))
        
    path += '/' + name
    
    npDofs = np.array(dofNames, dtype=np.str_)
    try:
      defType = self._unknownStr2Int[definedOn]
    except:
      raise ValueError('Unknown definition type: %s' % definedOn)
    try:
      entType = self._entryStr2Int[entryType]
    except:
      raise ValueError('Unknown entry type: %s' % entryType)
    
    self._h5[path].create_dataset('DOFNames', dtype=self._vlenStr, data=dofNames)
    self._h5[path].create_dataset('NumDOFs', shape=(1,), dtype=np.uint32, data=npDofs.shape[0])
    self._h5[path].create_dataset('DefinedOn', shape=(1,), dtype=np.uint32, data=defType)
    self._h5[path].create_dataset('EntryType', shape=(1,), dtype=np.uint32, data=entType)
    self._h5[path].create_dataset('Unit', shape=(1,), dtype=self._vlenStr, data=unit)
    self._h5[path].create_dataset('EntityNames', shape=(0,), dtype=self._vlenStr, chunks=True, maxshape=(None,))
    self._h5[path].create_dataset('StepNumbers', shape=(0,), dtype=np.uint32, chunks=True, maxshape=(None,))
    self._h5[path].create_dataset('StepValues', shape=(0,), dtype=np.float64, chunks=True, maxshape=(None,))
    
  def writeMeshResult(self, msStep, stepNum, stepValue, result, entity, data,
                      indices=None, overwrite=False):
    '''
    This method will write a mesh result for a given multisequence
    step in an existing HDF5 File

    Parameters
    ----------
    msStep: int
        the number of the multisequence step
    stepNum: int
        the step number that has to be appended or created
    stepValue: float
        the step value that has to be appended or created
    result: str
        the result name
    entity: str
        the name of a result region or group
    data:  ndarray
        the data which should be written. Data has to be in real/imag
        format
    indices: list, optional 
        node/element indices for whom to obtain the result
    overwrite: bool, optional
        determines if the result should be overwritten
    
    Examples
    ----------
    TODO:
    '''
    if not self.isWriteable:
      raise IOError('File was not opened in write mode!')
    
    if msStep not in self.getMultisequenceSteps():
      raise RuntimeError('Multisequence step %d does not exist.' % msStep)
    
    if result not in self.getMeshResultsForMultisequenceStep(msStep):
      raise RuntimeError('Result "%s" does not exist!' % result)

    msPath = '/Results/Mesh/' + self._getMultiStepStr(msStep)
    
    stepGroup = self._h5[msPath].require_group(self._getStepStr(stepNum))
    if 'StepValue' not in stepGroup.attrs.keys():
      stepGroup.attrs.create('StepValue', stepValue, shape=(), dtype=np.float64)
    elif stepGroup.attrs['StepValue'] != stepValue:
      raise RuntimeError('Step %d (multisequence step %d) already exists with a different step value' % (stepNum, msStep))
    
    extFile = None
    if self.useExtFiles:
      extFile = h5py.File(self._getExternalFilename(msStep, stepNum), mode='a')
      stepGroup = extFile['/']
    
    definedOn = self.getDefinedOnForResult(result)
    
    entityNumbers = []
    if definedOn == 'Nodes':
      if entity in self.regions:
        entityNumbers = self.getNodesForRegion(entity)
      elif entity in self.groups:
        entityNumbers = self.getNodesForGroup(entity)
      else:
        raise RuntimeError('Entity "%s" is neither a region nor a group!' % entity)
    elif definedOn == 'Elements':
      if entity in self.regions:
        entityNumbers = self.getElementsForRegion(entity)
      elif entity in self.groups:
        entityNumbers = self.getElementsForGroup(entity)
      else:
        raise RuntimeError('Entity "%s" is neither a region nor a group!' % entity)
    else:
      raise RuntimeError('Writing of results defined on %s is not implemented yet.' % definedOn)
    numEntities = entityNumbers.size
    
    idxArray = None
    if not indices is None:
      idxArray = []
      for idx in indices:
        try:
          idxArray.append(entityNumbers.index(idx))
        except:
          raise RuntimeError('%s with index %d is not available for region %s' %
                             (definedOn, idx, entity))
      idxArray = np.array(idxArray)
    
    resultGroup = stepGroup.require_group(result)
    entityGroup = resultGroup.require_group(entity)
    entryGroup = entityGroup.require_group(definedOn)
    
    isComplex = np.issubdtype(data.dtype, np.complex)
    analType = self.getAnalysisType(msStep)
    if analType == 'harmonic' or analType == 'eigenFrequency':
      if not isComplex:
        raise TypeError('Complex data is expected for a %s analysis.' % analType)
    elif isComplex:
      raise TypeError('Real data is expected for a %s analysis.' % analType)
    
    numDofs = self.getDofNamesForResult(result).size
    
    if idxArray is None:
      if data.shape[0] != numEntities:
        raise ValueError('Data shape mismatch: %d entities required, %d given' %
                         (numEntities, data.shape[0]))
    else:
      if data.shape[0] != idxArray.shape[0]:
        raise ValueError('Data shape mismatch: %d entities and %d indices given' %
                         (data.shape[0], idxArray.shape[0]))
    if data.shape[1] != numDofs:
      raise ValueError('Data shape mismatch: %d DoFs required, %d given' %
                       (numDofs, data.shape[1]))      
    
    if not overwrite and 'Real' in entryGroup:
      raise RuntimeError('Result %s already exists on %s in step %d and cannot be overwritten.' % (result, entity, stepNum))

    dsReal = entryGroup.require_dataset('Real', shape=(numEntities, numDofs),
                                        dtype=np.float64, exact=True)
    if not isComplex:
      if idxArray:
        dsReal[idxArray] = data[:]
      else:
        dsReal[:] = data[:]
    else:
      if not overwrite and 'Imag' in entryGroup:
        raise RuntimeError('Result %s already exists on %s in step %d and cannot be overwritten.' % (result, entity, stepNum))
      dsImag = entryGroup.require_dataset('Imag', shape=(numEntities, numDofs),
                                          dtype=np.float64, exact=True)
      if idxArray:
        dsReal[idxArray] = data[:].real
        dsImag[idxArray] = data[:].imag
      else:
        dsReal[:] = data[:].real
        dsImag[:] = data[:].imag
    
    if extFile:
      extFile.close()
    
    rdPath = msPath + '/ResultDescription/' + result + '/'
    
    entityNames = self._h5[rdPath+'EntityNames']
    if entity not in entityNames:
      entityNames.resize(entityNames.shape[0]+1, 0)
      entityNames[-1] = entity
    
    stepNumbers = self._h5[rdPath+'StepNumbers']
    if stepNum not in stepNumbers:
      stepNumbers.resize(stepNumbers.shape[0]+1, 0)
      stepNumbers[-1] = stepNum
    
    stepValues = self._h5[rdPath+'StepValues']
    if stepValue not in stepValues:
      stepValues.resize(stepValues.shape[0]+1, 0)
      stepValues[-1] = stepValue
    
    self._updateLastStepInfo(msStep, stepNum, stepValue)
    
  def writeMeshResultOverTime(self, msStep, result, entity, indices=None,
                              data=None, verbose=False, overwrite=True,
                              stepNumbers=None, stepValues=None):
    """
    The method writes the data for a mesh result for all steps of a given 
    multisequence step, result name and optional a name of a group 
    or region and/or for which result entites the data should be returned
    The data is returned as amplitude/phase or real/imag values
    
    Parameters
    ----------
    msStep:  int
        the number of the multisequence step
    result: str
        the result name for which to obtain step information
    entity: string, optional
        the name of a result region or group
    indices: list, optional 
        node/element indices for whom to obtain the result
    data:  ndarray
        the data which should be written. Data has to be in real/imag
        format

    Examples
    ----------
    TODO:
    """

    # File must be writeable
    if not self.isWriteable:
      raise IOError('File was not opened in write mode!')

    # Multisequence step must exist
    if msStep not in self.getMultisequenceSteps():
      raise RuntimeError('Multisequence step %d does not exist.' % msStep)
    
    # result description must exist
    if result not in self.getMeshResultsForMultisequenceStep(msStep):
      raise RuntimeError('Result "%s" does not exist!')

    # cast all arrays into numpy.ndarrays if necessary
    if data is not None: data = np.asarray(data)
    if indices is not None: indices = np.asarray(indices)
    if stepNumbers is not None: stepNumbers = np.asarray(stepNumbers)
    if stepValues is not None: stepValues = np.asarray(stepValues)
    
    # Try to determine steps from given parameters.
    # We use an OrderedDict, because we assume the steps in the data array have
    # the same order as the step numbers / values. 
    if stepNumbers is not None and stepValues is not None:
      if stepNumbers.shape != stepValues.shape:
        raise ValueError('Dimensions of stepNumbers and stepValues must match.')
      # Both step numbers and values given => just create mapping
      steps = OrderedDict(zip(stepNumbers, stepValues))
    else:
      availSteps = self.getStepsForMeshResult(msStep, result)
      steps = OrderedDict()
      if stepNumbers is not None:
        # Only step numbers are given => try to find matching step values
        for stepNum in stepNumbers:
          if stepNum in availSteps:
            steps[stepNum] = availSteps[stepNum]
          else:
            raise RuntimeError('Cannot create new step %d without a step value.' % stepNum)
      elif stepValues is not None:
        # Only step values are given => try to find matching step numbers
        revSteps = {(v, k) for (k, v) in availSteps.items()}
        for stepVal in stepValues:
          if stepVal in revSteps:
            steps[revSteps[stepVal]] = stepVal
          else:
            raise RuntimeError('Cannot create new step %d without a step value.' % stepNum)
      else:
        # Nothing given => use existing steps
        steps = OrderedDict(availSteps)
    
    # Make sure that steps could be determined
    nSteps = len(steps)
    if nSteps == 0:
      raise ValueError('stepNumbers and stepValues must be given when creating new steps.')
    
    # Is the result defined on nodes or elements?
    definedOn = self.getDefinedOnForResult(result)
    # Get number of DOFs for this result
    nDofs = len(self.getDofNamesForResult(result))

    # Obtain entity numbers of region or group
    entityNumbers = []
    if definedOn == 'Nodes':
      if entity in self.regions:
        entityNumbers = self.getNodesForRegion(entity)
      elif entity in self.groups:
        entityNumbers = self.getNodesForGroup(entity)
      else:
        raise RuntimeError('Entity "%s" is neither a region nor a group!' % entity)
    elif definedOn == 'Elements':
      if entity in self.regions:
        entityNumbers = self.getElementsForRegion(entity)
      elif entity in self.groups:
        entityNumbers = self.getElementsForGroup(entity)
      else:
        raise RuntimeError('Entity "%s" is neither a region nor a group!' % entity)
    else:
      raise RuntimeError('Writing of results defined on %s is not implemented yet.' % definedOn)
    numEntities = entityNumbers.size
    
    # get indices in result group
    idxArray = None
    if not indices is None:
      idxArray = np.zeros_like(indices)
      it = np.nditer(indices, flags=['c_index'])
      while not it.finished:
        try:
          idx = np.nonzero(indices == it[0])[0]
          idxArray[it.index] = idx
          it.iternext()
        except:
          raise RuntimeError('%s with index %d is not available for region %s' %
                             (definedOn, it[0], entity))
    del entityNumbers

    # check result array
    if idxArray is None:
      if data.shape[0] != numEntities:
        raise Exception('ERROR: data shape mismatch:\n\tnItems required\t%g\n\tnItems provided\t%g' %
                        (numEntities, data.shape[0]))
    else:
      if data.shape[0] != idxArray.shape[0]:
        raise ValueError('Data shape mismatch: %d entities and %d indices given' %
                         (data.shape[0], idxArray.shape[0]))
    if data.shape[1] != nSteps:
      raise Exception('ERROR: data shape mismatch:\n\tnSteps required\t%g\n\tnSteps provided\t%g' %
        (nSteps, data.shape[1]))
    if data.shape[2] != nDofs:
      raise Exception('ERROR: data shape mismatch:\n\tnDofs required\t%g\n\tnDofs provided\t%g' %
        (nDofs, data.shape[2]))
    
    # check data type (real / complex)
    isComplex = np.issubdtype(data.dtype, np.complex)
    aType = self.getAnalysisType(msStep)
    if isComplex:
      if aType != 'harmonic' and aType != 'eigenFrequency':
        raise Exception('ERROR: %s analysis cannot store complex values' % aType)
    else:
      if aType != 'transient' and aType != 'static':
        raise Exception('ERROR: %s analysis cannot store real values' % aType)

    msPath = '/Results/Mesh/' + self._getMultiStepStr(msStep)
      
    # iterate over results
    for idx, stepIdx in enumerate(steps.keys()):
      if verbose:
        print('INFO: step %g/%g' % (idx+1,nSteps))

      # create step group if necessary
      stepGroup = self._h5[msPath].require_group(self._getStepStr(stepIdx))
      if 'StepValue' not in stepGroup.attrs.keys():
        stepGroup.attrs.create('StepValue', steps[stepIdx], shape=(), dtype=np.float64)
      elif stepGroup.attrs['StepValue'] != steps[stepIdx]:
        raise RuntimeError('Step %d (multisequence step %d) already exists with a different step value' % (stepIdx, msStep))

      # if this file uses external files, open or create the file for the current step
      extFile = None
      if self.useExtFiles:
        extFile = h5py.File(self._getExternalFilename(msStep, stepIdx), mode='a')
        stepGroup = extFile['/']

      resultGroup = stepGroup.require_group(result)
      entityGroup = resultGroup.require_group(entity)
      entryGroup = entityGroup.require_group(definedOn)
      
      if not overwrite and 'Real' in entryGroup:
        raise RuntimeError('Result %s already exists on %s in step %d and cannot be overwritten.' % (result, entity, stepIdx))
      
      dsReal = entryGroup.require_dataset('Real', shape=(numEntities, nDofs),
                                          dtype=np.float64, exact=True)

      if not isComplex:
        if idxArray is not None:
          dsReal[idxArray] = data[:,idx,:]
        else:
          dsReal[:] = data[:]
      else:
        if not overwrite and 'Imag' in entryGroup:
          raise RuntimeError('Result %s already exists on %s in step %d and cannot be overwritten.' % (result, entity, stepIdx))
        dsImag = entryGroup.require_dataset('Imag', shape=(numEntities, nDofs),
                                            dtype=np.float64, exact=True)
        if idxArray is not None:
          dsReal[idxArray,:] = data[:,idx,:].real
          dsImag[idxArray,:] = data[:,idx,:].imag
        else:
          dsReal[...] = data[:,idx,:].real
          dsImag[...] = data[:,idx,:].imag
      
      if extFile:
        extFile.close()

    rdPath = msPath + '/ResultDescription/' + result + '/'
    
    entityNames = self._h5[rdPath+'EntityNames']
    if entity not in entityNames:
      entityNames.resize(entityNames.shape[0]+1, 0)
      entityNames[-1] = entity
    
    dsStepNums = self._h5[rdPath+'StepNumbers']
    stepNumbers = np.sort(np.union1d(dsStepNums[:], stepNumbers))
    if stepNumbers.shape != dsStepNums.shape:
      dsStepNums.resize(stepNumbers.shape)
      dsStepNums[:] = stepNumbers[:]
    
    dsStepValues = self._h5[rdPath+'StepValues']
    stepValues = np.sort(np.union1d(dsStepValues, stepValues))
    if stepValues.shape != dsStepValues.shape:
      dsStepValues.resize(stepValues.shape)
      dsStepValues[:] = stepValues[:]

    self._updateLastStepInfo(msStep, stepNumbers[-1], stepValues[-1])
    
    if verbose:
      print('INFO: finished!')


  def createGroup(self, name, nodeIdxes, elemIdxes = None, dim=0,
                   overwrite = False):
    """
    This method will create a new nodal or element group in an existing
    HDF5 File. The indexes of the mesh entities have to be provided.

    Parameters
    ----------
    name: str
        the name of the new group. It is not allowed to use the name
        of an existing region or group.
    nodeIdxes: np.array
        contains the nodes that define this group
    elemIdxes: np.array, optional
        contains the element indexes that define this group
        
            default: None
    dim: int, optional
        defines the dimension of the group. In case of a nodal group
        dim=0, otherwise 0 < dim < modelDim
        
            default: 0
    overwrite: bool, optional
        denotes whether an existing group will be overwritten
        
            default: False

    Examples
    ----------
    f.createGroup('sensors', np.array([5,6]))
    """

    if not self.isWriteable:
      raise IOError('H5 File is not writeable! Open file using writable flag!')

    if name in self.groups and not overwrite:
      raise ValueError('A group with name %s already exists!' % name)

    if name in self.regions:
      raise ValueError('A region with name %s already exists!' % name)

    if len(nodeIdxes) == 0:
      raise ValueError('nodeIdxes does not contain any indexes!')

    addr = '/Mesh/Groups/%s' % name
    if name in self.groups and overwrite:
      del self._h5[addr]
    self._h5.create_group(addr)

    attrs = self._h5[addr].attrs
    attrs.create('Dimension',dim, dtype=np.int32)

    nodeAddr = '/Mesh/Groups/%s/%s' % (name,'Nodes')
    self._h5.create_dataset(nodeAddr, data=nodeIdxes, dtype=np.int32)



if __name__ == "__main__":

    # for testing in the testsuite, this file is run with the TESTSUIT_DIR
    # as working directory
    # $ python nrfreader.py -v
    # Thus, you must cd to the correct location for interactive testing.
   
    # should load "example.cfs" from TESTSUITE_DIR
    f = NrfReader(fileName='PYTHON/doctest/nrfreader/exampleNrfReader.cfs',
                       writeable=False,
                       template=None,
                       externalFiles=False)
    
    # perform doctest
    import doctest
    from sys import exit
    result = doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    exit(result.failed)
