# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 09:35:12 2015
Last edited on Mar 08 2017

@author: martin, ivan lazarov, florian toth

Use Python 2.*, because Python 3.* is currently buggier!
"""

from datetime import datetime
import h5py
import numpy as np
import os
from nrfcommon import NrfCommon, ResultInfo

class NrfWriter(NrfCommon):
  
  
  def __init__(self, fileName=None, overwrite=False):
    super(NrfWriter, self).__init__()

    #set fileName    
    self.setFileName(fileName, overwrite)
    self._elemAttrs = dict()
    self._elemAttrs['Num1DElems'] = range(2,4)
    self._elemAttrs['Num2DElems'] = range(4,9)
    self._elemAttrs['Num3DElems'] = range(9,20) #20/21 not implemented
    self._elemAttrs['NumElems'] = range(2,20)   #20/21 not implemented
    self._elemAttrs['QuadraticElems'] = [5,7,8,10,12,13,15,17,18,19]
    self._elemAttrs[0] = 'UNDEF'
    self._elemAttrs[1] = 'POINT'
    self._elemAttrs[2] = 'LINE2'
    self._elemAttrs[3] = 'LINE2'
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
    self._elemAttrs[18] = 'WEDGE18'
    self._elemAttrs[19] = 'PYRA14'
    self._elemAttrs[20] = 'POLYGON'    # still not implemented in openCFS
    self._elemAttrs[21] = 'POLYHEDRON' # still not implemented in openCFS
    
  
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
  
  
  def setFileName(self, fileName, overwrite):
    """
    This method will create a new NRF file or open and clear an existing
    NRF file
    """

    if os.path.exists(fileName) and not overwrite:
      self._fileName = None
      self._h5 = None
      raise IOError('File %s already exists at that location!' % os.path.basename(fileName))
    elif fileName is None:
      raise IOError('A valid Filename has to provided in order to write a NRF file!')
      
    self._fileName = fileName
    self._h5 = h5py.File(self._fileName, 'w')
    self._initNrfFile()
    
  
  def _initNrfFile(self):
    """
    This method will create the basic groups for a NRF File
    """
    self._fileInfo = self._h5.create_group('FileInfo')
    
    content = np.array([1])
    self._fileInfo.create_dataset('Content', shape=(1,), dtype=np.uint32,
                                 data=content, chunks=True, maxshape=(5,))
    self._fileInfo.create_dataset('Creator', shape=(1,), dtype=self._vlenStr,
                                 data=os.path.basename(__file__)+' $Revision$')
    self._fileInfo.create_dataset('Date', shape=(1,), dtype=self._vlenStr,
                                 data=str(datetime.now()))
    self._fileInfo.create_dataset('Version', shape=(1,), dtype=self._vlenStr,
                                 data='0.9')
    
    self._mesh = self._h5.create_group('Mesh')
    
    self._nodes = self._mesh.create_group('Nodes')
    self._elements = self._mesh.create_group('Elements')
    self._groups = self._mesh.create_group('Groups')
    self._regions = self._mesh.create_group('Regions')
    self._results = self._h5.create_group('Results')
  
  def registerResult(self, resultName, resultDefinedOn):
    """
    This method will register new results.
    
    Parameters
    ----------
    
    resultName:  str
        describes the name of the new result
    resultDefinedOn: int
        describes the defined by id which may be one of
        
            NrfWriter.NODES
            
            NrfWriter.EDGES
            
            NrfWriter.FACES
            
            NrfWriter.ELEMENTS
            
            NrfWriter.REGIONS
            
            NrfWriter.ELEMENTGROUP
            
            NrfWriter.NODEGROUP
            
            NrfWriter.COILS
            
            NrfWriter.UNKNOWN
    """
    if resultName in self._resultInfo.keys():
      raise ValueError('resultName \'%s\' already exists!')
    self._resultInfo[resultName] = ResultInfo(typeId=resultDefinedOn)
  
  def setNodeCoordinates(self, nodeCoordinates):
    """
    This method will set and write the node coordinates
    """
    self._nodes.create_dataset('Coordinates', 
                               data=nodeCoordinates,
                               dtype='float64',
                               compression='gzip',
                               compression_opts=1)
    
    self._nodes.attrs.create('NumNodes',
                             data=len(nodeCoordinates),
                             shape=None,
                             dtype='uint32')

  def setElementConnectivity(self, elemConnectivity):
    """
    This method will set and write the element connectivity
    """
    self._elements.create_dataset('Connectivity', 
                                  data=elemConnectivity,
                                  dtype='uint32',
                                  compression='gzip',
                                  compression_opts=1)
                                  
  
  def setElementTypes(self, elemTypes):
    """
    This method will set and write the element types
    """
    self._elements.create_dataset('Types', 
                                  data=elemTypes,
                                  dtype='uint32',
                                  compression='gzip',
                                  compression_opts=1)
                                  
    if np.max(elemTypes) > 8:
      self._mesh.attrs.create('Dimension', 3 , dtype='uint32')
    else:
      self._mesh.attrs.create('Dimension', 2, dtype='uint32')
      
    for k, v in self._elemAttrs.items():
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

    newRegion.create_dataset('Nodes',
                             data=nodes,
                             dtype='int32',
                             compression='gzip',
                             compression_opts=1)
                             
    newRegion.create_dataset('Elements',
                             data=elements,
                             dtype='int32',
                             compression='gzip',
                             compression_opts=1)
    
    newRegion.attrs.create('Dimension', dim, dtype='uint32')
  
  def addGroup(self, name, nodes, elements=None, dim=2):
    """
    This method will create a new group data set
    """
    newGroup = self._groups.create_group(name)

    newGroup.create_dataset('Nodes',
                            data=nodes,
                            dtype='int32',
                            compression='gzip',
                            compression_opts=1)
                             
    if elements is not None and elements is not []:
      newGroup.create_dataset('Elements',
                              data=elements,
                              dtype='int32',
                              compression='gzip',
                              compression_opts=1)
    
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
    
    Parameters
    ----------
    msStep: int
        multisequence step number (value has to be >= 1)
    step: int
        step number
    result: str
        name of the result
    region: str
        name of the region for which the result is written

    realData: numpy.ndarray MxN
        real part of result data 

            M: nodes
            
            N: result components (e.g. [Rx, Ry, Rz])
    imagData: numpy.ndarray MxN
        imaginary part of result data
        
            M: nodes
            
            N: result components (e.g. [Rx, Ry, Rz])

            default: None
    isNodeResult: bool
        denotes whether the result is defined by nodes or elements
        
            default: True
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
    if len(realData.shape) == 1:
      realData=realData.reshape((realData.shape[0],1))
    entityGroup.create_dataset('Real',
                               data=realData,
                               dtype='float64',
                               compression='gzip',
                               compression_opts=1)
    if imagData is not None and len(imagData) != 0:
      if len(imagData.shape) == 1:
        imagData=imagData.reshape((imagData.shape[0],1))
      entityGroup.create_dataset('Imag',
                               data=imagData,
                               dtype='float64',
                               compression='gzip',
                               compression_opts=1)
    self._updateFileCapability(2)
                               
  def addHistoryResult(self, msStep, result, entity, 
                       realData, imagData=None,
                       definedOn=None,
                       isNodeResult=True):
    """
    This method will create a new history result data set.
    The result will be added for the provided group or region. 
    
    Parameters
    ----------
    msStep: int
        multisequence step nr. value has to be >= 1
    result: str
        name of the result
    entity: str
        name of the region or group for which the result is written
        
            Note: this region or group needs to exist prior writing the
            history result data
    realData: numpy.ndarray MxN
        real part of result data 
        
            M: nr. of history result steps
        
            N: nr. of history result entities (nodes/elements)
    imagData: numpy.ndarray MxN, optional
        imaginary part of result data
        
            M: nr. of history result steps
            
            N: nr. of history result entities (nodes/elements)
            
            default: None
    isNodeResult: bool
        denotes whether the result is defined by nodes or elements
        
            default: True
    """

    histResultGroup = self._results.require_group('History')
    
    
    definedOn = self._resultInfo[result].definedOnStr
    if definedOn == 'Nodes':
      if self._hasRegion(entity):
        entities = self.getNodesForRegion(entity)
      elif self._hasGroup(entity):
        entities = self.getNodesForGroup(entity)
      else:
        print('Entities with name %s do not exist!' % entity)
        return
    elif definedOn == 'Elements':
      if self._hasRegion(entity):
        entities = self.getElementsForRegion(entity)
      elif self._hasGroup(entity):
        entities = self.getElementsForGroup(entity)
      else:
        print('Entities with name %s do not exist!' % entity)
        return
    else:
      entities = [entity]
    
    msStepGroup = histResultGroup.require_group(self._getMultiStepStr(msStep))
    resultGroup = msStepGroup.require_group(result)
    
    entityGroup = resultGroup.require_group(definedOn)
    
    nDofs = self.getDofNamesForResult(result)
    
    for itemIdx, item in enumerate(entities):
      itemGroup = entityGroup.require_group('%s' % str(item))

      # real part
      if len(realData.shape) == 1:
        realData=realData.reshape((realData.shape[0],nDofs))
      itemGroup.create_dataset('Real',
                               data=realData[:,itemIdx*nDofs:(itemIdx+1)*nDofs].reshape(realData.shape[0],nDofs),
                               dtype='float64',
                               compression='gzip',
                               compression_opts=1)
      if imagData is not None and len(imagData) != 0:
        if len(imagData.shape) == 1:
          imagData=imagData.reshape((imagData.shape[0],nDofs))
        itemGroup.create_dataset('Imag',
                                 data=imagData[:,itemIdx*nDofs:(itemIdx+1)*nDofs].reshape(imagData.shape[0],nDofs),
                                 dtype='float64',
                                 compression='gzip',
                                 compression_opts=1)
    self._updateFileCapability(3)

  def addMeshResultDescription(self, msStep, result, dofNames,
                               entityNames, stepNumbers, stepValues, unit,
                               isNodeResult=True,
                               isHarmonic=True):
    """
    This method will add the result decription for a mesh result
    """
    self._addResultDescription(msStep, result, dofNames,
                               entityNames, stepNumbers, stepValues, unit,
                               isNodeResult=isNodeResult, isHarmonic=isHarmonic,
                               isMeshResult=True)
                               
  def addHistoryResultDescription(self, msStep, result, dofNames,
                               entityNames, stepNumbers, stepValues, unit,
                               isNodeResult=True,
                               isHarmonic=True):
    """
    This method will add the history result decription for a mesh result
    """
    self._addResultDescription(msStep, result, dofNames,
                               entityNames, stepNumbers, stepValues, unit,
                               isNodeResult=isNodeResult, isHarmonic=isHarmonic,
                               isMeshResult=False)
    
  def _addResultDescription(self, msStep, result, dofNames,
                           entityNames, stepNumbers, stepValues, unit,
                           isNodeResult=True, isHarmonic=True,
                           isMeshResult=True):
    """
    This method will add a result decription either for a history result
    or for a mesh result.
    """
    if isMeshResult:
      resultTypeGroup = self._results.require_group('Mesh')
    else:
      resultTypeGroup = self._results.require_group('History')
      
    msStepGroup = resultTypeGroup.require_group(self._getMultiStepStr(msStep))
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
    
    definedOnId = np.array([self._resultInfo[result].definedOnId])
    resultGroup.create_dataset('DefinedOn',
                               data=definedOnId,
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=1)
    nDofs = len(dofNames)
    if nDofs < 1:
      nDofs = 1
    resultGroup.create_dataset('NumDOFs',
                               data=np.array([nDofs]),
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=1)
    if nDofs == 1:
      entryType = 1
    elif 2 <= nDofs <= 3:
      entryType = 3
    else:
      entryType = 6
    resultGroup.create_dataset('EntryType',
                               data=np.array([entryType]),
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=1)
    resultGroup.create_dataset('StepNumbers',
                               data=stepNumbers,
                               dtype='uint32',
                               compression='gzip',
                               compression_opts=1)
    resultGroup.create_dataset('StepValues',
                               data=stepValues,
                               dtype='float64',
                               compression='gzip',
                               compression_opts=1)
    dt = h5py.special_dtype(vlen=unicode)
    #dt = h5py.special_dtype(vlen=bytes)
    resultGroup.create_dataset('EntityNames',
                               data=np.array(entityNames),
                               dtype=dt,
                               chunks=True,
                               compression='gzip',
                               compression_opts=1)
    resultGroup.create_dataset('DOFNames',
                               data=np.array(dofNames),
                               dtype=dt,
                               compression='gzip',
                               compression_opts=1)
    resultGroup.create_dataset('Unit',
                               data=np.array(unit),
                               dtype=dt,
                               compression='gzip',
                               compression_opts=1) 
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
                               
    if isMeshResult:
      for stepNr, stepValue in zip(stepNumbers, stepValues):
        stepStr = self._getStepStr(stepNr)
        stepGroup = msStepGroup.require_group(stepStr)
        if not 'StepValue' in stepGroup.attrs.keys():
          stepGroup.attrs.create('StepValue',
                                 data=stepValue,
                                 dtype='float64')
                               
if __name__=='__main__':
    pass
