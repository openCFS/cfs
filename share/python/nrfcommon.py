# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 09:35:12 2015
Last edited on Mar 08 2017

@author: martin, ivan lazarov, florian toth
"""

from __future__ import division, print_function

import os
import sys
from collections import OrderedDict
from datetime import datetime
import h5py
import numpy as np

class ResultInfoContainer():
  """
  This is a helper class for result information
  """
  def __init__(self):
    self.mesh = ResultInfo()
    self.history = ResultInfo()

class ResultInfo():
  """
  This is a helper class for result information
  """
  
  def __init__(self, typeId=None, typeStr=None):
    
    self._initParams()    
    
    if typeId is not None:
      self.definedOnId = typeId
      self.definedOnStr = self._unknownInt2Str[typeId]
    elif typeStr is not None:
      self.definedOnStr = typeStr
      self.definedOnId = self._unknownStr2Int[typeStr]
    else:
      self.definedOnId = None
      self.definedOnStr = None
      
  def _initParams(self):
    self._unknownStr2Int = {'Nodes': 1, 'Edges': 2, 'Faces': 3, 'Elements': 4,
                     'Regions': 7, 'ElementGroup': 8, 'NodeGroup': 9,
                     'Coils': 10, 'Unknowns': 11}
    self._unknownInt2Str = {i: s for (s, i) in self._unknownStr2Int.items()}
    # "5" = surface elements, but ParaView PlugIn handles only "4" at the moment        
    self._unknownInt2Str[5] = 'Elements'  
    

class NrfCommon(object):
  """
  This is a class contains common functions required by NrfReader and
  NrfWriter class.
  """

  _vlenStr = h5py.special_dtype(vlen=bytes)
  
  _unknownStr2Int = {'Nodes': 1, 'Edges': 2, 'Faces': 3, 'Elements': 4,
                     'Regions': 7, 'ElementGroup': 8, 'NodeGroup': 9,
                     'Coils': 10, 'Unknowns': 11}
  _unknownInt2Str = {i: s for (s, i) in _unknownStr2Int.items()}
  _unknownInt2Str[5] = 'Elements'
  
  _entryStr2Int = {'Unknown': 0, 'Scalar': 1, 'Vector': 3, 'Tensor': 6, 'String': 32}
  
  _entryInt2Str = {i: s for (s, i) in _entryStr2Int.items()}
  NODES = 1
  EDGES = 2
  FACES = 3
  ELEMENTS = 4
  REGIONS = 7
  ELEMENTGROUP = 8
  NODEGROUP = 9
  COILS = 10
  UNKNOWN = 11
  
  def __init__(self):
    self._initResultDict()
  
  def _initResultDict(self):
    """
    Initializes the result dictionary.
    """
    self._resultInfo = dict()
    acouResults = [('acouPressure',1),
                   ('acouAcceleration',1),
                   ('acouPotential',1),
                   ('acouVelocity',4),
                   ('acouNormalVelocity', 4),
                   ('acouPressureD1',1),
                   ('acouPressureD2',1),
                   ('acouForce', 7),
                   ('acouPotentialD1',1),
                   ('acouPotentialD2',1),
                   ('acouRhsLoad',1),
                   ('acouRhsLoadDensity',4),
                   ('acouDivLighthillTensor', 4),
                   ('acouRHSval', 1),
                   ('acouSurfRHSval', 4),
                   ('acouSpeedOfSound', 4),
                   ('acouPowerDensity', 4),
                   ('acouPower', 8),
                   ('acouIntensity', 4),
                   ('acouNormalIntensity', 4),
                   ('acouSurfIntensity', 4),
                   ('acouEnergy', 7),
                   ('acouPmlAuxVec', 1),
                   ('acouPmlAuxScalar', 1),
                   ('acouPseudoDensity', 4),
                   ('meanFluidVelocity', 1)]    
    elecResults = [('elecPotential', 1),
                   ('elecPotentialD1', 1),
                   ('elecFieldIntensity',4),
                   ('elecPolarisation', 0),
                   ('elecPseudoPolarization', 0),
                   ('elecForceVWP', 7),
                   ('elecCharge', 8),
                   ('elecChargeDensity',4),
                   ('elecFluxDensity', 4),
                   ('elecEnergy', 7),
                   ('elecEnergyDensity', 4),
                   ('elecForceDensity',4),
                   ('elecRhsLoad', 1)]
    magResults = [('magPotential', 1),
                  ('magPotentialD1', 1),
                  ('magTotalPotential', 1),
                  ('magReducedPotential', 1),
                  ('magRhsLoad', 1),
                  ('magFluxDensity', 4),
                  ('magFlux', 7),
                  ('magNormalFluxDensity', 4),
                  ('magFieldIntensity', 4),
                  ('magEddyCurrentDensity', 4),
                  ('magTotalCurrentDensity', 4),
                  ('magPotentialDiv', 0),
                  ('magForceLorentzDensity', 4),
                  ('magEddyPowerDensity', 4),
                  ('magEnergyDensity', 4),
                  ('magForceVWP', 7),
                  ('magForceLorentz', 7),
                  ('magEnergy', 7),
                  ('magEddyPower', 7),
                  ('magEddyCurrent', 7),
                  ('magCoilCurrentDensity', 4),
                  ('coiCurrent', 10),
                  ('coilInducedVoltage', 10),
                  ('coilInductance', 10),
                  ('coilLinkedFlux', 10),
                  ('magElemPermeability', 4)]
    mechResults = [('mechDisplacement', 1),
                   ('mechVelocity', 1),
                   ('mechAcceleration',1),
                   ('mechRhsLoad', 1),
                   ('mechStress', 4),
                   ('mechStrain', 4),
                   ('mechStructIntensity', 4),
                   ('mechNormalStructIntensity', 4),
                   ('vonMisesStress', 4),
                   ('vonMisesStrain', 4),
                   ('mechKinEnergyDensity', 4),
                   ('mechDeformationEnergyDensity', 4),
                   ('mechTotalEnergyDensity', 4),
                   ('mechKinEnergy', 7),
                   ('mechDeformEnergy', 7),
                   ('mechTotalEnergy', 7),
                   ('mechPower', 4),
                   ('mechWeight', 7),
                   ('mechDisplacedSurfVolume', 8),
                   ('mechNodalReactionForce', 1),
                   ('mechForce', 7),
                   ('mechPilesIrrStrain', 4),
                   ('mechPilesOverlap', 7),
                   ('mechPilesOverlapTensionElem', 4),
                   ('mechPilesOverlapBendingElem', 4)]
  
    for result, definedOn in elecResults + mechResults + magResults + acouResults:
      if definedOn != 0:
        self._resultInfo[result] = ResultInfo(typeId=definedOn)
    
  
  def _hasHistory(self):
    """
    Returns whether a history result section is available
    
    Examples
    --------
    >>> f._hasHistory()
    False
    """
    
    return 'History' in self._h5['Results'].keys()
    
  def _hasMesh(self):
    """
    Returns whether a mesh result section is available
    
    Examples
    --------
    >>> f._hasMesh()
    True
    """
    return 'Mesh' in self._h5['Results'].keys()
  
  def _hasRegion(self, region):
    """
    Returns whether the given region exists in the mesh
    
    Examples
    --------
    >>> f._hasRegion('domain')
    True
    """
    return region in self._h5['Mesh/Regions'].keys()
    
  def _hasGroup(self, group):
    """
    Returns whether the given group exists in the mesh
    
    Examples
    --------
    >>> f._hasGroup('mic1')
    False
    """
    return group in self._h5['Mesh/Groups'].keys()
    
  def _checkRegionExists(self, region):
    """
    Checks for the existance of a region and raises an exception
    if it does not exist!
    
    Returns whether the given group exists in the mesh
    
    Examples
    --------
    >>> f._checkRegionExists('domain')
    True
    """
    if not self._hasRegion(region):
      print("ERROR: region %s is not contained in the file!" % region)
      return False
    
    return True
    
  def _checkGroupExists(self, group):
    """
    Checks for the existance of a group and raises an exception
    if it does not exist!
    
    Examples
    --------
    >>> f._checkGroupExists('mic1')
    ERROR: group mic1 is not contained in the file!
    False
    """
    if not self._hasGroup(group):
      print("ERROR: group %s is not contained in the file!" % group)
      return False
    
    return True
    
  def _getRegions(self):
    """
    This method returns all region names available in the current file
    
    Examples
    --------
    >>> f._getRegions()
    ['domain']
    """
    return list(self._h5['Mesh/Regions'].keys())
    
  def _getGroups(self):
    """
    This method returns all group names available in the current file
    
    Examples
    --------
    >>> [key for key in f._getGroups()]
    []
    """
    return list(self._h5['Mesh/Groups'].keys())
    
  def _getDefinedOnStr(self, result):
    """
    This method returns the string by which a result is defined
    """
    #TODO:

  def _updateFileCapability(self, cap):
    """
    This method updates the capability list for this nrf file
    
    Parameters
    ----------
    cap   int
        denotes the capability as follows
        
            1: Mesh
            
            2: Mesh Results
            
            3: History Results
            
            4: Userdata
            
            5: Database
           
    Examples
    --------
    f._updateFileCapability(2)
    """
    caps = self._h5['/FileInfo/Content']
    if cap not in caps:
      caps.resize(caps.shape[0]+1, 0)
      caps[-1] = cap    
    
  def getBBoxForMesh(self):
    """
    Returns the bounding box of all coordinates
    
    Examples
    --------
    >>> f.getBBoxForMesh()
    array([[0., 1., 0., 1., 0., 0.]])
    """
    coords = self._h5['Mesh/Nodes/Coordinates'][:,:]
    bbox = np.zeros((1,6))
    bbox[0,:5:2] = coords.min(axis=0)
    bbox[0,1::2] = coords.max(axis=0)
    
    return bbox
  
  def getBBoxForRegion(self, region):
    """
    Returns the bounding box of the region coordinates
    as list: xmin, xmax, ymin, ymax, zmin, zmax
    
    Examples
    --------
    >>> f.getBBoxForRegion('domain')
    array([[0., 1., 0., 1., 0., 0.]])
    """    
    bbox = np.zeros((1,6))
    regionNodes = self._h5['Mesh/Regions/%s/Nodes' % region][:]-1
    regionCoords = self._h5['Mesh/Nodes/Coordinates'][regionNodes,:]
    bbox[0,:5:2] = regionCoords.min(axis=0)
    bbox[0,1::2] = regionCoords.max(axis=0)
    
    return bbox
    
  def getBBoxForGroup(self, group):
    """
    Returns the bounding box of the group coordinates
    as list: xmin, xmax, ymin, ymax, zmin, zmax
    
    Examples
    --------
    >>> f.getBBoxForGroup('mic1')
    ERROR: group mic1 is not contained in the file!
    array([[0., 0., 0., 0., 0., 0.]])
    """    
    
    bbox = np.zeros((1,6))
    if not self._checkGroupExists(group):
      return bbox
    groupNodes = self._h5['Mesh/Groups/%s/Nodes' % group][:]-1
    groupCoords = self._h5['Mesh/Nodes/Coordinates'][groupNodes,:]
    bbox[0,:5:2] = groupCoords.min(axis=0)
    bbox[0,1::2] = groupCoords.max(axis=0)
    
    return bbox

  def getElementsForGroup(self, group):
    """
    Returns a list of elements for this given group
    If the group does not exist or the group does not contain elements
    an empty list is returned
    
    Examples
    --------
    >>> f.getElementsForGroup('mic1')
    ERROR: group mic1 is not contained in the file!
    []
    """
    if not self._checkGroupExists(group):
      return []

    addr='Mesh/Groups/%s' % group
    if "Elements" in self._h5[addr].keys():
      return self._h5[addr+'/Elements'][:]
    else:
      return []
      
  def getElementsForRegion(self, region):
    """
    Returns a list of elements for this given region
    If the region does not exist or the region does not contain elements
    an empty list is returned
    
    Examples
    --------
    >>> f.getElementsForRegion('domain')
    array([1, 2, 3, 4], dtype=int32)
    """
    if not self._checkRegionExists(region):
      return []

    addr='Mesh/Regions/%s' % region
    if "Elements" in self._h5[addr].keys():
      return self._h5[addr+'/Elements'][:]
    else:
      return []
      
  def getElementsForResult(self, step, result):
    """
    Returns a list of elements that are available for the given result in the
    given multisequence step.
    If no elements are available an empty list is returned
    
    Examples
    --------
    >>> f.getElementsForResult(1, 'dummyResScalar')
    [1, 2, 3, 4]
    """
    elements = list()

    # check whether multisequence step exists    
    if not step in self.getMultisequenceSteps():
      return elements
    
    # check whether result exists in multisequence step
    if not result in self.getMeshResultsForMultisequenceStep(step):
      return elements
    
    # get result regions
    regions = self.getEntitiesForMeshResult(step, result)
    
    for region in regions:
      elements += list(self.getElementsForRegion(region))
      
    return elements
      
  def getNodesForGroup(self, group):
    """
    Returns a list of nodes for this given group
    If the group does not exist or the group does not contain nodes
    an empty list is returned
    
    Examples
    --------
    >>> f.getNodesForGroup('mic1')
    ERROR: group mic1 is not contained in the file!
    []
    """
    if not self._checkGroupExists(group):
      return []

    addr='Mesh/Groups/%s' % group
    if "Nodes" in self._h5[addr].keys():
      return self._h5[addr+'/Nodes'][:]
    else:
      return []
      
  def getNodesForRegion(self, region):
    """
    Returns a list of nodes for this given group
    If the group does not exist or the group does not contain nodes
    an empty list is returned
    
    Examples
    --------
    >>> f.getNodesForRegion('domain')
    array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=int32)
    """
    if not self._checkRegionExists(region):
      return []

    addr='Mesh/Regions/%s' % region
    if "Nodes" in self._h5[addr].keys():
      return self._h5[addr+'/Nodes'][:]
    else:
      return []
      
  def getNodesForResult(self, step, result):
    """
    Returns a list of nodes that are available for the given result in the
    given multisequence step.
    If no nodes are available an empty list is returned
    
    Examples
    --------
    >>> f.getNodesForResult(1, 'dummyResVector')
    [1, 2, 3, 4, 5, 6, 7, 8, 9]
    """
    nodes = list()

    # check whether multisequence step exists    
    if not step in self.getMultisequenceSteps():
      return nodes
    
    # check whether result exists in multisequence step
    if not result in self.getMeshResultsForMultisequenceStep(step):
      return nodes
    
    # get result regions
    regions = self.getEntitiesForMeshResult(step, result)
    
    for region in regions:
      nodes += list(self.getNodesForRegion(region))
      
    return nodes
    

  def getNearestNodeToCoord(self, x, y, z=0, group=None, region=None):
    '''
    Returns a int node number of the nearest neighbour to a given
    coordinates. 
    
    Parameters
    ----------
    x,y,z: in
        the coordinates of interest
    region: string
        the name of a region in the result file
    group: string
        the name of a group in the result file

    Returns
    -------
    out: int
        node number of the nearest neighbour

    Examples
    --------
    >>> f.getNearestNodeToCoord(0.6,0.6)
    5
    '''
    if group==None and region==None:
      coords = self._h5['Mesh/Nodes/Coordinates'][:,:]
      nodeIdx = np.array(range(coords.shape[0]))
    elif region != None:
      nodeIdx = self._h5['Mesh/Regions/%s/Nodes' % region][:]-1
      coords = self._h5['Mesh/Nodes/Coordinates'][nodeIdx,:]
    elif group != None:
      nodeIdx = self._h5['Mesh/Groups/%s/Nodes' % group][:]-1
      coords = self._h5['Mesh/Nodes/Coordinates'][nodeIdx,:]
    else:
      print("Unsure situation: Operating on all nodes!")      
      coords = self._h5['Mesh/Nodes/Coordinates'][:,:]
      nodeIdx = np.array(range(coords.shape[0]))
    minIdx = ((coords[:,0]-x)**2+(coords[:,1]-y)**2+(coords[:,2]-z)**2).argsort()[0]
    return nodeIdx[minIdx]+1
  
  def getRegionInfo(self, region):
    """
    This method returns the number of nodes and elements for a given region
    
    Parameters
    ----------
    region: string
            the name of a region in the result file
            
    Returns
    -------
    nNodes:    int
                the number of nodes for the given region
    nElements: int
                the number of elements for the given region
                
    Examples
    --------
    >>> f.getRegionInfo('domain')
    (9, 4)
    """
    if not self._checkRegionExists(region):
      return 0,0
    nNodes = self._h5['Mesh/Regions/%s/Nodes' % region][:].shape[0]
    nElems = self._h5['Mesh/Regions/%s/Elements' % region][:].shape[0]
    return nNodes, nElems
    
  def getGroupInfo(self, group):
    """
    This method returns the number of nodes and elements for a given group
    
    Parameters
    ----------
    group: string
        the name of a group in the result file
            
    Returns
    -------
    nNodes: int
        the number of nodes for the given group
    nElements: int
        the number of elements for the given group
        
    Examples
    --------
    >>> f.getGroupInfo('mic1')
    ERROR: group mic1 is not contained in the file!
    (0, 0)
    """
    if not self._checkGroupExists(group):
      return 0,0
    if 'Nodes' in self._h5['Mesh/Groups/%s' % group].keys():
      nNodes = self._h5['Mesh/Groups/%s/Nodes' % group][:].shape[0]
    if 'Elements' in self._h5['Mesh/Groups/%s' % group].keys():
      nElems = self._h5['Mesh/Groups/%s/Elements' % group][:].shape[0]
    return nNodes,nElems
  
  def getElementConnectivity(self):
    """
    This method returns all element connectivity information for all elements
    """
    return self._h5['Mesh/Elements/Connectivity'][:]
    
  def getElementTypes(self):
    """
    This method returns all element type information for all elements
    
    Examples
    --------
    >>> f.getElementConnectivity()
    array([[1, 2, 5, 4],
           [2, 3, 6, 5],
           [4, 5, 8, 7],
           [5, 6, 9, 8]], dtype=uint32)
    """
    return self._h5['Mesh/Elements/Types'][:]
  
  def getNodeCoordinates(self):
    """
    This method returns all node coordinates
    
    Examples
    --------
    >>> f.getNodeCoordinates()
    array([[0. , 0. , 0. ],
           [0.5, 0. , 0. ],
           [1. , 0. , 0. ],
           [0. , 0.5, 0. ],
           [0.5, 0.5, 0. ],
           [1. , 0.5, 0. ],
           [0. , 1. , 0. ],
           [0.5, 1. , 0. ],
           [1. , 1. , 0. ]])
    """
    return self._h5['Mesh/Nodes/Coordinates'][:]
  
  def getCoordForNode(self, node):
    """
    The method returns the coordinates of a node
    
    Examples
    --------
    >>> f.getCoordForNode(5)
    array([0.5, 0.5, 0. ])
    """
    return self._h5['Mesh/Nodes/Coordinates'][node-1,:]

  def getMultisequenceSteps(self):
    """
    The method returns the available multisequence step numbers
    
    Examples
    --------
    >>> f.getMultisequenceSteps()
    [1]
    """
    msSteps = []
    if self._hasHistory():
      for item in self._h5['Results/History']:
        msSteps.append(int(item[10:]))
    if self._hasMesh():
      for item in self._h5['Results/Mesh']:
        msIdx = int(item[10:])
        if not msIdx in msSteps:
          msSteps.append(msIdx)
    return msSteps
  
  def _getMultiStepStr(self, step):
    """
    Internal method returning the string for the given multisequence step
    
    Examples
    --------
    >>> f._getMultiStepStr(1)
    'MultiStep_1'
    """
    return 'MultiStep_%g' % step
  
  def _getStepStr(self, step):
    """
    Internal method returning the string for the given step 
    
    Examples
    --------
    >>> f._getStepStr(2)
    'Step_2'
    """
    return 'Step_%g' % step
  
  def _getExternalFilename(self, msStep, stepNum):
    """
    Chechs if the result file uses an externam file and return its name if so
    
    Examples
    --------
    f._getExternalFilename(1,2)
    """
    if not self.useExtFiles:
      raise RuntimeWarning('%s does not use external files.' % self.fileName)
    if msStep not in self.getMultisequenceSteps():
      raise RuntimeWarning('Multisequence step %d does not exist.' % msStep)
    
    stepGroup = None
    try:
      stepGroup = self._h5['/Results/Mesh/%s/%s' % (self._getMultiStepStr(msStep),
                                                   self._getStepStr(stepNum))] 
      fn = str(stepGroup.attrs['ExtHDF5FileName'], 'utf-8')
      fn = os.path.normpath(os.path.join(os.path.dirname(self.fileName), fn))
    except:
      nameTuple = os.path.splitext(self.fileName)
      fn = '%s_ms%d_step%d%s' % (nameTuple[0], msStep, stepNum, nameTuple[1])
      relFn = os.path.relpath(fn, os.path.dirname(fn))
      if stepGroup is not None:
        stepGroup.attrs.create('ExtHDF5FileName', relFn, shape=(1,), dtype=self._vlenStr)
    
    return fn
  
  def getAnalysisType(self, step=1):
    """
    The method returns the analysis type for the given sequence step 
    
    Examples
    --------
    >>> f.getAnalysisType(step=1)
    'transient'
    """
    stepStr = self._getMultiStepStr(step)

    msGroup = None    
    if self._hasMesh():
      meshRes = self._h5['Results/Mesh']
      if stepStr in meshRes:
        msGroup = meshRes[stepStr]
    elif self._hasHistory():
      histRes = self._h5['Results/History']
      if stepStr in histRes:
        msGroup = histRes[stepStr]
        
    if msGroup is None:
      return None
    
    analStr = ''
    if 'AnalysisType' in msGroup.attrs.keys():
      analStr = msGroup.attrs['AnalysisType']
    
    if analStr == 'static' or analStr == 'transient' \
        or analStr == 'harmonic' or analStr == 'eigenFrequency':
      return analStr
    else:
      return None
  
  def getHistoryResultsForMultisequenceStep(self, step):
    """
    This method returns the available history results for a given 
    multisequence step
    
    Examples
    --------
    >>> f.getHistoryResultsForMultisequenceStep(1)
    []
    """
    if not self._hasHistory():
      return []
    stepStr = self._getMultiStepStr(step)
    if not stepStr in self._h5['Results/History']:
      return []
    else:
      return list(self._h5['Results/History/%s/ResultDescription' % stepStr].keys())
      
  def getMeshResultsForMultisequenceStep(self, step):
    """
    This method returns the available mesh results for a given 
    multisequence step
    
    Examples
    --------
    >>> f.getMeshResultsForMultisequenceStep(1)
    ['dummyResScalar', 'dummyResVector']
    """
    if not self._hasMesh():
      return []
    stepStr = self._getMultiStepStr(step)
    if not stepStr in self._h5['Results/Mesh']:
      return []
    else:
      return list(self._h5['Results/Mesh/%s/ResultDescription' % stepStr].keys())
      
  def getStepsForHistoryResult(self, step, result):
    """
    This method returns step information for the given multisequence step and
    the refering history result
    
    Parameters
    ----------
    step:  int
        the number of the multisequence step
    result: str
        the result name for which to obtain step information
            
    Returns
    -------
    stepInfo: dict
        dict-keys:   step numbers
        dict-values: step values
        
    Examples
    --------
    >>> f.getStepsForHistoryResult(1, 'dummyResScalar')
    {}
    >>> f.getStepsForHistoryResult(1, 'dummyResVector')
    {}
    """
    if not result in self.getHistoryResultsForMultisequenceStep(step):
      return dict()
    addrStr = 'Results/History/%s/ResultDescription/%s' % (self._getMultiStepStr(step), result)
    stepIdxes = self._h5['%s/StepNumbers' % addrStr][:]
    stepVals = self._h5['%s/StepValues' % addrStr][:]
    return OrderedDict(zip(stepIdxes, stepVals))
    
  def getStepsForMeshResult(self, step, result):
    """
    This method returns step information for the given multisequence step and
    the refering mesh result
    
    Parameters
    ----------
    step:  int
        the number of the multisequence step
    result: str
        the result name for which to obtain step information
            
    Returns
    -------
    stepInfo: dict
        dict-keys:   step numbers
        dict-values: step values
        
    Examples
    --------
    >>> f.getStepsForMeshResult(1, 'dummyResScalar')
    OrderedDict([(1, 0.0), (2, 1.0)])
    >>> f.getStepsForMeshResult(1, 'dummyResVector')
    OrderedDict([(1, 0.0), (2, 1.0)])
    """
    if not result in self.getMeshResultsForMultisequenceStep(step):
      return dict()
    addrStr = 'Results/Mesh/%s/ResultDescription/%s' % (self._getMultiStepStr(step), result)
    stepIdxes = self._h5['%s/StepNumbers' % addrStr][:]
    stepVals = self._h5['%s/StepValues' % addrStr][:]
    return OrderedDict(zip(stepIdxes, stepVals))

  def getDefinedOnForResult(self, result):
    """
    Internal method returning the entity on which the given result is defined.
    This might be Nodes or Elements
    
    Examples
    --------
    >>> f.getDefinedOnForResult('dummyResScalar')
    'Elements'
    >>> f.getDefinedOnForResult('dummyResVector')
    'Nodes'
    """
    
    # initialize addrStr    
    addrStr=''    
    
    # iterate over multisequence steps and look for result
    for msStep in self.getMultisequenceSteps():
      if result in self.getMeshResultsForMultisequenceStep(msStep):
        addrStr = 'Results/Mesh/%s/ResultDescription/%s/DefinedOn' % (self._getMultiStepStr(msStep), 
                                                            result)
        break
      elif result in self.getHistoryResultsForMultisequenceStep(msStep):
        addrStr = 'Results/History/%s/ResultDescription/%s/DefinedOn' % (self._getMultiStepStr(msStep), 
                                                                         result)
        break
    
    if addrStr == '':
      print("ERROR: result %s not contained in result file!" % result)
      return ''

    definedOn = self._h5[addrStr][:][0]
    
    try:
      return self._unknownInt2Str[definedOn]
    except:
      print("ERROR: unknown ID %g for definedOn for result %s!" % (definedOn, result))
      return ''
      
  def getDofNamesForResult(self, result):
    """
    The method returns available information about the DOF names for a 
    given result.
    
    Examples
    --------
    >>> f.getDofNamesForResult('dummyResScalar')
    array([''], dtype=object)
    >>> f.getDofNamesForResult('dummyResVector')
    array(['x', 'y'], dtype=object)
    """
    
    # initialize addrStr    
    addrStr=''    
    
    # iterate over multisequence steps and look for result
    for msStep in self.getMultisequenceSteps():
      if result in self.getMeshResultsForMultisequenceStep(msStep):
        addrStr = 'Results/Mesh/%s/ResultDescription/%s/DOFNames' % (self._getMultiStepStr(msStep), 
                                                                     result)
        break
      elif result in self.getHistoryResultsForMultisequenceStep(msStep):
        addrStr = 'Results/History/%s/ResultDescription/%s/DOFNames' % (self._getMultiStepStr(msStep), 
                                                                        result)
        break
    
    if addrStr == '':
      print("ERROR: result %s not contained in result file!" % result)
      return ''

    return np.array([val.decode('UTF-8') for val in self._h5[addrStr][:]], dtype=object)


  def getEntitiesForMeshResult(self, step, result):
    """
    The method returns the entity numbers on which a result is defined
    for a given multisequence step
    
    Parameters
    ----------
    step:  int
        the number of the multisequence step
    result: str
        the result name for which to obtain step information
    
    Examples
    --------
    >>> f.getEntitiesForMeshResult(1, 'dummyResScalar')
    array(['domain'], dtype=object)
    >>> f.getEntitiesForMeshResult(1, 'dummyResVector')
    array(['domain'], dtype=object)
    """
    if not result in self.getMeshResultsForMultisequenceStep(step):
      return list()
    addrStr = 'Results/Mesh/%s/ResultDescription/%s' % (self._getMultiStepStr(step), result)
    return np.array([val.decode('UTF-8') for val in self._h5['%s/EntityNames' % addrStr][:]], dtype=object)

  def getEntitiesForHistoryResult(self, step, result):
    """
    The method returns the entity numbers on which a history result is
    defined for a given multisequence step
    
    Parameters
    ----------
    step:  int
        the number of the multisequence step
    result: str
        the result name for which to obtain step information
        
    Examples
    --------
    >>> f.getEntitiesForHistoryResult(1, 'dummyResScalar')
    []
    >>> f.getEntitiesForHistoryResult(1, 'dummyResVector')
    []
    """
    if not result in self.getHistoryResultsForMultisequenceStep(step):
      return list()
    addrStr = 'Results/History/%s/ResultDescription/%s' % (self._getMultiStepStr(step), result)
    return self._h5['%s/EntityNames' % addrStr][:]
    
  def getHistoryResult(self, step, result, entity, amplPhase=False):
    """
    The method returns the data for a history result defined through
    multisequence step, result name and the entity name of a group or region.
    The data is returned as numpy.ndarray which can be complex
    
    Parameters
    ----------
    step:  int
        the number of the multisequence step
    result: str
        the result name for which to obtain step information
    entity: str
        the name of a result region or group
            
    Returns
    -------
    data:   np.ndarray
        data array containing the history results without the step values
        
    Examples
    --------
    >>> f.getHistoryResult(1, 'dummyResScalar', 'domain', amplPhase=False)
    array([], dtype=float64)
    >>> f.getHistoryResult(1, 'dummyResVector', 'domain', amplPhase=False)
    array([], dtype=float64)
    """
    if not result in self.getHistoryResultsForMultisequenceStep(step):
      return np.array([])
    addrStr = 'Results/History/%s/%s' % (self._getMultiStepStr(step), result)
    thisType= list(self._h5[addrStr].keys())[0]
    addrStr += '/%s' % thisType
    if thisType == "Nodes" :
      theseItems = self.getNodesForGroup(entity)
    elif thisType == "Elements":
      theseItems = self.getElementsForGroup(entity)
    else:
      theseItems = [entity]
    nDofs = len(self.getDofNamesForResult(result))
    nItems = len(theseItems)
    nSteps = len(self.getStepsForHistoryResult(step, result))
    isComplex = 'Imag' in self._h5['%s/%s' % (addrStr,str(theseItems[0]))].keys()
    if isComplex:
      histResult=np.zeros((nSteps,nDofs*nItems),dtype=np.complex128)
    else:
      histResult=np.zeros((nSteps,nDofs*nItems))
    for itemIdx, item in enumerate(theseItems):
      realVals = self._h5['%s/%s/Real' % (addrStr,str(item))][:]
      if len(realVals.shape) == 1:
        realVals=realVals.reshape(realVals.shape[0],1)
      if isComplex:
        imagVals = self._h5['%s/%s/Imag' % (addrStr,str(item))][:]
        if len(imagVals.shape) == 1:
          imagVals=imagVals.reshape(imagVals.shape[0],1)
        histResult[:,nDofs*itemIdx:nDofs*(itemIdx+1)] = (realVals+1j*imagVals)
      else:
        histResult[:,nDofs*itemIdx:nDofs*(itemIdx+1)] = realVals
    
    return histResult


  def _getIndexMapping(self, msStep, result, regions=None, indices=None,
                       verbose=False):
    """
    This method will return an index mapping that for a given result
    with a given set of regions and/or entity indices.
    
    Parameters
    ---------
    msStep: int
        multisequence step index, 1-based
    result: str
        nacs result string, i.e. mechDisplacement
    regions: list of str
        list containing a set of regions to extract results indices from
        optional
        
            default=None
    indices: list of int
        list containint a set of indices to extract result indices from
        optional
        
            default=None
              
    Returns
    -------
    indexDict: dict
        keys: str
            region names
        values: numpy array, n x 2
            column 1: row indices in region result matrix
            column 2: row indices in global result matrix
                   
    Examples
    --------
    >>> np.testing.assert_equal(\
    f._getIndexMapping(1, 'dummyResScalar', regions=None, indices=None, verbose=False),\
    {'domain': np.array([[0, 0], [1, 1], [2, 2], [3, 3]]),\
            0: np.array([1, 2, 3, 4])})
    >>> np.testing.assert_equal(\
    f._getIndexMapping(1, 'dummyResVector', regions=None, indices=None, verbose=False),\
    {'domain': np.array([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 6], [7, 7], [8, 8]]),\
            0: np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])})
    """

    indexDict = dict()

    if not result in self.getMeshResultsForMultisequenceStep(msStep):
      return indexDict
    
    # check or create entity list
    resultRegions = self.getEntitiesForMeshResult(msStep, result)
    
    if regions is None:
      regions = resultRegions
    else:
      removeRegions = []
      for item in regions:
        if not item in resultRegions:
          print('Region with name %s is not available for result %s') % (item, result)
          removeRegions.append(item)
      for item in removeRegions:
        regions.remove(item)
        
    # determine list of result entities
    entityType = self.getDefinedOnForResult(result)
    
    # get list of indices
    if indices is None:
      indices = []
      for region in regions:
        if entityType == 'Nodes':
          indices += list(self.getNodesForRegion(region))
        elif entityType == 'Elements':
          indices += list(self.getElementsForRegion(region))
        else:
          raise ValueError('Unknown Result Entity Type found!')

    indices = np.array(indices)
    
    # remove doublets
    # here we obtain two return values. The resulting indices array with unique
    # index values, however sorted AND the idx values of the original array that
    # we will use to extract information and keep the original order of the indices
    # array
    tmp, idx = np.unique(indices, return_index=True)
    if len(idx) != len(indices):
      if verbose:
        print('%.0f doublets where removed from the indices array' % (len(indices)-len(idx)))
    # sort index array that leads to the unique indices array
    idx.sort()
    indices = indices[idx]
    if verbose:
      print('\ttotal nr. of indices:\t%.0f' % (len(indices)))

    unassignedIdx = np.arange(len(indices))
    foundIndices = []
    # now iterate over the regions and create the indexDict
    for region in regions:
      if entityType == 'Nodes':
        regionEntities = self.getNodesForRegion(region)
      elif entityType == 'Elements':
        regionEntities = self.getElementsForRegion(region)
      
      # make sure regionEntities are in ascending order      
      sortIdx = np.argsort(regionEntities)
      
      sortedIndex = np.searchsorted(regionEntities[sortIdx], indices[unassignedIdx])
      
      # find values of regionEntities that are available in indices
      indicesIndex = np.take(sortIdx, sortedIndex, mode='clip')
      
      # update unassigned index
      found = regionEntities[indicesIndex] == indices[unassignedIdx]
      idxArray = indicesIndex[found.nonzero()[0]]
      orgIdxArray = unassignedIdx[found.nonzero()[0]]
      
      mask = regionEntities[indicesIndex] != indices[unassignedIdx]
      if len(mask.nonzero()[0]) > 0:
        unassignedIdx = unassignedIdx[mask.nonzero()[0]]
      
      # now add information to the result dictionary
      if len(idxArray) > 0:
        indexDict[region] = np.array([idxArray,orgIdxArray]).transpose()
        foundIndices += orgIdxArray.tolist()
        if verbose:
          print('\tregion %s\tnr. of indices:\t%.0f' % (region,len(idxArray)))
      
        # remove values
        #sunassignedIdx = np.delete(unassignedIdx, assignedArray)

    # finally store the resulting index array
    indexDict[0] = indices
    
    # check whether all indices could be assigned
    
    nRemaining = len(indices)-len(foundIndices)
    if nRemaining > 0:
      print('\t%.0f indices could not be assigned to regions!\n%s' % \
            (nRemaining,repr(indexDict[0][np.delete(np.array(indices),np.array(foundIndices))]))) 
        
    return indexDict


  def getMeshResult(self, msStep, stepIndex, result, regions=None, indices=None, 
                    verbose = False):
    """
    The method returns the data for a mesh result defined through
    multisequence step, step number, result name and optional a name of a group 
    or region and/or for which result entites the data should be returned
    The data is returned as amplitude/phase or real/imag values
    
    Parameters
    ----------
    msStep:  int
        the number of the multisequence step
    stepIndex:  int
        the number of the analysis step
    result: str
        the result name for which to obtain step information
    regions: list, optional
        the name of a result regions or groups
    indices: list, optional 
        node/element indices for whom to obtain the result
            
    Returns
    -------
    data:   np.ndarray
        data array containing the mesh result without the step values
        
    Examples
    --------
    >>> f.getMeshResult(1, 2, 'dummyResScalar', regions=None, indices=None, verbose = False)
    array([[1.],
           [2.],
           [3.],
           [4.]])
    """
    
    steps = self.getStepsForMeshResult(msStep, result)
    
    if not stepIndex in steps.keys():
      return np.array([])
    
    # determine list of result entities
    entityType = self.getDefinedOnForResult(result)    
    
    # get index mapping for the current setup
    indexDict = self._getIndexMapping(msStep, result, regions, indices)
    
    # external data storage for subSteps?
    addrStr = 'Results/Mesh/%s/%s' % (self._getMultiStepStr(msStep),
                                      self._getStepStr(list(steps.keys())[0]))
    if 'ExtHDF5FileName' in self._h5[addrStr].attrs.keys():
      extData = True
    else:
      extData = False

    # nr. of entities
    nEntities = len(indexDict[0])
    if nEntities == 0:
      raise Exception('Could not find any matching entities for result %s'
                      % result)
    
    # get nr of DOFs for this result
    nDofs = len(self.getDofNamesForResult(result))
    
    # get analysis type
    aType = self.getAnalysisType(msStep)
    if aType == None:
      raise Exception('ERROR: unknown analysis type')
    isComplex = (aType == 'harmonic') or (aType == 'eigenfrequency')
    
    # initialize result array
    # format:
    #         item1-dof1 item1-dof2 item1-dof3 item2-dof1 item2-dof2 item2-dof3 ...
    # node 1
    # node 2
    # node 3
    
    if verbose:
      print('INFO: initializing data array for step %g with %g entities and %g DOFs for each step' %
        (nEntities, nDofs))
    resArray = np.zeros((nEntities,nDofs),
                        np.complex128 if isComplex else np.float64)

    
    # iterate over result regions as included in indexDict
    for region, idxMap in indexDict.items():
      # skip key with all index values      
      if region == 0:
        continue
    
      if extData:
        addrStr = 'Results/Mesh/%s/%s' % (self._getMultiStepStr(msStep),
                                          self._getStepStr(stepIndex))
        fn = str(self._h5[addrStr].attrs['ExtHDF5FileName'], 'utf-8')
        extFile = os.path.join(os.path.dirname(self.fileName), fn)
        thisFile = h5py.File(extFile,'r')
        addrStr = '%s/%s/%s/' % (result,
                                 region,
                                 entityType)
                                     
        if isComplex:
          resArray[idxMap[:,1],:] = (np.array(thisFile[addrStr+'Real'])[idxMap[:,0],:]) \
                                    + 1j * (np.array(thisFile[addrStr+'Imag'])[idxMap[:,0],:])
        else:
          resArray[idxMap[:,1],:] = np.array(thisFile[addrStr+'Real'])[idxMap[:,0],:]
  
        thisFile.close()
      
      else:
        addrStr = 'Results/Mesh/%s/%s/%s/%s/%s/' % (self._getMultiStepStr(msStep),
                                                        self._getStepStr(stepIndex),
                                                        result,
                                                        region,
                                                        entityType)
        if isComplex:
          resArray[idxMap[:,1],:] = (np.array(self._h5[addrStr+'Real'])[idxMap[:,0],:]) \
                                    + 1j * (np.array(self._h5[addrStr+'Imag'])[idxMap[:,0],:])
        else:
          resArray[idxMap[:,1],:] = (np.array(self._h5[addrStr+'Real'])[idxMap[:,0],:])
      
      
    return resArray
  
  def getMeshResultOverTime(self, msStep, result, regions=None, indices=None, 
                            verbose = False):
    """
    The method returns the data for a mesh result for all steps of a given 
    multisequence step, result name and optional a name of a group 
    or region and/or for which result entites the data should be returned
    The data is returned as amplitude/phase or real/imag values
    
    Parameters
    ----------
    msStep:  int
            the number of the multisequence step
    result: str
            the result name for which to obtain step information
    regions: list, optional
            the name of a result regions or groups
    indices: list, optional 
            node/element indices for whom to obtain the result
            
    Returns
    -------
    data:   np.ndarray
            data array containing the mesh result without the step values
    
    Examples
    --------
    >>> np.testing.assert_equal(\
    f.getMeshResultOverTime(1, 'dummyResScalar', regions=None, indices=None, verbose = False),\
    np.array([[[0.],[1.]], [[0.],[2.]], [[0.],[3.]], [[0.],[4.]]]))
    """
    
    steps = self.getStepsForMeshResult(msStep, result)
    nSteps = len(steps.keys())
    
    # determine list of result entities
    entityType = self.getDefinedOnForResult(result)    
    
    # get index mapping for the current setup
    if not isinstance(regions, list) and regions is not None:
      regions = [regions]
    indexDict = self._getIndexMapping(msStep, result, regions, indices)
    
    # external data storage for subSteps?
    addrStr = 'Results/Mesh/%s/%s' % (self._getMultiStepStr(msStep),
                                      self._getStepStr(list(steps.keys())[0]))
    if 'ExtHDF5FileName' in self._h5[addrStr].attrs.keys():
      extData = True
    else:
      extData = False

    # nr. of entities
    nEntities = len(indexDict[0])
    if nEntities == 0:
      raise Exception('Could not find any matching entities for result %s'
                      % result)

    # get nr of DOFs for this result
    nDofs = len(self.getDofNamesForResult(result))
    
    # get analysis type
    aType = self.getAnalysisType(msStep)
    if aType == None:
      raise Exception('ERROR: unknown analysis type')
    isComplex = (aType == 'harmonic') or (aType == 'eigenfrequency')
    
    # initialize result array
    # format:
    #         item1-dof1 item1-dof2 item1-dof3 item2-dof1 item2-dof2 item2-dof3 ...
    # node 1
    # node 2
    # node 3
    
    if verbose:
      print('INFO: initializing data array for %g steps with %g entities and %g DOFs for each step' %
        (nSteps, nEntities, nDofs))
    resArray = np.zeros((nEntities,nSteps,nDofs),
                        np.complex128 if isComplex else np.float64)

    
    # iterate over results
    # get step
    stepIdxes = list(steps.keys())
    stepIdxes.sort()
    for idx, stepIndex in enumerate(stepIdxes):
      if verbose:
        print('INFO: step %g/%g' % (idx+1,nSteps))
        
      # iterate over result regions as included in indexDict
      for region, idxMap in indexDict.items():
        
        # skip key with all index values      
        if region == 0:
          continue        
        
        if extData:
          addrStr = 'Results/Mesh/%s/%s' % (self._getMultiStepStr(msStep),
                                            self._getStepStr(stepIndex))
          fn = str(self._h5[addrStr].attrs['ExtHDF5FileName'], 'utf-8')
          extFile = os.path.join(os.path.dirname(self.fileName), fn)
          thisFile = h5py.File(extFile,'r')
          addrStr = '%s/%s/%s/' % (result,
                                   region,
                                   entityType)
                                       
          if isComplex:
            resArray[idxMap[:,1],idx,:] = (np.array(thisFile[addrStr+'Real']).reshape[idxMap[:,0],:]) \
                                        + 1j * (np.array(thisFile[addrStr+'Imag'])[idxMap[:,0],:])
          else:
            resArray[idxMap[:,1],idx,:] = np.array(thisFile[addrStr+'Real'])[idxMap[:,0],:]
    
          thisFile.close()
      
        else:
          addrStr = 'Results/Mesh/%s/%s/%s/%s/%s/' % (self._getMultiStepStr(msStep),
                                                      self._getStepStr(stepIndex),
                                                      result,
                                                      region,
                                                      entityType)
          if isComplex:
            resArray[idxMap[:,1],idx,:] = (np.array(self._h5[addrStr+'Real'])[idxMap[:,0],:]) \
                                        + 1j * (np.array(self._h5[addrStr+'Imag'])[idxMap[:,0],:])
          else:
            resArray[idxMap[:,1],idx,:] = (np.array(self._h5[addrStr+'Real'])[idxMap[:,0],:])
        
    return resArray  
    
  def printGridInformation(self):
    """
    The method prints information about the grid, the contained regions
    and groups
    """
    regionMaxLen=max(len(max(self.regions,key=len)),6)
    groupMaxLen=max(len(max(self.groups,key=len)),5)
    elemLen=max(len(str(self.nElements)),8)
    fieldDist=4
    print("+-----------------------------------------------------------+")
    print("|  Grid Information                                         |")
    print("+-----------------------------------------------------------+")
    print("|  File: %s" % os.path.basename(self.fileName))
    print("+-----------------------------------------------------------+")
    print("|  Elements: %.0f" % self.nElements)
    print("|  Nodes:    %.0f" % self.nNodes)
    print("+-----------------------------------------------------------+")
    print("| Region Information                                        |")
    print("| Region"+" "*(regionMaxLen-6+fieldDist)\
         +"Elements"+" "*(elemLen-8+fieldDist)\
         +"Nodes")
    for region in self.regions:
      nodes, elems = self.getRegionInfo(region)
      print("| %s%s%.0f%s%.0f" % (region," "*(regionMaxLen-len(region)+fieldDist+elemLen-len(str(elems))),
                              elems," "*(fieldDist),
                              nodes))
    print("+-----------------------------------------------------------+")
    print("| Group Information                                         |")
    print("| Group"+" "*(groupMaxLen-5+fieldDist)\
         +"Elements"+" "*(elemLen-8+fieldDist)\
         +"Nodes")
    for group in self.groups:
      nodes, elems = self.getGroupInfo(group)
      print("| %s%s%.0f%s%.0f" % (group," "*(groupMaxLen-len(group)+fieldDist+elemLen-len(str(elems))),
                              elems," "*(fieldDist),
                              nodes))
    print("+-----------------------------------------------------------+")

    
if __name__ == "__main__":

    # for testing in the testsuite, this file is run with the TESTSUIT_DIR
    # as working directory
    # $ python nrfcommon.py -v
    # Thus, you must cd to the correct location for interactive testing.

    # Load nrf reader
    import nrfreader as nrfR
    
    # should load "example.cfs" from TESTSUITE_DIR
    f = nrfR.NrfReader(fileName='PYTHON/doctest/nrfcommon/exampleNrfCommon.cfs',
                       writeable=False,
                       template=None,
                       externalFiles=False)

    # perform doctest
    import doctest
    result = doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    sys.exit(result.failed)
