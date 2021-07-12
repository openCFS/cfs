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
#  Revision: $Revision: 3503 $
#  Last Modified: $Date: 2016-11-25 14:32:11 +0100 (Fr, 25. Nov 2016) $
#  Last Author:   $Author: mmeiler $
#
"""
This file contains the NACS interface for Trelis.
It will modify the Trelis path and try to use h5py in order to be
able to create NACS Mesh Files.

The following conventions are given and are not further explained in the
interface code
Trelis Name  |     NACS Name
---------------------------
Blocks       |     Regions 
Side Sets    |     Groups (with element information)
Node Sets    |     Groups (only node information)
"""

import os
import sys

def initNacsEnv():
  # Is a NACS 2 installation available? 
  # If not the component can not be loaded and used
  if 'NACS2_ROOT_DIR' in os.environ:
   nacsDir=os.environ['NACS2_ROOT_DIR']
  else:
    nacsDir=r'C:\Program Files\SIMetris\NACS-2.2'
  
  # Binary directory of NACS
  nacsBinDir=os.path.normpath(os.path.join(nacsDir,'bin'))
  
  # Site-Packages of NACS containing numpy and h5py
  nacsPyDir=os.path.normpath(os.path.join(nacsDir,'lib','site-packages'))
  
  import site
  site.addsitedir(nacsPyDir)
  
  # append NACS binary directory to path variable
  if not nacsBinDir in os.environ['PATH']:
    os.environ['PATH'] += ';'+nacsBinDir
  
  # append NACS Site-Packages to PYTHONPATH
  if not nacsPyDir in sys.path:
    sys.path.append(nacsPyDir)
  
  sys.path.append('/usr/lib64/python2.7/site-packages')
  
# GUI stuff
pyqtVersion = None

# For Trelis 16 GUI does not work any more as there is no PyQt but it is also not really requiried
useGUI = False
if useGUI:
  try:
    import PyQt4
    from PyQt4.QtGui import QFileDialog, QMessageBox
    pyqtVersion = 4
  except Exception as e:
    print('PyQt4 module could not be loaded!\nError message:\n%s\n\n' % str(e))
    print('Now trying to load NACS PyQt4')
  
  if pyqtVersion is None:
    try:
      # note this is implemented for Windows only
      initNacsEnv()
      from nacs.widgets.qt import QFileDialog, QMessageBox
      pyqtVersion = 4
    except Exception as e:
      print('NACS PyQt4 module could not be loaded!\nError message:\n%s\n\n' % str(e))
      print('Now trying to load NACS PyQt5')
    
  if pyqtVersion is None:
    try:  
      import PyQt5
      from PyQt5.QtWidgets import QFileDialog, QMessageBox
    except Exception as e:
      print('PyQt5 module could not be loaded!\nError message:\n%s\n\n' % str(e))
      print('Now trying to load PyQt5')
  
  try:
    import widgets as wt
    reload(wt)
  except:
    print('Interface widgets module could not be loaded. Mesh can only be exported using python!')
else:    
  print('Did not even attempt to load GUI. Use writeNacsMeshFile(<filename>)')
  
try:
  # functionality of cubit interface (geo/mesh)
  import cubit
except Exception as e:
  print('Error loading cubit module\n%s' % str(e))
# framework parts that allow to register in the cubit framework
try:
  # framework part that allow to register in the cubit framework
  import emclaro
except Exception as e:
  print('Error loading emclaro module\n%s' % str(e))
try:
  # framework part that allow to register in the cubit framework
  import broker
except Exception as e:
  print('Error loading broker module\n%s' % str(e))

from datetime import datetime 

class NacsComponent():
  def __init__(self):
    pass

  
  def start_up(self, with_gui):
    """
    This callback is required to work together with trelis
    It handles the startup of our component
    """
    if with_gui > 0:
      self.ni = NacsInterface()
      
  def clean_up(self):
    """
    This callback is required to work together with trelis
    It handles the clean up process of our component
    """
    # Remove menus associated with this module
    emclaro.remove_menu_items("nacsinterface")
    
    # let broker knoiw that this module is cleaned up
    broker.clean_up_complete('nacsinterface')
    
  def interrupt_progress(self):
    """
    This callback is required to work together with trelis
    It handles interrupting a running progress
    """
    pass
  
class NacsInterface():
  def __init__(self):
    self._initMenus()
  
  def _initMenus(self):
    
    # write NACS mesh action    
    self.writeMeshAction = emclaro.PyAction()
    self.writeMeshAction.setText('Write NACS Meshfile ...')
    self.writeMeshAction.setMenuText('Write NACS Meshfile')
    self.writeMeshAction.setToolTip('Will allow to export a NACS Meshfile.')
    # connect the activated signal to the file write
    self.writeMeshAction.setActivateMethod('nacsinterface.writeNacsMeshFile')
    
    
    # add the menu items for this module (nacsinterface)
    emclaro.add_to_menu('&NACS', [self.writeMeshAction], 'nacsinterface')
  
def getNacsElemMap():
  """
  This method returns the element type map that maps Trelis to
  NACS element types.
  """
  elemMap=dict()
  elemMap['edge2'] = 2
  elemMap['edge3'] = 3
  elemMap['tri3'] = 4
  elemMap['tri6'] = 5
  elemMap['face4'] = 6
  elemMap['face8'] = 7
  elemMap['face9'] = 8
  elemMap['tet4'] = 9
  elemMap['tet10'] = 10
  elemMap['hex8'] = 11
  elemMap['hex20'] = 12
  elemMap['hex27'] = 13
  elemMap['pyramid5'] = 14
  elemMap['pyramid13'] = 15
  elemMap['wedge6'] = 16
  elemMap['wedge15'] = 17
  return elemMap

def getNacsElemDefinition():
  """
  This method returns the index arrays that map Trelis to NACS elements.
  Since the definition of elements differs between Trelis and NACS
  the connectivity has to be adapted here.
  In order to get the correct strings we index arrays for each
  element type
  """
  import numpy as np
  elemDef=dict()
  elemDef['edge2'] = np.array(range(2))
  elemDef['edge3'] = np.array(range(3))
  elemDef['tri3'] =  np.array(range(3))
  elemDef['tri6'] =  np.array(range(6))
  elemDef['face4'] = np.array(range(4))
  elemDef['face8'] = np.array(range(8))
  elemDef['face9'] = np.array(range(9))
  elemDef['tet4'] =  np.array(range(4))
  elemDef['tet10'] = np.array(range(10))
  elemDef['hex8'] =  np.array(range(8))
  elemDef['hex20'] = np.array(list(range(12))+list(range(16,20))+list(range(12,16)))
  elemDef['hex27'] = np.array(range(27))
  elemDef['pyramid5'] = np.array(range(5))
  elemDef['pyramid13'] = np.array(range(13))
  elemDef['wedge6'] = np.array(range(6))
  elemDef['wedge15'] = np.array(list(range(9))+list(range(12,15))+list(range(9,12)))
  return elemDef

def writeNacsMeshFile(nmfFilePath=None, scaleFactor=1.0, 
                      useCompression=True, level=1, 
                      verbose=True):
  
  initNacsEnv()
  import numpy as np
  from nacsfilewriter import NacsFileWriter
  
  # get element map
  elemMap = getNacsElemMap()
  
  # get element definitions (modified order of nodes for connectivity)
  elemDef = getNacsElemDefinition()
  
  # first get geometry information
  nVolumes = cubit.get_volume_count()
  nSurfaces = cubit.get_surface_count()
  nCurves = cubit.get_curve_count()
  nVertices = cubit.get_vertex_count()
  
  # second get element information
  nHex = cubit.get_hex_count()
  nTet = cubit.get_tet_count()
  nWedge = len(cubit.parse_cubit_list('wedge','in volume all'))
  nPyra = cubit.get_pyramid_count()
  nTri = cubit.get_tri_count()
  nQuad = cubit.get_quad_count()
  
  if nmfFilePath is None:
    saveDialog = wt.SaveDialog()
    retVal = saveDialog.exec_()
    if retVal ==QFileDialog.Rejected:
      if verbose:
        print('NACS Mesh Export was canceled by user.')
      return
    
    nmfFilePath = saveDialog.getFileName()
    compression, level = saveDialog.getCompressionSettings()
    useCompression = compression is not None
    
    scaleFactor = saveDialog.getModelScaleFactor()

  t0 = datetime.now()
  if verbose:
    print('-------- NACS Mesh Export ------------')
    print(t0)    
    
  if (abs(scaleFactor) < 1e-20):
    print('Scale factor is unreasonable low and assumed to be zero!')
    return
   
  
  # dimension
  if nHex + nTet + nPyra + nWedge > 0:
    dim = 3
  elif nTri + nQuad > 0:
    dim = 2
  elif nVolumes > 0:
    dim = 3
  elif nSurfaces > 0:
    dim = 2
  else:
    dim = 0
  
  if verbose:
    print('Creating NACS Mesh File Writer ...')
  
  out = NacsFileWriter(fileName=nmfFilePath,
                       useCompression=useCompression,
                       compressionLevel=level)
  
  # we don't need any gaps - therefore compress all Ids
  
  localt0 = datetime.now()
  if verbose:
    print('Compressing Ids ...')
  #cubit.cmd('compress all')
  cubit.cmd('compress element')
  cubit.cmd('compress node')
  if verbose:
    print('\ttook %.3f s' % ((datetime.now()-localt0).total_seconds()))
  
  # get node information and store it in a coordinate array
  localt0 = datetime.now()
  if verbose:
    print('Processing Nodes ...')
  nodeIds = cubit.get_entities('node')
  nodeCoords = np.zeros((len(nodeIds),3)) 
  for idx, node in enumerate(nodeIds):
    nodeCoords[idx,:] = cubit.get_nodal_coordinates(node)
  out.setNodeCoordinates(nodeCoords*scaleFactor)
  # free memory
  del nodeCoords
  if verbose:
    print('\ttook %.3f s' % ((datetime.now()-localt0).total_seconds()))
  
  # get region information
  blockIds = cubit.get_entities('block')
  
  # group information
  sideSetIds = cubit.get_entities('sideset')
  nodeSetIds = cubit.get_entities('nodeset')
  
  regionNames = dict()
  
  # get nr. of elements
  nElems = 0
  nRegionElems = 0
  nGroupElems = 0
  for blockId in blockIds:
    nRegionElems += cubit.get_exodus_element_count(blockId,'block')
    regionName = cubit.get_exodus_entity_name('block', blockId) 
  for sideSetId in sideSetIds:
    nGroupElems += cubit.get_exodus_element_count(sideSetId,'sideset')
    groupName = cubit.get_exodus_entity_name('sideset', sideSetId)
    nElems += nGroupElems
  for nodeSetId in nodeSetIds:
    groupName = cubit.get_exodus_entity_name('nodeset', nodeSetId)
  
  nElems = nRegionElems + nGroupElems
  
  if nElems == 0:
    msg = 'Seems as if there are no elements in this model!\n'
    msg += 'Did you set up blocks, side sets and node sets?'
    QMessageBox.error(None,'No Elements present!', msg)

  localt0 = datetime.now()
  if verbose:
    print('Processing Elements ...')
  elemConn = np.zeros((nElems,3))
  elemType = np.zeros((nElems,1))
  
  # now we start to iterate over all element types and write the connectivity
  for elemNr in range(1,nRegionElems+1):
    eTypeStr = cubit.get_element_type(elemNr)
    eTypeId = cubit.get_element_type_id(elemNr)
    conn=np.array(cubit.get_expanded_connectivity(eTypeStr, eTypeId))
    nElemNodes = len(conn)
    elemStr = eTypeStr + '%.0f' % nElemNodes
    # if elemConn has not enough columns, add columns
    if nElemNodes > elemConn.shape[1]:
      newColumns = nElemNodes - elemConn.shape[1]
      elemConn = np.hstack((elemConn,np.zeros((nElems, newColumns))))
    # special element handling for 2nd order elements
    elemConn[elemNr-1,:nElemNodes]=conn[elemDef[elemStr]]

    try:
      elemType[elemNr-1,0]=elemMap[elemStr]
    except Exception as e:
      print('%s\nelemNr:\t%.0f\nnElemConn:\t%.0f\nelemConn:\t%s' % (str(e), 
                                                                    elemNr, 
                                                                    nElemNodes, 
                                                                    str(conn)))
  if verbose:
    print('\ttook %.3f s' % ((datetime.now()-localt0).total_seconds()))
    
  # now write regions
  localt0 = datetime.now()
  if verbose:
    print('Processing Regions ...')
  
  nRegionElemsProcessed = 0
  for blockId in blockIds:
    regionName = cubit.get_exodus_entity_name('block', blockId)
    if verbose:
      print('\t%s' % regionName)
    b_nodes= cubit.vectori()
    b_unknown = cubit.vectori()
    b_edges= cubit.vectori()
    b_tris= cubit.vectori()
    b_faces= cubit.vectori()
    b_pyramids= cubit.vectori()
    b_wedges= cubit.vectori()
    b_tets= cubit.vectori()
    b_hexes= cubit.vectori()
  
    cubit.get_block_elements_and_nodes(blockId,
                                       b_nodes,
                                       b_unknown,
                                       b_edges,b_tris,b_faces,
                                       b_pyramids,b_wedges,b_tets,b_hexes)
    regionElems = []
    if dim == 3:
      for hexId in list(b_hexes):
        regionElems.append(cubit.get_global_element_id('hex', hexId))
      for tetId in list(b_tets):
        regionElems.append(cubit.get_global_element_id('tet', tetId))
      for pyramidId in list(b_pyramids):
        regionElems.append(cubit.get_global_element_id('pyramid', pyramidId))
      for wedgeId in list(b_wedges):
        regionElems.append(cubit.get_global_element_id('wedge', wedgeId))
    else:
      for triId in list(b_tris):
        regionElems.append(cubit.get_global_element_id('tri', triId))
      for faceId in list(b_faces):
        regionElems.append(cubit.get_global_element_id('face', faceId))
    
    regionNodes = np.unique(b_nodes)
    regionElems = np.unique(regionElems)
    nRegionElemsProcessed += regionElems.shape[0]
    out.addRegion(regionName, regionNodes, regionElems, dim)
    if verbose:
      print('\t\t%.0f elements' % regionElems.shape[0])

  print('%.0f - %.0f' % (nRegionElems, nRegionElemsProcessed))

  # free memory
  del b_nodes
  del b_unknown
  del b_edges
  del b_tris
  del b_faces
  del b_pyramids
  del b_wedges
  del b_tets
  del b_hexes
  
  if verbose:
    print('\ttook %.3f s' % ((datetime.now()-localt0).total_seconds()))

  # Trelis python interface has a short coming in terms of 
  # global to local element type mapping between sideset elements
  # Therefore we have to take care about it using a dictionary
  sideSetElemIdx = nRegionElems
  sideSetElementDict = dict()
  
  localt0 = datetime.now()
  if verbose:
    print('Processing Groups ...')
  # now write groups with elements and nodes
  for sideSetId in sideSetIds:
    groupName = cubit.get_exodus_entity_name('sideset', sideSetId)
    if verbose:
      print('\t%s - element group' % groupName)
    groupNodes = []
    groupElems = []
    if dim == 3:
      for surfaceId in cubit.get_sideset_surfaces(sideSetId):
        for triId in cubit.get_surface_tris(surfaceId):
          eTypeStr = 'tri'
          elemId = '%s%.0f' % (eTypeStr, triId)
          if not elemId in sideSetElementDict:
            sideSetElementDict[elemId] = sideSetElemIdx+1
            conn = cubit.get_expanded_connectivity(eTypeStr, triId)
            elemConn[sideSetElemIdx,:len(conn)] = conn
            elemType[sideSetElemIdx,0]=elemMap[eTypeStr+'%.0f' % len(conn)]
            sideSetElemIdx += 1
          groupElems.append(sideSetElementDict[elemId])
          groupNodes+=cubit.get_expanded_connectivity(eTypeStr, triId)
        for faceId in cubit.get_surface_quads(surfaceId):
          eTypeStr = 'face'
          elemId = '%s%.0f' % (eTypeStr, faceId)
          if not elemId in sideSetElementDict:
            sideSetElementDict[elemId] = sideSetElemIdx+1
            conn = cubit.get_expanded_connectivity(eTypeStr, faceId)
            elemConn[sideSetElemIdx,:len(conn)] = conn
            elemType[sideSetElemIdx,0]=elemMap[eTypeStr+'%.0f' % len(conn)]
            sideSetElemIdx += 1
          groupElems.append(sideSetElementDict[elemId])
          groupNodes+=cubit.get_expanded_connectivity(eTypeStr, faceId)
    
    # if there are any elements in the list it has to be a 2d elem 
    if len(groupElems) > 0:
      groupDim = 2
    else:
      groupDim = 1

    for curveId in cubit.get_sideset_curves(sideSetId):
      for edgeId in cubit.get_curve_edges(curveId):
        eTypeStr = 'edge'
        elemId = '%s%.0f' % (eTypeStr, edgeId)
        if not elemId in sideSetElementDict:
          sideSetElementDict[elemId] = sideSetElemIdx+1
          conn = cubit.get_expanded_connectivity(eTypeStr, edgeId)
          elemConn[sideSetElemIdx,:len(conn)] = conn
          elemType[sideSetElemIdx,0]=elemMap[eTypeStr+'%.0f' % len(conn)]
          sideSetElemIdx += 1
        groupElems.append(sideSetElementDict[elemId])
        groupNodes+=cubit.get_expanded_connectivity(eTypeStr, edgeId)
    groupNodes = np.unique(groupNodes)
    if len(groupElems) > 0:
      groupElems = np.unique(groupElems)
    else:
      groupElems = None
    out.addGroup(groupName, groupNodes, groupElems, groupDim)

  # we have finished the element connectivity data collection. 
  out.setElementConnectivity(elemConn)
  out.setElementTypes(elemType)
  
  # free memory
  del elemConn
  del elemType
  
  # now write groups with only nodes
  for nodeSetId in nodeSetIds:
    groupName = cubit.get_exodus_entity_name('nodeset', nodeSetId)
    if verbose:
      print('\t%s - node group' % groupName)
    groupNodes = cubit.get_nodeset_nodes_inclusive(nodeSetId)
    groupElems = None
    groupDim = 0
    out.addGroup(groupName, groupNodes, groupElems, groupDim)
  
  if verbose:
    print('\ttook %.3f s' % ((datetime.now()-localt0).total_seconds()))
  
  localt0 = datetime.now()
  if verbose:
    print('Closing Mesh File')
  out.close() 
  if verbose:
    print('\ttook %.3f' % ((datetime.now()-localt0).total_seconds()))
  
  if verbose:
    print('Mesh Export took %.3f s' % ((datetime.now()-t0).total_seconds()))
    print(datetime.now())
    print('-------- NACS Mesh Export Finished ----------')
