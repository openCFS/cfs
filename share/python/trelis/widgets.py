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
#  Revision: $Revision: 3546 $
#  Last Modified: $Date: 2017-01-19 14:39:51 +0100 (Do, 19. Jan 2017) $
#  Last Author:   $Author: mmeiler $
#

# functionality of cubit interface (geo/mesh)
try:
  import cubit
except:
  print('Can not import cubit module!')
# framework parts that allow to register in the cubit framework

import os
import re
import sip
import sys
import traceback

pyQt4Loaded = False
try:
  from PyQt4.QtGui import QDialog, QTableWidget, QTableWidgetItem, QHeaderView, \
                          QFormLayout, QLineEdit, QLabel, QComboBox, QTextEdit, \
                          QPushButton, QFileDialog, QBrush, QColor, QGroupBox, \
                          QDoubleSpinBox, QHBoxLayout, QVBoxLayout, QValidator
  pyQt4Loaded = True
except:
  from nacsinterface import initNacsEnv
  initNacsEnv()
  from nacs.widgets.qt import QDialog, QTableWidget, QTableWidgetItem, QHeaderView, \
                          QFormLayout, QLineEdit, QLabel, QComboBox, QTextEdit, \
                          QPushButton, QFileDialog, QBrush, QColor, QGroupBox, \
                          QDoubleSpinBox, QHBoxLayout, QVBoxLayout, QValidator
  pyQt4Loaded = True
if not pyQt4Loaded:
  try:
    from PyQt5.QWidgets import QDialog, QTableWidget, QTableWidgetItem, QHeaderView, \
                              QFormLayout, QLineEdit, QLabel, QComboBox, QTextEdit, \
                              QPushButton, QFileDialog, QGroupBox, \
                              QDoubleSpinBox, QHBoxLayout, QVBoxLayout
    from PyQt5.QtGui import QBrush, QColor, QValidator
  except:
    print('Can not load PyQt4 or PyQt5!')
  
def getElementTypeMap():
  eTypeMap = dict()
  eTypeMap[1] = ['quad4','shell4','tri3','trishell3','tetra4','hex8','pyramid5','wedge6']
  eTypeMap[2] = ['quad8','shell8','tri6','trishell6','tetra10','hex20','pyramid13','wedge15']
  eTypeMap['tetra4'] = 'tetra10'
  eTypeMap['tetra10'] = 'tetra4'
  eTypeMap['hex8'] = 'hex20'
  eTypeMap['hex20'] = 'hex8'
  eTypeMap['pyramid5'] = 'pyramid13'
  eTypeMap['pyramid13'] = 'pyramid5'
  eTypeMap['wedge6'] = 'wedge15'
  eTypeMap['wedge15'] = 'wedge6'
  eTypeMap['quad4'] = 'quad8'
  eTypeMap['shell4'] = 'shell8'
  eTypeMap['tri3'] = 'tri6'
  eTypeMap['trishell3'] = 'trishell6'
  eTypeMap['quad8'] = 'quad4'
  eTypeMap['shell8'] = 'shell4'
  eTypeMap['tri6'] = 'tri3'
  eTypeMap['trishell6'] = 'trishell3'
  return eTypeMap
  
  
class SaveDialog(QDialog):
  """
  This class realizes a widget that will show information about
  the current project, defined NACS regions and groups, nr. of
  nodes and elements and allows a precheck before writing the
  meshfile
  """
  def __init__(self, parent=None):
    super(SaveDialog, self).__init__(parent)
    
    self._fileName = None
    
    self._initUI()
    self._initSignals()
  
  def _initUI(self):
    
    # --- list widgets for region and group information ---    
    self._regionList = QTableWidget()
    self._regionList.setColumnCount(4)
    self._regionList.setRowCount(0)
    self._regionList.setHorizontalHeaderLabels(['ID','Region Name',
                                                'Nodes','Elements'])
    self._regionList.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
    
    self._groupList  = QTableWidget()
    self._groupList.setColumnCount(5)
    self._groupList.setRowCount(0)
    self._groupList.setHorizontalHeaderLabels(['ID','Group Name','Type',
                                                'Nodes','Elements'])
    self._groupList.horizontalHeader().setResizeMode(QHeaderView.ResizeToContents)
    
    # --- create the mesh info groupbox ---
    meshInfo_grpbx = QGroupBox()
    meshInfo_lform = QFormLayout()
    
    nNodes_lbl = QLabel('Nr. of Nodes')
    self.nNodes_txt = QLineEdit()
    self.nNodes_txt.setEnabled(False)
    self.nNodes_txt.setReadOnly(True)
    
    nElems_lbl = QLabel('Nr. of Elements')
    self.nElems_txt = QLineEdit()
    self.nElems_txt.setReadOnly(True)
    
    dim_lbl = QLabel('Dimension')
    self.dim_txt = QLineEdit()
    self.dim_txt.setReadOnly(True)
    
    meshOrder_lbl = QLabel('Mesh Order')
    self.meshOrder_cbx = QComboBox()
    self.meshOrder_cbx.addItems(['1','2'])
    self.meshOrder_cbx.setCurrentIndex(0)

    # a text edit for displaying status messages 
    self.status_txt = QTextEdit()
    self.status_txt.setReadOnly(True)
    
    
    # mesh file label, line edit and push button
    meshFilePath_lbl = QLabel('Meshfile')
    self.meshFilePath_txt = QLineEdit()
    self.meshFilePath_txt.setReadOnly(True)
    self.selectFile_pbtn = QPushButton('Select File ...')
    
    scaleFactor_lbl = QLabel('Model Scale Factor')
    self.scaleFactor_spbx = ScientificDoubleSpinBox()
    #self.scaleFactor_spbx = QDoubleSpinBox()
    #self.scaleFactor_spbx.setDecimals(15)
    self.scaleFactor_spbx.setMinimum(1e-12)
    self.scaleFactor_spbx.setValue(1.0)
    
    self.compression_cbx = QComboBox()
    self.compression_cbx.addItems(['None', 'Gzip'])
    self.compression_cbx.setCurrentIndex(1)
    self.compression_cbx.setVisible(False)
    
    self.compressionLevel_cbx = QComboBox()
    self.compressionLevel_cbx.addItems([str(x) for x in range(1,10)])
    self.compressionLevel_cbx.setVisible(False)
    
    # 
    self._accept_pbtn = QPushButton('Save')
    self._accept_pbtn.setEnabled(False)
    self._reject_pbtn = QPushButton('Abort') 
    
    meshInfo_lform.addRow(nNodes_lbl, self.nNodes_txt)
    meshInfo_lform.addRow(nElems_lbl, self.nElems_txt)
    meshInfo_lform.addRow(dim_lbl, self.dim_txt)
    meshInfo_lform.addRow(meshOrder_lbl, self.meshOrder_cbx)
    
    meshInfo_grpbx.setLayout(meshInfo_lform)
    
    main_lvbox = QVBoxLayout()
    line1_lhbox = QHBoxLayout()
    line3_lhbox = QHBoxLayout()
    line4_lhbox = QHBoxLayout()
    line5_lhbox = QHBoxLayout()
    
    line1_lhbox.addWidget(self._regionList)
    line1_lhbox.addWidget(self._groupList)
    
    line3_lhbox.addWidget(meshFilePath_lbl)
    line3_lhbox.addWidget(self.meshFilePath_txt)
    line3_lhbox.addWidget(self.selectFile_pbtn)
    line3_lhbox.setStretchFactor(self.meshFilePath_txt, 1)
    
    #line4_lhbox.addStretch(1)
    line4_lhbox.addWidget(scaleFactor_lbl)
    line4_lhbox.addWidget(self.scaleFactor_spbx)
    line4_lhbox.addWidget(self.compression_cbx)
    line4_lhbox.addWidget(self.compressionLevel_cbx)
    
    line5_lhbox.addStretch(1)
    line5_lhbox.addWidget(self._accept_pbtn)
    line5_lhbox.addWidget(self._reject_pbtn)
    
    main_lvbox.addLayout(line1_lhbox)
    main_lvbox.addWidget  (meshInfo_grpbx)
    main_lvbox.addWidget(self.status_txt)
    main_lvbox.addLayout(line3_lhbox)
    main_lvbox.addLayout(line4_lhbox)
    main_lvbox.addLayout(line5_lhbox)
    
    self.setWindowTitle('NACS Mesh Export - Precheck')
    
    self.setLayout(main_lvbox)
    
  def _initSignals(self):
    self._accept_pbtn.clicked.connect(self._accept)
    self._reject_pbtn.clicked.connect(self._reject)
    self.selectFile_pbtn.clicked.connect(self._selectFile)
    self.compression_cbx.currentIndexChanged.connect(self._changeCompression)
    self.meshOrder_cbx.currentIndexChanged.connect(self._changeMeshOrder)
    
  def _selectFile(self):
    fileName = unicode(QFileDialog.getSaveFileName(None,
                                                         
                       caption='NACS Mesh Export - File Selection',
                       filter='NACS Meshfile (*.nmf)'))
    if fileName == '':
      return
    self._fileName = fileName
    self.meshFilePath_txt.setText(self._fileName)
    self._accept_pbtn.setEnabled(True)
  
  def showEvent(self, event):
    if not event.spontaneous():
      self._initData()
  
  def _changeCompression(self):
    useCompression = str(self.compression_cbx.currentText()).lower() != 'none' 
    self.compressionLevel_cbx.setEnabled(useCompression)
    
  def _changeMeshOrder(self):
    meshOrder = int(str(self.meshOrder_cbx.currentText()))
    eTypeMap = getElementTypeMap()
    for blockId in cubit.get_entities('block'):
      eType = self._getBlockElementType(blockId)
      if not eType.lower() in eTypeMap[meshOrder]:
        eTypeNew = eTypeMap[eType.lower()]
        cubit.cmd('block %s element type %s' % (str(blockId), eTypeNew.upper()))
  
  def getFileName(self):
    return self._fileName
  
  def getCompressionSettings(self):
    compression = str(self.compression_cbx.currentText()).lower()
    if compression == 'none':
      compression = None
    level = int(str(self.compressionLevel_cbx.currentText()))
    
    return compression, level
  
  def getModelScaleFactor(self):
    return float(self.scaleFactor_spbx.value())
  
  def _accept(self):
    self.accept()
  
  def _reject(self):
    self.reject()
  
  def _getBlockElementType(self, blockId):
    # the simplest way to determine the grid order
    blockElementType = cubit.get_block_element_type(blockId)

    eTypeMap = getElementTypeMap()
        
    if blockElementType.lower() not in eTypeMap[1] and \
       blockElementType.lower() not in eTypeMap[2]:
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
      dim, msg  = self._getDimension()
      
      if dim == 3:
        if len(list(b_hexes)) > 0:
          #item = cubit.get_global_element_id('hex', list(b_hexes)[0])
          item = b_hexes[0]
          nNodes = len(cubit.get_expanded_connectivity('hex', item))
          blockElementType = 'hex%d' % nNodes
        elif len(list(b_tets)) > 0:
          #item = cubit.get_global_element_id('tet', list(b_tets)[0])
          item = b_tets[0]
          nNodes = len(cubit.get_expanded_connectivity('tet', item))
          blockElementType = 'tetra%d' % nNodes
        elif len(list(b_pyramids)) > 0:
          #item = cubit.get_global_element_id('pyramid', list(b_pyramids)[0])
          item = b_pyramids[0]
          nNodes = len(cubit.get_expanded_connectivity('pyramid', item))
          blockElementType = 'pyramid%d' % nNodes
        elif len(list(b_wedges)) > 0:
          #item = cubit.get_global_element_id('wedge', list(b_wedges)[0])
          item = b_wedges[0]
          nNodes = len(cubit.get_expanded_connectivity('wedge', item))
          blockElementType = 'wedge%d' % nNodes
      else:
        if len(list(b_tris)) > 0:
          #item = cubit.get_global_element_id('tri', list(b_tris)[0])
          item = b_tris[0]
          nNodes = len(cubit.get_expanded_connectivity('tri', item))
          blockElementType = 'tri%d' % nNodes
        elif len(list(b_faces)) > 0:
          #item = cubit.get_global_element_id('face', list(b_faces)[0])
          item = b_faces[0]
          nNodes = len(cubit.get_expanded_connectivity('face', item))
          blockElementType = 'quad%d' % nNodes
    
    return blockElementType
  
  def _getBlockGridOder(self, blockId):
    """
    This method returns the block grid order based on the blockElementType
    of a block given by its blockId
    """
    
    # get block element type
    blockElementType = self._getBlockElementType(blockId)
    
    # initialize grid order of block
    blockGridOrder = None
    
    eTypeMap = getElementTypeMap()
    if blockElementType.lower() in eTypeMap[1]:
      return 1
    elif blockElementType.lower() in eTypeMap[2]:
      return 2
    else:
      return blockGridOrder
  
  def _getDimension(self):
    """
    This method returns the model dimension. Therefore the 
    number of volumes or surfaces or the number of 3d or 2d elements
    will be taken into account
    """
    
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
    
    msg = ''
    
    # dimension
    if (nVolumes + nSurfaces + nCurves + nVertices) == 0:
      msg = 'Dimension could not be determined.\n' \
            '\tThere seems to be no geometry present.'
      dim = 0
    elif nHex + nTet + nPyra + nWedge > 0:
      dim = 3
    elif nTri + nQuad > 0:
      dim = 2
    elif nVolumes > 0:
      dim = 3
    elif nSurfaces > 0:
      dim = 2
    elif nCurves > 0:
      dim = 1
      msg = 'NACS requires at least 2d geometry items!\n' \
            '\tPlease set up your model prior NACS Mesh Export!'
    else:
      dim = 0
      msg = 'NACS requires at least 2d geometry items!\n' \
            '\tPlease set up your model prior NACS Mesh Export!'
    
    return dim, msg
  
  def _initData(self):
    
    try:
      # clear tables
      self._regionList.clearContents()
      self._groupList.clearContents()
      self.status_txt.clear()
      
      errorBrush = QBrush(QColor('red'))
      hasErrors = False
      
      # dimension
      dim, msg = self._getDimension()
      
      if msg != '':
        self.status_txt.append(msg)
      
      if dim < 2:
        hasErrors = True
      
      self.dim_txt.setText(str(dim))
     
      self.nNodes_txt.setText(str(cubit.get_node_count()))
      self.nElems_txt.setText(str(cubit.get_element_count()))
      
      # get region information
      blockIds = cubit.get_entities('block')
      
      # group information
      sideSetIds = cubit.get_entities('sideset')
      nodeSetIds = cubit.get_entities('nodeset')
      
      if len(blockIds) == 0:
        msg = 'No Blocks were created therefore no regions would be available in the mesh file.\n'\
              'Please specify Blocks, Side Sets and Node Sets prior NACS Mesh Export!'
        self.status_txt.append(msg)
        hasErrors = True
        
      if len(sideSetIds)+len(nodeSetIds) == 0:
        msg = 'No Side Sets and Node Sets were created therefore no groups would be available in the mesh file.\n'\
              'Typically you forgot to define those! Please check your setup'
        self.status_txt.append(msg)
      # get nr. of elements
      nElems = 0
      names = []
      self._regionList.setRowCount(len(blockIds))
      gridOrder = None
      updateGridOrder = False
      for idx, blockId in enumerate(blockIds):
        # collect data
        nRegionElems = cubit.get_exodus_element_count(blockId,'block')
        regionName = cubit.get_exodus_entity_name('block', blockId)
        
        # create QTableWidgetItems
        blockIdItem = QTableWidgetItem(str(blockId))
        regionNameItem = QTableWidgetItem(regionName)
        elementsItem = QTableWidgetItem(str(nRegionElems))
        
        if not regionName in names:
          names.append(regionName)
        else:
          msg = 'Block with ID %.0f and name %s was already used!\n' % (blockId, regionName) \
              +  '\tNACS region and group names allow no duplicates!\n' \
              +  '\tPlease modify the Block Setup'
          self.status_txt.append(msg)
          regionNameItem.setForeground(errorBrush)
          hasErrors = True
        nElems += nRegionElems
        
        # send data to QTableWidget
        self._regionList.setItem(idx,0,blockIdItem)
        self._regionList.setItem(idx,1,regionNameItem)
        self._regionList.setItem(idx,3,elementsItem)
        
        blockGridOrder = self._getBlockGridOder(blockId)
        

        if blockGridOrder is None:
          msg = 'Element type of block %s can not be determined!\n' % (regionName) \
              + '\tSwitching grid order will not be available!\n'
          self.status_txt.append(msg)
          hasErrors = True
          
        self.meshOrder_cbx.setEnabled(blockGridOrder is None)
        if blockGridOrder is not None:
          if gridOrder is None:
            gridOrder = blockGridOrder
          elif gridOrder != blockGridOrder:
            gridOrder = min(gridOrder, blockGridOrder)
            updateGridOrder = True
          self.meshOrder_cbx.setEnabled(True)
          self.blockSignals(True)
          self.meshOrder_cbx.setCurrentIndex(gridOrder-1)
          self.blockSignals(False)
          if updateGridOrder:
            self._changeMeshOrder()
        else:
          self.meshOrder_cbx.setEnabled(False)
              
      
      self._groupList.setRowCount(len(sideSetIds)+len(nodeSetIds))  
      for idx, sideSetId in enumerate(sideSetIds):
        # collect data
        nGroupElems = cubit.get_exodus_element_count(sideSetId,'sideset')
        groupName = cubit.get_exodus_entity_name('sideset', sideSetId)
        
        # create QTableWidgetItems
        sideSetIdItem = QTableWidgetItem(str(sideSetId))
        groupNameItem = QTableWidgetItem(groupName)
        typeItem = QTableWidgetItem('ELEMENT')
        elementsItem = QTableWidgetItem(str(nGroupElems))
        
        if not groupName in names:
          names.append(groupName)
        else:
          msg = 'Side Set with ID %.0f and name %s was already used!\n' % (sideSetId, groupName) \
              + '\tNACS region and group names allow no duplicates!\n' \
              + '\tPlease modify the Side Set Setup'
          self.status_txt.append(msg)
          groupNameItem.setForeground(errorBrush)
          hasErrors = True
        nElems += nGroupElems
        
        # send data to QTableWidget
        self._groupList.setItem(idx,0,sideSetIdItem)
        self._groupList.setItem(idx,1,groupNameItem)
        self._groupList.setItem(idx,2,typeItem)
        self._groupList.setItem(idx,4,elementsItem)
        
      for idx, nodeSetId in enumerate(nodeSetIds):
        # collect data
        groupName = cubit.get_exodus_entity_name('nodeset', nodeSetId)
        
        # create QTableWidgetItems
        nodeSetIdItem = QTableWidgetItem(str(nodeSetId))
        groupNameItem = QTableWidgetItem(groupName)
        typeItem = QTableWidgetItem('NODE')
        
        if not groupName in names:
          names.append(groupName)
        else:
          msg = 'Node Set with ID %.0f and name %s was already used!\n' % (nodeSetId, groupName) \
              + '\tNACS region and group names allow no duplicates!\n' \
              + '\tPlease modify the Node Set Setup'
          self.status_txt.append(msg)
          groupNameItem.setForeground(errorBrush)
          hasErrors = True
  
        # send data to QTableWidget
        self._groupList.setItem(idx+len(sideSetIds),0,nodeSetIdItem)
        self._groupList.setItem(idx+len(sideSetIds),1,groupNameItem)
        self._groupList.setItem(idx+len(sideSetIds),2,typeItem)
        
      self._regionList.resizeColumnsToContents()
      self._groupList.resizeColumnsToContents()
      
      self._accept_pbtn.setDisabled(hasErrors)
    except Exception:
      print('Error\n%s' % traceback.format_exc())

class FloatValidator(QValidator):

  def validate(self, string, position):
    if valid_float_string(string):
      if sip.getapi('QString') == 1:
        return (self.Acceptable, position)
      else:
        return (self.Acceptable, string, position)
    if sip.getapi('QString') == 1:
      string = str(string) 
    if string == "" or string[position-1] in 'e.-+':
      if sip.getapi('QString') == 1:
        return (self.Intermediate, position)
      else:
        return (self.Intermediate, string, position)
    if sip.getapi('QString') == 1:
      return (self.Invalid, position)
    else:
      return (self.Invalid, string, position)

  def fixup(self, text):
    match = _float_re.search(text)
    if match:
      return match.groups()[0]
    else:
      return ""
    

class ScientificDoubleSpinBox(QDoubleSpinBox):

  def __init__(self, *args, **kwargs):
    super(ScientificDoubleSpinBox, self).__init__(*args, **kwargs)
    import numpy as np
    self.setMinimum(-np.inf)
    self.setMaximum(np.inf)
    self.validator = FloatValidator()
    self.setDecimals(1000)

  def validate(self, text, position):
    return self.validator.validate(text, position)

  def fixup(self, text):
    return self.validator.fixup(text)

  def valueFromText(self, text):
    return float(text)

  def textFromValue(self, value):
    return format_float(value)

  def stepBy(self, steps):
    text = self.cleanText()
    groups = _float_re.search(text).groups()
    decimal = float(groups[1])
    decimal += steps
    new_string = "{:g}".format(decimal) + (groups[3] if groups[3] else "")
    self.lineEdit().setText(new_string)

def format_float(value):
  """Modified form of the 'g' format specifier."""
  string = "{:g}".format(value).replace("e+", "e")
  string = re.sub("e(-?)0*(\d+)", r"e\1\2", string)
  return string
  
# Regular expression to find floats. Match groups are the whole string, the
# whole coefficient, the decimal part of the coefficient, and the exponent
# part.
_float_re = re.compile(r'(([+-]?\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)')

def valid_float_string(string):
  match = _float_re.search(string)
  if match:
    return match.groups()[0] == string
  else:
    return False