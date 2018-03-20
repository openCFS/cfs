import paraview
import numpy

from paraview import vtk

from paraview.vtk.util import numpy_support

from paraview.vtk import vtkPVCatalyst as catalyst
from paraview.vtk import vtkPVPythonCatalystPython as pythoncatalyst
from paraview.vtk import vtkParallelCorePython
from paraview.vtk import vtkPVClientServerCoreCorePython as CorePython
from paraview.vtk import vtkPVServerManagerApplicationPython as ApplicationPython

import paraview.simple

import sys, ntpath
import types

COPROCESSOR_MAP = {}

def delete_coprocessor(key):
  if not key in COPROCESSOR_MAP:
    return
  tmp_co = COPROCESSOR_MAP.pop(key)
  del tmp_co

# initialize the coprocessor
# this is mainly taken from https://gitlab.kitware.com/paraview/paraview/blob/master/Examples/Catalyst/PythonFullExample/coprocessor.py
# and adapted 
def coprocessor_initialize():
    
    print('initialize coprocessor!')

    paraview.options.batch = True
    paraview.options.symmetric = True
    if not CorePython.vtkProcessModule.GetProcessModule():
        pvoptions = None
        if paraview.options.batch:
            pvoptions = CorePython.vtkPVOptions();
            pvoptions.SetProcessType(CorePython.vtkPVOptions.PVBATCH)
            if paraview.options.symmetric:
                pvoptions.SetSymmetricMPIMode(True)
        ApplicationPython.vtkInitializationHelper.Initialize(sys.executable, CorePython.vtkProcessModule.PROCESS_BATCH, pvoptions)

    import paraview.servermanager as pvsm
    # we need ParaView 4.2 since ParaView 4.1 doesn't properly wrap
    # vtkPVPythonCatalystPython
    if pvsm.vtkSMProxyManager.GetVersionMajor() < 4 or (pvsm.vtkSMProxyManager.GetVersionMajor() == 4 and pvsm.vtkSMProxyManager.GetVersionMinor() < 2):
        print('Must use ParaView v4.2 or greater')
        return None

    paraview.options.batch = True
    paraview.options.symmetric = True

    coprocessor_coProcessor = catalyst.vtkCPProcessor()
    pm = paraview.servermanager.vtkProcessModule.GetProcessModule()

    pipeline = pythoncatalyst.vtkCPPythonScriptPipeline()
    pipeline.Initialize("cpscript.py")
    coprocessor_coProcessor.AddPipeline(pipeline)
  
    return coprocessor_coProcessor

# this sends the data
# if no coprocessor exists, one will be created
# PLEASE NOTE: NODE DATA NOT YET IMPLEMENTED
def send_data(key, xml, host, port):
  
  global COPROCESSOR_MAP
  
  catalyst_receive_key = host + ':' + str(port)
  
  if not key in COPROCESSOR_MAP:
    print("Setting up new coprocessor map for " + key)
    COPROCESSOR_MAP[key] = {}
  
  if not catalyst_receive_key in COPROCESSOR_MAP[key]:
    print("Setting up new coprocessor map for " + key + ": " + catalyst_receive_key)
    COPROCESSOR_MAP[key][catalyst_receive_key] = coprocessor_initialize();
    
  coprocessor_coProcessor = COPROCESSOR_MAP[key][catalyst_receive_key]

  if not coprocessor_coProcessor: # is NoneObject if error with python
    return

  time_steps = xml.xpath('//calculation/process/sequence/result/item/@step')
  time = 1
  if len(time_steps) > 0:
    time = int(time_steps[-1]) # get the latest value of the step attr which will be the total timesteps
  timeStep = time

  # create the datapackage we are sending
  dataDescription = catalyst.vtkCPDataDescription()
  dataDescription.SetTimeData(float(time), timeStep)
  dataDescription.AddInput(key)

  if coprocessor_coProcessor.RequestDataDescription(dataDescription):
    # only send data if a catalyst is connected
    
    # build get the grid!
    polyData = get_grid_obj(xml) # <-- vtkPolyData

    # get the node and element data
    node_data_arr, element_data_arr = get_data_arrays(xml)

    # and attach it to the data
    for data_arr in element_data_arr:
      polyData.GetCellData().AddArray(element_data_arr[data_arr])
    #for data_arr in DATA_ARRAY_MAP[key].node_data_arr:
    #  polyData.GetCellData().AddArray(data_arr) ///TODO do this for nodes
    
    # attack the data to the data that needs to be send
    dataDescription.GetInputDescriptionByName(key).SetGrid(polyData)
    
    dataDescription.host = host
    dataDescription.port = port
    dataDescription.dataset_key = key
    
    # and send the data to the correct ip:port
    coprocessor_coProcessor.CoProcess(dataDescription)

def get_grid_obj(xml):

  DIMENSIONS = xml.xpath('//domain/grids/grid/@dimensions')[0]

  # first, reserve space so all node ids will fit in. potentially having a lot of nodes from 0 to x empty
  # aka we map the node and element ids 1:1
  ELEMENT_COUNT = max(list(map(int, xml.xpath('//grid/regionList/region/element/@id')))) + 1
  NODE_COUNT = max(list(map(int, xml.xpath('//grid/nodeList/node/@id')))) + 1
  
  node_id_str = xml.xpath('//grid/nodeList/node/@id')
  node_x_str = xml.xpath('//grid/nodeList/node/@x')
  node_y_str = xml.xpath('//grid/nodeList/node/@y')
  node_z_str = xml.xpath('//grid/nodeList/node/@z')
  
  # initialize all cells at x, y, z = 0
  node_x = [0] * NODE_COUNT
  node_y = [0] * NODE_COUNT
  node_z = [0] * NODE_COUNT

  print("element count: " + str(ELEMENT_COUNT))
  print("node count: " + str(NODE_COUNT))

  # iterate over all xells in the xml
  for idx in range(len(node_id_str)):
    this_node_id = int(node_id_str[idx])
    # and set the coordinates
    node_x[this_node_id] = float(node_x_str[idx])
    node_y[this_node_id] = float(node_y_str[idx])
    node_z[this_node_id] = float(node_z_str[idx])

  pts = vtk.vtkPoints()
  
  # insert all points
  for idx in range(NODE_COUNT):
    pts.InsertNextPoint(node_x[idx], node_y[idx], node_z[idx])

  # we now have all points from the xml as vtkPoints
  # the mapping of the id is 1:1, points not defined in the
  # xml have x=y=z=0

  cell_list = [0] * ELEMENT_COUNT # initialize all cells as 0

  # same procedure with the elements
  for region_name in xml.xpath('//grid/regionList/region/@name'):
    element_arr = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element')
    this_region_element_count = len(element_arr)
    
    print('region_name: '+ region_name)
    
    for element_idx in range(this_region_element_count):
      
      element = element_arr[element_idx]
      type = element.attrib['type']
      
      element_id = int(element.attrib['id'])
      
      if type == 'UNDEF':
        print('warning: undefined element')
      elif type == 'POINT':
        point = vtk.vtkVertex()
        point.GetPointIds().SetId(0,int(element.attrib['node_0']))
        cell_list[element_id] = point
      elif type == 'LINE2':
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0,int(element.attrib['node_0']))
        line.GetPointIds().SetId(1,int(element.attrib['node_1']))
        cell_list[element_id] = line
      elif type == 'LINE3':
        quadraticEdge = vtk.vtkQuadraticEdge()
        for i in range(3):
          quadraticEdge.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quadraticEdge
      elif type == 'TRIA3':
        polygon = vtk.vtkPolygon()
        for i in range(4):
          polygon.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = polygon
      elif type == 'TRIA6':
        quadraticTriangle = vtk.vtkQuadraticTriangle()
        for i in range(6):
          quadraticTriangle.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quadraticTriangle
      elif type == 'QUAD4':
        quad = vtk.vtkQuad()
        for i in range(4):
          quad.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quad
      elif type == 'QUAD8':
        quadraticQuad = vtk.vtkQuadraticQuad()
        for i in range(8):
          quadraticQuad.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quadraticQuad
      elif type == 'QUAD9':
        polygon = vtk.vtkPolygon() # node 9 in the middle, then a circle around
        polygon.GetPointIds().SetId(0,int(element.attrib['node_8']))
        polygon.GetPointIds().SetId(1,int(element.attrib['node_0']))
        polygon.GetPointIds().SetId(2,int(element.attrib['node_4']))
        polygon.GetPointIds().SetId(3,int(element.attrib['node_1']))
        polygon.GetPointIds().SetId(4,int(element.attrib['node_5']))
        polygon.GetPointIds().SetId(5,int(element.attrib['node_2']))
        polygon.GetPointIds().SetId(6,int(element.attrib['node_6']))
        polygon.GetPointIds().SetId(7,int(element.attrib['node_3']))
        polygon.GetPointIds().SetId(8,int(element.attrib['node_7']))
        cell_list[element_id] = polygon
      elif type == 'TET4':
        tetra = vtk.vtkTetra()
        for i in range(4):
          tetra.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = tetra
      elif type == 'TET10':
        quadraticTetra = vtk.vtkQuadraticTetra()
        for i in range(10):
          quadraticTetra.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quadraticTetra
      elif type == 'HEXA8':
        node_id = [0, 0, 0, 0, 0, 0, 0, 0]
        
        # See  https://datascience.lanl.gov/data/ParaViewCatalyst4_2Tutorial.pdf
        node_id[0] = int(element.attrib['node_4'])
        node_id[1] = int(element.attrib['node_5'])
        node_id[2] = int(element.attrib['node_1'])
        node_id[3] = int(element.attrib['node_0'])
        
        node_id[4] = int(element.attrib['node_7'])
        node_id[5] = int(element.attrib['node_6'])
        node_id[6] = int(element.attrib['node_2'])
        node_id[7] = int(element.attrib['node_3'])

        hexa8 = vtk.vtkHexahedron()
        for i in range(8):
          hexa8.GetPointIds().SetId(i,node_id[i])

        cell_list[element_id] = hexa8
      elif type == 'HEXA20':
        quadraticHexahedron = vtk.vtkQuadraticHexahedron()
        for i in range(20):
          quadraticHexahedron.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quadraticHexahedron
      elif type == 'HEXA27': # same as 20; TODO: complete
        quadraticHexahedron = vtk.vtkQuadraticHexahedron()
        for i in range(20):
          quadraticHexahedron.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = quadraticHexahedron
      elif type == 'PYRA5':
        pyramid = vtk.vtkPyramid()
        for i in range(5):
          pyramid.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = pyramid
      elif type == 'PYRA13': # same as 5; TODO: complete
        pyramid = vtk.vtkPyramid()
        for i in range(5):
          pyramid.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = pyramid
      elif type == 'PYRA14': # same as 5; TODO: complete
        pyramid = vtk.vtkPyramid()
        for i in range(5):
          pyramid.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = pyramid
      elif type == 'WEDGE6':
        wedge = vtk.vtkWedge()
        for i in range(6):
          wedge.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = wedge
      elif type == 'WEDGE15': # same as 6; TODO: complete
        wedge = vtk.vtkWedge()
        for i in range(6):
          wedge.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = wedge
      elif type == 'WEDGE18': # same as 6; TODO: complete
        wedge = vtk.vtkWedge()
        for i in range(6):
          wedge.GetPointIds().SetId(i,int(element.attrib['node_'+str(i)]))
        cell_list[element_id] = wedge
      elif type == 'POLYGON':
        polygon = vtk.vtkPolygon() # node 9 in the middle, then a circle around
        
        for node_id in range(100): # 100 poly limit
          if not ('node_'+str(i)) in element.attrib:
            break;
          polygon.GetPointIds().SetId(node_id,int(element.attrib['node_'+str(node_id)]))

        cell_list[element_id] = polygon
      elif type == 'POLYHEDRON':
        print('WARNING: type "' + type + '" not supported (yet)!')
      else:
        print('WARNING: type "' + type + '" not supported (yet)!')
  
  pdo = vtk.vtkUnstructuredGrid()

  pdo.SetPoints(pts)
  
  # Allocate memory for elements
  pdo.Allocate(ELEMENT_COUNT)

  for cell in cell_list:
    # if cell not present in xml, add a single point element using point #0
    if cell == 0:
      # use node #0 which is just at 0, 0, 0:
      dummy_vertex = vtk.vtkVertex()
      dummy_vertex.GetPointIds().SetId(0,0)
      
      pdo.InsertNextCell(dummy_vertex.GetCellType(), dummy_vertex.GetPointIds())
    else:
      pdo.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

  
  pdo.tmp_bounds = pts.GetBounds()
  pdo.tmp_pts = pts

  return pdo

def get_data_arrays(xml):

  NODE_COUNT = max(map(int, xml.xpath('//domain/grids/grid/@nodes')))+1
  ELEMENT_COUNT = max(map(int, xml.xpath('//domain/grids/grid/@elements')))+1

  node_data_arr = {}
  element_data_arr = {}
  
  # the result tags within sequence are unique
  # where the result tags within "results" are appended.
  # since we only want to send the latest values,
  # we take the latest tag <results><result>
  # with the correct name and region

  for sequence_result_element in xml.xpath('//sequence/result'):
    seq_data_name = sequence_result_element.attrib['data']
    seq_region_name = sequence_result_element.attrib['location']
    
    result_element = None
    
    try:
      result_element = xml.xpath('//results/result[@name="' + seq_data_name + '"][@region="' + seq_region_name +'"][last()]')[0]
      break
    except:
      continue
    
    data_name = result_element.attrib['name']
    data_defOn = result_element.attrib['solution']
    region_name = result_element.attrib['region']
    dofs = result_element.attrib['dofs']
  
    #print("data_name: " + data_name)
    #print("data_defOn: " + data_defOn)
    #print("region_name: " + region_name)
    #print("dofs: " + dofs)
    
    if data_defOn == "element":
       
      dofs = int(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/@dofs')[0])

      print("element result: " + data_name + " has dofs=" + str(dofs))
      
      #initialize value 2D array with zeros
      tmp_element_data_arr = numpy.zeros((ELEMENT_COUNT, dofs)).astype(numpy.float)
       
      value_data = numpy_support.numpy_to_vtk(tmp_element_data_arr, deep=True)
      value_data.SetName(region_name + '/' + data_name)

      for tmp_element in list(result_element):
        
        # set all tuples from xml into vtkData Array
        tuple_arr = numpy.zeros(dofs, dtype=float)
         
        for tuple_idx in range(dofs):
          tuple_arr[tuple_idx] = tmp_element.attrib['v_' + str(tuple_idx)]

        value_data.SetTuple(int(tmp_element.attrib['id']), tuple_arr)

      # add that array with a recognizeable key
      element_data_arr[region_name + '/' + data_name] = value_data

  # return node and element vtkDataArray
  return node_data_arr, element_data_arr
