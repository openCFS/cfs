import paraview
import numpy

from paraview import vtk

from paraview.vtk import vtkPVCatalyst as catalyst
from paraview.vtk.util import numpy_support

import vtkPVPythonCatalystPython as pythoncatalyst
import paraview.simple

import vtkParallelCorePython
import sys, ntpath
import types

COPROCESSOR_MAP = {}
DATA_ARRAY_MAP = {}

def coprocessor_initialize():
    
    print('initialize coprocessor!')

    paraview.options.batch = True
    paraview.options.symmetric = True
    import vtkPVClientServerCoreCorePython as CorePython
    try:
        import vtkPVServerManagerApplicationPython as ApplicationPython
    except:
        paraview.print_error("Error: Cannot import vtkPVServerManagerApplicationPython")

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
    
    import vtkPVPythonCatalystPython as pythoncatalyst
    pipeline = pythoncatalyst.vtkCPPythonScriptPipeline()
    pipeline.Initialize("cpscript.py")
    coprocessor_coProcessor.AddPipeline(pipeline)
  
    return coprocessor_coProcessor

def send_data(key, xml, host, port):
  
  global COPROCESSOR_MAP
  
  catalyst_receive_key = host + ':' + str(port)
  
  if not key in COPROCESSOR_MAP:
    print("Setting up new coprocessor map for " + key)
    COPROCESSOR_MAP[key] = {}
    DATA_ARRAY_MAP[key] = types.SimpleNamespace()
    DATA_ARRAY_MAP[key].node_data_arr, DATA_ARRAY_MAP[key].element_data_arr = get_init_data_arrays(xml)
    
  
  if not catalyst_receive_key in COPROCESSOR_MAP[key]:
    print("Setting up new coprocessor map for " + key + ": " + catalyst_receive_key)
    
    dataDescription = catalyst.vtkCPDataDescription()
    dataDescription.SetTimeData(0, 0)
    dataDescription.AddInput(key)
  
    polyData = get_grid_obj(xml) # <-- vtkPolyData
        
    dataDescription.GetInputDescriptionByName(key).SetGrid(polyData)

    dataDescription.host = host
    dataDescription.port = port
    dataDescription.dataset_key = key
    
    COPROCESSOR_MAP[key][catalyst_receive_key] = coprocessor_initialize();
    COPROCESSOR_MAP[key][catalyst_receive_key].CoProcess(dataDescription)
    
  coprocessor_coProcessor = COPROCESSOR_MAP[key][catalyst_receive_key]

  if not coprocessor_coProcessor: # is NoneObject if error with python
    return

  time_steps = xml.xpath('//calculation/process/sequence/result/item/@step')
  time = 1
  if len(time_steps) > 0:
    time = int(time_steps[-1]) # get the latest value of the step attr which will be the total timesteps
  timeStep = time

  dataDescription = catalyst.vtkCPDataDescription()
  dataDescription.SetTimeData(float(time), timeStep)
  dataDescription.AddInput(key)

  if coprocessor_coProcessor.RequestDataDescription(dataDescription):
      polyData = get_grid_obj(xml) # <-- vtkPolyData
      
      for data_arr in DATA_ARRAY_MAP[key].element_data_arr:
        polyData.GetCellData().AddArray(DATA_ARRAY_MAP[key].element_data_arr[data_arr])
      #for data_arr in DATA_ARRAY_MAP[key].node_data_arr:
      #  polyData.GetCellData().AddArray(data_arr) ///TODO do this for nodes
      
      dataDescription.GetInputDescriptionByName(key).SetGrid(polyData)
      
      set_node_element_data(xml, DATA_ARRAY_MAP[key].node_data_arr, DATA_ARRAY_MAP[key].element_data_arr)
      
      dataDescription.host = host
      dataDescription.port = port
      dataDescription.dataset_key = key
      coprocessor_coProcessor.CoProcess(dataDescription)

  # the following 2 comments are presented by the catalyst example
  # if we are running through Python we need to finalize extra stuff
  # to avoid memory leak messages.
  #if ntpath.basename(sys.executable) == 'python':
  #    import vtkPVServerManagerApplicationPython as ApplicationPython
  #    ApplicationPython.vtkInitializationHelper.Finalize()

  # note to self: help(pipeline) is an interesting workaround to deleting the pipeline effectively.
  #     del or removeAll apparently does not work
  #pipeline = coprocessor_coProcessor.GetPipeline(0)
  #help(pipeline)
  #del pipeline

  #coprocessor_coProcessor.RemoveAllPipelines()
  #coprocessor_coProcessor.RemoveAllObservers()
  #del coprocessor_coProcessor

def get_grid_obj(xml):

  DIMENSIONS = xml.xpath('//domain/grids/grid/@dimensions')[0]
  
  #ELEMENT_COUNT = int(xml.xpath('//domain/grids/grid/@elements')[0])
  #NODE_COUNT = int(xml.xpath('//domain/grids/grid/@nodes')[0])
  
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

  for idx in range(len(node_id_str)):
    this_node_id = int(node_id_str[idx])
    node_x[this_node_id] = float(node_x_str[idx])
    node_y[this_node_id] = float(node_y_str[idx])
    node_z[this_node_id] = float(node_z_str[idx])

  pts = vtk.vtkPoints()
  
  for idx in range(NODE_COUNT):
    pts.InsertNextPoint(node_x[idx], node_y[idx], node_z[idx])

  cell_list = [0] * ELEMENT_COUNT # initialize all cells as 0

  for region_name in xml.xpath('//grid/regionList/region/@name'):
    this_region_element_count = len(xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element'))
    
    print('region_name: '+ region_name)
    
    type_arr = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element/@type')
    
    element_arr = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element')
    
    for element_idx in range(this_region_element_count):
      type = type_arr[element_idx]
      
      element = element_arr[element_idx]
      
      element_id = int(element.attrib['id'])
      
      if type == 'QUAD4':
        node_0_id = int(element.xpath('@node_0')[0])
        node_1_id = int(element.xpath('@node_1')[0])
        node_2_id = int(element.xpath('@node_2')[0])
        node_3_id = int(element.xpath('@node_3')[0])
        
        
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0,node_0_id)
        quad.GetPointIds().SetId(1,node_1_id)
        quad.GetPointIds().SetId(2,node_2_id)
        quad.GetPointIds().SetId(3,node_3_id)
        
        #poly_elements.InsertNextCell(quad)
        cell_list[element_id] = quad
      elif type == 'HEXA8':
        node_id = [0, 0, 0, 0, 0, 0, 0, 0]
        
        # See  https://www.evl.uic.edu/aej/524/lecture05.html
        node_id[0] = int(element.xpath('@node_' + str(4))[0])
        node_id[1] = int(element.xpath('@node_' + str(5))[0])
        node_id[2] = int(element.xpath('@node_' + str(1))[0])
        node_id[3] = int(element.xpath('@node_' + str(0))[0])
        
        node_id[4] = int(element.xpath('@node_' + str(7))[0])
        node_id[5] = int(element.xpath('@node_' + str(6))[0])
        node_id[6] = int(element.xpath('@node_' + str(2))[0])
        node_id[7] = int(element.xpath('@node_' + str(3))[0])
        
        #for tmp_indx in range(8):
        #  print(str(tmp_indx) + ': x=' + str(pts.GetPoint(node_id[tmp_indx])[0]) + ': y=' + str(pts.GetPoint(node_id[tmp_indx])[1]) + ': z=' + str(pts.GetPoint(node_id[tmp_indx])[2]))
        
        #return
        
        hexa8 = vtk.vtkHexahedron()
        for i in range(8):
          hexa8.GetPointIds().SetId(i,node_id[i])
        
        #poly_elements.InsertNextCell(hexa8)
        cell_list[element_id] = hexa8
      elif type == 'asdf1':
        print('WARNING: type "' + type + '" not supported (yet)!')
      elif type == 'asdf2':
        print('WARNING: type "' + type + '" not supported (yet)!')
      else:
        print('WARNING: type "' + type + '" not supported (yet)!')
  
  poly_elements = vtk.vtkCellArray()

  for cell in cell_list:
    if cell == 0:
      # use node #0 which is just at 0, 0, 0:
      dummy_quad = vtk.vtkQuad()
      dummy_quad.GetPointIds().SetId(0,0)
      dummy_quad.GetPointIds().SetId(1,0)
      dummy_quad.GetPointIds().SetId(2,0)
      dummy_quad.GetPointIds().SetId(3,0)
      poly_elements.InsertNextCell(dummy_quad)
    else:
      poly_elements.InsertNextCell(cell)

  pdo = vtk.vtkPolyData()

  pdo.SetPoints(pts)
  
  # Allocate memory for elements
  pdo.Allocate(ELEMENT_COUNT)
  
  pdo.SetPolys(poly_elements)
  
  pdo.tmp_bounds = pts.GetBounds()
  pdo.tmp_pts = pts
  
  return pdo

def get_init_data_arrays(xml):

  NODE_COUNT = max(map(int, xml.xpath('//domain/grids/grid/@nodes')))+1
  ELEMENT_COUNT = max(map(int, xml.xpath('//domain/grids/grid/@elements')))+1

  node_data_arr = {}
  element_data_arr = {}

  for region_name in xml.xpath('//grid/regionList/region/@name'):
    this_region_element_count = len(xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element'))
    
    print('region_name: '+ region_name)
    
    data_name_arr = xml.xpath('//streaming/process/sequence/result[@location="' + region_name + '"]/@data')
    data_defOn_arr = xml.xpath('//streaming/process/sequence/result[@location="' + region_name + '"]/@definedOn')
    
    print(data_name_arr)
    print(data_defOn_arr)
    
    type_arr = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element/@type')

    for result_index in range(len(data_name_arr)):
      
      data_name = data_name_arr[result_index]
      data_defOn = data_defOn_arr[result_index]

      if data_defOn == "element":
        
        dofs = int(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/@dofs')[0])
        
        print("element result: " + data_name + " has dofs=" + str(dofs))
        tmp_element_data_arr = numpy.zeros((ELEMENT_COUNT, dofs)).astype(numpy.double)
        
        value_data = numpy_support.numpy_to_vtk(tmp_element_data_arr)
        value_data.SetName(region_name + '/' + data_name)
  
        element_data_arr[region_name + '/' + data_name] = value_data
      elif data_defOn == "node":
        dofs = int(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/@dofs')[0])
        node_data_zeros = numpy.zeros((NODE_COUNT, dofs))
        
        
        value_data = numpy_support.numpy_to_vtk(node_data_zeros)
        value_data.SetName(region_name + '/' + data_name)
        
        node_data_arr[region_name + '/' + data_name] = value_data
      else:
        print("unknown defon")

  
  return node_data_arr, element_data_arr

def set_node_element_data(xml, node_data_arr, element_data_arr):
  NODE_COUNT = max(map(int, xml.xpath('//domain/grids/grid/@nodes')))+1
  ELEMENT_COUNT = max(map(int, xml.xpath('//domain/grids/grid/@elements')))+1
  
  for region_name in xml.xpath('//grid/regionList/region/@name'):
    this_region_element_count = len(xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element'))
    
    print('region_name: '+ region_name)
    
    data_name_arr = xml.xpath('//streaming/process/sequence/result[@location="' + region_name + '"]/@data')
    data_defOn_arr = xml.xpath('//streaming/process/sequence/result[@location="' + region_name + '"]/@definedOn')
    
    print(data_name_arr)
    print(data_defOn_arr)
    
    type_arr = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element/@type')

    for result_index in range(len(data_name_arr)):
      
      data_name = data_name_arr[result_index]
      data_defOn = data_defOn_arr[result_index]
      
      if data_defOn == "element":
        
        dofs = int(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/@dofs')[0])
        
        print("element result: " + data_name + " has dofs=" + str(dofs))
        
        this_data = element_data_arr[region_name + '/' + data_name]
        
        for tmp_element in xml.xpath('(//results/result[@name="' + data_name + '"][@region="' + region_name +'"])[last()]/item'):
          
          tuple_arr = numpy.zeros(dofs)
          
          for tuple_idx in range(dofs):
            tuple_arr[tuple_idx] = tmp_element.attrib['v_' + str(tuple_idx)]

          this_data.SetTuple(int(tmp_element.attrib['id']), tuple_arr)
        
      #elif data_defOn == "node":
        #dofs = int(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/@dofs')[0])
        #
        #print("node result: " + data_name + " has dofs=" + str(dofs))
        #node_data_arr = numpy.zeros((NODE_COUNT, dofs))
        
        #value_data = numpy_support.numpy_to_vtk(node_data_arr)
        #value_data.SetName(region_name + '/' + data_name)
        #pdo.GetCellData().SetScalars(value_data) TODO: fix
      else:
        print("unknows result")
  
  return
