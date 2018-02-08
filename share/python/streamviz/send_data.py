import paraview
import numpy

from paraview import vtk

from paraview.vtk import vtkPVCatalyst as catalyst
from paraview.vtk.util import numpy_support

import vtkPVPythonCatalystPython as pythoncatalyst
import paraview.simple

import vtkParallelCorePython
import sys, ntpath


COPROCESSOR_MAP = {}

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
    COPROCESSOR_MAP[key] = {}
  
  if not catalyst_receive_key in COPROCESSOR_MAP[key]:
    COPROCESSOR_MAP[key][catalyst_receive_key] = coprocessor_initialize();
    
  coprocessor_coProcessor = COPROCESSOR_MAP[key][catalyst_receive_key]

  if not coprocessor_coProcessor: # is NoneObject if error with python
    return

  time = xml.xpath('//calculation/process/sequence/result/item/@step')[-1] # get the latest value of the step attr which will be the total timesteps
  timeStep = time

  dataDescription = catalyst.vtkCPDataDescription()
  dataDescription.SetTimeData(time, timeStep)
  dataDescription.AddInput(key)

  print("crocess time", time, " timeStep", timeStep)

  if coprocessor_coProcessor.RequestDataDescription(dataDescription):
      imageData = get_grid_obj(xml) # <-- vtkPolyData
      
      dataDescription.GetInputDescriptionByName(key).SetGrid(imageData)
      
      for data_arr in get_data_arrays(xml):
        imageData.GetCellData().AddArray(data_arr)
      
      dataDescription.host = host
      dataDescription.port = port
      dataDescription.dataset_key = key
      coprocessor_coProcessor.CoProcess(dataDescription)

  # the following 2 comments are presented by the catalyst example
  # if we are running through Python we need to finalize extra stuff
  # to avoid memory leak messages.
  if ntpath.basename(sys.executable) == 'python':
      import vtkPVServerManagerApplicationPython as ApplicationPython
      ApplicationPython.vtkInitializationHelper.Finalize()

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
  
  ELEMENT_COUNT = int(xml.xpath('//domain/grids/grid/@elements')[0])
  NODE_COUNT = int(xml.xpath('//domain/grids/grid/@nodes')[0])
  
  node_id_str = xml.xpath('//grid/nodeList/node/@id')
  node_x_str = xml.xpath('//grid/nodeList/node/@x')
  node_y_str = xml.xpath('//grid/nodeList/node/@y')
  node_z_str = xml.xpath('//grid/nodeList/node/@z')
  
  node_id = []
  node_x = []
  node_y = []
  node_z = []

  print("element count: " + str(NODE_COUNT))
  for idx in range(NODE_COUNT):
     node_id.append(int(node_id_str[idx]))
     node_x.append(float(node_x_str[idx]))
     node_y.append(float(node_y_str[idx]))
     node_z.append(float(node_z_str[idx]))

  pts = vtk.vtkPoints()
  
  for idx in range(NODE_COUNT):
    pts.InsertNextPoint(node_x[idx], node_y[idx], node_z[idx])

  poly_elements = vtk.vtkCellArray()

  pdo = vtk.vtkPolyData()
  
  for region_name in xml.xpath('//grid/regionList/region/@name'):
    this_region_element_count = len(xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element'))
    
    print('region_name: '+ region_name)
    
    type_arr = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element/@type')
    
    for element_idx in range(this_region_element_count):
      type = type_arr[element_idx]
      
      element = xml.xpath('//grid/regionList/region[@name="' + region_name + '"]/element')[element_idx]
      
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
        
        poly_elements.InsertNextCell(quad)
      elif type == 'asdf':
        print('WARNING: type "' + type + '" not supported (yet)!')
      else:
        print('WARNING: type "' + type + '" not supported (yet)!')
  
  pdo.SetPoints(pts)
  
  # Allocate memory for elements
  pdo.Allocate(ELEMENT_COUNT)
  
  pdo.SetPolys(poly_elements)
  
  return pdo

def get_data_arrays(xml):
  
  DIMENSIONS = xml.xpath('//domain/grids/grid/@dimensions')[0]
  
  ELEMENT_COUNT = int(xml.xpath('//domain/grids/grid/@elements')[0])
  NODE_COUNT = int(xml.xpath('//domain/grids/grid/@nodes')[0])
  
  data_arr = []
  
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
        element_data_arr = numpy.zeros((ELEMENT_COUNT, dofs)).astype(numpy.float)
        
        for dimension in range(dofs):
          element_data_arr[:,dimension] = numpy.array(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/item/@v_' + str(dimension))).astype(numpy.float)
        
        value_data = numpy_support.numpy_to_vtk(element_data_arr)
        value_data.SetName(region_name + '/' + data_name)

        data_arr.append(value_data)
      elif data_defOn == "node":
        dofs = int(xml.xpath('//results/result[@name="' + data_name + '"][@region="' + region_name +'"]/@dofs')[0])
        
        print("node result: " + data_name + " has dofs=" + str(dofs))
        #node_data_arr = numpy.zeros((NODE_COUNT, dofs))
        
        #value_data = paraview.numpy_support.numpy_to_vtk(node_data_arr)
        #value_data.SetName(region_name + '/' + data_name)
        #pts.SetData(value_data)
      else:
        print("unknows result")
  
  return data_arr
