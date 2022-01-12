import numpy as np
import matplotlib
matplotlib.use('Agg') # importing matplotlib and this line fixes segfaults when spamming plot too often
import matplotlib.pyplot as plt
from io import StringIO
import threading
import json
import svgwrite
import time
import traceback
from lxml import etree

# get the latest values
def get_values(GLOBAL_RAW_DATA_DICT, UPDATE_EVENTS, key, iteration_num):
  data = {}
  
  xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[key])
  
  # iteration_num is the latest iteration that was shown by the GUI
  # we will wait for a new data to arrive f 
  xp = xml.xpath('//process/iteration[last()]/@number')
  latest_received_iteration = int(xp[0]) if len(xp) > 0 else -1
  
  if latest_received_iteration > 0:
    if iteration_num >= latest_received_iteration: # wait 9 seconds to receive new xml
      
      e = None
      
      if key in UPDATE_EVENTS:
        e = UPDATE_EVENTS[key]
      else:
        e = threading.Event()
        UPDATE_EVENTS[key] = e
      
      print('request too early, waiting')
  
      # client waits 10 seconds to give the http connection a 1 second buffer to avoid errors
      e.wait(9)
      
      # remove the event to prevent cluttering
      if key in UPDATE_EVENTS:
        del UPDATE_EVENTS[key]
    
    if iteration_num >= latest_received_iteration: # still nothing new, tell the client javascript that there is nothing new
      print('no new data!')
      return '"no_new_data"'
  else:
    time.sleep(2)

  xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[key]) # load xml again in case there is a new one available

  data['results'] = {}

  for result in xml.xpath('//calculation/process/sequence/result'):
    dat = result.xpath('item[last()]')[0]
    if 'value' in dat:
      name = result.get('data') + '_' + result.get('location')
      data['results'][name] = name + ': ' + dat.get('value') + dat.get('unit')
  
  data['iterations'] = {}
  
  xp = xml.xpath('//process/iteration[last()]')
  if len(xp) > 0:
    value_map = xml.xpath('//process/iteration[last()]')[0]
    for name, value in value_map.items():
      data['iterations'][name] = name + ': ' + value
  
  return json.dumps(data)


# helper for error string when lengths are not matching
def length_error(values, name, ref_values):
  return 'eror: ' + name + ' has ' + len(values) + ' values, but ' + len(ref_values) + ' are exected' 

# helper for plot which plots only the 1d plots
def plot1d(data, xml, y1_it_names, y2_it_names, y1_res_names, y2_res_names, logscale_y1, logscale_y2):
  fig, ax1 = plt.subplots()
  fig.set_size_inches(10,8)
  
  # the x-data is either the iterations or the region result sequence steps
  x = [int(v) for v in xml.xpath('//process/iteration/@number')]
  ax1.set_xlabel("iteration")
  if len(x) == 0:
    x = [float(v) for v in xml.xpath('/cfsStreaming/cfsInfo/calculation/process/sequence/result[last()]/item/@step_val')]
    ax1.set_xlabel("step_val")
  
  y1_label = ""
  
  # set the plot mode to lines or dottet lines when we are not larger 20 iterations
  plot_mode = '-' if len(x) > 30 else '-o'
  
  for y1_it_name in y1_it_names:
    y1_it = [float(i) for i in xml.xpath('//process/iteration/@' + y1_it_name)]
  
    if len(x) != len(y1_it):
      data += length_error(y1_it, y1_it_name, x)
      continue
  
    y1_label += y1_it_name + ', '
    ax1.plot(x, y1_it, plot_mode, label=y1_it_name)
    
  for y1_res_name in y1_res_names:
    name, region = y1_res_name.split(':')
    # convert to list of floats
    y1_res = [float(i) for i in xml.xpath('/cfsStreaming/cfsInfo/calculation/process/sequence/result[@data="' + name + '"][@location="' + region + '"]/item/@value')]
  
    if len(x) != len(y1_res):
      data += length_error(y1_res, y1_res_name, x)
      continue

    y1_label += y1_res_name + ', '
    ax1.plot(x, y1_res, plot_mode, label=y1_res_name)
    
  ax1.set_ylabel(y1_label[:-2]) # remove trailing ', '
  
  if logscale_y1:
    ax1.set_yscale('log')
  
  ax1.legend(loc='lower left', bbox_to_anchor=(0, 1))
  ax2 = ax1.twinx()

  # "advance the colors" aka prevent color resetting when plotting for other axis
  ax2._get_lines.prop_cycler = ax1._get_lines.prop_cycler

  y2_label = ""
  
  for y2_it_name in y2_it_names:
    y2_it = [float(i) for i in xml.xpath('//process/iteration/@' + y2_it_name)]
  
    if len(x) != len(y2_it):
      data += length_error(y2_it, name, x)
      continue
  
    y2_label += y2_it_name + ', '
    ax2.plot(x, y2_it, plot_mode, label=y2_it_name)
    
  for y2_res_name in y2_res_names:
    name, region = y2_res_name.split(':')
    y2_res = [float(i) for i in xml.xpath('/cfsStreaming/cfsInfo/calculation/process/sequence/result[@data="' + name + '"][@location="' + region + '"]/item/@value')]

    if len(x) != len(y2_res):
      data += length_error(y2_res, y2_res_name, x)
      continue

    y2_label += y2_res_name + ', '
    ax2.plot(x, y2_res, plot_mode, label=y2_res_name)
    
  ax2.set_ylabel(y2_label[:-2])
  
  if logscale_y2:
    ax2.set_yscale('log')

  ax2.legend(loc='lower right', bbox_to_anchor=(1, 1))
  
  fig.tight_layout()

  imgdata = StringIO()
  fig.savefig(imgdata,format = 'svg')
  
  # this 'all' might cause some probolems if this tool
  # is used by multiple people at the same time
  plt.close('all')
  
  imgdata.seek(0)  # rewind the data

  svg_dta = imgdata.readlines()  # this is svg data

  # just to be sure:
  del imgdata

  # transfer everything from svg_data to our "data" output variable
  for l in svg_dta:
    data += l
  
  del svg_dta

  return data


# returns html with embedded svg OR error code 
def plot(key, UPDATE_EVENTS, GLOBAL_RAW_DATA_DICT, y1_it_names, y2_it_names, y1_res_names, y2_res_names, \
         iteration_num, logscale_y1, logscale_y2, result_view_arr, show_bloch):
  data = ""
  
  xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[key])
  
  # iteration_num is the latest iteration that was shown by the GUI
  # we will wait for a new data to arrive f
  xp = xml.xpath('//process/iteration[last()]/@number')
  latest_received_iteration = int(xp[0]) if len(xp) > 0 else -1
  
  if latest_received_iteration > 0:
    if iteration_num >= latest_received_iteration: # wait 9 seconds to receive new xml
      
      e = None
      
      if key in UPDATE_EVENTS:
        e = UPDATE_EVENTS[key]
      else:
        e = threading.Event()
        UPDATE_EVENTS[key] = e
      
      print('request too early, waiting')
  
      # client waits 10 seconds to give the http connection a 1 second buffer to avoid errors
      e.wait(9)
      
      # remove the event to prevent cluttering
      if key in UPDATE_EVENTS:
        del UPDATE_EVENTS[key]
    
    if iteration_num >= latest_received_iteration: # still nothing new, tell the client javascript that there is nothing new
      print('no new data!')
      return 'no_new_data'

  else:
    time.sleep(2)

  xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[key]) # load xml again in case there is a new one available

  # bloch plot
  if show_bloch:
    
    x_axis_values = [] # 1D array with 1, 2, 3, <latest step>, <latest step +1>
    
    y_axis_values_array = [] # 2D array. first Index is the mode, second the step. Latest step is first repeated
    
    last_result_elem = xml.xpath('//sequenceStep/eigenFrequency/result[last()]')[0]
    
    fig = plt.figure()

    plt.xlabel("direction")
    plt.ylabel("frequency [Hz]")

    # plot stuff
    
    first_wavevector = list(last_result_elem)[0]
    
    for tmp in first_wavevector:
      y_axis_values_array.append([])
    
    for wave_vector_elm in last_result_elem:
      this_step = int(wave_vector_elm.attrib['step'])
      
      x_axis_values.append(this_step)
      
      for mode in wave_vector_elm:
        y_axis_values_array[int(mode.attrib['nr'])-1].append(float(mode.attrib['frequency']))
    
    # append the first vector as last element
    for mode in first_wavevector:
        y_axis_values_array[int(mode.attrib['nr'])-1].append(float(mode.attrib['frequency']))
    x_axis_values.append(len(last_result_elem))
    
    for mode_nr in range(len(first_wavevector)):
      plt.plot(x_axis_values, y_axis_values_array[mode_nr], '-o')
    
    fig.tight_layout()
    plt.gca().set_ylim(bottom=0) # set bottom of plot to be always 0
    
    imgdata = StringIO()
    fig.savefig(imgdata, format = 'svg')
    
    # this 'all' might cause some probolems if this tool
    # is used by multiple people at the same time
    plt.close('all')
    
    imgdata.seek(0)  # rewind the data
  
    svg_dta = imgdata.readlines()  # this is svg data

    # transfer everything from svg_data to our "data" output variable
    for l in svg_dta:
      data += l
  
    # just to be sure:
    del imgdata
    del svg_dta
  
  # 2D plots
  if len(result_view_arr) > 0:
    domain_elem = xml.xpath('//header/domain')[0]
    
    _2d_orientation = "error"
    
    if int(domain_elem.attrib['nx']) == 1:
      _2d_orientation = 'yz'
    elif int(domain_elem.attrib['ny']) == 1:
      _2d_orientation = 'xz'
    elif int(domain_elem.attrib['nz']) == 1:
      _2d_orientation = 'xy'
    else:
      print('Error: not a 2D grid!')
      data += '<br/>Error: not a 2D grid!'
      # return data to avoid complete failure
      return data
  
    for this_result_name in result_view_arr:
      
      result_elm = xml.xpath('//results/result[@name="' + this_result_name + '"][last()]')[0]
      
      ax_x_descr = "x"
      ax_y_descr = "x"
      
      x_coords = []
      y_coords = []
      values = []
      
      values = np.array(result_elm.xpath('item/@v_0')).astype(np.float)

      if _2d_orientation == 'yz':
        ax_x_descr = 'y'
        ax_y_descr = 'z'
      elif _2d_orientation == 'xz':
        ax_x_descr = 'x'
        ax_y_descr = 'z'
      elif _2d_orientation == 'xy':
        ax_x_descr = 'x'
        ax_y_descr = 'y'
        
      x_node_coords = []
      y_node_coords = []
      if _2d_orientation == 'yz':
        x_node_coords = np.array(xml.xpath('//grid/nodeList/node/@y')).astype(np.float)
        tmparr = np.array(xml.xpath('//grid/nodeList/node/@z')).astype(np.float)
      elif _2d_orientation == 'xz':
        x_node_coords = np.array(xml.xpath('//grid/nodeList/node/@x')).astype(np.float)
        tmparr = np.array(xml.xpath('//grid/nodeList/node/@z')).astype(np.float)
      elif _2d_orientation == 'xy':
        x_node_coords = np.array(xml.xpath('//grid/nodeList/node/@x')).astype(np.float)
        tmparr = np.array(xml.xpath('//grid/nodeList/node/@y')).astype(np.float)

      y_node_coords = np.full(len(tmparr), max(tmparr)) - tmparr

      imgdata = StringIO()

      value_pos = 0
      if result_elm.attrib['solution'] == 'element':
        
        dwg = svgwrite.Drawing(profile='tiny')
        
        max_value = max(values)
        min_value = min(values)
        
        value_factor = 255/((max_value-min_value) if max_value != min_value else 1.0)
        
        for this_poly_element in xml.xpath('//regionList/region[@name="' + result_elm.attrib['region'] + '"]/element'):
          node_count = int(this_poly_element.attrib['nodes'])
          
          points = []
        
          for node_id in range(node_count):
            node_number = int(this_poly_element.attrib['node_'+str(node_id)])
            
            # we need to subtract one since the node start counting at 1, but we start counting at 0
            coord_tuple = (x_node_coords[node_number-1], y_node_coords[node_number-1])
            points.append(coord_tuple)

          color_value = (values[value_pos]-min_value) * value_factor

          #def tri_func(value, center):
          #  if value < center:
          #    return max((center-value)*2, 0)
          #  return max((value-center)*2, 0)

          #tmp_color = svgwrite.rgb(tri_func(color_value, 0), tri_func(color_value, 255/2), tri_func(color_value, 255))
          tmp_color = svgwrite.rgb(color_value, 0, 255-color_value)
          
          polygon = dwg.polygon(points=points, fill=tmp_color, stroke='black', stroke_width=0)
          dwg.add(polygon)
          
          value_pos += 1
          
        dwg.viewbox(min(x_node_coords), min(y_node_coords), max(x_node_coords), max(y_node_coords))
        dwg.write(fileobj=imgdata)
          
      elif result_elm.attrib['solution'] == 'node':
        fig = plt.figure()
        
        plt.colorbar()
        
        fig.tight_layout()
        
        fig.savefig(imgdata, format = 'svg')
        
        # this 'all' might cause some probolems if this tool
        # is used by multiple people at the same time
        plt.close('all')
      
        plt.tripcolor(x_node_coords, y_node_coords, values)

      else:
        data += "<br>unknown element solution type<br>"
        continue

      imgdata.seek(0)  # rewind the data
    
      svg_dta = imgdata.readlines()  # this is svg data
    
      # just to be sure:
      del imgdata

      data += "<br/>"
      data += this_result_name + ":<br/>"
      
      # transfer everything from svg_data to our "data" output variable
      for l in svg_dta:
        data += l
  
  fig = plt.figure()

  # 1D plot only if we have no 2D plots
  if len(result_view_arr) == 0:
    try:
      data = plot1d(data, xml, y1_it_names, y2_it_names, y1_res_names, y2_res_names, logscale_y1, logscale_y2)
    except Exception as e:
      data += 'exception occured: ' + str(e) + '<br>'
      print("".join(traceback.TracebackException.from_exception(e).format()))
      pass

  # data iteration num is used by javascript to determine the latest updated version
  # this is transmitted when requesting new data
  data += '<div id="iteration_num">' + str(latest_received_iteration) + '</div>'
  data += '<div id ="status">' + xml.xpath('//cfsInfo/@status')[0] + '</div>'
  
  return data
