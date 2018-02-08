import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import threading
import json

def get_values(GLOBAL_DATA_DICT, UPDATE_EVENTS, key, iteration_num):
  data = {}
  
  xml = GLOBAL_DATA_DICT[key]
  
  # iteration_num is the latest iteration that was shown by the GUI
  # we will wait for a new data to arrive f
  latest_received_iteration = int(xml.xpath('//process/iteration[last()]/@number')[0])

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


  xml = GLOBAL_DATA_DICT[key] # load xml again in case there is a new one available

  data['results'] = {}

  for result in xml.xpath('//calculation/process/sequence/result'):
    name = result.xpath('@data')[0]
    
    data['results'][name] = name + ': ' + result.xpath('item[last()]')[0].xpath('@value')[0] + \
                              ' ' + result.xpath('item[last()]')[0].xpath('@unit')[0]
  
  data['iterations'] = {}
  
  value_map = xml.xpath('//process/iteration[last()]')[0]
  for this_tuple in value_map.items():
    name = this_tuple[0]
    value = this_tuple[1]
    
    data['iterations'][name] = name + ': ' + value
  
  return json.dumps(data)


# returns html with embedded svg OR error code 
def plot(key, UPDATE_EVENTS, GLOBAL_DATA_DICT, x_name, y1_it_names, y2_it_names, y1_res_names, y2_res_names, iteration_num, logscale_y1, logscale_y2):
  
  data = ""
  
  xml = GLOBAL_DATA_DICT[key]
  
  # iteration_num is the latest iteration that was shown by the GUI
  # we will wait for a new data to arrive f
  latest_received_iteration = int(xml.xpath('//process/iteration[last()]/@number')[0])
  
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


  xml = GLOBAL_DATA_DICT[key] # load xml again in case there is a new one available

  fig = plt.figure()

  x = xml.xpath('//process/iteration/@' + x_name)
  
  fig, ax1 = plt.subplots()
  t = x
  
  y1_label = ""
  
  y1_min = []
  y1_max = []
  
  for y1_it_name in y1_it_names:
    # convert to list of floats
    y1_it = [float(i) for i in xml.xpath('//process/iteration/@' + y1_it_name)]
  
    if len(t) != len(y1_it):
      data += 'error len(' + y1_it_name + ') != len(' + x_name + ')!<br>'
      data += str(len(y1_it_name)) + '!=' + str(len(x_name)) + '<br><br>'
      continue
  
    y1_label += y1_it_name + ', '
    ax1.plot(t, y1_it, '-')
    
    y1_min.append(min(y1_it))
    y1_max.append(max(y1_it))

  for y1_res_name in y1_res_names:
    # convert to list of floats
    y1_res = [float(i) for i in xml.xpath('//calculation/process/sequence/result[@data="' + y1_res_name + '"]/item/@value')]
  
    if len(t) != len(y1_res):
      data += 'error len(' + y1_res_names + ') != len(' + x_name + '):<br>'
      data += str(len(y1_res_names)) + '!=' + str(len(x_name)) + '<br><br>'
      continue

    y1_label += y1_res_name + ', '
    ax1.plot(t, y1_res, '-')
    
    y1_min.append(min(y1_res))
    y1_max.append(max(y1_res))

  if x_name == "name":
    ax1.set_xlabel("iteration")
  else:
    ax1.set_xlabel(x_name)
  
  ax1.set_ylabel(y1_label[:-2])
  
  if logscale_y1:
    ax1.set_yscale('log')
    if y1_min:
      ax1.set_ylim(min(y1_min), max(y1_max))
  else:
    ax1.autoscale(True, 'both')
  
  ax2 = ax1.twinx()
  
  y2_label = ""
  
  y2_min = []
  y2_max = []
  
  for y2_it_name in y2_it_names:
    # convert to list of floats
    y2_it = [float(i) for i in xml.xpath('//process/iteration/@' + y2_it_name)]
  
    if len(t) != len(y2_it):
      data += 'error len(' + y2_it_name + ') != len(' + x_name + ')!<br>'
      data += str(len(y2_it_name)) + '!=' + str(len(x_name)) + '<br><br>'
      continue
  
    y2_label += y2_it_name + ', '
    ax2.plot(t, y2_it, '-')
    
    y2_min.append(min(y2_it))
    y2_max.append(max(y2_it))

  for y2_res_name in y2_res_names:
    # convert to list of floats
    y2_res = [float(i) for i in xml.xpath('//calculation/process/sequence/result[@data="' + y2_res_name + '"]/item/@value')]
  
    if len(t) != len(y2_res):
      data += 'error len(' + y2_res_names + ') != len(' + x_name + ')!<br>'
      data += str(len(y2_res_names)) + '!=' + str(len(x_name)) + '<br><br>'
      continue

    y2_label += y2_res_name + ', '
    ax2.plot(t, y2_res, '-')
    
    y2_min.append(min(y2_res))
    y2_max.append(max(y2_res))

  ax2.set_ylabel(y2_label[:-2])
  
  if logscale_y2:
    ax2.set_yscale('log')
    if y2_min:
      ax2.set_ylim(min(y2_min), max(y2_max))
  else:
    ax2.autoscale(True, 'both')
  
  fig.tight_layout()
  
  imgdata = StringIO()
  fig.savefig(imgdata, format = 'svg')
  
  # this 'all' might cause some probolems if this tool
  # is used by multiple people at the same time
  # TODO: make sure it doesn't break the program
  plt.close('all')
  
  imgdata.seek(0)  # rewind the data

  svg_dta = imgdata.readlines()  # this is svg data

  # just to be sure:
  del imgdata

  # data iteration num is used by javascript to determine the latest updated version
  # this is transmitted when requesting new data
  data += '<div id="iteration_num">' + str(latest_received_iteration) + '</div>'
  data += '<div id ="status">' + xml.xpath('//cfsInfo/@status')[0] + '</div>'
  
  # transfer everything from svg_data to our "data" output variable
  for l in svg_dta:
    data += l
  
  return data
