import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
import threading

# returns html with embedded svg OR error code 
def plot(key, UPDATE_EVENTS, GLOBAL_DATA_DICT, x_name, y1_it_names, y2_it_names, y1_res_names, y2_res_names, iteration_num):
  
  data = ""
  
  xml = GLOBAL_DATA_DICT[key]
  
  # iteration_num is the latest iteration that was shown by the GUI
  # we will wait for a new data to arrive f
  latest_received_iteration = xml.xpath('//process/iteration[last()]/@number')[0]
  
  if iteration_num >= latest_received_iteration: # wait 9 seconds to receive new xml
    e = threading.Event()
    
    print('request too early, waiting')

    if not key in UPDATE_EVENTS:
      UPDATE_EVENTS[key] = []
    
    UPDATE_EVENTS[key].append(e)
    
    # client waits 10 seconds to give the http connection a 1 second buffer to avoid errors
    e.wait(9)
    
    # remove the event to prevent cluttering
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
  
  for y1_it_name in y1_it_names: 
    y1_it = xml.xpath('//process/iteration/@' + y1_it_name)
  
    if len(t) != len(y1_it):
      data += 'error len(' + y1_it_name + ') != len(' + x_name + ')!<br>'
      data += str(len(y1_it_name)) + '!=' + str(len(x_name)) + '<br><br>'
      continue
  
    y1_label += y1_it_name + ', '
    ax1.plot(t, y1_it, '-')

  for y1_res_name in y1_res_names:
    y1_res = xml.xpath('//calculation/process/sequence/result[@data="' + y1_res_name + '"]/item/@value')
  
    if len(t) != len(y1_res):
      data += 'error len(' + y1_res_names + ') != len(' + x_name + '):<br>'
      data += str(len(y1_res_names)) + '!=' + str(len(x_name)) + '<br><br>'
      continue

    y1_label += y1_res_name + ', '
    ax1.plot(t, y1_res, '-')

  ax1.set_xlabel(x_name)
  
  ax1.set_ylabel(y1_label[:-2])
  
  ax1.autoscale(True, 'both')
  
  ax2 = ax1.twinx()
  
  y2_label = ""
  
  for y2_it_name in y2_it_names: 
    y2_it = xml.xpath('//process/iteration/@' + y2_it_name)
  
    if len(t) != len(y2_it):
      data += 'error len(' + y2_it_name + ') != len(' + x_name + ')!<br>'
      data += str(len(y2_it_name)) + '!=' + str(len(x_name)) + '<br><br>'
      continue
  
    y2_label += y2_it_name + ', '
    ax2.plot(t, y2_it, '-')

  for y2_res_name in y2_res_names:
    y2_res = xml.xpath('//calculation/process/sequence/result[@data="' + y2_res_name + '"]/item/@value')
  
    if len(t) != len(y2_res):
      data += 'error len(' + y2_res_names + ') != len(' + x_name + ')!<br>'
      data += str(len(y2_res_names)) + '!=' + str(len(x_name)) + '<br><br>'
      continue

    y2_label += y2_res_name + ', '
    ax2.plot(t, y2_res, '-')

  ax2.set_ylabel(y2_label[:-2])
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
  data += '<div id="iteration_num">' + latest_received_iteration + '</div>'
  
  for l in svg_dta:
    data += l
  
  return data

  return "asdf ldjsfiojse"
