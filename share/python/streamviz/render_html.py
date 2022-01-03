from flask.globals import request

from strings import *

import urllib
from audioop import reverse

import datetime
import os
import psutil
import objgraph
import random, string
import base64
import html
from lxml import etree

#The html file with marker; see strings.py for all marker definitions 
html_raw_data = "html file not loaded"
with open('template_html/index.html', 'r') as myfile:
  html_raw_data = myfile.read()

def randomword(length):
   letters = string.ascii_lowercase
   return ''.join(random.choice(letters) for i in range(length))

def get_human_readable_bytes(memory_in_bytes):
  if memory_in_bytes > 1024*1024*1024:
    return ("%.3f" % (memory_in_bytes/(1024*1024*1024))) + ' GByte'
  elif memory_in_bytes > 1024*1024:
    return ("%.3f" % (memory_in_bytes/(1024*1024))) + ' MByte'
  elif memory_in_bytes > 1024:
    return ("%.3f" % (memory_in_bytes/1024)) + ' KByte'
  else:
    return str(memory_in_bytes) + ' Byte'

#reder the selection menu
def render_menu(GLOBAL_RAW_DATA_DICT, current_site):
  ret_string  = '<div class="btn-group">'
  if current_site == 'index':
    ret_string += '<a class="nav-link text-white bg-primary" href="/">overview</a>'
  else:
    ret_string += '<a class="nav-link" href="/">overview</a>'
  if current_site == 'status':
    ret_string += '<a class="nav-link text-white bg-primary" href="/status">status</a>'
  else:
    ret_string += '<a class="nav-link" href="/status">status</a>'
  if current_site == 'status_memory_pics':
    ret_string += '<a class="nav-link text-white bg-primary" href="/status_memory_pics">memory status</a>'
  if current_site == 'status_log':
    ret_string += '<a class="nav-link text-white bg-primary" href="/status_log">log</a>'

  ret_string += '</div>'
  return ret_string # menu is too slow to render and useless anyway

  hosts = {}
  for key in GLOBAL_RAW_DATA_DICT:
    this_host = key[:key.index('/')]
    
    if not this_host in hosts:
      hosts[this_host] = {}
    
    project_time_rest = key[key.index('/')+1:]
    
    this_problem = project_time_rest[:project_time_rest.index('/')]
    
    if not this_problem in hosts[this_host]:
      hosts[this_host][this_problem] = {}

    this_timekey = project_time_rest[project_time_rest.index('/')+1:]
    
    hosts[this_host][this_problem][this_timekey] = key
  
  for this_host in hosts:
    ret_string += '<div class="dropdown">' + "\n"
    ret_string += '<button class="btn btn-secondary dropdown-toggle" type="button" id="dropdownMenuButton_' + this_host
    ret_string += '" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">'
    ret_string += this_host + '</button> <div class="dropdown-menu small" aria-labelledby="dropdownMenuButton_' + this_host + '">' + "\n"
    
    for this_problem in hosts[this_host]:
      if len(hosts[this_host][this_problem]) == 1:
        this_timekey = list(hosts[this_host][this_problem].keys())[0]
        xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[hosts[this_host][this_problem][this_timekey]])
        ret_string += '<a class="dropdown-item" href="/view/' + hosts[this_host][this_problem][this_timekey] + '">'
        ret_string += '[' + xml.xpath('//cfsInfo/@status')[0] + '] '
        ret_string += this_problem + '/' + this_timekey + '</a>' + "\n"
        
      else:
        ret_string += '<div class="dropdown-submenu dropright">' + "\n"
        ret_string += '<button class="btn btn-secondary dropdown-toggle" type="button" id="dropdownMenuButton_' + this_host + '_' + this_problem
        ret_string += '" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">'
        ret_string += this_problem + '</button> <div class="dropdown-menu" aria-labelledby="dropdownMenuButton_' + this_host + '_' + this_problem + '">' + "\n"
      
        for this_timekey in hosts[this_host][this_problem]:
          xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[hosts[this_host][this_problem][this_timekey]])
          ret_string += '<a class="dropdown-item" href="/view/' + hosts[this_host][this_problem][this_timekey] + '">'
          ret_string += '[' + xml.xpath('//cfsInfo/@status')[0] + '] '
          ret_string += this_timekey + '</a>' + "\n"
      
        ret_string += '</div></div>' + "\n"
    
    ret_string += '</div></div>' + "\n"

  ret_string += '</div>'
  return ret_string

def render_status(GLOBAL_RAW_DATA_DICT, max_memory_in_bytes, MEMORY_BYTE_RATIO, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, GLOBAL_KEY_DEQUEUE, app, GLOBAL_STAT_VARS):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_RAW_DATA_DICT, 'status'))
  
  process = psutil.Process(os.getpid())
  
  memory_in_bytes = process.memory_info().rss
  body_data = '<div>'
  body_data += '<a class="nav-link" href="/status_memory_pics">memory status</a>'
  body_data += '<a class="nav-link" href="/status_log">Receive log</a>'
  body_data += 'used memory: '
  
  body_data += get_human_readable_bytes(memory_in_bytes) + '<br />'
  
  body_data += 'maximum memory: ' + get_human_readable_bytes(max_memory_in_bytes) + '<br />'
  body_data += 'byte ratio: ' + str(MEMORY_BYTE_RATIO) + '<br />'

  body_data += 'total_problem_count: ' + str(GLOBAL_STAT_VARS['total_problem_count']) + '<br />'
  body_data += 'total_iteration_count: ' + str(GLOBAL_STAT_VARS['total_iteration_count']) + '<br />'
  body_data += 'total_bytes_received: ' + get_human_readable_bytes(GLOBAL_STAT_VARS['total_bytes_received']) + '<br />'
  body_data += 'last_deleted: ' + GLOBAL_STAT_VARS['last_deleted'] + '<br />'
  body_data += 'last_deleted_date: ' + GLOBAL_STAT_VARS['last_deleted_date'] + '<br />'
  
  body_data += '<br />size per problem:<br />' + "\n"
  body_data += '<table class="table table-sm table-bordered"><thead><tr><th>key</th><th>raw xml size</th><th>estimated ram size</th><th>delete</th><th>download</th></tr></thead><tbody>' + "\n"
  
  for tmp_key in GLOBAL_DATA_SIZE_DICT:
    body_data += '<tr><td>' + tmp_key + '</td><td>' + get_human_readable_bytes(GLOBAL_DATA_SIZE_DICT[tmp_key]) + '</td><td>' + get_human_readable_bytes(MEMORY_BYTE_RATIO*GLOBAL_DATA_SIZE_DICT[tmp_key]) + '</td>'
    if tmp_key in GLOBAL_RAW_DATA_DICT:
      body_data += '<td><a href="/status?delete=' + html.escape(tmp_key) + '">x</a></td>'
    else:
      body_data += '<td></td>'
    
    body_data += '<td><a href="/download_xml/' + tmp_key + '">download</a></td></tr>' + "\n"
  
  body_data += '</tbody></table>' + "\n"

  body_data += '</div>'
  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  return retdata


def render_status_log(GLOBAL_RAW_DATA_DICT, MEMLOG_RECEIVELOG):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_RAW_DATA_DICT, 'status_log'))
  
  process = psutil.Process(os.getpid())
  
  memory_in_bytes = process.memory_info().rss
  body_data = '<div><table class="table table-sm table-bordered"><thead><tr><th>key</th><th>last_updated</th><th>xml_size_in_byte</th>'
  body_data += '<th>estimated_xml_size_in_byte</th><th>memory_before</th><th>memory_after</th><th>deleted_count</th></tr></thead><tbody>'

  for log_entry in reversed(MEMLOG_RECEIVELOG):
    body_data += '<tr>'
    body_data += '<td>' + log_entry['key'] + '</td>'
    body_data += '<td>' + log_entry['last_updated'] + '</td>'
    body_data += '<td>' + get_human_readable_bytes(log_entry['xml_size_in_byte']) + '</td>'
    body_data += '<td>' + get_human_readable_bytes(log_entry['estimated_xml_size_in_byte']) + '</td>'
    body_data += '<td>' + get_human_readable_bytes(log_entry['memory_before']) + '</td>'
    body_data += '<td>' + get_human_readable_bytes(log_entry['memory_after']) + '</td>'
    body_data += '<td>' + str(log_entry['deleted_count']) + '</td>'
    body_data += '</tr>'

  body_data += '</tbody></table></div>'
  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  return retdata

def render_status_memory_pics(GLOBAL_RAW_DATA_DICT, max_memory_in_bytes, MEMORY_BYTE_RATIO, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, GLOBAL_KEY_DEQUEUE, app):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_RAW_DATA_DICT, 'status_memory_pics'))
  
  seed = str(randomword(16))

  objgraph.show_refs([GLOBAL_RAW_DATA_DICT], filename=(seed + "-GLOBAL_RAW_DATA_DICT.png"))
  objgraph.show_refs([GLOBAL_DATA_SIZE_DICT], filename=(seed + "-GLOBAL_DATA_SIZE_DICT.png"))
  objgraph.show_refs([GLOBAL_UPDATED_DICT], filename=(seed + "-GLOBAL_UPDATED_DICT.png"))
  objgraph.show_refs([GLOBAL_KEY_DEQUEUE], filename=(seed + "-GLOBAL_KEY_DEQUEUE.png"))
  objgraph.show_refs([app], filename=(seed + "-app.png"))
  
  body_data = 'memory graphs for streamviz: <br />'
  
  data_uri = base64.b64encode(open(seed + "-GLOBAL_RAW_DATA_DICT.png", 'rb').read()).decode('utf-8')
  body_data += '<img src="data:image/png;base64,{0}">'.format(data_uri) + '<br />'

  data_uri = base64.b64encode(open(seed + "-GLOBAL_DATA_SIZE_DICT.png", 'rb').read()).decode('utf-8')
  body_data += '<img src="data:image/png;base64,{0}">'.format(data_uri) + '<br />'

  data_uri = base64.b64encode(open(seed + "-GLOBAL_UPDATED_DICT.png", 'rb').read()).decode('utf-8')
  body_data += '<img src="data:image/png;base64,{0}">'.format(data_uri) + '<br />'

  data_uri = base64.b64encode(open(seed + "-GLOBAL_KEY_DEQUEUE.png", 'rb').read()).decode('utf-8')
  body_data += '<img src="data:image/png;base64,{0}">'.format(data_uri) + '<br />'

  data_uri = base64.b64encode(open(seed + "-app.png", 'rb').read()).decode('utf-8')
  body_data += '<img src="data:image/png;base64,{0}">'.format(data_uri) + '<br />'

  os.remove(seed + "-GLOBAL_RAW_DATA_DICT.png")
  os.remove(seed + "-GLOBAL_DATA_SIZE_DICT.png")
  os.remove(seed + "-GLOBAL_UPDATED_DICT.png")
  os.remove(seed + "-GLOBAL_KEY_DEQUEUE.png")
  os.remove(seed + "-app.png")

  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  return retdata

#render main page
def render_index(GLOBAL_METADATA_DICT, request):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(0, 'index'))
  
  TABLE_DATA = []
  
  columns = ['host', 'status', 'problem', 'started', 'updated', 'objective', 'iterations']
  
  restricted_conditions = {}

  # get all conditions, e.g. if host should be eam080, then
  # there will be a get variable restrict_host
  for key in request.args:
    if key.find("restrict_") != -1:
      restricted_conditions[key[9:]] = request.args[key] 

  for key in GLOBAL_METADATA_DICT:
    this_table_data = GLOBAL_METADATA_DICT[key]

    add_this_element = True # assume we can add this element
    
    for key in restricted_conditions:
      if this_table_data[key] != restricted_conditions[key]:
        add_this_element = False
        break # get out of this loop, already excluded by condition
    
    if add_this_element:
      TABLE_DATA.append(this_table_data)
  
  #sort by column request.args['sort']
  if len(request.args.getlist("sort")) > 0:
    sorted_with = request.args.getlist('sort')[0]
    if sorted_with == 'iterations':
      TABLE_DATA = sorted(TABLE_DATA, key=lambda x: int(x[sorted_with]),  reverse=(request.args.get("reverse", 'False') == 'True'))
    else:
      TABLE_DATA = sorted(TABLE_DATA, key=lambda x: x[sorted_with],  reverse=(request.args.get("reverse", 'False') == 'True'))
  else:
    TABLE_DATA = sorted(TABLE_DATA, key=lambda x: x["updated"],  reverse=True)
  
  body_data  = '<table class="table table-sm" id="table_overview">'
  body_data += ' <thead><td></td>'
  for s in columns:
    #remove this restriction by clicking on the header of the table
    if 'restrict_' + s in request.args:
      #keep all other restrictions!
      new_args = request.args.copy()
      del new_args['restrict_' + s]
      body_data += '<td><a href="?' + urllib.parse.urlencode(new_args) + '">' + s + '</a></td>'
    else:
      body_data += '<td>' + s + '</td>'
  body_data += '</thead><tbody><tr><td>sort:</td>'
  
  for s in columns:
    body_data += '<td>'
    
    new_args = request.args.copy()
    
    is_reverse = (request.args.get("reverse", 'False') == 'True')
    
    if 'reverse' in new_args:
      del new_args['reverse']
    
    # if sort ist set, render a button to sort the other way or to cancel sorting.
    if request.args.get("sort") == s:
      #invert the arguments and add the opposite sorting button
      if is_reverse:
        new_args['reverse'] = 'False'
        body_data += '<a href="?' + urllib.parse.urlencode(new_args) + '">asc</a>|'
      else:
        new_args['reverse'] = 'True'
        body_data += '<a href="?' + urllib.parse.urlencode(new_args) + '">desc</a>|'
        
      #add a button to cancel sorting
      if 'sort' in new_args:
        del new_args['sort']
      if 'reverse' in new_args:
        del new_args['reverse']
      body_data += '<a href="?' + urllib.parse.urlencode(new_args) + '">x</a>'
      
    # if sort ist not set, render a button to sort ascending and descending
    else:
      new_args['sort'] = s
      body_data += '<a href="?' + urllib.parse.urlencode(new_args) + '">asc</a>|'
      new_args['reverse'] = 'True'
      body_data += '<a href="?' + urllib.parse.urlencode(new_args) + '">desc</a>'
    
    body_data += '</td>'
  
  body_data += '</tr>'
  
  for dataset in TABLE_DATA:
    body_data += '<tr>'
    body_data += '<td><a href="/view/' + dataset['key'] + '">view</a></td>'
    for key in columns: # take columns and NOT dataset because of sorting
      if not 'restrict_' + key in request.args and key in ['host', 'status', 'problem']:
        # we do a little value polishing
        value = str(dataset[key]) if key != 'host' else dataset[key][:dataset[key].find('.')] 
        #keep all other restrictions!
        new_args = request.args.copy()
        new_args['restrict_' + key] = str(dataset[key])
        body_data += '<td><a href="?' + urllib.parse.urlencode(new_args) + '">' + value + '</a></td>'
      else:
        body_data += '<td>' + str(dataset[key]) + '</td>'
    body_data += '</tr>\n'
  
  body_data += '</tbody></table>'

  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  return retdata

#creates a human readable time format
def get_dd_hh_mm_ss_fromsecs(td):
  
    ret_string = ""
    
    if td.days > 0:
      ret_string += str(td.days) + "d "
    
    td_hours = td.seconds//3600
    if td_hours > 0:
      ret_string += str(td_hours) + "h "
    
    td_mins = (td.seconds//60)%60
    if td_mins > 0:
      ret_string += str(td_mins) + "m "
      
    td_secs = td.seconds%60
    
    ret_string += str(td_secs) + "s"
    
    return ret_string

#render view of one simulation
def render_view(GLOBAL_RAW_DATA_DICT, GLOBAL_USE_PARAVIEW, key, client_ip):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_RAW_DATA_DICT, 'view'))
  
  if not key in GLOBAL_RAW_DATA_DICT:
    retdata = retdata.replace(settings['html_template']['key_content'], "simulation not found!")
    return retdata
  
  xml = etree.fromstring(GLOBAL_RAW_DATA_DICT[key])
  
  settings_data = '<div class="col-sm-4 table-borderless" id="settings">'
  settings_data += '<div id="simulation_id">' + key + '</div>'
  settings_data += '<ul class="list-group"><li class="list-group-item  border border-white">'
  settings_data += '    <ul class="list-group">'

  # status information  
  settings_data += '        <li><h5>global information</h5></li>'
  settings_data += '        <li class="list-group-item  border border-white">status: ' + xml.xpath('//cfsInfo/@status')[0] + '</li>'
  if xml.xpath('//cfsInfo/@status')[0] == 'finished':
    wall_in_secs = float(str(xml.xpath('//timer/@wall')[0]))
    wall_td = datetime.timedelta(0, wall_in_secs)
    cpu_in_secs = float(str(xml.xpath('//timer/@cpu')[0]))
    cpu_td = datetime.timedelta(0, cpu_in_secs)
    settings_data += '        <li class="list-group-item border border-white">time: wall=' + get_dd_hh_mm_ss_fromsecs(wall_td) + ' cpu=' + get_dd_hh_mm_ss_fromsecs(cpu_td) + '</li>'
    # not all systems have memory
    mp = xml.xpath('//memory/@peak')
    if len(mp) > 0:
      settings_data += '        <li class="list-group-item border border-white">peak memory: ' + str(mp[0]) + ' MB</li>'  
  settings_data += '    </ul></li>'
  
  # when streamviz.py could not import paraview(.vtk) we switch catalyst off
  if GLOBAL_USE_PARAVIEW:
    settings_data += '    <li class="list-group-item border border-white">'
    settings_data += '    send data to catalyst (v5.5*):<br/>'
    settings_data += '    ip: <input type="text" id="catalyst_ip" value="' + client_ip + '"/> port: <input type="text" id="catalyst_port" value="22222"/><br/>'
    settings_data += '    <button class="btn btn-primary" id="catalyst_send_button">send data</button>'
    settings_data += '    auto update: <input type="checkbox" id="catalyst_send_auto_update"/>'
    settings_data += '    </li>'
  
  
  # plot settings with log scale, regions results and optimization iter
  settings_data += '    <li class="list-group-item border border-white"><h5>scalar results</h5>'
  settings_data += '    <table class="table table-sm" id="iteration">'

  settings_data += '        <thead><td>y1</td><td>y2</td></thead>'
  settings_data += '        <tbody>'

  # add the log scale selectors
  settings_data += '         <tr><td><input type="checkbox" id="logscale_y1" />'
  settings_data += '             <td><input type="checkbox" id="logscale_y2" />'
  settings_data += '             <td>' + html.escape('<logarithmic scaling>') + '</td></tr>'
  
  # process region results (mechTotalEnergy, ...)
  for result in xml.xpath('//calculation/process/sequence/result'):
    dat = result.xpath('item[last()]')[0]
    value = dat.get('value')     
    if value is not None: # not every result has a value
      name = result.get('data') + ':' + result.get('location')
      settings_data += '        <tr><td><input class="form-control" type="checkbox" name="result_selector_y1" value="' + name + '" /></td>'
      settings_data +=             '<td><input class="form-control" type="checkbox" name="result_selector_y2" value="' + name + '" /></td>'
      settings_data += '            <td id="td_seqres_container_' + name + '">' + name + ': ' + dat.get('value') + ' ' + dat.get('unit') + '</td></tr>'

  # process optimzation iterations
  last_iteration_elem = xml.xpath('//process/iteration[last()]')
  if len(last_iteration_elem) > 0:
    for value in last_iteration_elem[0].items():
      if value[0] != 'name':
        settings_data +=         '<tr><td><input class="form-control" type="checkbox" name="iteration_selector_y1" value="' + value[0] + '" /></td>'
        settings_data +=             '<td><input class="form-control" type="checkbox" name="iteration_selector_y2" value="' + value[0] + '" /></td>'
        settings_data +=             '<td id="td_seqit_container_' + value[0] + '">' + value[0] + '</td></tr>'
  settings_data += '    </tbody></table></li>'
  
  
  # offer plot for element 
  close_result_list = False
  
  has_grid = len(xml.xpath('//grid/nodeList/node')) > 0
  
  if not has_grid:
    settings_data += '    <li class="list-group-item border border-white"><h5>2D results</h5> not available, no mesh sent</li>'
  
  if int(xml.xpath('//grids/grid/@dimensions')[0]) == 2 and has_grid:
    close_result_list = True
    # we have a 2-Dimensional grid here. Therefore we offer to render it
    settings_data += '    <li class="list-group-item  border border-white"><h5>2D results</h5> <button id="disable_results">hide 2d results</button>'
    settings_data += '    <table class="table table-sm" id="result">'
    settings_data += '        <tbody>'
    
    for result in xml.xpath('//results/result'):
      # currently we only display 1D data as color image - could be extendend
      # same is for nodal data, which would be needed to be interpreted
      if int(result.attrib['dofs']) == 1:
        name = result.attrib['name']
        
        settings_data += '        <tr><td><input class="form-control" type="radio" name="result_selector_view" value="' + name + '" /></td>'
        settings_data +=             '<td>' + name + '</td></tr>'
  
  if len(xml.xpath("//eigenFrequency/result/wave_vector"))>0:
    if not close_result_list:
      settings_data += '    <li class="list-group-item border border-white">results:'
      settings_data += '    <table class="table table-sm" id="result">'
      settings_data += '        <thead><td>view</td><td></td></thead><tbody>'
    
    close_result_list = True
    settings_data += '        <tr><td><input class="form-control" type="checkbox" name="result_selector_view_bloch" value="true" /></td>'
    settings_data +=             '<td>[bloch] eigenfrequency</td></tr>'
  
  if close_result_list:
    settings_data += '    </tbody></table></li>'
    
  settings_data += '</ul></div>'

  body_data = '<div class="col-sm-8" id="content">' + key + '<br>'
  body_data += '<div id="iteration_plot"><div id="iteration_num">-1</div>'
  if len(xml.xpath('//objective/@type')) > 0:
    body_data += '<div id ="objective">' + xml.xpath('//objective/@type')[0] + '</div>'
  body_data += '<div id ="status">' + xml.xpath('//cfsInfo/@status')[0] + '</div>'
  body_data += '</div></div>'
  
  retdata = retdata.replace(settings['html_template']['key_content'], settings_data + body_data)
  return retdata
