from flask.globals import request

from strings import *

import urllib
from audioop import reverse

import datetime

#The html file with marker; see strings.py for all marker definitions 
html_raw_data = "html file not loaded"
with open('template_html/index.html', 'r') as myfile:
  html_raw_data = myfile.read().replace('\n', '')


#reder the selection menu
def render_menu(GLOBAL_DATA_DICT, current_site):
  ret_string  = '<div class="btn-group">'
  ret_string += '<a class="nav-link" href="/">overview</a>'
  
  hosts = {}
  for key in GLOBAL_DATA_DICT:
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
        ret_string += '<a class="dropdown-item" href="/view/' + hosts[this_host][this_problem][this_timekey] + '">'
        ret_string += '[' + GLOBAL_DATA_DICT[hosts[this_host][this_problem][this_timekey]].xpath('//cfsInfo/@status')[0] + '] '
        ret_string += this_problem + '/' + this_timekey + '</a>' + "\n"
        
      else:
        ret_string += '<div class="dropdown-submenu dropright">' + "\n"
        ret_string += '<button class="btn btn-secondary dropdown-toggle" type="button" id="dropdownMenuButton_' + this_host + '_' + this_problem
        ret_string += '" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">'
        ret_string += this_problem + '</button> <div class="dropdown-menu" aria-labelledby="dropdownMenuButton_' + this_host + '_' + this_problem + '">' + "\n"
      
        for this_timekey in hosts[this_host][this_problem]:
          ret_string += '<a class="dropdown-item" href="/view/' + hosts[this_host][this_problem][this_timekey] + '">'
          ret_string += '[' + GLOBAL_DATA_DICT[hosts[this_host][this_problem][this_timekey]].xpath('//cfsInfo/@status')[0] + '] '
          ret_string += this_timekey + '</a>' + "\n"
      
        ret_string += '</div></div>' + "\n"
    
    ret_string += '</div></div>' + "\n"

  ret_string += '</div>'
  return ret_string

#render main page
def render_index(GLOBAL_DATA_DICT, GLOBAL_UPDATED_DICT, request):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_DATA_DICT, 'index'))
  
  TABLE_DATA = []
  
  columns = ['host', 'status', 'problem', 'started', 'updated', 'iterations']
  
  restricted_conditions = {}

  # get all conditions, e.g. if host should be eam080, then
  # there will be a get variable restrict_host
  for key in request.args:
    if key.find("restrict_") != -1:
      restricted_conditions[key[9:]] = request.args[key] 

  for key in GLOBAL_DATA_DICT:
    this_table_data = {}
    
    xml = GLOBAL_DATA_DICT[key]
    
    this_host = key[:key.index('/')]    
    project_time_rest = key[key.index('/')+1:]
    this_problem = project_time_rest[:project_time_rest.index('/')]
    this_started = project_time_rest[project_time_rest.index('/')+1:]
    
    this_table_data['key'] = key
    this_table_data['host'] = this_host
    this_table_data['status'] = xml.xpath('//cfsInfo/@status')[0]
    this_table_data['problem'] = this_problem
    this_table_data['started'] = this_started
    this_table_data['updated'] = GLOBAL_UPDATED_DICT[key]
    this_table_data['iterations'] = int(xml.xpath('//process/iteration[last()]/@number')[0])
    
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
      if not 'restrict_' + key in request.args:
        #keep all other restrictions!
        new_args = request.args.copy()
        new_args['restrict_' + key] = str(dataset[key]) 
        body_data += '<td><a href="?' + urllib.parse.urlencode(new_args) + '">' + str(dataset[key]) + '</a></td>'
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
def render_view(GLOBAL_DATA_DICT, key, client_ip):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_DATA_DICT, 'view'))
  
  if not key in GLOBAL_DATA_DICT:
    retdata = retdata.replace(settings['html_template']['key_content'], "simulation not found!")
    return retdata
  
  xml = GLOBAL_DATA_DICT[key]
  
  settings_data = '<div class="col-sm-4" id="settings"><h4></h4>'
  settings_data += '<div id="simulation_id">' + key + '</div>'
  settings_data += '<ul class="list-group"><li class="list-group-item">'
  settings_data += '    <ul class="list-group">'
  
  settings_data += '        <li class="list-group-item">status: [' + xml.xpath('//cfsInfo/@status')[0] + ']</li>'
  if xml.xpath('//cfsInfo/@status')[0] == 'finished':
    wall_in_secs = float(str(xml.xpath('//timer/@wall')[0]))
    
    wall_td = datetime.timedelta(0, wall_in_secs)
    
    settings_data += '        <li class="list-group-item">wall: ' + get_dd_hh_mm_ss_fromsecs(wall_td) + '</li>'
    
    cpu_in_secs = float(str(xml.xpath('//timer/@cpu')[0]))
    
    cpu_td = datetime.timedelta(0, cpu_in_secs)
    
    settings_data += '        <li class="list-group-item">CPU time: ' + get_dd_hh_mm_ss_fromsecs(cpu_td) + '</li>'
    settings_data += '        <li class="list-group-item">peak memory: ' + xml.xpath('//memory/@peak')[0] + ' MB</li>'
  
  settings_data += '    </ul></li>'
  
  settings_data += '    <li class="list-group-item">'
  settings_data += '    send data to catalyst (v5.5*):<br/>'
  settings_data += '    ip: <input type="text" id="catalyst_ip" value="' + client_ip + '"/> port: <input type="text" id="catalyst_port" value="22222"/><br/>'
  settings_data += '    <button class="btn btn-primary" id="catalyst_send_button">send data</button>'
  settings_data += '    auto update: <input type="checkbox" id="catalyst_send_auto_update"/>'
  settings_data += '    </li>'
  
  if int(xml.xpath('//grids/grid/@dimensions')[0]) == 2:
    # we have a 2-Dimensional grid here. Therefore we offer to render it
    settings_data += '    <li class="list-group-item">results:'
    settings_data += '    <table class="table table-sm" id="result">'
    settings_data += '        <thead><td>view</td><td></td></thead><tbody>'
    
    for result in xml.xpath('//results/result'):
      if int(result.attrib['dofs']) == 1:
        # we can only display 1D data as color image!
        name = result.attrib['name']
        
        settings_data += '        <tr><td><input class="form-control" type="checkbox" name="result_selector_view" value="' + name + '" /></td>'
        settings_data +=             '<td>' + name + '</td></tr>'
  
    settings_data += '    </tbody></table></li>'
  
  
  settings_data += '    <li class="list-group-item">iterations and scalar results:'
  settings_data += '    <table class="table table-sm" id="iteration">'
  #hide the x selector for now
  settings_data += '        <thead><td style="display: none">x</td><td>y1</td><td>y2</td></thead><tbody>'
  
  for result in xml.xpath('//calculation/process/sequence/result'):
    name = result.xpath('@data')[0]
    
    #hide the x selector for now
    settings_data += '        <tr><td style="display: none"></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="result_selector_y1" value="' + name + '" /></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="result_selector_y2" value="' + name + '" /></td>'
    settings_data += '            <td id="td_seqres_container_' + name + '">' + name + ': ' + result.xpath('item[last()]')[0].xpath('@value')[0]
    settings_data += ' ' + result.xpath('item[last()]')[0].xpath('@unit')[0] + '</td></tr>'

  for value in xml.xpath('//process/iteration[last()]')[0].items():
    #hide the x selector for now
    if value[0] == 'number':
      # don't show the number in y selection:
      settings_data += '        <tr style="display:none">'
    else:
      settings_data += '        <tr>'
    
    settings_data += '<td style="display: none"><input class="form-control" type="radio" name="iteration_selector_x" value="' + value[0] + '" /></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="iteration_selector_y1" value="' + value[0] + '" /></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="iteration_selector_y2" value="' + value[0] + '" /></td>'
    settings_data +=             '<td id="td_seqit_container_' + value[0] + '">' + value[0] + '</td></tr>'

  settings_data += '    </tbody></table></li>'
  settings_data += '    <li class="list-group-item">use logscale:<br />'
  settings_data += '        y1: <input type="checkbox" id="logscale_y1" /><br />'
  settings_data += '        y2: <input type="checkbox" id="logscale_y2" />'
  settings_data += '    </li>'
  settings_data += '</ul></div>'

  body_data  = '<div class="col-sm-8" id="content"><h5 class="card-title">' + key + ':</h5>'
  body_data += '<div id="iteration_plot"><div id="iteration_num">-1</div>'
  if len(xml.xpath('//objective/@type')) > 0:
    body_data += '<div id ="objective">' + xml.xpath('//objective/@type')[0] + '</div>'
  body_data += '<div id ="status">' + xml.xpath('//cfsInfo/@status')[0] + '</div>'
  body_data += '</div></div>'
  
  retdata = retdata.replace(settings['html_template']['key_content'], settings_data + body_data)
  return retdata
