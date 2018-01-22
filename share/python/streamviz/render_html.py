from flask.globals import request

from strings import *

#The html file with marker; see strings.py for all marker definitions 
html_raw_data = "html file not loaded"
with open('template_html/index.html', 'r') as myfile:
  html_raw_data = myfile.read().replace('\n', '')

#reder the selection menu
def render_menu(GLOBAL_DATA_DICT, current_site):
  ret_string  = '<div class="btn-group">'
  ret_string += '<a class="nav-link" href="/">Overview</a>'
  
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
    ret_string += '<button class="btn btn-secondary dropdown-toggle" type="button" id="dropdownMenuButton_' + this_host + "\n"
    ret_string += '" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">'
    ret_string += this_host + '</button> <div class="dropdown-menu" aria-labelledby="dropdownMenuButton_' + this_host + '">' + "\n"
    
    for this_problem in hosts[this_host]:
      ret_string += '<div class="dropdown-submenu dropright">' + "\n"
      ret_string += '<button class="btn btn-secondary dropdown-toggle" type="button" id="dropdownMenuButton_' + this_host + '_' + this_problem + "\n"
      ret_string += '" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">'
      ret_string += this_problem + '</button> <div class="dropdown-menu" aria-labelledby="dropdownMenuButton_' + this_host + '_' + this_problem + '">' + "\n"
      
      for this_timekey in hosts[this_host][this_problem]:
        ret_string += '<a class="dropdown-item" href="/view/' + hosts[this_host][this_problem][this_timekey] + '">' + this_timekey + '</a>' + "\n"
      
      ret_string += '</div></div>' + "\n"
    
    ret_string += '</div></div>' + "\n"

  ret_string += '</div>'
  return ret_string

#render main page
def render_index(GLOBAL_DATA_DICT, request):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_DATA_DICT, 'index'))
  
  body_data = 'body'
  
  
  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  retdata = retdata.replace(settings['html_template']['key_settings'], '')
  return retdata

#render view of one simulation
def render_view(GLOBAL_DATA_DICT, key):
  retdata = html_raw_data
  retdata = retdata.replace(settings['html_template']['key_menu'], render_menu(GLOBAL_DATA_DICT, 'view'))
  
  if not key in GLOBAL_DATA_DICT:
    retdata = retdata.replace(settings['html_template']['key_content'], "simulation not found!")
    return retdata
  
  settings_data = '<h4>Settings</h4>'
  settings_data = '<div id="simulation_id">' + key + '</div>'
  settings_data += '<ul class="list-group"><li class="list-group-item">Simulation Data:'
  settings_data += '    <ul class="list-group">'
  settings_data += '        <li class="list-group-item">1</li>'
  settings_data += '        <li class="list-group-item">2</li>'
  settings_data += '    </ul></li>'
  settings_data += '    <li class="list-group-item">Results:'
  settings_data += '    <table class="table table-sm" id="result">'
  settings_data += '        <thead><td>view</td><td></td></thead><tbody>'
  
  for result in GLOBAL_DATA_DICT[key].xpath('//streaming/process/sequence/result'):
    name = result.xpath('@data')[0]
    
    settings_data += '        <tr><td><input class="form-control" type="checkbox" name="result_selector_view" value="' + name + '" /></td>'
    settings_data +=             '<td>' + name + '</td></tr>'

  settings_data += '    </tbody></table></li>'
  settings_data += '    <li class="list-group-item">Iterations and scalar Results:'
  settings_data += '    <table class="table table-sm" id="iteration">'
  settings_data += '        <thead><td>x</td><td>y1</td><td>y2</td></thead><tbody>'
  
  for result in GLOBAL_DATA_DICT[key].xpath('//calculation/process/sequence/result'):
    name = result.xpath('@data')[0]
    
    settings_data += '        <tr><td></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="result_selector_y1" value="' + name + '" /></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="result_selector_y2" value="' + name + '" /></td>'
    settings_data += '            <td>' + name + ': ' + result.xpath('item[last()]')[0].xpath('@value')[0]
    settings_data += ' ' + result.xpath('item[last()]')[0].xpath('@unit')[0] + '</td></tr>'

  for value in GLOBAL_DATA_DICT[key].xpath('//process/iteration[last()]')[0].items():
    settings_data += '        <tr><td><input class="form-control" type="radio" name="iteration_selector_x" value="' + value[0] + '" /></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="iteration_selector_y1" value="' + value[0] + '" /></td>'
    settings_data +=             '<td><input class="form-control" type="checkbox" name="iteration_selector_y2" value="' + value[0] + '" /></td>'
    settings_data +=             '<td>' + value[0] + '</td></tr>'

  settings_data += '    </tbody></table></li>'
  settings_data += '<!-- <li class="list-group-item"><button class="btn btn-primary" onclick="trigger_update()" >Update Data</button><br>'
  settings_data += 'auto refresh: <input type="checkbox" name="enable_autorefresh" checked="checked" /></li> -->'
  settings_data += '</ul>'
  
  retdata = retdata.replace(settings['html_template']['key_settings'], settings_data)
  
  body_data  = '<h5 class="card-title">View current simulation "' + key + '":</h5>'
  body_data += '<div id="iteration_plot"><div id="iteration_num">-1</div></div>'
  
  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  return retdata
