from flask.globals import request

from strings import *

#The html file with marker; see strings.py for all marker definitions 
html_raw_data = "html file not loaded"
with open('template_html/index.html', 'r') as myfile:
  html_raw_data = myfile.read().replace('\n', '')

#reder the selection menu
def render_menu(GLOBAL_DATA_DICT, current_site):
  ret_string = '<ul>'
  
  ret_string += '<li class="'
  if current_site == 'index':
    ret_string += 'menu_active'
  ret_string += '"><a href="/">Overview</a></li>'
  
  ret_string += '<li class="'
  if current_site == 'view':
    ret_string += 'menu_active'
  ret_string += '">View<ul>'
  for key in GLOBAL_DATA_DICT:
    ret_string += '<li>'
    ret_string += '<a href="' + settings["api"]["view_url"] + '/' + key + '">['
    ret_string += GLOBAL_DATA_DICT[key].xpath('//cfsStreaming/cfsInfo/@status')[0] + '] ' + key + '</a>'
    ret_string += '</li>'

  ret_string += '</ul></li>'

  ret_string += '</ul>'
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
  settings_data += '<ul><li>Simulation Data:'
  settings_data += '    <ul>'
  settings_data += '        <li>1</li>'
  settings_data += '        <li>2</li>'
  settings_data += '    </ul></li>'
  settings_data += '    <li>Results:'
  settings_data += '    <table id="result">'
  settings_data += '        <thead><td>view</td><td></td></thead><tbody>'
  
  for result in GLOBAL_DATA_DICT[key].xpath('//streaming/process/sequence/result'):
    name = result.xpath('@data')[0]
    
    settings_data += '        <tr><td><input type="checkbox" name="result_selector_view" value="' + name + '" /></td>'
    settings_data +=             '<td>' + name + '</td></tr>'

  settings_data += '    </tbody></table></li>'
  settings_data += '    <li>Iterations and scalar Results:'
  settings_data += '    <table id="iteration">'
  settings_data += '        <thead><td>x</td><td>y1</td><td>y2</td></thead><tbody>'
  
  for result in GLOBAL_DATA_DICT[key].xpath('//calculation/process/sequence/result'):
    name = result.xpath('@data')[0]
    
    settings_data += '        <tr><td></td>'
    settings_data +=             '<td><input type="checkbox" name="result_selector_y1" value="' + name + '" /></td>'
    settings_data +=             '<td><input type="checkbox" name="result_selector_y2" value="' + name + '" /></td>'
    settings_data += '            <td>' + name + ': ' + result.xpath('item[last()]')[0].xpath('@value')[0]
    settings_data += ' ' + result.xpath('item[last()]')[0].xpath('@unit')[0] + '</td></tr>'

  for value in GLOBAL_DATA_DICT[key].xpath('//process/iteration[last()]')[0].items():
    settings_data += '        <tr><td><input type="radio" name="iteration_selector_x" value="' + value[0] + '" /></td>'
    settings_data +=             '<td><input type="checkbox" name="iteration_selector_y1" value="' + value[0] + '" /></td>'
    settings_data +=             '<td><input type="checkbox" name="iteration_selector_y2" value="' + value[0] + '" /></td>'
    settings_data +=             '<td>' + value[0] + '</td></tr>'

  settings_data += '    </tbody></table></li>'
  settings_data += '</ul>'
  
  settings_data += '<button onclick="trigger_update()" >Update Data</button><br>'
  settings_data += 'auto refresh: <input type="checkbox" name="enable_autorefresh" checked="checked" />'
  
  retdata = retdata.replace(settings['html_template']['key_settings'], settings_data)
  
  body_data = '<h1>View current simulation "' + key + '"</h1><div id="iteration_plot"><div id="iteration_num">-1</div></div>'
  
  retdata = retdata.replace(settings['html_template']['key_content'], body_data)
  return retdata
