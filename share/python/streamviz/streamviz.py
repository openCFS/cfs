#!/usr/bin/env python3

# webserver:
from flask import Flask
from flask import request

# xml parser
from lxml import etree

# settings
from strings import *

# for forced garbage collection
import gc

# html render module
import render_html

# plot module
import api_plot

# for global update dict:
import datetime

# for keeping track of all keys
from collections import deque

# for geting ram usage
import os
import psutil

import sys
import argparse

# for sending data to catalyst
GLOBAL_USE_PARAVIEW = None
try:
  import paraview
  import paraview.vtk
  import send_data
  GLOBAL_USE_PARAVIEW = True
except ImportError as ie:
  print('Tolerable import error (no paraview catalyst service):', ie)
  GLOBAL_USE_PARAVIEW = False

#for threading.lock
import threading

from flask.wrappers import Request

app = Flask(__name__)

MEMLOG_FILENAME = "streamviz_memlog_" + str(datetime.datetime.now()) + ".log"
MEMLOG_LOCK = threading.Lock()
#with open(MEMLOG_FILENAME, "a") as myfile:
#  myfile.write("#key last_updated xml_size_in_byte(6) estimated_xml_size_in_byte(7) memory_before memory_after deleted_count\n")

MEMLOG_RECEIVELOG = []

# global data dict always contains the latest xml file (in parsed version)
GLOBAL_RAW_DATA_DICT = {}
GLOBAL_DATA_SIZE_DICT = {}

# contains the last update
GLOBAL_UPDATED_DICT = {}
GLOBAL_OBJECTIVE_DICT = {}

GLOBAL_METADATA_DICT = {} # will be copied into array before sorting

GLOBAL_KEY_DEQUEUE = deque([])

# ajax functions will wait for an even aka new data
# in order to do that they will put and event into this dictionary: key => Array(Event)
# on new cfs_recieve, the events will be set() and deleted from the dict, releasing all
# events and firing the http return.
# current timeout for events should be set to 9 seconds server side and 10 seconds 
# client side within javascript. This will ensure performance and prevent useless
# data transfer and html updating
UPDATE_EVENTS = {}

# This header should be set by the reverse proxy.
# Make sure to have the real IP here
# Don't use the X-Forward-For from the request
PROXY_REAL_IP_HEADER = "X-Real-IP"

# ration <bytes used reported by os> divided by <sum of sizes of saved xml strings>
MEMORY_BYTE_RATIO = -1

# status vars
GLOBAL_STAT_VARS = {}
GLOBAL_STAT_VARS['total_problem_count'] = 0
GLOBAL_STAT_VARS['total_iteration_count'] = 0
GLOBAL_STAT_VARS['last_deleted'] = ""
GLOBAL_STAT_VARS['last_deleted_date'] = ""
GLOBAL_STAT_VARS['total_bytes_received'] = 0

@app.route('/', methods = ['GET'])
def index():
  if request.method == 'GET':
    return render_html.render_index(GLOBAL_METADATA_DICT, request)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route('/status', methods = ['GET'])
def status():
  key_to_delete = request.args.get('delete', '')
  if key_to_delete != '':
    if key_to_delete in GLOBAL_RAW_DATA_DICT:
      oldest_data = GLOBAL_RAW_DATA_DICT.pop(key_to_delete) # delete oldest xml keys
      del oldest_data
    if key_to_delete in GLOBAL_UPDATED_DICT:
      oldest_updated = GLOBAL_UPDATED_DICT.pop(key_to_delete) 
      del oldest_updated
    if key_to_delete in UPDATE_EVENTS:
      oldest_update_event = UPDATE_EVENTS.pop(key_to_delete)
      oldest_update_event.set()
      del oldest_update_event
    if key_to_delete in GLOBAL_OBJECTIVE_DICT:
      oldest_updated = GLOBAL_OBJECTIVE_DICT.pop(key_to_delete) 
      del oldest_updated
    if key_to_delete in GLOBAL_METADATA_DICT:
      oldest_entry = GLOBAL_METADATA_DICT.pop(key_to_delete)
      del oldest_entry

    send_data.delete_coprocessor(key_to_delete)
    gc.collect()
  
  return render_html.render_status(GLOBAL_RAW_DATA_DICT, MAX_MEMORY, MEMORY_BYTE_RATIO, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, GLOBAL_KEY_DEQUEUE, app, GLOBAL_STAT_VARS)

@app.route('/status_log', methods = ['GET'])
def status_log():
    return render_html.render_status_log(GLOBAL_RAW_DATA_DICT, MEMLOG_RECEIVELOG)

@app.route('/status_memory_pics', methods = ['GET'])
def status_memory_pics():
    return render_html.render_status_memory_pics(GLOBAL_RAW_DATA_DICT, MAX_MEMORY, MEMORY_BYTE_RATIO, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, GLOBAL_KEY_DEQUEUE, app)

@app.route(settings["api"]["values"] + '/<path:key>', methods = ['GET', 'POST'])
def values(key):
  if request.method == 'GET':
    return api_plot.get_values(GLOBAL_RAW_DATA_DICT, UPDATE_EVENTS, key, int(request.args.get('iteration_num')))
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route(settings["api"]["view_url"] + '/<path:key>', methods = ['GET', 'POST'])
def view(key):
  if request.method == 'GET':
    client_ip = request.remote_addr
    if PROXY_REAL_IP_HEADER in request.headers:
      client_ip = request.headers.get(PROXY_REAL_IP_HEADER, "127.0.0.1")
    return render_html.render_view(GLOBAL_RAW_DATA_DICT, GLOBAL_USE_PARAVIEW, key, client_ip)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route(settings["api"]["plot_url"] + '/<path:key>', methods = ['GET', 'POST'])
def plot(key):
  return api_plot.plot(key, UPDATE_EVENTS, GLOBAL_RAW_DATA_DICT, \
                       request.args.getlist('y1_it'), request.args.getlist('y2_it'), \
                       request.args.getlist('y1_res'), request.args.getlist('y2_res'), \
                       int(request.args.get('iteration_num')), \
                       (request.args.get('logscale_y1') == 'true'), (request.args.get('logscale_y2') == 'true'), \
                       request.args.getlist('view_results'), \
                       (request.args.get('view_result_bloch', 'false') == 'true'))


@app.route(settings["api"]["catalyst_send"] + '/<path:key>', methods = ['GET', 'POST'])
def send_data_func(key):
  send_data.send_data(key, etree.fromstring(GLOBAL_RAW_DATA_DICT[key]), request.args.get('ip'), request.args.get('port'))
  return ""

@app.route('/', methods = ['POST'])
@app.route(settings["api"]["recieve_url"], methods = ['GET', 'POST'])
def cfs_recieve_blank():
  return cfs_recieve("")

@app.route('/download_xml/<path:key>')
def download_xml(key):
  return GLOBAL_RAW_DATA_DICT[key]

@app.route(settings["api"]["recieve_url"] + '/', methods = ['GET', 'POST'])
@app.route(settings["api"]["recieve_url"] + '/<path:url_key>', methods = ['GET', 'POST'])
def cfs_recieve(url_key = ""):
  global MEMORY_BYTE_RATIO, GLOBAL_RAW_DATA_DICT, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, UPDATE_EVENTS, GLOBAL_STAT_VARS
  process = psutil.Process(os.getpid())
  
  memory_in_bytes = process.memory_info().rss
  
  total_xml_size = 0
  deleted_count = 0
  for tmp_key in GLOBAL_RAW_DATA_DICT:
    total_xml_size += GLOBAL_DATA_SIZE_DICT[tmp_key]

  GLOBAL_STAT_VARS['total_bytes_received'] += total_xml_size

  if memory_in_bytes > MAX_MEMORY:
    if MEMORY_BYTE_RATIO < 0:
      # we need to calculate the ratio:
      MEMORY_BYTE_RATIO = memory_in_bytes / total_xml_size

    xml_size_to_reduce = total_xml_size - (MAX_MEMORY / MEMORY_BYTE_RATIO) # size in bytes of xml

    while xml_size_to_reduce > 0:
      oldest_key = GLOBAL_KEY_DEQUEUE.pop() # get oldest xml key
      if oldest_key in GLOBAL_RAW_DATA_DICT:
        oldest_data = GLOBAL_RAW_DATA_DICT.pop(oldest_key) # delete oldest xml keys
        xml_size_to_reduce -= GLOBAL_DATA_SIZE_DICT.pop(oldest_key)
        del oldest_data
      if oldest_key in GLOBAL_UPDATED_DICT:
        oldest_updated = GLOBAL_UPDATED_DICT.pop(oldest_key) 
        del oldest_updated
      if oldest_key in UPDATE_EVENTS:
        oldest_update_event = UPDATE_EVENTS.pop(oldest_key)
        oldest_update_event.set()
        del oldest_update_event
      if oldest_key in GLOBAL_OBJECTIVE_DICT:
        oldest_updated = GLOBAL_OBJECTIVE_DICT.pop(oldest_key) 
        del oldest_updated
      if oldest_key in GLOBAL_METADATA_DICT:
        oldest_entry = GLOBAL_METADATA_DICT.pop(oldest_keys)
        del oldest_entry

      deleted_count += 1
      send_data.delete_coprocessor(oldest_key)

      GLOBAL_STAT_VARS['last_deleted'] = oldest_key
      GLOBAL_STAT_VARS['last_deleted_date'] = str(datetime.datetime.now())
      
    # collect garbage
    gc.collect()


  data = ""
  
  # read the data from input stream directly
  if request.method == 'POST':
    chunk_size = 4096
    
    while True:
      chunk = request.stream.read(chunk_size)
      if len(chunk) == 0:
        print("received data")
        break
      
      data += chunk.decode("utf-8")

    xml = etree.fromstring(data)
    
    key = xml.xpath('//environment/@host')[0] + '/' + xml.xpath('//progOpts/@problem')[0]
    
    if len(xml.xpath('//cfsInfo/header/@id')) > 0:
      if not url_key == '':
        key += url_key + '/' + xml.xpath('//cfsInfo/header/@id')[0]
    else:
      if not url_key == '':
        key += url_key

    key += '/' + xml.xpath('//header/environment/@started')[0][5:]

    if not key in GLOBAL_RAW_DATA_DICT: # only add key if not already in
      GLOBAL_KEY_DEQUEUE.appendleft(key)
      GLOBAL_STAT_VARS['total_problem_count'] += 1

    GLOBAL_RAW_DATA_DICT[key] = data
    GLOBAL_DATA_SIZE_DICT[key] = len(data)
    
    GLOBAL_UPDATED_DICT[key] = str(datetime.datetime.now())[5:-7]

    objective_string = "" # empty on default, if simulation then it will be "compliance: <value>"
    
    try:
      objective_name = str(xml.xpath("//header/objective/@name")[0])
      objective_string = objective_name + ": " + str(xml.xpath("(//process/iteration/@" + objective_name + ")[last()]")[0])
    except Exception as e:
      pass
    
    GLOBAL_OBJECTIVE_DICT[key] = objective_string
    
    metadata_dict = {}
    
    this_host = key[:key.index('/')]    
    project_time_rest = key[key.index('/')+1:]
    this_problem = project_time_rest[:project_time_rest.index('/')]
    this_started = project_time_rest[project_time_rest.index('/')+1:]
    
    this_table_data = {}
    this_table_data['key'] = key
    this_table_data['host'] = this_host
    this_table_data['status'] = xml.xpath('//cfsInfo/@status')[0]
    this_table_data['problem'] = this_problem
    this_table_data['started'] = this_started
    this_table_data['updated'] = GLOBAL_UPDATED_DICT[key]
    this_table_data['objective'] = GLOBAL_OBJECTIVE_DICT[key]
    
    iteration_num_array = xml.xpath('//process/iteration[last()]/@number')
    
    if len(iteration_num_array) > 0:
      this_table_data['iterations'] = int(iteration_num_array[0])
    else:
      this_table_data['iterations'] = -1
    
    GLOBAL_METADATA_DICT[key] = this_table_data
    
    # set the event to trigger sending of the new xml data
    # make sure to set the event AFTER putting the new xml
    # into GLOBAL_RAW_DATA_DICT
    if key in UPDATE_EVENTS:
      UPDATE_EVENTS[key].set()

    memory_in_bytes_after_adding = psutil.Process(os.getpid()).memory_info().rss

    with MEMLOG_LOCK:
      #with open(MEMLOG_FILENAME, "a") as myfile:
      #  myfile.write(key + " " + GLOBAL_UPDATED_DICT[key] + " " + str(total_xml_size) + " " + str(MEMORY_BYTE_RATIO*total_xml_size) + " " + str(memory_in_bytes) + " " + str(memory_in_bytes_after_adding) + " " + str(deleted_count) + "\n")

      log_entry = {}
      log_entry['key'] = key
      log_entry['last_updated'] = GLOBAL_UPDATED_DICT[key]
      log_entry['xml_size_in_byte'] = total_xml_size
      log_entry['estimated_xml_size_in_byte'] = MEMORY_BYTE_RATIO*total_xml_size
      log_entry['memory_before'] = memory_in_bytes
      log_entry['memory_after'] = memory_in_bytes_after_adding
      log_entry['deleted_count'] = deleted_count
      
      MEMLOG_RECEIVELOG.append(log_entry)

    GLOBAL_STAT_VARS['total_iteration_count'] += 1
    return 'data recieved!\n'
  else:
    return 'expected a POST'



if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='streamviz web server')
  parser.add_argument('host', nargs='?', help='host to assign the server', default='127.0.0.1')
  parser.add_argument("-port", help='port to assign the server', type=int, default=5000)
  parser.add_argument("-mem", help='maximal memory in giga bytes', type=float, default=2.0)
  args = parser.parse_args()
  # in bytes
  MAX_MEMORY = args.mem*1024**3 # = 2 GByte
    
  # need multithreading to server multiple clients.
  # Otherwise pushing xml while ajaxing it with the event
  # synchronization will not be possible
  app.run(threaded=True, port=args.port, host=args.host)
