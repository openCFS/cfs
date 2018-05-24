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

# for sending data to catalyst
import send_data

#for threading.lock
import threading

from flask.wrappers import Request

app = Flask(__name__)

MEMLOG_FILENAME = "streamviz_memlog_" + str(datetime.datetime.now()) + ".log"
MEMLOG_LOCK = threading.Lock()
with open(MEMLOG_FILENAME, "a") as myfile:
  myfile.write("#key last_updated xml_size_in_byte(6) estimated_xml_size_in_byte(7) memory_before memory_after deleted_count\n")

MEMLOG_RECEIVELOG = []

# global data dict always contains the latest xml file (in parsed version)
GLOBAL_DATA_DICT = {}
GLOBAL_RAW_DATA_DICT = {}
GLOBAL_DATA_SIZE_DICT = {}

# contains the last update
GLOBAL_UPDATED_DICT = {}

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

# in bytes
MAX_MEMORY = 2*1024*1024*1024 # = 2 GByte
#MAX_MEMORY = 250*1024*1024 # = 250 MByte

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
    return render_html.render_index(GLOBAL_DATA_DICT, GLOBAL_UPDATED_DICT, request)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route('/status', methods = ['GET'])
def status():
  key_to_delete = request.args.get('delete', '')
  if key_to_delete != '':
    if key_to_delete in GLOBAL_DATA_DICT:
      oldest_data = GLOBAL_DATA_DICT.pop(key_to_delete) # delete oldest xml keys
      del oldest_data
    if key_to_delete in GLOBAL_UPDATED_DICT:
      oldest_updated = GLOBAL_UPDATED_DICT.pop(key_to_delete) 
      del oldest_updated
    if key_to_delete in UPDATE_EVENTS:
      oldest_update_event = UPDATE_EVENTS.pop(key_to_delete)
      oldest_update_event.set()
      del oldest_update_event

    send_data.delete_coprocessor(key_to_delete)
    gc.collect()
  
  return render_html.render_status(GLOBAL_DATA_DICT, MAX_MEMORY, MEMORY_BYTE_RATIO, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, GLOBAL_KEY_DEQUEUE, app, GLOBAL_STAT_VARS)

@app.route('/status_log', methods = ['GET'])
def status_log():
    return render_html.render_status_log(GLOBAL_DATA_DICT, MEMLOG_RECEIVELOG)

@app.route('/status_memory_pics', methods = ['GET'])
def status_memory_pics():
    return render_html.render_status_memory_pics(GLOBAL_DATA_DICT, MAX_MEMORY, MEMORY_BYTE_RATIO, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, GLOBAL_KEY_DEQUEUE, app)

@app.route(settings["api"]["values"] + '/<path:key>', methods = ['GET', 'POST'])
def values(key):
  if request.method == 'GET':
    return api_plot.get_values(GLOBAL_DATA_DICT, UPDATE_EVENTS, key, int(request.args.get('iteration_num')))
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route(settings["api"]["view_url"] + '/<path:key>', methods = ['GET', 'POST'])
def view(key):
  if request.method == 'GET':
    client_ip = request.remote_addr
    if request.headers.has_key(PROXY_REAL_IP_HEADER):
      client_ip = request.headers.get(PROXY_REAL_IP_HEADER, "127.0.0.1")
    return render_html.render_view(GLOBAL_DATA_DICT, key, client_ip)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route(settings["api"]["plot_url"] + '/<path:key>', methods = ['GET', 'POST'])
def plot(key):
  return api_plot.plot(key, UPDATE_EVENTS, GLOBAL_DATA_DICT, request.args.get('x'), \
                       request.args.getlist('y1_it'), request.args.getlist('y2_it'), \
                       request.args.getlist('y1_res'), request.args.getlist('y2_res'), \
                       int(request.args.get('iteration_num')), \
                       (request.args.get('logscale_y1') == 'true'), (request.args.get('logscale_y2') == 'true'), \
                       request.args.getlist('view_results'), \
                       (request.args.get('view_result_bloch', 'false') == 'true'))


@app.route(settings["api"]["catalyst_send"] + '/<path:key>', methods = ['GET', 'POST'])
def send_data_func(key):
  send_data.send_data(key, GLOBAL_DATA_DICT[key], request.args.get('ip'), request.args.get('port'))
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
  global MEMORY_BYTE_RATIO, GLOBAL_DATA_DICT, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, UPDATE_EVENTS, GLOBAL_STAT_VARS
  process = psutil.Process(os.getpid())
  
  memory_in_bytes = process.memory_info().rss
  
  total_xml_size = 0
  deleted_count = 0
  for tmp_key in GLOBAL_DATA_DICT:
    total_xml_size += GLOBAL_DATA_SIZE_DICT[tmp_key]

  GLOBAL_STAT_VARS['total_bytes_received'] += total_xml_size

  if memory_in_bytes > MAX_MEMORY:
    if MEMORY_BYTE_RATIO < 0:
      # we need to calculate the ratio:
      MEMORY_BYTE_RATIO = memory_in_bytes / total_xml_size

    xml_size_to_reduce = total_xml_size - (MAX_MEMORY / MEMORY_BYTE_RATIO) # size in bytes of xml

    while xml_size_to_reduce > 0:
      oldest_key = GLOBAL_KEY_DEQUEUE.pop() # get oldest xml key
      if oldest_key in GLOBAL_DATA_DICT:
        oldest_data = GLOBAL_DATA_DICT.pop(oldest_key) # delete oldest xml keys
        xml_size_to_reduce -= GLOBAL_DATA_SIZE_DICT.pop(oldest_key)
        del oldest_data
      if oldest_key in GLOBAL_UPDATED_DICT:
        oldest_updated = GLOBAL_UPDATED_DICT.pop(oldest_key) 
        del oldest_updated
      if oldest_key in UPDATE_EVENTS:
        oldest_update_event = UPDATE_EVENTS.pop(oldest_key)
        oldest_update_event.set()
        del oldest_update_event

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
        print("new file")
        break
      
      data += chunk.decode("utf-8")

    xml = etree.fromstring(data)
    
    print("status: " + xml.xpath('//cfsInfo/@status')[0])
    
    key = xml.xpath('//environment/@host')[0] + '/' + xml.xpath('//progOpts/@problem')[0]
    
    if len(xml.xpath('//cfsInfo/header/@id')) > 0:
      if not url_key == '':
        key += url_key + '/' + xml.xpath('//cfsInfo/header/@id')[0]
    else:
      if not url_key == '':
        key += url_key

    key += '/' + xml.xpath('//header/environment/@started')[0]

    if not key in GLOBAL_DATA_DICT: # only add key if not already in
      GLOBAL_KEY_DEQUEUE.appendleft(key)
      GLOBAL_STAT_VARS['total_problem_count'] += 1

    GLOBAL_DATA_DICT[key] = xml
    GLOBAL_RAW_DATA_DICT[key] = data
    GLOBAL_DATA_SIZE_DICT[key] = len(data)
    
    GLOBAL_UPDATED_DICT[key] = str(datetime.datetime.now())
    
    # set the event to trigger sending of the new xml data
    # make sure to set the event AFTER putting the new xml
    # into GLOBAL_DATA_DICT
    if key in UPDATE_EVENTS:
      UPDATE_EVENTS[key].set()

    memory_in_bytes_after_adding = psutil.Process(os.getpid()).memory_info().rss

    with MEMLOG_LOCK:
      with open(MEMLOG_FILENAME, "a") as myfile:
        myfile.write(key + " " + GLOBAL_UPDATED_DICT[key] + " " + str(total_xml_size) + " " + str(MEMORY_BYTE_RATIO*total_xml_size) + " " + str(memory_in_bytes) + " " + str(memory_in_bytes_after_adding) + " " + str(deleted_count) + "\n")

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

    print("recieved data!\n")
    return 'data recieved!\n'
  else:
    return 'expected a POST'



if __name__ == "__main__":
  # need multithreading to server multiple clients.
  # Otherwise pushing xml while ajaxing it with the event
  # synchronization will not be possible
  app.run(threaded=True, port=5000, host='127.0.0.1')
