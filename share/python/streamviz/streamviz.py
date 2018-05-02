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

from flask.wrappers import Request

app = Flask(__name__)

# global data dict always contains the latest xml file (in parsed version)
GLOBAL_DATA_DICT = {}
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

@app.route('/', methods = ['GET'])
def index():
  if request.method == 'GET':
    return render_html.render_index(GLOBAL_DATA_DICT, GLOBAL_UPDATED_DICT, request)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route('/status', methods = ['GET'])
def status():
    return render_html.render_status(GLOBAL_DATA_DICT, MAX_MEMORY, MEMORY_BYTE_RATIO)

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
  
@app.route(settings["api"]["recieve_url"] + '/', methods = ['GET', 'POST'])
@app.route(settings["api"]["recieve_url"] + '/<path:url_key>', methods = ['GET', 'POST'])
def cfs_recieve(url_key = ""):
  global MEMORY_BYTE_RATIO, GLOBAL_DATA_DICT, GLOBAL_DATA_SIZE_DICT, GLOBAL_UPDATED_DICT, UPDATE_EVENTS
  process = psutil.Process(os.getpid())
  
  memory_in_bytes = process.memory_info().rss
  
  if memory_in_bytes > MAX_MEMORY:
    total_xml_size = 0
    for tmp_key in GLOBAL_DATA_DICT:
      total_xml_size += GLOBAL_DATA_SIZE_DICT[tmp_key]

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
      send_data.delete_coprocessor(oldest_key)
      
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

    GLOBAL_DATA_DICT[key] = xml
    GLOBAL_DATA_SIZE_DICT[key] = len(data)
    
    GLOBAL_UPDATED_DICT[key] = str(datetime.datetime.now())
    
    # set the event to trigger sending of the new xml data
    # make sure to set the event AFTER putting the new xml
    # into GLOBAL_DATA_DICT
    if key in UPDATE_EVENTS:
      UPDATE_EVENTS[key].set()

    print("recieved data!\n")
    return 'data recieved!\n'
  else:
    return 'expected a POST'



if __name__ == "__main__":
  # need multithreading to server multiple clients.
  # Otherwise pushing xml while ajaxing it with the event
  # synchronization will not be possible
  app.run(threaded=True, port=5000, host='127.0.0.1')
