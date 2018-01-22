#!/usr/bin/env python2

# webserver:
from flask import Flask
from flask import request

# xml parser
from lxml import etree

# settings
from strings import *

# html render module
import render_html

# plot module
import api_plot

app = Flask(__name__)

# global data dict always contains the latest xml file (in parsed version)
GLOBAL_DATA_DICT = {}

# ajax functions will wait for an even aka new data
# in order to do that they will put and event into this dictionary: key => Array(Event)
# on new cfs_recieve, the events will be set() and deleted from the dict, releasing all
# events and firing the http return.
# current timeout for events should be set to 9 seconds server side and 10 seconds 
# client side within javascript. This will ensure performance and prevent useless
# data transfer and html updating
UPDATE_EVENTS = {}

@app.route('/', methods = ['GET', 'POST'])
def index():
  if request.method == 'GET':
    return render_html.render_index(GLOBAL_DATA_DICT, request)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'


@app.route(settings["api"]["view_url"] + '/<path:key>', methods = ['GET', 'POST'])
def view(key):
  if request.method == 'GET':
    return render_html.render_view(GLOBAL_DATA_DICT, key)
  else:
    return 'expected a GET, use "' + settings["api"]["recieve_url"] + '" to send data'

@app.route(settings["api"]["plot_url"] + '/<path:key>', methods = ['GET', 'POST'])
def plot(key):
  return api_plot.plot(key, UPDATE_EVENTS, GLOBAL_DATA_DICT, request.args.get('x'), \
                       request.args.getlist('y1_it'), request.args.getlist('y2_it'), \
                       request.args.getlist('y1_res'), request.args.getlist('y2_res'), \
                       request.args.get('iteration_num'))

@app.route(settings["api"]["recieve_url"], methods = ['GET', 'POST'])
def cfs_recieve_blank():
  return cfs_recieve("")
  
@app.route(settings["api"]["recieve_url"] + '/', methods = ['GET', 'POST'])
@app.route(settings["api"]["recieve_url"] + '/<path:key>', methods = ['GET', 'POST'])
def cfs_recieve(key = ""):
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
    
    if len(xml.xpath('//cfsStreaming/@id')) > 0:
      if key == '':
        key = xml.xpath('//cfsStreaming/@id')[0] + '/' + xml.xpath('//header/environment/@started')[0]
      else:
        key += '/' + xml.xpath('//cfsStreaming/@id')[0] + '/' + xml.xpath('//header/environment/@started')[0]
    else:
      if key == '':
        key = xml.xpath('//header/environment/@started')[0]
      else:
        key += '/' + xml.xpath('//header/environment/@started')[0]

    GLOBAL_DATA_DICT[key] = xml
    
    # set the event to trigger sending of the new xml data
    # make sure to set the event AFTER putting the new xml
    # into GLOBAL_DATA_DICT
    if key in UPDATE_EVENTS:
      for e in UPDATE_EVENTS[key]:
        e.set()

    print("recieved data!\n")
    return 'data recieved!\n'
  else:
    return 'expected a POST'



if __name__ == "__main__":
  # need multithreading to server multiple clients.
  # Otherwise pushing xml while ajaxing it with the event
  # synchronization will not be possible
  app.run(threaded=True)
