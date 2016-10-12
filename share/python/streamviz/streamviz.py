# start for development with the following commands in a terminal
# export FLASK_APP=streamviz.py
# flask run

from flask import Flask
from flask import request 
app = Flask(__name__)

@app.route('/receive_cfs_streaming', methods=['GET','POST']) 
def login():
  if request.method == 'POST': 
    print request.data
    return '' # anything more informative will be printed by cfs
  else: 
    return 'expected a POST'
  
