#!/usr/bin/env python

# the purpose of this tool is to process snopt output files. There are three purposes
# 1) extract the convergence information as a nice tabular file, e.g. for gnuplot or plotviz.py
# 2) when a second file is given, filter from the second file the lines corresponding to major iterations of the snopt file
# 3) provide the information of 1) as a function to be called directly from plotviz.py

import os.path
import sys
from cfs_utils import *


# concenience function for plotviz.py
def process(file):
  header, data = parsesnopt(file)
  return tostring(header, data)

# make nice strings out of lists
def tostring(header, data):
  h = []
  h.append('# ' + ' '.join(header))

  d = [' '.join(l) for l in data]

  return h, d 


# parses the snopt file and gives back header and data
# data is filled where necessary but everything is strings
# @param file a snopt file
# @return header is single array of strings and data is vector of vectors of strings
def parsesnopt(file):
  f = open(file, 'r')
  lines = f.readlines()

  # single array of strings
  header = []
 
  # array of array with string content
  data = []
  
  mode = 0 # 0 is before or within pivoting, 1 is within iterations
 
  for o in lines:
    l = o.strip()
    ls = l.split()
    #print(mode,l)
    if len(ls) > 1 and (ls[0] == 'Itns' or ls[0] == 'Itn'): # it for Pivot stuff
      if mode == 0 and not 'Pivot' in l:
        # we need to skip the last columns. The last contains symbols we remove in the body
        header = l.split()
        mode = 1
      elif 'Pivot' in l:
        mode = 0 # pivot switches it off
      continue # just skip further Itns lines   
    if mode == 0:
      continue # we are before the data
    
    if l.startswith('SNOPTA'):
      #assert len(data) > 2
      break # finished reading file

    # skip too short lines with 'line numbers'?!
    if len(l) < 3:
      continue
    
    # now process the content lines
    # remove not consistent and not plotable symbols
    for s in ['_','r', 's', 'n', 'M', '(', ')', 'm', 'R', 'T','l','i']:
      l = l.replace(s,' ')
      
    data.append(l.split())       
   
  # no try to fix the stuff
  hdl = len(header)
  mil, mal = minmax(data)

  # the zero iteration has no step value but also further iterations might miss a step (e.g. for mutiple majors!!!)
  idx = header.index('Step')
  for l in data:
    if not '.' in l[idx]:
      assert len(l) < hdl
      l.insert(idx, '0.0')

  # when no constraint is given, the feasible column is empty
  if not 'MeritFunction' in header:
    assert mal < hdl 
    assert header.index('Feasible') == 5
    for i in range(0,len(data)): 
      data[i].insert(5,'0.0')

  mil, mal = minmax(data)

  # when there is a penalty it might not be set fror the first iteration.
  if 'Penalty' in header:
    assert header[-1] == 'Penalty' # assume last, can easily be generalized
    assert header[-2] == 'condHz'  # also float
    assert header[-3] == 'nS'      # no float!
    for i, l in enumerate(data):
      if not ('.'  in l[-1] and '.' in l[-2]): # we don't have both, condHz and Penalty.
        if '.' in l[-1]: # what shall be condHz is ns, hence Penalty is missing
          data[i].append('0.0')# Penalty
        else:
          # both is missing
          data[i].append('0.0') # condHz
          data[i].append('0.0') # Penalty
  else:
    assert header[-1] == 'condHz'  # float
    assert header[-2] == 'nS'      # no float!
    for i, l in enumerate(data):
      if not '.'  in l[-1]:
        data[i].append('0.0')# condHz

  mil, mal = minmax(data)
  
  # now BSwap and possibly ns might be missing a
  assert mil == mal == hdl or (hdl <= mal+2 and hdl <= mil+2) # we assume that missing BSwap and nS is the only issue
  assert 'BSwap' in header and 'nS' in header
  for i, l in enumerate(data):
    if len(l) < hdl:
      assert len(l) >= hdl-2
      if len(l) == hdl-1:
        data[i].insert(header.index('BSwap'),'0')
      else: # note that an insert also affects the l!
        data[i].insert(header.index('BSwap'),'0')
        data[i].insert(header.index('nS'),'0')

  mil, mal = minmax(data)

  # we validate, that MerritFunction or Optimal is a float
  idx = header.index('MeritFunction') if 'MeritFunction' in header else header.index('Optimal') 
  for l in data:
    assert '.' in l[idx]
    assert len(l) == hdl

  
  # we might have the strange case, that majors are doubled. In that case later might have no step (which is inserted as 0.0 further above)
  # we choose to take the last one
  idx = header.index('Major')
  majors = [int(l[idx]) for l in data]
  remove = [] # indices to be removed
  for i in range(0,len(majors)-1):
    if majors[i] == majors[i+1]:
      remove.append(i) 

  # now remove the stuff. We dont't want to alter the list during removing items of it, so do it in two parts
  base = 0 # index for removed elements
  for i in remove:
    del data[base + i]
    base -= 1 # data got smaller, hence decrement

  #print(header, len(header))
  #for l in data:
  #  print(l,len(l))
  #print(hdl,mil,mal)


 
  return header, data
  
# give min / max for data
def minmax(data):
  mil = min([len(t) for t in data[1:]]) 
  mal = max([len(t) for t in data[1:]])
  return mil, mal

# parse dat file and extend relevant header line by snopt_major columns
# return list of header string (one is extended!) and list of content strings
def parseother_and_extend_snopt_major(file):   
  # copy and paste from plotviz.py process()
  f = open(file, 'r')
  lines = f.readlines()
 
  header = []
  body = [] 
 
  for l in lines:
    h = l.strip()
    if h.startswith('#'): 
      header.append(h)
    else:
      body.append(l.strip())  
  
  if header:
    cols = len(body[0].split()) # number of columns we expect

    # loop through header and extend the first which has cols splittings
    for idx, line in enumerate(header):
      l = line.strip()[1:].strip() # skip #
      for test in [':', '(', '|', ',', '\t']:
        if l.count(test) == cols:
          header[idx] += ' \t'
          if test == ':':
            header[idx] += str(cols+1) + ':'
          elif test == '(':                  
            header[idx] += '(' + str(cols+1) + ') '
          else: 
            header[idx] +=  test
          header[idx] += 'snopt_major'
          return header, body
  
  return header, body
   
# plotviz.py is imported by plotviz.py, so guard argparse
if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(description="process a .snopt output file. When a second file is given, it is filtered by snopt file's majors")
  parser.add_argument("input", nargs='+', help="first file shall be a snopt file, optioanlly a second .plot.dat or .density.xml")
  
  args = parser.parse_args()

  for input in args.input:
    if not os.path.exists(input):
      print('Error: given file not found', input)
      sys.exit(1)

  if len(args.input) > 2: # 0 case handled by arparse 
    print('Error: give on or two files as arguments')
    sys.exit(1)       

  header, data = parsesnopt(args.input[0])
  
  if len(args.input) == 1:
    print('#', '  '.join(header))
    for l in data:
      print('  '.join(l))
  else:
    idx = header.index('nObj' if 'nObj' in header else 'nCon')
    neval = [int(l[idx]) for l in data]

    if args.input[1].endswith('.dat'):
      # we want to filter the content which belongs to majors. 
      comments, body = parseother_and_extend_snopt_major(args.input[1]) # list of strings
      
      # e.g. neval[-1] = 31 --> len(body) shall be 31. 
      if neval[-1] != len(body):
        print('# Warning: expect',neval[-1],'data lines in ', args.input[1],'but have',len(body), 'possible missmatch with',args.input[0])
      
      # print the comments where one line shall be enriched by major_snopt
      for l in comments: 
        print(l) 
      # fault tolerant printing all data  
      for idx in range(min(len(neval), len(body))): # we could add somehow the major
        major = neval[idx] -1
        print(body[major] + ' \t' + str(idx)) # snopt is one-based.
    elif args.input[1].endswith('.density.xml'):
      xml = open_xml(args.input[1])
      for i in range(neval[-1]):
        if i == (neval[0]-1):
          neval.pop(0)
        else:
        #  pathname = '//set[@id="' + str(i) + '"]'
        #  print(pathname)
          remove(xml, '//set[@id="' + str(i) + '"]', unique = True)
      out = args.input[1][0:-12] + '.snopt.density.xml'    
      xml.write(out, xml_declaration=True, encoding='UTF-8')
      print('created file',out)
    else:
      print('Warning: File type of second argument not permissible!\n EXIT')
# here could be an else case for the import into plotviz part   
