#!/usr/bin/env python

# the purpose of the tool is to provide easy matplotlib visualization for dat files like .plot.dat, .grad.dat and .snopt
# it extracts iformation from a header line with like #(1) iter \t(2) compliance ....
import os.path
import sys
import re
import glob
import datetime
import numpy as np
from itertools import cycle

known_functions = ('abs', 'sqrt', 'pi', 'sin', 'cos')

# when we use plotviz from postproc.py, we don't need this stuff
if __name__ == '__main__':
  import argparse
  import matplotlib
  import matplotlib.pyplot as plt
  from matplotlib.ticker import MaxNLocator
  import snopt # our snopt.py helper fticklabel_formator process

  # in case we have --dashed we use c_cms_y and c_cms_y2
  # https://stackoverflow.com/questions/7358118/matplotlib-black-white-colormap-with-dashes-dots-etc
  # probably to be combined with --black
  from cycler import cycler
  color_c = cycler('color', ['k'])
  markr_c = cycler('marker', ['', '.', 'o'])
  style_c_y = cycler('linestyle', ['-', '--', ':', '-.'])
  style_c_y2 = cycler('linestyle', ['-.', ':','--', '-'])
  c_cms_y = color_c * markr_c * style_c_y
  c_cms_y2 = color_c * markr_c * style_c_y2

# having two y2-axis we need to handle colors manually, otherwise they repeat
# https://matplotlib.org/stable/gallery/color/named_colors.html
# 'gold' has index 13 for y2 axis
# in case of --black will be replaced below by all 'black'
colors = ['tab:green','tab:red','tab:purple','tab:blue','tab:orange','black','tab:brown','tab:gray','tab:olive','blue','tab:cyan','tab:pink','cornflowerblue', 
          'gold','peru','blueviolet', 'coral','yellowgreen'] 
colors += colors # repeat such that it should be really enough

# this will cycle through different marker/line styles
markers = [".","x","+","v","^","<",">","s","*"]
lines = ["-","--",":","-."]
markercycler = cycle(markers)
linecycler = cycle(lines)

# parses the header lines for the first hint on column names.
# Tries to be smart!!
# @param data the processed content() result. To fill up missing hints in the comment
# @param comments a list of lines which start with a hashtag
# @return meta array
def header(data, comments):

  meta = []
  for l in reversed(comments):
    # we assume the last comment contains the labels. Special case are ignored
    if l.startswith('##') or l.startswith('--') or '----' in l:
      continue
  
    if l.startswith('#'):
      l = l.strip()[1:].strip() # remove trailing single '#'
    assert not l.startswith('#')

    # now the cases
    if l.count(':') > 1: # #1:iter  2:compliance  3:duration
      ll = l.split(':') # ['1', 'iter\t2', 'compliance\t3', 'duration']
      for i, tl in enumerate(ll):
        # ['1'] / ['iter', '2'] / ['compliance', '3'] / ['duration']
        # or in the grad.dat case stuff like ['d_compliance', '/', 'd_s0_px', '9']
        t = tl.split()
        if i == 0:
          assert len(t) == 1 # skip this the first number
        elif i == len(ll)-1:
          meta.append(''.join(t))
        else:
          meta.append(''.join(t[:-1])) # all but the last: ['d_compliance', '/', 'd_s0_px', '9'] -> 'd_compliance/d_s0_px'

    elif l.startswith('(1)'):
      # we cut the '(x)' as we have our own numbering
      for t in l.split('('): #  ['', '1)No ', '2)Date ', '3)Time ', '4)CO2 ', '5)Temp ', '6)Humi']
        if t.find(')') >= 0:
          meta.append(t[t.find(')')+1:].strip())

    elif l.count('\t') >= 1:
      t = l.split('\t')
      meta = [s.strip() for s in t]

    elif l.count('|') >= 1:
      t = l.split('|')
      meta = [s.strip() for s in t]

    elif l.count(',') >= 1:
      t = l.split(',')
      meta = [s.strip() for s in t]

    else:
      # can be the postroc.py case with (1) in separate line and possibly empty 'id'
      
      # lables with '/' are from sensitivity like 'd_compliance / d_s0_px' -> compress to have a single value
      l = l.replace(' / ', '/')
      meta = l.split()
      if meta[-1] == 'id' and len(meta) == len(data[0])+1:
        meta.remove('id')
        
    # skip the loop, we've done it    
    break      

  # for CO2 we combine the columns Date and Tine to one datetime -> leave only time
  if 'Time' in meta:
    meta.remove('Date')
  if 'time' in meta and 'date' in meta: # fan_control.dat
    meta.remove('date') 
  if 'Zeit' in meta:
    meta.remove('Zeit') 


  # fill missing columns (or all)
  for _ in range(len(data[0]) - len(meta)):
    meta.append('anonymous')

  # handle openCFS history case # t (s) (K) -> ['t', '(s)', '(K)']
  if len(meta) > len(data[0]) and len(meta) >= 2 and '(' not in meta[0] and '(' in meta[1]:
    meta[0] += meta[1]
    del meta[1]

  if len(meta) > len(data[0]): # even the above case might not fix all
    print('Error: found too much meta data ', meta, '=', len(meta), ' for data ', data[0], '=',len(data[0]))
    sys.exit()     

  return meta

# print header in a nice way
def print_header(meta,data,inputs):
  assert len(meta) == len(data) == len(inputs)
  
  ml = max([len(max(m, key = len)) for m in meta]) # find largest meta key
  
  print('key: ' + 'label'.ljust(ml) + ' : first value         : file')
  cnt = 1
  for mi, m in enumerate(meta):
    for li, l in enumerate(m):
      # at time of printing we did not smooth or grad yet
      val = '*' if 'grad_' in l or 'smooth_' in l else str(data[mi][0][li]) 
      # file = inputs[mi-1] if mi > 0 else '-' # first pseudo input is index
      file = inputs[mi]
      print('{:3d}'.format(cnt) + ': ' + l.ljust(ml) + ' : ' + str(val).ljust(19) + ' : ' + file) 
      cnt += 1


# parse the content and give back an array. Combines date and time to datetime
# data is either float or datetime
def content(body):
  if len(body) == 0:
    print('Error: no content lines in the given file')
    sys.exit()
  
  assert not body[0].strip().startswith('#')
  
  data = []

  # for what do we split
  key = None # default
  if body[0].count('|') > 1:
    key = '|'
  elif body[0].count(',') > 1:
    key = ','
  elif body[0].count('\t') > 1:
    key = '\t'

  # no tabs mean spaces, or comma and we need to check for date and time splitted
  s = body[0].split(key)

  # we expect a date first or second (after counter for co2 data)
  dt0 = len(s) > 1 and check('%Y-%m-%d',s[0]) and check('%H:%M:%S', s[1])
  dt1 = len(s) > 2 and check('%Y-%m-%d',s[1]) and check('%H:%M:%S', s[2])
  # german style
  dt0 = dt0 or (len(s) > 1 and check('%d.%m.%Y',s[0]) and check('%H:%M', s[1]))
  for l in body:
    if len(l.strip()) == 0:
      continue # skip empty lines)
    # handle nasty data where numbers are with , within quotes: 2023/03/11 17:03:15,"18,06","56,3",""
    if key == ',' and l.count('"') > 0 and l.count('"') % 2 == 0: # we have an even number of quotes
      s = ''
      inner = False
      for c in l: 
        if c == '"':  
          inner = True if not inner else False 
        else: 
          s += '.' if c == ',' and inner else c 
      l = s 
    t = [s.strip() for s in l.split(key)]
    if dt0:
      data.append([t[0] + ' ' + t[1]] + t[2:]) # combine date and time to datetime in case
    elif dt1:    
      data.append([t[0]] + [t[1] + ' ' + t[2]] + t[3:]) 
    else:
      data.append(t) 

    if len(data[0]) != len(data[-1]):
      print("Error: inconsistent number of entities with line", data[-1] , len(data[-1]),'vs. first line',data[0],len(data[0]))
      sys.exit()

  # check for datetime in the first columns and replace in case
  for c in range(min(len(data[0]),3)):
     # ..., 28.03.21,  28.03.2021, 2021-01-04
    for frmt in ['%Y-%m-%d %H:%M:%S', '%d.%m.%Y %H:%M', '%Y/%m/%d %H:%M:%S','%d.%m.%y', '%d.%m.%Y', '%m.%Y', '%Y-%m-%d']:
      if check(frmt, str(data[0][c])):
        for l in data:
          l[c] = datetime.datetime.strptime(l[c], frmt)


  # convert str to float
  for i, c in enumerate(data[0]):
    if type(c) == str:
      for r in data:
        try:
          if r[i].endswith('%'): # remove trailing % from bluetooth sensor for rel humidity 
            r[i] = r[i][:-1]
          r[i] = 0 if r[i] == '' else float(r[i])
        except ValueError as ve:
          pass
  # print(data)
  return data

# validate datetime by try and except
def check(format, test):
  try:
    time = datetime.datetime.strptime(test, format)
    #  8.9819 -> 9819-08-01 00:00:00 9819
    if time.year == 0 or (time.year >= 2010 and time.year <= 2030):
      return True
    else:
      return False
  except ValueError:
    return False
  
# gives back the file index and the corresponding 0-based index.
# @param key if string, search for uniqueness, if number, make double 0-based
# @return 0-based file index, 0-based column within file, label
def find_index(meta, key):
  fi = -1  # file index
  idx = -1 # relative within file index
  if all(map(str.isdigit, key)):
    k = int(key)-1 # key from print_header() as 0-based
    base = 0
    for f, m in enumerate(meta): # traverse files
      if f > 0:
        base += len(meta[f-1])
      if k < base + len(m) : # is this our matching file?
        idx = k-base
        fi = f
        break
    if idx == -1:
      print('Error: given key', k+1, 'out of range')
      sys.exit()
  else:
    for file, m in enumerate(meta):
      for i, t in enumerate(m):
        contained = t.startswith(key) if not noautocomplete else t == key  
        if contained:
          if idx > -1:
            print("Error: key not unique '", key, "' consider using --noautocomplete")
            sys.exit()
          idx = i
          fi = file
    if idx == -1:
      print("key not found '", key, "'")
      sys.exit()
  return fi, idx, meta[fi][idx]

# transform list of list of keys to list of tuples (1-based-id, key)
def flatten_meta(meta_list):
  res = []
  cnt = 1
  for array in meta_list:
    for key in array:
      res.append((cnt, key))
      cnt += 1
 
  return res 

# contrary to find_index. A label in arg (e.g. arg.x) is checked for multiple occurence in meta.
# Only if so, the key in arg is replaced by all indices of meta. 
def resolve_multiple(meta, args):
  if args is None:
    return None 
  ret = []
  flat = flatten_meta(meta)
  for key in args:
    ids = [str(e[0]) for e in flat if key == e[1]]
    if len(ids) > 1:
      ret.extend(ids)
    else:
      ret.append(key)
  return ret   

# gets back the column by index descrition
# file-index is 1-based when the key is Not Note (-x) and encodes bar
# @param key if None return range, can be a list
# @return arrays of file-index, data column, label
def column(meta, data, key, bars):
  if key == None:
    # for default if -x is not given. Shall not be callend for y2 is None
    n = len(data)
    return [i for i in range(n)] , [range(len(data[i])) for i in range(n)], [''] * n
  else:
    if type(key) != list:
      if any(op in key for op in ['+','-','*','/']):
        return eval_expression(meta, data, key, bars)
      fi, idx, label = find_index(meta, key)
      return (-(fi + 1) if bars and key in bars else fi+1), [d[idx] for d in data[fi]], label
    else:
      fl = []
      dl = []
      ll = []
      for k in key:
        f, d, l = column(meta,data,k,bars)
        fl.append(f)
        dl.append(d)
        ll.append(l)
      return fl, dl, ll  

# evaluate a expression like '2+compliance-$5'
def eval_expression(meta, data, key, bars):
  # remove whitespaces and split expression at operators
  exp = re.split(r'([\+\-\*\/\(\)])', key.replace(' ',''))
  # remove empty strings (non-printable characters get inserted before ( and ) by bash)
  exp = [kk for kk in exp if kk]
  # extract column numbers/labels which we want the data for
  k = []
  for kk in exp:
    if kk in known_functions:
      continue
    if kk[0] == '$':
      k.append(kk[1:])
    if re.search('[a-zA-Z]', kk):
      k.append(kk)
  # get data
  f, d, l = column(meta, data, k, bars)
  # assert all file indices are equal. else d might belong to different x values
  assert(f[:-1] == f[1:])
  # convert to numpy to be able to operate on arrays
  d = np.array(d)
  # replace '$*' in expression by 'd[i]' and in label by l[i]
  label = exp.copy()
  exp_idx = 0
  for i in range(len(exp)):
    if exp[i][0] == '$' or re.search('[a-zA-Z]', exp[i]):
      if exp[i] in known_functions:
        exp[i] = 'np.' + exp[i] # use numpy functions
        continue
      exp[i] = 'd[{:d}]'.format(exp_idx)
      label[i] = l[exp_idx]
      exp_idx += 1
  # evaluate expression with data from d
  d = eval(''.join(exp))
  return f[0], d, ''.join(label)


# extract kwnown extension of filename
def filename_base(filename):
  if filename.endswith('.plot.dat'):
    return filename[:-9]
  name, _ = os.path.splitext(filename)
  return name


# modify multiple occurence of lables by augmenting them with the filename
def fix_labels(labels, fileindices, filenames):
  names = [fn for fn in filenames]
  # reduce filenames by cutting common tailing parts
  for _ in range(int(filenames[0].count('.'))):  
    # tails of all filenames
    tails = [f[f.rfind('.'):] for f in names] # .dat from killme_scpip.plot.dat
    # reduce names if all extensions are the same. If not the case we repeat the for loop with with same data (no change)
    if tails.count(tails[0]) == len(tails):
      names = [f[:f.rfind('.')] for f in names] # killme_scpip.plot.dat -> killme_scpip.plot 
  res = []
  for fi, l in enumerate(labels):
    if labels.count(l) == 1:
      res.append(l)
    else:
      res.append(l + ' (' + names[fileindices[fi]-1] + ')') # fileindices are 1-based
  return res    
        
## find idx within column for the given start and end value which can be datetime
# @param x single column of datetime
# @return start and end idx within x range. end index is exclusive      
def bounds(x, start, end):
  
  if start >= x[-1] or end <= x[0]:
    return 0,0
  
  sidx = 0
  for i in range(len(x)-1): # loop until we are not too early
    #print(i,'cmp',x[i],start,x[i] > start)
    if x[i] > start:
      sidx = i
      break
         
  eidx = len(x)
  for i in range(len(x)-1,0,-1): # list(range(4,0,-1)) -> [4, 3, 2, 1]
    if x[i-1] > end:
      eidx = i-1
    else:    
      break      
      
  return sidx, eidx   
   
# process and input file, returns meta and data
# also used by postproc.py
def process(input):
  file = open(input, 'r')
  lines = file.readlines()
 
  comments = []
  body = [] 
 
  # we assume first comments (and assume the last comment to be the header descriotion)
  # then the body. For comment/body/comment we ignore comments after body
  for l in lines:
    h = l.strip()
    # it seems excel creates utf-8 bom at file start, simply skip it
    if ord(h[0]) == 0xfeff:
      h = h[1:]
    if h.startswith('#') or h.startswith('iter') or h.startswith('Temp') or h.startswith('---'):
      if len(body) == 0: # ignore comments after body 
        comments.append(h)
    else:
      body.append(l)  
  data = content(body)
  meta = header(data,comments)
  return meta, data        

# process results from a .info.xml - not optimization iterations!
# adds pairs of key/value for any result > 1 element
def process_info_xml_results(input):
  import lxml
  import lxml.etree
  xml = lxml.etree.parse(input, lxml.etree.XMLParser(remove_comments=True))
  results = xml.xpath('//result')

  meta = [] # list of headers  
  tmp_data = [] # list of data corresponding to meta but probably inconsistent length
  max_len = 0 # to fill up for all results
  
  for r in results:
    items = r.xpath('item')
    if len(items) > 1:
      type  = r.attrib['data'] # heatTemperature
      loc   = r.attrib['location'] # outlet_nodes
      defon = r.attrib['definedOn'] # 'node', 'element'
      key = type + '_' + loc
      if 'step_val' in items[0].attrib and items[0].attrib['step_val'] != items[-1].attrib['step_val']:
        meta.append(key + '-step')
        tmp_data.append([float(x.attrib['step_val']) for x in items])
      if 'id' in items[0].attrib:
        meta.append(key + '-' + defon)
        tmp_data.append([int(x.attrib['id']) for x in items])
      if 'value' in items[0].attrib:
        meta.append(key + '-value')
        tmp_data.append([float(x.attrib['value']) for x in items])
      if 'x' in items[0].attrib:
        meta.append(key + '-x')
        tmp_data.append([float(x.attrib['x']) for x in items])
      if 'y' in items[0].attrib:
        meta.append(key + '-y')
        tmp_data.append([float(x.attrib['y']) for x in items])
      if 'z' in items[0].attrib:
        meta.append(key + '-z')
        tmp_data.append([float(x.attrib['z']) for x in items])
      max_len = max(max_len, len(items))  
  
  # we need to fill up data which might not be constant when defined on different entitites
  data = []
  for r in range(max_len):
    row = []
    for c in tmp_data:
      row.append(c[r] if r < len(c) else np.nan) # pyplot skips nan which is very nice
    data.append(row)  
  return meta, data        

# convenience function which gives a y(2)-axis label. 
# @param args either args.ylabel or args.y2label. If not None this is returned
# @oaram legend array of legends. Returns a list of unique entries
def label(args, legends):
  if args is not None:
    return args
  else:
    # we make manually unique without list(set(legends)) to keep order
    res = []
    for l in legends:
      if not l in res:
        res.append(l)
    
    return ', '.join(res) 

# reserve artificial data for smoothing or gradient if given in key
# the real numerics can only be later when we have the x-key and the range, but 
# for this we need the reservation first. The values are the original ones to be processed later
# @param label smooth or grad
def reserve_artificial(input, meta, keys, label):
  if keys is not None and len(keys) > 0:
    # data is list of lines with all columns, convert first for our selected colums
    fix, cols, labels = column(meta,data,keys,[])
    assert len(fix) == len(cols) == len(labels)
    
    for i, c in enumerate(cols):
      meta[fix[i]-1].append(label + '_' + labels[i])
      for j, v in enumerate(c):
        data[fix[i]-1][j].append(v)
        
  return meta, data        
  
# apply smoothing for the already existing (original) data  
def apply_smooth(data, labels, window, poly):
  if sum('smooth_' in l for l in labels) > 0:
    import scipy.signal
  for i, l in enumerate(labels):
    if 'smooth_' in l: # label might be grad_smooth_compliance.
      c = data[i]
      w = min(window, 2*int(len(c)/2) +1) # window nees to be odd
      if not 'grad_smooth_' in l: # we smooth data for grad but don't repeat the smooth_ output of the original data
        print("smooth " + labels[i] + "' with Savitzky–Golay filter with poly", poly, "and window",w)
      sc = scipy.signal.savgol_filter(c, w, poly)

      data[i] = sc # replace data by smoothed data
  return data 

# helper which gives the difference val2 - val1- in hours if values are datetime
def diff(val1, val2):
  d = val2 - val1
  if type(d) == datetime.timedelta:
    return d.seconds * 3600.0
  else:
    return d

# apply smoothing for the already existing (original) data  
# @param fi file index is 1-based and negative if bar data
def apply_grad(fi, data, labels, all_x):
  for i, l in enumerate(labels):
    if 'grad_' in l: # label might be grad_smooth_compliance.
      #print(fi, i, abs(fi[i])-1)
      x = all_x[abs(fi[i])-1]
      c = data[i]
      dc = np.zeros(len(c))
      assert len(x) == len(c)
      assert len(c) > 1
      #print(c)
      #print(x)
      for vi, v in enumerate(c):
        if vi == 0: # forward difference
          vn = c[vi+1]
          dc[vi] = (vn - v) / diff(x[vi], x[vi+1])
        elif vi == len(c) -1: # backward difference
          vp = c[vi-1]
          dc[vi] = (v - vp) / diff(x[vi-1], x[vi])
        else: # central difference
          vn = c[vi+1]
          vp = c[vi-1]
          dc[vi] = (vn - vp) / diff(x[vi-1], x[vi+1])
      data[i] = dc # replace data 
  return data 

# plotviz.py is imported by postproc.py, so guard argparse
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Simple gnuplot replacement for standard plots. Needs a header comment')
  parser.add_argument("input", nargs='+', help="one or more .plot.dat or similar tabular text files. Also .info.xml (result) or .snopt")
  parser.add_argument("-x", nargs='*',  help="index or label for the abscissa (optional). Space separated list for multiple inputs")
  parser.add_argument("-y", nargs='+',  help="indices or labels for the ordinate. Multiple separated by space.\
                                              Expressions are possible in single quotes (e.g. '0.5*compliance -$3/ 2') with references to columns prefixed with $.")
  parser.add_argument("-y2", nargs='*', help="optional indices or labels for the secondary ordinate")
  parser.add_argument("-z", nargs='*',  help="trigger 3D plots wich requires a single -x and -y component")
  parser.add_argument("--range", help='range (day fractions for datetime) to be shown. negative for final range', type=float)
  parser.add_argument("--irange", help='index-based range to be shown. negative for final range', type=int)
  parser.add_argument("--shift", help='shift range (day fractions for datetime)', type=float)
  parser.add_argument("--ishift", help='index-based shift range', type=int)
  parser.add_argument("--ylim", nargs=2, help='range for y-axis', type=float)
  parser.add_argument("--y2lim", nargs=2, help='range for y2-axis', type=float)
  parser.add_argument("--xlabel", help='optional label for the abscissa')
  parser.add_argument("--ylabel", help='optional label for the primary ordinate')
  parser.add_argument("--y2label", help='optional label for the secondary ordinate')
  parser.add_argument("--zlabel", help='optional label for the 3d data')
  parser.add_argument("--marker", help="optional matplotlib marker: e.g. . , o or ' ' to disable marker")
  parser.add_argument("--linestyle", help='optional matplotlib linestyle: e.g. - dashed')
  parser.add_argument("--legend", nargs='*', help="(partially) overwrite labels in the legend as space separated list of strings")
  parser.add_argument("--legend_loc", help="string for matplotlib.legend(loc)")
  parser.add_argument("--legend_ncol", help="number of columns of legend", type=int, default=1)
  parser.add_argument("--title", help='optional title for the plot')
  parser.add_argument("--xscale", help="scaling type from choice, google matplotlib xscale", choices=["linear", "log", "symlog", "logit"],default='linear')
  parser.add_argument("--yscale", help="scaling type from choice, google matplotlib yscale", choices=["linear", "log", "symlog", "logit", "logrel", "percentage"],default='linear')
  parser.add_argument("--y2scale", help="like --yscale but for y2 axis", choices=["linear", "log", "symlog", "logit", "logrel", "percentage"],default='linear')
  parser.add_argument("--bar", nargs='*', help="indices from y or y2 which are to displayed as bars instead of plots")
  parser.add_argument("--barwidth", help="barplots for datetime need manual adjustment", type=float, default=.8)
  parser.add_argument("--smooth", nargs='*', help="create new smoothed data for given fields")
  parser.add_argument("--smooth_window", help="window size of Savitzky–Golay filter", type=int, default = 7)
  parser.add_argument("--smooth_poly", help="polynomial order of Savitzky–Golay filter", type=int, default = 3)
  parser.add_argument("--grad", nargs='*', help="give finite difference gradients, you might want to smooth first")
  parser.add_argument("--dashed", help="cylce through different line styles. Use with --black for b/w",action='store_true')
  parser.add_argument("--black", help="change all line colors to black. Use with --dashed",action='store_true')
  parser.add_argument("--save", help='write to given filename using the extension')
  parser.add_argument("--noautocomplete", help='supress searching only for beginning of key',action='store_true')
  parser.add_argument("--noshow", help='supress popping up the image window', action='store_true')
    
  args = parser.parse_args()
  
  global noautocomplete 
  noautocomplete = args.noautocomplete
  
  if args.black:
    colors = ['black'] * 20
  
  # array of headers per file
  meta = [] 
  # matrix of data per file
  data = []
  
  # handle Windows and macOS debugging
  if len(args.input) == 1:
    args.input = glob.glob(args.input[0]) # replace with more content in case there are Wildcards
  for input in args.input:
    if not os.path.exists(input):
      print('Error: no valid .dat or .snopt file given', input)
      sys.exit()
    m = None
    d = None  
    if input.endswith('.snopt'):
      comments, body = snopt.process(input)
      d = content(body)
      m = header(d,comments)
    elif input.endswith('.info.xml'):
      m, d = process_info_xml_results(input)  
    else:    
      m, d = process(input)
    meta.append(m)
    data.append(d)  
  
  # insert artificial index for linspace x-axis as 0th file
  #meta.insert(0,['index']) # for each file a list of header names
  # the index 'file' are row column lists
  #index_file = []
  #for row in range(len(data[0])): # for each file a list of column vector
  #  index_file.append([row+1]) # 1-based index
  #data.insert(0,index_file)              

  # if one of the arguments is a multiple key in meta, it is replaced by multiple ids (argument becomes larger)
  # we first do this for artificial data
  # for real smoothing we wait for the restrictions
  args.smooth = resolve_multiple(meta, args.smooth)
  meta, data  = reserve_artificial(input, meta, args.smooth, 'smooth')
  
  # for real grad we need the corresponding x and for that we need to wait for print
  args.grad   = resolve_multiple(meta, args.grad)
  meta, data  = reserve_artificial(input, meta, args.grad, 'grad')
  
  args.x = resolve_multiple(meta, args.x)
  args.y = resolve_multiple(meta, args.y)
  args.y2 = resolve_multiple(meta, args.y2)
  args.z = resolve_multiple(meta, args.z)

  # after printing meta is fixed
  print_header(meta, data, args.input)

  # sanity checks   
  if not args.y:
    print('Usage: provide at least -y and possibly -x and -y2. Key/label may be space separated list')
    sys.exit()
  
  if args.x != None and len(args.x) != len(args.input):
    print('Error: on multiple inputs either have no -x or -x with keys for all input files')
    sys.exit()
  
  # number of input arrays of file index, column data, label
  fix, x, xlabel = column(meta,data,args.x,[])
  assert len(x) == len(args.input)
  
  # fiy is 1-base with positive idx for plot and negative for bar
  fiy,  y, ylabel = column(meta,data,args.y,args.bar)

  # fiy2 is also 1-based and encodes bars
  fiy2, y2, y2lbl = column(meta,data,args.y2,args.bar) if args.y2 else ([],[],[])
  
  fiz, z, zlabel = column(meta,data,args.z,args.bar) if args.z else ([],[],[])
  
  for i in range(1,len(x)):
    if type(x[i][0]) != type(x[0][0]):
      print('Error: inconsisten data type for your x axis choice',fix,xlabel)
      sys.exit()      
  
  has_dt = type(x[0][0]) == datetime.datetime # we validated common type for all x before
  
  # now do restrictions
  # for restrictions, this is the start index and end index for the x array of columns
  # currently multiple input needs to be datetime
  start_idx = [0] * len(x)
  end_idx   = [len(t) for t in x]
  
  delta = None # either float or timedelta in the datetime case. Used in the datetime case for adjustment of time axis when plotting

  # no datetime means two options: range/shift by index or by given -x. 
  # when multiple data is given and the files have different length, range/shit makes only sense if -x is given.
  # -x is either given neither or for all files!
  assert not (not has_dt and len(x) > 1 and args.x != None and len(args.x) >= 1 and len(args.x) != len(x)) 
  
  mil = min([len(t) for t in x]) 
  mal = max([len(t) for t in x])
  if (mil != mal) and args.x is None:
    print('Error: data size for mupltiple input files varies',mil,'...',mal,' and -x is not given')
    sys.exit() 

  # when -x is not given, x is range, otherwise it can be anything       
  miv = min([t[0] for t in x])  # min value is earliest for datetime
  mav = max([t[-1] for t in x]) # max value is latest for datetime
  delta  = mav-miv # diff value is a timedifference for datetime 
  if has_dt:
    print('earliest datatapoint:',miv,'latest datapoint:',mav,'days:',round(((delta.days * 86400 + delta.seconds)/86400),3))
  else:
    print('smallest x val:',miv,'largest x val:',mav,'diff:',delta)
    
  if args.range or args.shift:
    r = abs(float(args.range)) if args.range else delta
    s = float(args.shift) if args.shift else 0
    if has_dt:
      r = datetime.timedelta(seconds=abs(args.range) * 86400) if args.range else delta # postive timedelta is easier to handle
      s = datetime.timedelta(seconds=(args.shift if args.shift else 0) * 86400)

    if r > delta:
      print('Warning: given range larger than data range')
    if s > delta:  
     print('Warning: given shift larger than data range')
    
    miv += s
    start = miv            if (args.range == None or args.range > 0) else max(mav-r,miv)
    end   = min(miv+r,mav) if (args.range == None or args.range > 0) else mav
    delta = end - start
        
    for i, t in enumerate(x):
      si, ei = bounds(t, start, end)
      # out of data range returns 0,0 and end index 0 cannot display data
      print('actual restriction', (t[si] if ei != 0 else '-'),'to',(t[ei-1] if ei != 0 else '-'),'which are',(ei-si),'/',len(x[i]),'datapoints:',args.input[i])
      start_idx[i] = si 
      end_idx[i] = ei
  elif args.irange or args.ishift:
    for i, t in enumerate(x):
      si = args.ishift if args.shift else 0
      ei = end_idx[i]
      start_idx[i] = si                      if (args.range == None or args.range > 0) else max(ei-args.range, si)
      end_idx[i]   = min(si+args.irange, ei) if (args.range == None or args.range > 0) else ei
  else:
    print('plot all',mal,'datapoints') 
        
  # restrict the data - actually meant for 2D plots, see what happens for 3D     
  for i in range(len(x)):
    x[i] = x[i][start_idx[i]:end_idx[i]]
    
  for i in range(len(y)):
    yfactor = 100/(y[i][0])
    idx = abs(fiy[i])-1 # 1-based and +/- to encode bar
    y[i] = y[i][start_idx[idx]:end_idx[idx]]
    if args.yscale == 'percentage':
      for j in range(len(y[i])):
        y[i][j] *= yfactor
    # adjust to value above final value for relative logarithmic scale 'logrel'
    if args.yscale == 'logrel':
      for j in range(len(y[i])):
        y[i][j] = y[i][j] - y[i][-1]
  if args.yscale == 'percentage':
    args.yscale = 'linear'
  if args.yscale == 'logrel':
    args.yscale = 'log'

  for i in range(len(y2)):
    yfactor = 100/(y2[i][0])
    idx = abs(fiy2[i])-1
    y2[i] = y2[i][start_idx[idx]:end_idx[idx]]
    if args.yscale == 'percentage':
      for j in range(len(y[i])):
        y[i][j] *= yfactor
    # adjust to value above final value for relative logarithmic scale 'logrel'
    if args.y2scale == 'logrel':
      for j in range(len(y[i])):
        y[i][j] = y[i][j] - y[i][-1]
  if args.yscale == 'percentage':
    args.yscale = 'linear'
  if args.y2scale == 'logrel':
    args.y2scale = 'log'
  # finished with restrictions



  # now we can do smooth and grad, identified by the predix smooth_ and grad_ in the label
  y  = apply_smooth(y, ylabel, args.smooth_window, args.smooth_poly) 
  y2 = apply_smooth(y2, y2lbl, args.smooth_window, args.smooth_poly)
  z  = apply_smooth(z, zlabel, args.smooth_window, args.smooth_poly) 
  y  = apply_grad(fiy, y, ylabel, x)
  y2 = apply_grad(fiy2, y2, y2lbl, x)
  z  = apply_grad(fiz, z, zlabel, x)
  
  if args.dashed:      
    plt.rc('axes', prop_cycle=c_cms_y)

  # now plot the stuff on potentially in 1D by the y-axis or in 2D/3D (3d=warped)
  fig = None
  ax = None
  if not args.z:
    fig, ax = plt.subplots() 
    lines = []

    for i in range(len(y)):
      # y is a list of data. if fiy we know the current file index and take the x with the proper file index for you y columns
      # fiy(2) is 1-based and encodes bar with a negative value
      if fiy[i] > 0:
        lines.append(ax.plot(x[fiy[i]-1],y[i], color=colors[i], marker=args.marker if args.marker else next(markercycler), linestyle=args.linestyle if args.linestyle else next(linecycler))[0]) # returns multiple results and we want only the first
      else:
        lines.append(ax.bar(x[-fiy[i]-1],y[i], width=args.barwidth, color=colors[i], linestyle=args.linestyle if args.linestyle else next(linecycler))) # has only one return
    if args.ylim:
      ax.set_ylim(args.ylim)

    if args.y2:
      if args.dashed:      
        plt.rc('axes', prop_cycle=c_cms_y2)
      ax2 = ax.twinx()
      for i in range(len(y2)):
        if fiy2[i] > 0:
          lines.append(ax2.plot(x[fiy2[i]-1],y2[i], color=colors[13+i], marker=args.marker if args.marker else next(markercycler), linestyle=args.linestyle)[0]) # start with gold
        else:
          lines.append(lines.append(ax2.bar(x[-fiy2[i]-1],y2[i], width=args.barwidth, color=colors[13+i], linestyle=args.linestyle)))
      if args.y2lim:
        ax2.set_ylim(args.y2lim)

    labels = fix_labels(ylabel, fiy, args.input) + fix_labels(y2lbl, fiy2, args.input)
    if args.legend:
      if len(args.legend) > len(labels):
        print('Error: more entries given with --legend', len(args.legend), ' than lines', len(labels))
        sys.exit()
      labels[0:len(args.legend)] = args.legend
    plt.legend(lines, labels, loc=args.legend_loc, ncol=args.legend_ncol)

  else: # here comes the z-case
    # https://towardsdatascience.com/an-easy-introduction-to-3d-plotting-with-matplotlib-801561999725
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    #print(len(x), len(y), len(y2), len(z))
    if not (len(x) == 1 and len(y) == 1 and len(y2) == 0 and len(z) > 0):
      print('Error: 3D plots require one -x and -y, no -y2 and at least one -z')
      sys.exit(1)
    # we span the spaces
    x0 = np.sort(np.unique(x[0]))
    y0 = np.sort(np.unique(y[0]))
    
    if len(x0) * len(y0) != len(z[0]):
      print('Warning: Data not fit for 3D: |x|=',len(x0),'(' + str(len(x[0])) + ')','|y|=',len(y0),'(' + str(len(x[0])) + ')','z=',len(z[0]),'should be',len(x0)*len(y0),'fill with 0.0')
 
    assert len(x[0]) == len(y[0]) == len(z[0])
    Z = np.zeros((len(y0),len(x0),len(z))) # often len(z) is 1  
    for i, yv in enumerate(y0):
      for j, xv in enumerate(x0):
        # quite expensive search, note that x and y are of size x0*y0
        zv = None
        for k in range(len(x[0])):
          if x[0][k] == xv and y[0][k] == yv:
            if zv is not None:
              print('Error: data pair x=',xv,'y=',yv,'not unique:',zv,z[0][k]) 
            zv = [zi[k] for zi in z]
            #zv = z[0][k]
        if zv == None:
          print('miss',xv,yv)
        else:
          Z[i,j] = zv
    
    X, Y = np.meshgrid(x0,y0)       
    if len(z) == 1: # closed surfase for one value
      ax.plot_surface(X, Y, Z[:,:,0], rstride=1, cstride=1, cmap='jet', edgecolor='black')
    else: # grid for more data
      for i in range(len(z)):
        ax.plot_wireframe(X, Y, Z[:,:,i],color=colors[i])
    
    if args.save and '.vtr' in args.save:
      from pyevtk.hl import gridToVTK 
      pd = {}
      for i in range(len(z)): 
        pd[zlabel[i]] = np.atleast_3d(Z[:,:,i]).copy() # we need to copy to prevent assert (data.flags['C_CONTIGUOUS'] or data.flags['F_CONTIGUOUS'])
      gridToVTK(args.save[:-4], x0, y0, np.zeros(1), pointData = pd)  
      print('wrote',args.save)
    
    
  # common stuff for 2D and 3D
  ax.set_xscale(args.xscale)  
  ax.set_yscale(args.yscale)
  # ax.ticklabel_format(useOffset=False)  causes AttributeError: This method only works with the ScalarFormatter
  if args.y2:
    ax2.set_yscale(args.y2scale)
    # ax2.ticklabel_format(useOffset=False) causes AttributeError: This method only works with the ScalarFormatter
 
  # when the timespan is too short, we skip the day information squeezed in by matplotlib
  if has_dt:
    if delta.days < 2: 
      ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
  else:
    # name integer format when we assue iterations or such
    if args.xscale == 'linear' and abs(x[0][-1]-x[0][0]) > 5:
      ax.xaxis.set_major_locator(MaxNLocator(integer=True))

  if args.xlabel:
    ax.set_xlabel(args.xlabel)
  elif xlabel[0] != '':
    ax.set_xlabel(xlabel[0]) 
  
  # beautify the x-labels -> write timestamps diagonal
  plt.gcf().autofmt_xdate()
   
  # save for None
  ax.set_ylabel(label(args.ylabel,ylabel))
  if args.y2:
    ax2.set_ylabel(label(args.y2label,y2lbl))

  if args.z:
    ax.set_zlabel(label(args.zlabel,zlabel))
    ax.set_zscale(args.zscale) 

  plt.title(label(args.title,args.input))
 
  if args.save and not '.vtr' in args.save:
    print("write image to '" + args.save + "'")
    plt.savefig(args.save)
    
  if not args.noshow:
    #print('show ' + str(len(x[0])) + ' of ' + str(len(data[0])) + ' datapoints')
    plt.show()

# here could be an else case for the import plotviz part   
