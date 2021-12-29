#!/usr/bin/env python

# the purpose of the tool is to provide easy matplotlib visualization for dat files like .plot.dat and .grad.dat
# it extracts iformation from a header line with like #(1) iter \t(2) compliance ....
import os.path
import sys

import datetime
import numpy as np

# when we use plotviz from postrpoc.py, we don't need this stuff
if __name__ == '__main__':
  import argparse
  import matplotlib
  import matplotlib.pyplot as plt
  from matplotlib.ticker import MaxNLocator
  import snopt # our snopt.py helper for process



# having two y2-axis we need to handle colors manually, otherwise they repeat
# https://matplotlib.org/stable/gallery/color/named_colors.html
# 'gold' has index 13 for y2 axis
colors = ['tab:green','tab:red','tab:purple','tab:blue','tab:orange','black','tab:brown','tab:gray','tab:olive','blue','tab:cyan','tab:pink','cornflowerblue', 
          'gold','peru','blueviolet', 'coral','yellowgreen'] 
colors += colors # repeat such that it should be really enough

# parses the header lines for the first hint on column names.
# Tries to be smart!!
# @param data the processed content() result. To fill up missing hints in the comment
# @param comments a list of lines which start with a hashtag
def header(data, comments):

  meta = []
  
  for l in reversed(comments):
    # we assume the last comment contains the labels. ## are ignored.
    if l.startswith('##') or l.startswith('--'):
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

  if len(meta) > len(data[0]):
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
      print('{:3d}'.format(cnt) + ': ' + l.ljust(ml) + ' : ' + str(data[mi][0][li]).ljust(19) + ' : ' + inputs[mi])
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
    for frmt in ['%Y-%m-%d %H:%M:%S', '%d.%m.%Y %H:%M', '%Y/%m/%d %H:%M:%S','%d.%m.%y', '%d.%m.%Y', '%Y-%m-%d']:
      if check(frmt, str(data[0][c])):
        for l in data:
          l[c] = datetime.datetime.strptime(l[c], frmt)


  # convert str to float
  for i, c in enumerate(data[0]):
    if type(c) == str:
      for r in range(len(data)):
        try:
          if data[r][i].endswith('%'): # remove trailing % from bluetooth sensor for rel humidity 
            data[r][i] = data[r][i][:-1]
          data[r][i] = 0 if data[r][i] == '' else float(data[r][i])
        except ValueError as ve:
          pass
  #print(data)
  return data

# validate datetime by try and except
def check(format, test):
  try:
    datetime.datetime.strptime(test, format)
    return True
  except ValueError:
    return False
  
# gives back the file index and the corresponding 0-based index.
# @param ky if string, search for uniqueness, if number, make double 0-based
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
    for f, m in enumerate(meta):
      for i, t in enumerate(m):
        if t.startswith(key):
          if idx > -1:
            print("Error: key not unique '", key, '"')
            sys.exit()
          idx = i
          fi = f
  
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
def column(meta, data, key,bars):
  if key == None:
    # for default if -x is not given. Shall not be callend for y2 is None
    n = len(data)
    return [i for i in range(n)] , [range(len(data[i])) for i in range(n)], [''] * n
  else:
    if type(key) != list:
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
 
  for l in lines:
    h = l.strip()
    if h.startswith('#') or h.startswith('iter') or h.startswith('Temp') or h.startswith('---'): 
      comments.append(h)
    else:
      body.append(l)  

  data = content(body)
  meta = header(data,comments)
  
  return meta, data        


# convenience function which gives a y(2)-axis label. 
# @param args either args.ylabel or args.y2label. If not None this is returned
# @oaram legend array of legends. Returns a list of unique entries
def label(args, legends):
  if args:
    return args
  else:
    # we make manually unique without list(set(legends)) to keep order
    res = []
    for l in legends:
      if not l in res:
        res.append(l)
    
    return ', '.join(res) 

# create artificial data for the numerical smooth keys if any
def append_smooth(input, meta, data, smooth, window, poly):
  if smooth is not None and len(smooth) > 0:
    # data is list of lines with all columns, convert first for our selected colums
    fix, cols, labels = column(meta,data,smooth,[])
    assert len(fix) == len(cols) == len(labels)
    
    import scipy.signal
    for i, c in enumerate(cols):
      w = min(window, 2*int(len(c)/2) +1) # window nees to be odd
      print("generate 'smooth_" + labels[i] + "' Savitzky–Golay filtered data with poly", poly, "and window",w)
      
      
      sc = scipy.signal.savgol_filter(c, w, poly)
      meta[fix[i]-1].append('smooth_' + labels[i])
      for vi, v in enumerate(sc):
        data[fix[i]-1][vi].append(v)

  return meta, data

# create artificial data for the finite difference data keys if any
# the gradient is are for inbetween data but mooved to the left and the last value is doubled
def append_grad(input, meta, data, grad):
  if grad is not None and len(grad) > 0:
    fix, cols, labels = column(meta,data,grad,[])
    assert len(fix) == len(cols) == len(labels)
    for i, c in enumerate(cols):
      dc = np.diff(c,n=1)
      dc = np.append(dc,dc[-1]) # repeat last value
      assert len(c) == len(dc)
      meta[fix[i]-1].append('grad_' + labels[i])
      for vi, v in enumerate(dc):
        data[fix[i]-1][vi].append(v)

  return meta, data



# plotviz.py is imported by postproc.py, so guard argparse
if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Simple gnuplot replacement for standard plots. Needs a header comment')
  parser.add_argument("input", nargs='+', help="one or two .plot.dat or similar tabular text files")
  parser.add_argument("-x", nargs='*',  help="index or label for the abscissa (optional). Space separated list for multiple inputs")
  parser.add_argument("-y", nargs='+',  help="indices or labels for the ordinate. Multiple separated by space")
  parser.add_argument("-y2", nargs='*', help="optional indices or labels for the secondary ordinate")
  parser.add_argument("-z", nargs='*',  help="trigger 3D plots wich requires a single -x and -y component")
  parser.add_argument("--range", help='range (day fractions for datetime) to be shown. negative for final range', type=float)
  parser.add_argument("--shift", help='shift range (day fractions for datetime)', type=float)
  parser.add_argument("--xlabel", help='optional label for the abscissa')
  parser.add_argument("--ylabel", help='optional label for the primary ordinate')
  parser.add_argument("--y2label", help='optional label for the secondary ordinate')
  parser.add_argument("--zlabel", help='optional label for the 3d data')
  parser.add_argument("--marker", help='optional matplotlib marker: e.g. . , o v')
  parser.add_argument("--legend", nargs='*', help="(partially) overwrite labels in the legend as space separated list of strings")
  parser.add_argument("--legend_loc", help="string for matplotlib.legend(loc)")
  parser.add_argument("--title", help='optional title for the plot')
  parser.add_argument("--yscale", help="scaling type from choice, google matplotlib yscale", choices=["linear", "log", "symlog", "logit"],default='linear')
  parser.add_argument("--y2scale", help="like --yscale but for y2 axis", choices=["linear", "log", "symlog", "logit"],default='linear')
  parser.add_argument("--bar", nargs='*', help="indices from y or y2 which are to displayed as bars instead of plots")
  parser.add_argument("--barwidth", help="barplots for datetime need manual adjustment", type=float, default=.8)
  parser.add_argument("--smooth", nargs='*', help="create new smoothed data for given fields")
  parser.add_argument("--smooth_window", help="window size of Savitzky–Golay filter", type=int, default = 7)
  parser.add_argument("--smooth_poly", help="polynomial order of Savitzky–Golay filter", type=int, default = 3)
  parser.add_argument("--grad", nargs='*', help="give finite difference gradients, you might want to smooth first")
  parser.add_argument("--save", help='write to given filename using the extension')
  parser.add_argument("--noshow", help='supress popping up the image window', action='store_true')
    
  args = parser.parse_args()
  
  # array of headers per file
  meta = [] 
  # matrix of data per file
  data = []
  
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
    else:    
      m, d = process(input)
    meta.append(m)
    data.append(d)  

  # if one of the arguments is a multiple key in meta, it is replaced by multiple ids (argument becomes larger)
  # we first do this for artificial data
  args.smooth = resolve_multiple(meta, args.smooth)
  meta, data = append_smooth(args.input, meta, data, args.smooth, args.smooth_window, args.smooth_poly)
  
  args.smooth = resolve_multiple(meta, args.grad)
  meta, data = append_grad(args.input, meta, data, args.grad)
  
  args.x = resolve_multiple(meta, args.x)
  args.y = resolve_multiple(meta, args.y)
  args.y2 = resolve_multiple(meta, args.y2)
  args.z = resolve_multiple(meta, args.z)



  print_header(meta, data, args.input)
   
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
  else:
    print('plot all',mal,'datapoints') 
        
  # restrict the data - actually meant for 2D plots, see what happens for 3D     
  for i in range(len(x)):
    x[i] = x[i][start_idx[i]:end_idx[i]]
    
  for i in range(len(y)):
    idx = abs(fiy[i])-1 # 1-based and +/- to encode bar
    y[i] = y[i][start_idx[idx]:end_idx[idx]]   

  for i in range(len(y2)):
    idx = abs(fiy2[i])-1
    y2[i] = y2[i][start_idx[idx]:end_idx[idx]]   

        
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
        lines.append(ax.plot(x[fiy[i]-1],y[i], color=colors[i], marker=args.marker)[0]) # returns multiple results and we want only the first
      else:
        lines.append(ax.bar(x[-fiy[i]-1],y[i], width=args.barwidth, color=colors[i])) # has only one return
      
    if args.y2:
      ax2 = ax.twinx()
      for i in range(len(y2)):
        if fiy2[i] > 0:
          lines.append(ax2.plot(x[fiy2[i]-1],y2[i], color=colors[13+i], marker=args.marker)[0]) # start with gold
        else:
          #print(-fiy2[i]-1,x[-fiy2[i]-1],y2[i])
          lines.append(lines.append(ax2.bar(x[-fiy2[i]-1],y2[i], width=args.barwidth, color=colors[13+i])))

    labels = fix_labels(ylabel, fiy, args.input) + fix_labels(y2lbl, fiy2, args.input)
    if args.legend:
      if len(args.legend) > len(labels):
        print('Error: more entries given with --legend', len(args.legend), ' than lines', len(labels))
        sys.exit()
      labels[0:len(args.legend)] = args.legend
    if args.legend_loc:
      if not len(args.legend_loc) == 2:
        print('Error: --legend_loc requires two values (e.g. 1.05 1) to form the bbox_to_anchor attribute to plt.legend()')
    plt.legend(lines, labels, loc=args.legend_loc)

  else: # here comes the z-case
    # https://towardsdatascience.com/an-easy-introduction-to-3d-plotting-with-matplotlib-801561999725
    fig = plt.figure()
    ax = plt.axes(projection="3d")
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
            zv = [z[i][k] for i in range(len(z))]
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
  ax.set_yscale(args.yscale)
  if args.y2:
    ax2.set_yscale(args.y2scale)
 
  # when the timespan is too short, we skip the day information squeezed in by matplotlib
  if has_dt:
    if delta.days < 2: 
      ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
  else:
    # name integer format when we assue iterations or such
    if abs(x[0][-1]-x[0][0]) > 5:
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

  plt.title(label(args.title,args.input))
 
  if args.save and not '.vtr' in args.save:
    print("write image to '" + args.save + "'")
    plt.savefig(args.save)
    
  if not args.noshow:
    #print('show ' + str(len(x[0])) + ' of ' + str(len(data[0])) + ' datapoints')
    plt.show()

# here could be an else case for the import plotviz part   
