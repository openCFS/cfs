#!/usr/bin/env python
import argparse
import sys
import os.path
import yaml
from cfs_utils import *


## this collects the information from a timer
# the data is label, wall and cpu. Extend for calls when needed
# the root Timer has the label 'total'. wall and cpu are float
class Timer:
  def __init__(self, label, wall, cpu, counter = -1, sub = False, id = None, parent_id = None):
    self.id = id
    self.label = label
    self.wall = wall
    self.cpu = cpu
    self.cnt = counter
    self.sub = sub
    self.parent_id = parent_id
    self.speed = cpu / wall if wall >= 0.02 else None
    self.children = [] # set in order() with parent
    self.parent = None
    self.wall_std = None 
    self.cpu_std = None

  # prefix is '' for independent, '* ' when no parent and '* -' or '* --' for children
  def prefix(self):
    if not self.sub:
      return ''
    if self.parent is None:
      return '* '
    if self.parent.parent is None:
      return '* -'
    return '* --'
  
# helper to get extract a timer from an xml node
# extracts the data from a node. If no id attribute is given tries to be smart
def set_timer(node):  
  # id and parent for modern cfs
  id = int(node.attrib['id']) if 'id' in node.attrib else -1
  parent_id = int(node.attrib['parent']) if 'parent' in node.attrib else None
  w = float(node.attrib['wall'])
  c = float(node.attrib['cpu'])
  i = int(node.attrib['calls'])
  # the label is more complicated. Older info.xml might not have label yet
  l = node.tag
  if l == 'timer':
    l = node.attrib['label'] if 'label' in node.attrib else None
    if l == None: # do our best
      l = node.getparent().tag if node.getparent().tag != "summary" else node.getparent().getparent().tag
      # is this the root node?
      if l == "cfsInfo":
        l = 'total'
        i = -1 # skip counter 

  # the optional attribute sub="true" indicates that this is sub-element 
  # and shall not be considered for missing_time
  # not any sub-element has parent set
  s = node.attrib['sub'] if 'sub' in node.attrib else None
  return Timer(l, w, c, i, s == 'true', id, parent_id)      
          
## structure timers to include sub-timers and sort as requested. Exclude timer with cnt==0
def order(timers, sort_by_id=False):
  if not sort_by_id:
    timers.sort(key=lambda t: t.wall, reverse=True)
  if sum([t.id == -1 for t in timers]) > 1: # root is allowed to have -1
    return timers

  # make a map by id for easy access
  map = {}
  for t in timers:
    map[t.id] = t

  main = [] # timer of first level = independent timers + parentless sub-timers
  for t in timers:
    if t.cnt == 0:
      continue
    if not t.sub or t.parent_id == None or t.parent_id not in map:
      main.append(t)
    else:
      # we are a sub-timer with known parent, add to tree structure
      map[t.parent_id].children.append(t)  
      t.parent = map[t.parent_id]

  if sort_by_id:  
    main.sort(key=lambda t: t.id)
     
  def add(lst, res):
    for t in lst:
      res.append(t)
      if t.children:
        add(t.children, res)

  # flatten timers
  flat = []
  add(main, flat)
  return flat
          
## extracts all timers from info.xml and give back as array of Timer objects
#@param gap shall a timer added as the last one which contains the gap from the first 
#           minus the sum of the rest
#@param what 'timers' (default) or 'memory': selects which metrics to extract
def read_info(xml, gap = False, sort_by_id=False, what='timers'):
  if what == 'memory':
    return read_memory(xml)
  res = []
  all = xml.xpath("//*[contains(local-name(),'timer')]") # some stuff is renamed like snopt_timer
  for node in all:
    res.append(set_timer(node))  
  res = order(res, sort_by_id)

  if gap:
    wall = res[0].wall  
    cpu  = res[0].cpu
    for t in res[1:]: # skip first
      if not t.sub:  
        wall -= t.wall # hope we stay non-negative :)
        cpu -= t.cpu
    res.append(Timer('not_measured', wall, cpu))
  return res


## extract memory metrics from an info.xml
# When cfs was run with -d the per-matrix memory breakdown is present under
# <setup><matrices> and is extracted as one entry per matrix type plus a grand
# 'total' (MB, using the per-type totalMemory attribute). Without -d (no per-matrix
# info) the process memory <memory final peak> is returned instead as a fallback
# (final, peak and total=peak).
# The entry labelled 'total' is kept first, as write_timers() takes the first entry
# as the total the percentage column relates to.
def read_memory(xml):
  res = []

  # pick the <matrices> node holding the per-matrix breakdown (marked by a child
  # element carrying totalMemory); other <matrices> nodes only contain timers
  chosen = [m for m in xml.xpath("//*[local-name()='matrices']")
            if m.xpath("./*[@totalMemory]")]
  if chosen:
    grand = 0.0
    for typeNode in chosen[0]:
      total = typeNode.get('totalMemory')
      if total is None:
        continue
      total = float(total) * 1.0 # already in MB
      grand += total
      res.append(Timer(typeNode.tag.rsplit('}', 1)[-1], 0.0, total))
    res.insert(0, Timer('total', 0.0, grand))
  else:
    # fallback for non-detailed runs: process memory final/peak
    mem = xml.xpath("//*[local-name()='memory']")
    if mem:
      peak  = mem[0].get('peak')
      final = mem[0].get('final')
      peak  = float(peak)  if peak  is not None else 0.0
      final = float(final) if final is not None else 0.0
      res.append(Timer('total', 0.0, peak if peak else final))
      res.append(Timer('final', 0.0, final))
      res.append(Timer('peak', 0.0, peak))
  return res

## aggregate the values for each entry with a list of a list of Timer objects
# @param mode how to aggregate over the runs: 'min' (default), 'mean' or 'max', see --mode
# @return a timer object with the aggregated values. Note that gap is not recomputed
#         only for 'mean' wall_std/cpu_std are set, otherwise they stay None
def aggregate_timer(timers, mode='min'):
  # we have not an object with an info.xml timer status but each info.xml is a list of Timer

  # check if all timers have same length
  assert(all(len(t) == len(timers[0]) for t in timers)), 'number of timers in info.xml files differ'

  pick = min if mode == 'min' else max if mode == 'max' else lambda v: sum(v) / len(v)

  # sample standard deviation, needs at least two runs to be defined
  def std(v):
    if mode != 'mean' or len(v) < 2:
      return None
    mean = sum(v) / len(v)
    return (sum((x - mean)**2 for x in v) / (len(v) - 1))**0.5

  res = []
  item = len(timers[0])
  for e in range(item):
    l = timers[0][e].label
    wall = [t[e].wall for t in timers]
    cpu  = [t[e].cpu  for t in timers]
    t = Timer(l, pick(wall), pick(cpu))
    t.wall_std = std(wall)
    t.cpu_std  = std(cpu)
    res.append(t)
  return res

## create gnuplot output from timer
# timer list of Timer objects
def gnuplot(timer, header=True, out=sys.stdout):
  if header:
    line = '#'
    for i in range(len(timer)):
      line += '(' + str(i*2+1) + '):' + timer[i].label + ' \t(' + str(i*2+2) + ')<-cpu \t'
    print(line, file=out)

  line = ''
  for t in timer:
    line += str(t.wall) + ' \t' + str(t.cpu) + ' \t'
  print(line, file=out)

## write timer standard analysis
# @param eps show the deviation from the reference only above this relative value - default is 10%
def write_timers(input_file_path, timers, format='txt', brief=False, wall=True, cpu=True,
                 cnt=False, ref=None, threshold=0.0, mode='min', eps=0.1, out=sys.stdout):
  ## shared preprocessing: aggregate multiple runs into a single timer set
  if type(timers[0]) == list and len(timers) > 1:
    timer = [aggregate_timer(timers, mode)]
    meta = len(timers)
  else:
    timer = timers
    meta = 0

  ## ---------------------------- yaml ----------------------------
  if format == 'yaml':
    # to write the timer structure into the yaml file
    struct = timers[0]
    # write a timer only above --threshold
    apply_threshold = [(t, s) for t, s in zip(timer[0], struct) if not (t.wall < threshold and t.wall >= 0)]
    # labels are not unique, so a list keeps all timers and preserves their order
    if mode == 'mean':
      data = {'runs': 1 if meta == 0 else meta,
              'files': [os.path.basename(f) for f in input_file_path],
              'analysis': {'mode': mode, 'threshold': threshold},
              'timers': [{'label': t.label, 'id': s.id, 'parent': s.parent_id, 'sub': s.sub, 'calls': s.cnt,
                          'wall': t.wall, 'wall_std': t.wall_std, 'cpu': t.cpu, 'cpu_std': t.cpu_std}
                         for t, s in apply_threshold]}
    else:
      data = {'runs': 1 if meta == 0 else meta,
              'files': [os.path.basename(f) for f in input_file_path],
              'analysis': {'mode': mode, 'threshold': threshold},
              'timers': [{'label': t.label, 'id': s.id, 'parent': s.parent_id, 'sub': s.sub, 'calls': s.cnt,
                          'wall': t.wall, 'cpu': t.cpu}
                         for t, s in apply_threshold]}

    if out is sys.stdout:
        yaml.dump(data, sys.stdout, default_flow_style=False, sort_keys=False)
    else:
        yaml_file = os.path.join(os.getcwd(), out + '.yaml')
        with open(yaml_file, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    return 0

  ## ---------------------------- txt ----------------------------
  if out is not sys.stdout:
    out = os.path.join(os.getcwd(), out + '.txt')
    f = open(out, 'w')
  else:
    f = out

  # header
  tag = {'min': 'MIN', 'mean': 'AVG', 'max': 'MAX'}[mode] # all three chars wide to keep the columns aligned
  title = 'TIMER (sec, {})'.format(tag) if meta else 'TIMER (sec)'
  max_label = max(max([len(t.label) for t in timer[0]]) + 3, len(title) + 1) # '*--' prefix; never shorter than the title
  speed_head = ' : C/W' if not meta else '' # cpu/wall, aggregated runs have none
  head = title.ljust(max_label) + ':'
  head += ' cnt :' if cnt else ''
  if brief:
    if not ref:
      head += '   WALL' if wall else ''
      head += ' ~ ' if wall and cpu else ''
      head += '  CPU  ' if cpu else ''
      head += speed_head
      for m in range(meta):
        head += '  | ' if m == 0 else '  : '
        head += 'WALL_{:d}'.format(m) if wall else ''
        head += ' ~ ' if wall and cpu else ''
        head += ' CPU_{:d}'.format(m) if cpu else ''
    else:
      head += tag + '_CPU' + speed_head + '  |     REF' + '    REL DIFF'
  else:
    if not ref:
      head += '______WALL______' if wall else ''
      head += ' ~ ' if wall and cpu else ''
      head += '______CPU_______' if cpu else ''
      head += speed_head
      for m in range(meta):
        head += '  |' if m == 0 else '  :'
        head += ' WALL_{:d} '.format(m) if wall else ''
        head += '   ~' if wall and cpu else ''
        head += '  CPU_{:d} '.format(m) if cpu else ''
        head += '  '
    else:
      head += '  ___' + tag + '_CPU____' + speed_head + '  |  ___CPU_REF____ ' + ' ___REL DIFF___'
  print(head, file=f)

  total_wall = max(timer[0][0].wall, 1e-3)
  total_cpu  = max(timer[0][0].cpu, 1e-3)
  if ref:
    total_cpu_ref  = max(ref[0].cpu, 1e-3)

  # timer
  for e,t in enumerate(timer[0]):
    if t.wall < threshold and t.wall >= 0:
      continue
    # format for time display
    format_wall = '{: 6.0f}' if t.wall >= 10000 else '{: 7.3f}'
    format_cpu = '{: 6.0f}' if t.cpu >= 10000 else '{: 7.3f}'

    prefix = t.prefix()
    l = prefix + t.label

    line = l.ljust(max_label)
    line += ':'
    if cnt:
      if t.cnt >= 0:
        line += '{:=4d} :'.format(t.cnt)
      else:
        line += '     :'.format(t.cnt)
    if wall and not ref:
      # wall time
      line += format_wall.format(t.wall)
      if not brief:
        line += ' [{:.1%}]'.format(t.wall/total_wall).rjust(9)
    line += ' ~ ' if wall and cpu else ''
    if cpu:
      # cpu time
      line += format_cpu.format(t.cpu)
      if not brief:
        line += '[{:.1%}]'.format(t.cpu/total_cpu).rjust(9)

    if not meta:
      line += ' : {:.1f}'.format(t.speed) if t.speed and t.speed >= 1 else ' ' * 6

    rel_eps = 0.1
    if not ref:
      # iterate over all timers
      for m in range(meta):
        line += '  | ' if m == 0 else '  : '
        if wall:
          line += format_wall.format(timers[m][e].wall)
          if not brief:
            # show change if influence on total time is > 2% (min only)
            if mode == 'min' and t.wall > total_wall * 0.02:
              diff = rel_diff(t.wall, timers[m][e].wall)
              line += '(+)' if diff > rel_eps else '(-)' if diff < -rel_eps else '   '
            else:
              line += '   '
        line += ' ~ ' if wall and cpu else ''
        if cpu:
          line += format_cpu.format(timers[m][e].cpu)
          if not brief:
            # show change if influence on total time is > 5% (min only)
            if mode == 'min' and t.cpu > total_cpu * 0.03:
              diff = rel_diff(t.cpu, timers[m][e].cpu)
              line += '(+)' if diff > rel_eps else '(-)' if diff < -rel_eps else '   '
            else:
              line += '   '
    else:
      line += '  | ' + format_cpu.format(ref[e].cpu)
      if not brief:
        line += ' [{:.1%}]'.format(ref[e].cpu/total_cpu_ref).rjust(9)
      # deviation from the reference, positive means slower - shown only above eps. Read it
      # against the 'total' row to tell a timer which really behaves differently from one
      # which just went along with a globally slower or faster run
      diff = rel_diff(ref[e].cpu, t.cpu)
      line += '{:+.1%}'.format(diff).rjust(9) if abs(diff) > eps else ' ' * 9
    print(line.rstrip(), file=f)

  if f is not sys.stdout:  # close written txt file
    f.close()
  return 0


## write memory metrics extracted from info.xml files.
# For txt each run is one column, with rows aligned by label (different runs
# may have different sets of matrix types); for yaml (which feeds 'compare' a single
# timers list) the first run is written and the number of runs is recorded.
# @param input_file_path the info.xml files the values were read from
# @param timers list of lists of Timer records, cpu holds the value in MB
def write_memory(input_file_path, timers, format='txt', out=sys.stdout):
  ## ---------------------------- yaml ----------------------------
  if format == 'yaml':
    struct = timers[0]
    data = {'runs': len(timers),
            'files': [os.path.basename(f) for f in input_file_path],
            'analysis': {'what': 'memory', 'unit': 'MB', 'mode': None, 'threshold': 0.0},
            'timers': [{'label': t.label, 'id': t.id, 'parent': t.parent_id, 'sub': t.sub,
                        'calls': t.cnt, 'wall': t.wall, 'cpu': t.cpu}
                       for t in struct]}
    if out is sys.stdout:
      yaml.dump(data, sys.stdout, default_flow_style=False, sort_keys=False)
    else:
      with open(os.path.join(os.getcwd(), out + '.yaml'), 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
    return 0

  ## ---------------------------- txt ----------------------------
  if out is not sys.stdout:
    out = os.path.join(os.getcwd(), out + '.txt')
  f = open(out, 'w') if out is not sys.stdout else out

  # per-run value map and the union of labels (run0 order first, then new ones)
  per_run = [{t.label: t.cpu for t in run} for run in timers]
  labels = []
  for run in timers:
    for t in run:
      if t.label not in labels:
        labels.append(t.label)

  total0 = max(per_run[0].get('total', 0.0), 1e-3)
  title = 'MEMORY (MB)'
  if len(timers) > 1:
    title += '  (one column per run)'
  max_label = max(max([len(l) for l in labels]) + 3, len(title) + 1)

  head = title.ljust(max_label)
  for m in range(len(timers)):
    head += '  | ' if m == 0 else '  : '
    head += 'MB_{:d}'.format(m)
  print(head, file=f)

  for lbl in labels:
    line = lbl.ljust(max_label) + ':'
    for m in range(len(timers)):
      v = per_run[m].get(lbl)
      line += '  | ' if m == 0 else '  : '
      line += ' {:8.3f}'.format(v) if v is not None else '       -'
    line += '  [{:.1%}]'.format(per_run[0].get(lbl, 0.0) / total0).rjust(9)
    print(line, file=f)

  if f is not sys.stdout:
    f.close()
  return 0

        
## read analysed timers in the yaml format written by write_timers,
## either from a yaml file or - if yaml_file is None - from stdin.
# @return (list of Timer objects, as read_info() gives it, analysis dict)
def read_timers_yaml(yaml_file=None):
  if yaml_file is None:
    data = yaml.safe_load(sys.stdin)
    name = '<stdin>'
  else:
    with open(yaml_file) as f:
      data = yaml.safe_load(f)
    name = yaml_file
  if not isinstance(data, dict) or 'timers' not in data:
    raise RuntimeError("no 'timers' in " + name)
 
  res = []
  for t in data['timers']:
    timer = Timer(t['label'], t['wall'], t['cpu'], t.get('calls', -1), t.get('sub', False), t.get('id'), t.get('parent'))
    timer.wall_std = t.get('wall_std')
    timer.cpu_std  = t.get('cpu_std')
    res.append(timer)
 
  # restore the tree which order() built, such that prefix() gives the sub-timer level again
  map = {t.id: t for t in res if t.id is not None}
  for t in res:
    if t.parent_id is not None and t.parent_id in map:
      t.parent = map[t.parent_id]
      t.parent.children.append(t)
  return res, data.get('analysis', {})


# searches in list of timers for timers which have not the name 'timer' and are as such ignored
def check_invalid_timers(infos):
  for info in infos:
    elem = open_xml(info).xpath("//*[@cpu]") # sum stuff is renamed like snopt_timer
    for e in elem:
      if e.tag != 'timer':
        print("Warning: ignored element '" + str(e.tag) + "': ",info)

# check the timers required to pass against the reference, the deviations themselves are
# shown by the DIFF column of write_timers()
# @param eps compares relative error against eps (optional, default 10%)
# @param labels only check the timers with these labels, None (default) checks all of them
# @return the tuple (slower, checked): slower is True if one of the checked timers is more
#         than eps slower than the reference, checked holds a (label, deviation) pair for
#         every timer which took part in the comparison
def compare_timers(timer_ref, timer, eps=0.1, skip_noise=None, labels=None):
  assert(len(timer_ref) == len(timer))

  if labels:
    # a typo in --timer must not compare nothing and silently report 'times are good'
    known = set([t.label for t in timer])
    unknown = [l for l in labels if l not in known]
    if unknown:
      print(' * error: no timer labelled ' + ', '.join(unknown) + " - known are: " + ', '.join(sorted(known)))
      return True, []

  slower = False
  checked = []
  # decides if comparison fails - if timers from --timer are slower
  for time_ref, time in zip(timer_ref, timer):
    if time.label == 'not_measured':
      continue

    if labels and time.label not in labels:
      continue

    if skip_noise and time_ref.cpu < skip_noise and time.cpu < skip_noise:
      continue

    diff = rel_diff(time_ref.cpu, time.cpu)
    checked.append((time.label, diff))
    if diff > eps:
      slower = True

  return slower, checked

      
## --wall shows only the wall times, --cpu only the cpu times, none of them both
# @return the tuple (wall, cpu) as print_timer() expects it
def wall_cpu(args):
  return (False if args.cpu else True, False if args.wall else True)

## read a list of info.xml files into a list of a list of Timer objects
def read_infos(files, sort_by_id, what='timers'):
  return [read_info(open_xml(f), gap=True, sort_by_id=sort_by_id, what=what) for f in files]
  
## check if a yaml file is inputed
def is_yaml(f):
  return os.path.basename(f).lower().endswith(('.yaml', '.yml'))

## relative deviation of a cpu time from the reference cpu time, positive means slower.
# This is what the DIFF column shows and what compare_timers() decides upon, hence a timer
# is reported as failed exactly when the table gives it a positive value above eps
def rel_diff(cpu_ref, cpu):
  return (cpu - cpu_ref)/cpu_ref if cpu_ref != 0.0 else 0.0

## 'analyse': read the timers (or memory) of the given info.xml files, aggregate and output them
def cmd_analyse(args):
  # gnuplot can only be written to stdout
  if args.format == 'gnuplot' and args.output != 'stdout':
      main.error("gnuplot format can only be written to stdout, not to a file")

  wall, cpu = wall_cpu(args)
  timers = read_infos(args.info, args.appearance, what=args.what)
  out = sys.stdout if args.output == 'stdout' else args.output

  if args.what == 'memory':
    if args.format == 'gnuplot':
      main.error("gnuplot format is not supported for memory metrics")
    return write_memory(args.info, timers, format=args.format, out=out)

  if args.format == 'gnuplot':
    gnuplot(aggregate_timer(timers, args.mode), out=sys.stdout)
  else:
    write_timers(args.info, timers, format=args.format, brief=args.brief,
                wall=wall, cpu=cpu, cnt=True, ref=None,
                threshold=args.threshold, mode=args.mode, out=out)
  if args.format != 'yaml' and out is not sys.stdout:
    check_invalid_timers(args.info)
  return 0


def cmd_compare(args):
  # read the test run: yaml file, info xml file, or - if nothing was given - yaml from stdin
  test = args.test  # nargs='?' gives a string or None
  analysis_vars = None
  if test is None:
    if sys.stdin.isatty():
      compare.error('no TEST_INPUT given and nothing piped to stdin')
    timer, analysis_vars = read_timers_yaml()  # reads stdin
    test = '<stdin>'
  elif is_yaml(test):
    timer, analysis_vars = read_timers_yaml(test)
  else:
    timer = read_info(open_xml(test), gap=True, sort_by_id=args.appearance)
 
  # read the reference: yaml file or info xml file
  analysis_vars_ref = None
  if is_yaml(args.ref):
    timer_ref, analysis_vars_ref = read_timers_yaml(args.ref)
  else:
    timer_ref = read_info(open_xml(args.ref), gap=True, sort_by_id=args.appearance)
 
  # check for same analysis types (only needed if both sides came from yaml)
  if analysis_vars is not None and analysis_vars_ref is not None:
    if analysis_vars_ref['mode'] != analysis_vars['mode']:
      raise Exception('The reference and test runs have not been analysed with the same mode: ref: {}; test: {}'.format(analysis_vars_ref['mode'], analysis_vars['mode']))
    elif analysis_vars_ref['threshold'] != analysis_vars['threshold']:
      raise Exception('The reference and test runs have used different thresholds during analysis: ref: {}; test: {}'.format(analysis_vars_ref['threshold'], analysis_vars['threshold']))
 
  # perform comparison - the table shows the deviation of every timer, compare_timers() only
  # decides whether the timers required to pass are within eps
  write_timers([test], [timer], brief=args.brief, wall=False, cpu=True, cnt=False, ref=timer_ref,
                 eps=args.eps)
  slower, checked = compare_timers(timer_ref, timer, eps=args.eps, skip_noise=args.skip_noise, labels=args.timer)

  if not checked:
    if slower: # an unknown --timer label, compare_timers() already said which one
      return 1
    print(' * warning: no timer was compared, all of them are below --skip-noise of {:g} sec'.format(args.skip_noise))
    return 0

  # name what was compared, --timer defaults to all timers
  if args.timer is None:
    which = 'all timers'
  elif len(args.timer) == 1:
    which = 'the timer ' + args.timer[0]
  else:
    which = 'the timers ' + ', '.join(args.timer)

  # report the timer which came closest to the limit, it is the one the verdict rests on.
  # Naming it adds nothing when a single label was requested, it can only be that one
  label, diff = max(checked, key=lambda c: c[1])
  most = '' if len(checked) == 1 else ' max.'
  worst = '' if args.timer is not None and len(args.timer) == 1 else ' ({:s})'.format(label)
  print(' * The cpu time of {:s} has changed by{:s} {:+.1%}{:s}'.format(which, most, diff, worst))

  if slower:
    # how many of the compared timers are responsible, the ERROR lines above name them
    failed = len([c for c in checked if c[1] > args.eps])
    count = '' if len(checked) == 1 else ', {:d} of {:d} compared timers are too slow'.format(failed, len(checked))
    print(' * This exceeds the limit of +{:.1%}{:s}'.format(args.eps, count))
    print(' * The comparison failed, a regression has been detected.')
    return 1

  print(' * This is below the limit of +{:.1%}'.format(args.eps))
  print(' * The comparison passed and no regression has been detected.')
  return 0

# ----------------------------------------------------------------
# analysing the info.xml files is the default and needs no sub-command, 'compare' is the
# only sub-command. Options which both share are collected in parsers of their own and
# handed over via parents=[]
display = argparse.ArgumentParser(add_help=False)
display.add_argument('--brief', help="brief analysis output, for cmd(stdout) and txt", action='store_true')
display.add_argument('--threshold', help="show only wall time above this seconds", type=float, default=0.0)
display.add_argument('--appearance', help="do not sort by wall time but by id", action='store_true')

main = argparse.ArgumentParser(prog='performance.py', parents=[display],
  formatter_class=argparse.RawDescriptionHelpFormatter,
  description="Analyse the timers of the info.xml files of a series of cfs performance runs.\n"
              "The timers are read, aggregated over the runs with --mode and are output. With more\n"
              "than one file the aggregated values come first, then one column per run.\n"
              "'not_measured' is the total minus the sum of the non-sub timers.",
  epilog="""examples:
  performance.py bracket-*.info.xml
  performance.py bracket-*.info.xml --mode mean --threshold 1
  performance.py bracket-*.info.xml --brief --cpu 
  performance.py bracket-*.info.xml --output bracket_analysis --format yaml --mode mean 
  performance.py bracket-*.info.xml --format yaml --mode mean 
  # memory requires a cfs run with -d for the per-matrix breakdown; without -d
  # only the process memory final/peak is available.
  performance.py detail.info.xml --what memory
  performance.py detail.info.xml --what memory -f yaml | performance.py compare --ref ref_mem.yaml

See 'performance.py compare -h' to compare a run or aggregated runs against a reference instead""")
main.add_argument('info', help="the info.xml file(s) to analyse, aggregated by --mode", nargs='+', metavar='INFO_XML')
main.add_argument('-o', '--output', help="default is stdout, otherwise write the desired filename", default='stdout')
main.add_argument('-f', '--format', help="define the output format, default is txt", choices=('txt', 'yaml', 'gnuplot'), default='txt')
main.add_argument('--what', help="which metrics to extract, default is timers", choices=('timers', 'memory'), default='timers')
times = main.add_mutually_exclusive_group() # --wall and --cpu exclude each other
times.add_argument('--wall', help="show only wall times", action='store_true')
times.add_argument('--cpu', help="show only cpu times", action='store_true')
main.add_argument('--mode', help="how to aggregate the timers of multiple runs - default is min", choices=('min', 'mean', 'max'), default='min')
main.set_defaults(func=cmd_analyse)

# the 'compare' sub-command. argparse cannot mix an optional sub-command with the INFO_XML
# positional of the default analysis, hence it is a parser of its own, dispatched below
compare = argparse.ArgumentParser(prog='performance.py compare', parents=[display],
  formatter_class=argparse.RawDescriptionHelpFormatter,
  description="Compare the timers of the given test runs against a reference.\n"
              "One can compare a specific timers or all timers\n"
              "Prints the table when a timer deviates by more than --eps.",
  epilog="""examples:
  performance.py compare --ref reference.info.xml bracket.info.xml
  performance.py compare --ref reference.yaml new-*.yaml --brief
  performance.py compare --ref reference.info.xml new-*.info.xml --eps 0.2 --skip-noise 0.5
  performance.py compare --ref reference.info.xml new-*.info.xml --timer total
  
  Analysis and compare can be piped, f.e.: 
  performance.py info1.xml info2.xml --format yaml --mode min | performance.py compare --ref agg_timers.yaml""")
compare.add_argument('test', help=" test run(s): takes either a single info xml or already aggregated results in the form of a single .yaml file. Optional, when info.xml files precede 'compare' or the analysis is piped in via stdin", nargs='?', metavar='TEST_INPUT')
compare.add_argument('--ref', help="reference run(s): takes either a single info xml or already aggregated results in the form of a single .yaml file", required=True)
compare.add_argument('--eps', help="allowed relative deviation from the reference - default is 0.1 = 10%%", type=float, default=0.1)
compare.add_argument('--skip-noise', '--skip_noise', help="suppress too small times", type=float, default=1e-1, dest='skip_noise')
compare.add_argument('--timer', help="labels of timers required to be below eps to pass - default: all", nargs='+', metavar='LABEL')
# compare always shows both times and aggregates with 'min', so it has no --wall/--cpu/--mode.
# cmd_compare() reads them nevertheless, hence they are set here
compare.set_defaults(func=cmd_compare, wall=False, cpu=False, mode='min')

if __name__ == '__main__':  # guard, such that measure.py can import from here
  argv = sys.argv[1:]
  if 'compare' in argv:
    if argv.index('compare') > 0:
      compare.error("analyse and compare are separate invocations, pipe them: "
                    "performance.py RUNS*.xml -f yaml | performance.py compare --ref REF")
    args = compare.parse_args(argv[1:])
  else:
    args = main.parse_args(argv)
  sys.exit(args.func(args))