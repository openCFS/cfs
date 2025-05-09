#!/usr/bin/env python
import argparse
import glob
import os
import numpy as np

# known lables we work with
known_labels = ['fast', 'slow', 'unstable', 'broken']

# parsed information form a log file by test
class Test:
  # @param testname from log e.g. Optimization_Stress_local_vts -> unclear which '_' is from name or '/'
  # @param path from dir like 'Technical/PythonMesher' 
  def __init__(self, testname, path, passed = None, time = None):
    # From log or properties: Optimization_Stress_local_vts where _ can be / or original 
    self.testname = testname 

    self.path = path  # Technical/PythonMesher find from the Testsuite dir
    assert path is None or len(path) > 2

    # properties from the log
    self.passed = passed 
    self.time = time

    # the CMakeLists.txt file name
    self.cmakelists = None
    # labels known_labels from the CMakeLists.txt. Makes not much sense if set_property > 1
    self.label = []
    # labels which are not within known_labels.
    self.unknown = []     
    
    # count in the path's CMakeLists.txt files. Set to >= 0 by read_cmake_lists() 
    self.add_test = -1
    self.set_property = -1
    
# parse the log file. The log file can contain almost everything, we filter for lines like
# 385/445 Test #385: Optimization_TwoScale_HomRect_C1_non_regular ..........................................................   Passed    1.70 sec
# @param logile filename 
# @param optionally all CMakeLists.txt is Testsuite, all found tests are removed from the list. Then path is set in tests
# @return a list of Test objects
def parse_log(files, logfile = None):
  tests = []
  if logfile:
    file = open(logfile,'r')
    for line in file.readlines():
      if '....' in line:
        assert line.count('/') == 1 and line.count('#') == 1 and line.count(': ') == 1, line
        name = line[line.find(': ')+2:line.find(' ...')] # e.g. 'Optimization_Stress_local_vts'
        
        # the name has no clear path info ('_' has double meaning), so search for it in files and reduce files
        idx = -1
        path = None # can be None e.g. for Coupledfield_MechSmooth_MechSmoothMWE2D_evalIter
        for i, f in enumerate(files):
          t = f.replace('_','/') # split for _ and /
          if name.split('_') == t.split('/'): 
            idx = i
            path = f
            break
        if idx != -1:
          del files[idx]
             
        end = line[line.rindex('...')+3:] # '   Passed    0.65 sec\n' or '***Failed    0.73 sec\n'
        if not ('Passed' in end or 'Failed' in end or 'Timeout' in end):
          print('error: no Passed/Failed/Timeout in line',line)
          os.sys.exit()
        assert end.count('ed') + end.count('Timeout') == 1 and end.count('sec') == 1
        if 'Timeout' in end: #162/477 Test #162: Singl[...]rch ...............***Timeout 300.06 sec
          time =   end[end.find('out')+3:end.find('sec')].strip()
        else:  
          time = end[end.find('ed')+2:end.find('sec')].strip()
        tests.append(Test(name, path, 'Passed' in end, float(time)))

  # add test which were not logged (not removed as found)
  for f in files:
    tests.append(Test(None, f))
    
  return tests     
        
# serches the TESTSUIT and PYTHON dir and gives the structure 
# we skip TESTSUIT but not PYTHON (as in ctest logs)
# e.g. Technical/PythonMesher (w/o TESTSUIT) or PYTHON/doctest/hdf5_tools
def parse_dir(dir): 
  org = dir        
    
  if not os.path.exists(dir + 'TESTSUIT/') or  not os.path.exists(dir + 'PYTHON/') or  not os.path.exists(dir + 'TestMacros.cmake'):
    print("dir argument '" + org + "' seems not to point to Testsuite including TESTSUITE, PYTHON, TestMacros.cmake, ...")
    os.sys.exit()   
    
  fg1 = glob.glob(dir + 'TESTSUIT/**/CMakeLists.txt' ,recursive = True)
  fg2 = glob.glob(dir + 'PYTHON/**/CMakeLists.txt' ,recursive = True)
  fg = fg1 + fg2
  assert len(fg) > 10 and 'TESTSUIT' in fg[0] # '../TESTSUIT/Technical/PythonMesher/CMakeLists.txt'

  # directories to be matched agains log. W/o TESTSUIT but with PYTHON    
  files = []  
  for f in fg: # Technical/PythonMesher
    if 'TESTSUIT' in f:
      tmp = f[f.find('TESTSUIT')+len('TESTSUIT/'):-len('/CMakeLists.txt')]
      if len(tmp) > 2: # skip the root CMakeLists which does not add tests 
        files.append(tmp)
    elif 'PYTHON' in f:
      tmp = f[f.find('PYTHON'):-len('/CMakeLists.txt')]
      if len(tmp) > len('PYTHON') + 2:
        files.append(tmp)    
 
  return files        

# helper to sort a label list to known and unknown lables
def sort_labels(test, known, unknown):
  for t in test:
    if t in known_labels:
      known.append(t)
    else:
      unknown.append(t)  

# count for every test the add_test and set_property in the CMakeLists.txt 
def read_cmake_lists(tests, dir):
  for t in tests:
    if t.path:
      t.add_test = 0
      t.set_property = 0
      t.cmakelists = dir + ('TESTSUIT/' if not t.path.startswith('PYTHON') else '') + t.path + '/CMakeLists.txt'
      file = open(t.cmakelists,'r')
      for l in file:
        s = l.strip().lower() # easy comparison by lower case
        if s.startswith('#'): # skip comments
          continue
        if 'add_test(' in s:  
          t.add_test += 1
        labels = parse_set_property(l)
        if labels:
          t.set_property += 1
          sort_labels(labels, t.label, t.unknown)

# helper which checks if the label has a proper timing
def check_timing(timing, time, label,fast,slow):
  assert not ('slow' in label and 'fast' in label),label
  assert timing == 'slow' or timing == 'fast'
  if timing == 'fast':
    if time <= fast and not 'fast' in label:
      return False
    if time > fast and 'fast' in label:
      return False
  if timing == 'slow':
    if time >= slow and not 'slow' in label:
      return False
    if time < slow and 'slow' in label:
      return False
  return True   

# helper wich finds the proper timing label. None for normal
def timing_label(timing, time,fast,slow):
  assert 'fast' in known_labels and 'slow' in known_labels
  assert timing == 'slow' or timing == 'fast'

  if timing == 'fast' and time <= fast:
    return 'fast'
  if timing == 'slow' and time >= slow:
    return 'slow'
  return None  

# rewrite CMakeLists.txt such that it matches the timings.
# assumes 'fast' and 'slow' in single set_property)
def write_timing_to_cmake_lists(timing, tests,fast,slow):
  assert timing == 'slow' or timing == 'fast'
  for t in tests:
    if t.path and t.time:
      if check_timing(timing, t.time, t.label, fast, slow):
        continue
      if t.add_test != 1:
        print('- skip writing ',t.cmakelists,'it has',t.add_test,'add_test()')
        continue
      print('- rewrite',t.cmakelists,'for', timing,'timing with',t.time,'sec')
      org = open(t.cmakelists,'r')
      input = org.readlines()
      org.close() 
      new = open(t.cmakelists,'w')
      for line in input:
        # we expect only a single timing in org
        assert not ('slow' in t.label and 'fast' in t.label)
        labels = parse_set_property(line) 
        if not labels or not (timing in labels): 
          new.write(line)
        else:
          print('  - remove',line.strip())  
      label = timing_label(timing,t.time, fast, slow)
      if label:
        s = set_property(label)
        new.write('\n' + s + '\n')
        print('  - add',s)
      new.close()      

def write_TEST_LABEL_property(tests):
  for t in tests:
    if t.path and '${TEST_LABELS}' not in t.unknown:
      print('- append',t.cmakelists,'for TESTSUITE_KEYWORDS guarded  ${TEST_LABELS} property')
      file = open(t.cmakelists,'a')
      file.write('\n')
      file.write('if(TESTSUITE_KEYWORDS)\n')
      file.write('   set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS ${TEST_LABELS})\n')
      file.write('endif()\n')
      file.close()


# read labels if the line is a set_property line.
# set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "fast")
# SET_PROPERTY(TEST ${TEST_NAME} APPEND PROPERTY LABELS ${TEST_LABELS})
# Else return []
def parse_set_property(line):
  s = line.strip().lower()
  if s.startswith('#'):
    return []
  # we assume single lines
  if not 'set_property(' in s:
    return []
  # extract testname
  assert line.index(' APPEND') > line.index('TEST '), line
  # t = l[l.index('TEST ') + len('TEST '):l.index(' APPEND')]
  # extract labels
  assert line.index(')') > line.index(' LABELS '), line
  t = line[line.index(' LABELS ') + len(' LABELS '):line.index(')')]
  t = t.replace('"','')
  labels = t.split(';') # cmake uses semicolons to define lists or ['fast']
  return labels
       

# helper which gives a set_property() string for the given label
def set_property(label):
  return 'set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "' + label + '")'
  


# conditionally print tests by list
def print_list(list, attr = 'testname'): 
  if args.list:
    for i,v in enumerate(list):
      if type(v) == str:
        s = v
      else:
        assert attr == 'testname' or attr == 'path'
        if attr == 'testname':
          s = v.testname
        else:
          s = v.path    
      print(i+1,s)
    print()    

# helper wich gives a nice average including empty list
def nice_avg(list):
  if len(list) == 0:
    return '-'
  else: 
    return round(np.mean(list),2)

# give message and to too slow tests for tests starting with type string as done for debug tests in .gitlab-ci.yml
def type_test(tests, type, slow):
  all_tests = []
  slow_tests = []
  non_slow = 0
  for t in tests:
    if t.time and t.testname.startswith(type):
      all_tests.append(t.time)
      if t.time >= slow:
        slow_tests.append(t.testname + ' time=' + str(t.time))
      else:
        non_slow += t.time  
  msg  = str(len(slow_tests)) + "/" + str(len(all_tests)) + " logged tests starting with '" + type + '" are too slow. ' 
  msg += str(round(non_slow,2)) + ' of ' + str(round(np.sum(all_tests),2)) + ' sec to be tested'   
  return msg, slow_tests

     
parser = argparse.ArgumentParser()
parser.add_argument("dir", help="path to Testsuite project directory")  
parser.add_argument("--log", help="log file of a ctest run. Usually a release log for write_timing fast and debug for slow")
parser.add_argument("--write_timing", help="writes specified timing, removes potential old other timing", choices = ['slow','fast'])
parser.add_argument("--write_missing_TEST_LABEL_property", help="append guarded TEST_LABEL property where missing", action='store_true')
parser.add_argument("--fast", help="tests below this time in seconds (release) are considered fast", type=float, default=1.0)
parser.add_argument("--slow", help="tests above this time in seconds (release) are considered slow", type=float, default=5.0)
parser.add_argument("--list", help="print individual items in status output", action='store_true')
parser.add_argument("--hist", help="plot histogram of testing time", action='store_true')
parser.add_argument("--range", help="highest displayable time for histogram in sec", type=float)
args = parser.parse_args()

if not args.dir.endswith('/'):
  args.dir += '/'

# test directorires 
files = parse_dir(args.dir)

# all test items with propers
tests = parse_log(files, args.log) # args.log might be None

read_cmake_lists(tests, args.dir)


no_path = [x.testname for x in tests if not x.path]
print('logged tests without a CMakeLists.txt:',len(no_path))
print_list(no_path)

no_add_test = [x for x in tests if x.add_test == 0]
print('CMakeLists.txt without add_test():',len(no_add_test))
print_list(no_add_test, 'path')

multiple_add_test = [x.path + ': ' + str(x.add_test) + ' add_test()' for x in tests if x.add_test > 1]
print('CMakeLists.txt with more than one add_test():',len(multiple_add_test))
print_list(multiple_add_test)

no_set_property = [x for x in tests if x.path and x.set_property == 0]
print('CMakeLists.txt without set_property():',len(no_set_property))
print_list(no_set_property, 'path')

fast_property = [x for x in tests if x.path and 'fast' in x.label]
print("CMakeLists.txt with property 'fast':",len(fast_property))
print_list(fast_property, 'path')


slow_property = [x.path + ' time=' + (str(x.time) if x.time else '-') for x in tests if x.path and 'slow' in x.label]
print("CMakeLists.txt with property 'slow':",len(slow_property))
print_list(slow_property)

unstable_property = [x for x in tests if x.path and 'unstable' in x.label]
print("CMakeLists.txt with property 'unstable':",len(unstable_property))
print_list(unstable_property, 'path')

broken_property = [x for x in tests if x.path and 'broken' in x.label]
print("CMakeLists.txt with property 'broken':",len(broken_property))
print_list(broken_property, 'path')

no_test_labels_property = [x for x in tests if x.path and '${TEST_LABELS}' not in x.unknown]
print("CMakeLists.txt with no property '${TEST_LABELS}':",len(no_test_labels_property))
print_list(no_test_labels_property, 'path')

if args.log:
  fast_label = [x.time for x in tests if x.time is not None and 'fast' in x.label]
  print("logged tests with property 'fast' (",len(fast_label),') sum up to', round(np.sum(fast_label),2),'sec with avg',nice_avg(fast_label),'sec')
  
  slow_label = [x.time for x in tests if x.time is not None and 'slow' in x.label]
  print("logged tests with property 'slow' (",len(slow_label),') sum up to', round(np.sum(slow_label),2),'sec with avg',nice_avg(slow_label),'sec')
  
  tested = [x.time for x in tests if x.passed != None]
  failed = [x for x in tests if x.passed == False] # not None!
  passed = [x.time for x in tests if x.passed == True] # not None!
  
  fast    = [x for x in tests if x.passed != None and x.time <= args.fast]
  tfast   = [x.time for x in fast]
  slow    = [x for x in tests if x.passed != None and x.time >= args.slow]
  tslow   = [x.time for x in slow]
  normal  = [x for x in tests if x.passed != None and x.time > args.fast and x.time < args.slow]
  tnormal = [x.time for x in normal]
  
  
  print('total time for all',len(tested),'tests in log:',np.sum(tested),'sec')
  print('according to log', len(fast),'fast tests <=',args.fast,'sec take', round(np.sum(tfast),2),'sec with avg',nice_avg(tfast),'sec')
  print_list(fast)
  
  
  print('according to log', len(normal),'normal tests',args.fast,'...',args.slow,'sec take', round(np.sum(tnormal),2),'sec with avg',nice_avg(tnormal),'sec')
  print_list(normal)
  
  normal_label = [x.time for x in tests if x.time is not None and ('fast' not in x.label and 'slow' not in x.label)]
  print("logged tests w/o time property (",len(normal_label),') sum up to', round(np.sum(normal_label),2),'sec with avg',nice_avg(normal_label),'sec')
  
  print('according to log', len(slow),'slow tests >=',args.slow,'sec take', round(np.sum(tslow),2),'sec with avg',nice_avg(tslow),'sec')
  if args.list:
    for i,f in enumerate(slow):
      print(i+1,f.testname,'failed' if f.testname in failed else 'passed')
    print()  
  
  
  # print timing accoriding to parallel debug run without slow
  msg, too_slow = type_test(tests, 'Solver', args.slow)
  print(msg)
  print_list(too_slow)
  
  msg, too_slow = type_test(tests, 'Singlefield', args.slow)
  print(msg)
  print_list(too_slow)
  
  msg, too_slow = type_test(tests, 'Coupledfield', args.slow)
  print(msg)
  print_list(too_slow)
  
  msg, too_slow = type_test(tests, 'Optimization', args.slow)
  print(msg)
  print_list(too_slow)
  
  msg, too_slow = type_test(tests, 'CFSdat', args.slow)
  print(msg)
  print_list(too_slow)
  
  msg, too_slow = type_test(tests, 'Technical', args.slow)
  print(msg)
  print_list(too_slow)
  
  
  print('slowest passed test in log:',np.max(passed),'sec')
  print('slowest failed test in log:',str(np.max([x.time for x in failed])) if len(failed) > 0 else '-','sec')
  
  if args.hist:
    import matplotlib.pyplot as plt
    range = args.range if args.range else np.max(tested) 
    if args.range:
      print('change history range to 0 ...',args.range,' instead of 0 ...',np.max(tested))
    else:
      print('-> show full history range 0 ...',np.max(tested),'use --range to reduce')  
    plt.hist(tested,bins=50, range=[0,range])
    if args.fast < range:
      plt.plot([args.fast,args.fast],[0,20],label='fast')
    if args.slow < range:  
      plt.plot([args.slow,args.slow],[0,20],label='slow')
    plt.legend()
    plt.xlabel('test time in sec')
    plt.ylabel('number of passed tests (total=' + str(len(passed)) + ')')
    plt.show()  
  
  print('passed tests:', len(passed))
  print('failed tests:',len(failed))
  if args.list:
    for i,x in enumerate(failed):
      print(i+1,x.testname,x.time,'sec')
    print()  
  
  print('test directories in Testsuite:',len([x for x in tests if x.path]))
  
  untested = [x for x in tests if not x.testname]
  print('Testsuite tests not in log:',len(untested))
  print_list(untested, 'path')
  
  #not_fast = [x.path + ' time=' + str(x.time) + ' labels=' + str(x.label) for x in fast and x.path and not 'fast' in x.label]
  not_fast = [x.path + ' time=' + str(x.time) + ' labels=' + str(x.label) for x in fast if x.path and not 'fast' in x.label]
  print("CMakeLists.txt which should be 'fast' but are not:",len(not_fast))
  print_list(not_fast)
  
  not_normal = [x.path + ' time=' + str(x.time) + ' labels=' + str(x.label) for x in normal if x.path and ('fast' in x.label or 'slow' in x.label)]
  print('CMakeLists.txt which have speed tag but are normal:',len(not_normal))
  print_list(not_normal)
  
  not_slow = [x.path + ' time=' + str(x.time) + ' labels=' + str(x.label) for x in slow if x.path and not 'slow' in x.label]
  print("CMakeLists.txt which should be 'slow' but are not:", len(not_slow))
  print_list(not_slow)
  
  if not args.hist:
    print('-> use --hist [with --range] to visual timing of log',args.log)
  
  if args.write_timing:
    write_timing_to_cmake_lists(args.write_timing, tests, args.fast, args.slow)
  else:
    print("-> use --write_timing to modify CMakeLists.txt according to log timing")

if not args.list:
  print('-> use --list to list the items')


if args.write_missing_TEST_LABEL_property:
  write_TEST_LABEL_property(tests)
