#!/usr/bin/env python
# we refrain from using cfs_utils.py and optimization_tools.py, hece we need to reimplement a lot :(
# use lxml as it can be used from python2 and python3. It is based on libxml2 but the libxml2 python interface is not available for python3
import lxml
import lxml.etree
import argparse
import os.path

# The purpose of this script is to compare parts of the info.xml against a reference info.xml.
# To be used in the testsuite where mesh results cannot be used, e.g. eigenvalue problems


### Functions ###

# calculates the relative error or compares and eventually prints against eps
# @param good anything which can be converted to a float or be an xml attribute or even a list 
# @param test see good
# @param eps if given compares relative error against eps
# @param testinfo for error message if testinfo is given
# @param skip_noise skip too small values which are below the given value (e.g. lower modes for Bloch waves)
# @return the relative error or 0.0 for skip_noise if eps is not given or boolean if eps is given. True is bad and False is good
def has_rel_error(good, test, eps = 0.0, testinfo = None, skip_noise = None, verbose = False):
  #print 'has_rel_error with ' + str(type(good)) 
  # called by list?
  if isinstance(good, list):
    assert(isinstance(test, list))
    for good_item, test_item in zip(good, test):
      if has_rel_error(good_item, test_item, eps, testinfo, skip_noise,verbose):
        return True
    return False # no fail
    
  good_val = float(good)
  test_val = float(test)
 
  diff = abs((test_val - good_val) / good_val) if good_val != 0.0 else abs(test_val)
  if skip_noise and abs(good_val) < skip_noise and abs(test_val) < skip_noise:
    diff = 0.0
  if abs(good_val) <= 180 and abs(good_val) >= 180 - eps and abs(test_val) <= 180 and abs(test_val) > 180 - eps:
    diff = 0.0       
  #print(skip_noise, good_val, test_val, diff,eps,verbose)  

  if diff > eps:
    print(' * error: difference of ' + (testinfo + ': ' if testinfo else '') + str(good) + ' (ref) against ' + str(test) + ' (test) is ' + str(diff) + ' > eps ' + str(eps))
    return True     
  else:
    if verbose:
      print(' - good:',good,'vs. test',test,'->',diff,'<',eps,(('skip_noise=' + str(skip_noise)) if skip_noise else ''))
    return False

# compare the whole data and gives the maximal difference. 
# @param testinfo if the data is given and testinfo is given print a message
# @param maxdiff start value for maxdiff
def max_diff(good_list, test_list, testinfo = None, maxdiff = None, skip_noise = None, eps=0.0):
  assert(len(good_list) == len(test_list))
  max_diff = maxdiff if maxdiff else 0.0

  for good, test in zip(good_list, test_list):
    max_diff = max(max_diff, has_rel_error(good, test, skip_noise=skip_noise, eps=eps))

  if good_list and testinfo:
    print(' - the maximal relative difference comparing ' + testinfo + ' is ' + str(max_diff))   

  return max_diff

# sorts eigenfrequencies such that we can compare them as arpack might given them unsorted
# @return list of sorted floats 
def sort_eigenfrequencies(freq_list):
  # extract the content from the xmlAttr as a list, thanks to google  
  unsorted = [float(x) for x in freq_list]
  return sorted(unsorted)
     
# checks if we need to sort a list of xmlAttr
def need_to_sort(attr_list):
  # extract the content from the xmlAttr as a list, thanks to google  
  unsorted = [float(x) for x in attr_list]
  
  # we cannot check the sorted index list because of multiuple values: index = sorted(range(len(unsorted)), key=lambda k: unsorted[k])
  for i in range(2, len(unsorted)):
    if unsorted[i-1] > unsorted[i]:
      return True

  return False    
    

# return 0 if nothing to check, 1 if tested and -100 if failed
def compare_eigenfrequencies(ref, tst, eps, skip_noise):
  freq_ref = ref.xpath('//sequenceStep/eigenFrequency/result/mode/@frequency') # xmlAttr
  freq_tst = tst.xpath('//sequenceStep/eigenFrequency/result/mode/@frequency')

  dmp_ref = ref.xpath('//sequenceStep/eigenFrequency/result/mode/@damping') 
  dmp_tst = tst.xpath('//sequenceStep/eigenFrequency/result/mode/@damping') 

  print('test for standard eigenfrequencies ... ' + ('found' if freq_ref else 'none')) 

  if not freq_ref:
    return 0

  print(' - compare for ' + str(len(freq_ref)) + ' eigenfrequencies')
  print(' - has damping: ' + str(dmp_ref != []))

  if len(freq_ref) != len(freq_tst):
    print('error: number of eigenfrequencies differ: ref=' + str(freq_ref) + ' test=' + str(freq_tst))
    return -99

  # we might need sorting as for the bloch_waves!
  for i in range(len(freq_ref)):
     if has_rel_error(freq_ref[i], freq_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       print('error for ref=' + str(freq_ref[i]) + ' test=' + str(freq_tst[i]))
       return -99
     if dmp_ref and has_rel_error(dmp_ref[i], dmp_tst[i], eps, 'damping', skip_noise=skip_noise):
       return -99

  # prints the report
  max_diff(freq_ref, freq_tst, 'frequency',eps=eps)
  # does/ prints nothing if no damping
  max_diff(dmp_ref, dmp_tst, 'damping',eps=eps)
  
  print(' - no significant difference found') 
  return 1

# @return see compare_eigenfrequencies
def compare_bloch_waves(ref, tst, eps, skip_noise):
  wave_ref = ref.xpath('//sequenceStep/eigenFrequency/result/wave_vector') # xmlNode  
  wave_tst = tst.xpath('//sequenceStep/eigenFrequency/result/wave_vector')

  print('test for Bloch wave vectors ... ' + ('found' if wave_ref else 'none')) 

  if not wave_ref:
    return 0
 
  modes = len(wave_ref[0].xpath('mode/@nr'))
  print(' - compare for ' + str(len(wave_ref)) + ' wave vector with ' + str(modes) + ' modes')

  if len(wave_ref) != len(wave_tst):
    print(' * error: test hast ' + str(len(wave_tst)) + ' wave vectors')
    return -99
  
  assert(modes > 0)

  maxdiff = 0.0
  
  for w in range(len(wave_ref)):
    freq_ref = wave_ref[w].xpath('mode/@frequency')
    freq_tst = wave_tst[w].xpath('mode/@frequency')
    
    if len(freq_ref) != len(freq_tst):
      print(' * error: for mode ' + str(w) + ' number of modes differs: ' + str(len(freq_tst)) + ' instead of ' + str(len(freq_ref)))
      return -99
    
    assert(len(freq_ref) == modes)
    
    # we need to sort as arpack can return another sorting with gfortran/openblas and intel/mkl
    if need_to_sort(freq_ref) or need_to_sort(freq_tst):
      print(' - we need to sort the frequencies for wave vector ' + str(w) + ': ref=' + str(need_to_sort(freq_ref)) + ' test=' + str(need_to_sort(freq_tst)))  
    
    srt_ref = sort_eigenfrequencies(freq_ref)
    srt_tst = sort_eigenfrequencies(freq_tst)
    # use the list feature
    if has_rel_error(srt_ref, srt_tst, eps, skip_noise=skip_noise): 
      print('error: for wave vector ' + str(w) + ': ref=' + str(srt_ref) + ' test=' + str(srt_tst))
      return -99
          
    # global stuff
    maxdiff = max_diff(srt_ref, srt_tst, maxdiff = maxdiff,eps=eps)  

  print(' - maximal relative frequency difference over all modes in all wave vectors: ' + str(maxdiff))
  return 1

# @return see compare_eigenfrequencies
def compare_output_hystOperator(ref, tst, skip_noise, eps):
  px_ref = ref.xpath('//sequenceStep/result/Teststep/@Px') # xmlNode  
  px_tst = tst.xpath('//sequenceStep/result/Teststep/@Px')

  py_ref = ref.xpath('//sequenceStep/result/Teststep/@Py') # xmlNode  
  py_tst = tst.xpath('//sequenceStep/result/Teststep/@Py')
  
  pz_ref = ref.xpath('//sequenceStep/result/Teststep/@Pz') # xmlNode  
  pz_tst = tst.xpath('//sequenceStep/result/Teststep/@Pz')

  pabs_ref = ref.xpath('//sequenceStep/result/Teststep/@Pabs') # xmlNode  
  pabs_tst = tst.xpath('//sequenceStep/result/Teststep/@Pabs')

  phi_ref = ref.xpath('//sequenceStep/result/Teststep/@phi') # xmlNode  
  phi_tst = tst.xpath('//sequenceStep/result/Teststep/@phi')
  
  theta_ref = ref.xpath('//sequenceStep/result/Teststep/@theta') # xmlNode  
  theta_tst = tst.xpath('//sequenceStep/result/Teststep/@theta')
  
  print('test for output of hystOperator ... ' + ('found' if px_ref else 'none')) 

  if not px_ref:
    return 0

  print(' - compare for ' + str(len(px_ref)) + ' test steps' + ' (skip_noise:'+str(skip_noise)+')')

  if len(px_ref) != len(px_tst):
    print('error: number of test steps differ: ref=' + str(px_ref) + ' test=' + str(px_tst))
    return -99  

  for i in range(len(px_ref)):
     if has_rel_error(px_ref[i], px_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       return -99

  for i in range(len(py_ref)):
     if has_rel_error(py_ref[i], py_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       return -99
     
  for i in range(len(pz_ref)):
     if has_rel_error(pz_ref[i], pz_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       return -99

  for i in range(len(pabs_ref)):
     if has_rel_error(pabs_ref[i], pabs_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       return -99

  for i in range(len(phi_ref)):
     if has_rel_error(phi_ref[i], phi_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       return -99
     
  for i in range(len(theta_ref)):
     if has_rel_error(theta_ref[i], theta_tst[i], eps, skip_noise=skip_noise): # parses attribute 'frequency="254.660604829521"'
       return -99

  # prints the report
  max_diff(px_ref, px_tst, 'Px',eps=eps, skip_noise=skip_noise)
  max_diff(py_ref, py_tst, 'Py',eps=eps, skip_noise=skip_noise)
  max_diff(pz_ref, pz_tst, 'Pz',eps=eps, skip_noise=skip_noise)
  max_diff(pabs_ref, pabs_tst, 'Pabs',eps=eps, skip_noise=skip_noise)
  max_diff(phi_ref, phi_tst, 'phi',eps=eps, skip_noise=skip_noise)
  max_diff(theta_ref, theta_tst, 'theta',eps=eps, skip_noise=skip_noise)
  
  print(' - no significant difference found') 
  return 1
   

# return 0 if nothing to check, 1 if tested and -100 if failed
def compare_homogenized_tensor(ref, tst, eps, skip_noise, last,verbose):
  tensor_ref = ref.xpath('//process/iteration/homogenizedTensor/tensor/real') # xmlAttr
  tensor_tst = tst.xpath('//process/iteration/homogenizedTensor/tensor/real')

  print('test for homogenized tensor ... ' + ('found' if tensor_ref else 'none')) 

  if not tensor_ref:
    return 0

  if len(tensor_ref) != len(tensor_tst):
    if not last:  
      print('error: ' + str(len(tensor_ref)) + ' homogenized tensors in reference and ' + str(len(tensor_tst)) + ' in test!')
      return -99
    else:
      print('warning: ' + str(len(tensor_ref)) + ' homogenized tensors in reference and ' + str(len(tensor_tst)) + ' in test!')  

  tensor_ref = tensor_ref[-1].text.split()
  tensor_tst = tensor_tst[-1].text.split()
  err = has_rel_error(tensor_ref, tensor_tst, eps, skip_noise=skip_noise,verbose=verbose)

  # prints the report
  max_diff(tensor_ref, tensor_tst, 'the homogenized tensor', skip_noise=skip_noise,eps=eps)
  if err:
    return -99
  else: 
    print(' - no significant difference found')
    return 1


## compare the number of steps for iterCoupledPDE
def compare_iter_coupled_pde(ref,tst,verbose = False):
  steps_ref = ref.xpath('//PDE/iterCoupledPDE/convergence/step') 
  steps_tst = tst.xpath('//PDE/iterCoupledPDE/convergence/step') 
  
  print('test for coupled pde iterations ...'  + ('found' if steps_tst else 'none'))
  
  # do we have //PDE/iterCoupledPDE/convergence/step at all?
  if not steps_ref: 
    return 0 

  assert len(steps_ref) > 0
  
  if len(steps_ref) != len(steps_tst):
    print('ref has ',len(steps_ref), '//PDE/iterCoupledPDE/convergence/step, but test only',len(steps_tst))
    return -99
  
  for i in range(len(steps_ref)):
    srv = int(steps_ref[i].attrib['numIters'])
    stv = int(steps_tst[i].attrib['numIters'])
    if verbose:
      print('compare_iter_coupled_pde: ref step: number:',steps_ref[i].attrib['number'],'numIters:',srv)
      print('compare_iter_coupled_pde: test step: number:',steps_tst[i].attrib['number'],'numIters:',stv)
    assert steps_ref[i].attrib['number'] == steps_ref[i].attrib['number']
    if srv != stv:
      print('iterCoupledPDE step',steps_ref[i].attrib['number'],'has',stv,'numIters but expected',srv)
      return -99
    # one could here also compare converged for all quantities for all iterations
    
  print('numIters for all iterCoupledPDE steps match')
  return 1

## compare the regions for matching names, ids, types, ...
def compare_regions(ref,tst,verbose = False):
  regions_ref = ref.xpath('//header/domain/regions/region') 
  regions_tst = tst.xpath('//header/domain/regions/region') 
  
  print('test for regions ...'  + ('found' if regions_tst else 'none'))
  
  # do we have //header/domain/regions at all?
  if not regions_ref: 
    return 0 

  assert len(regions_ref) > 0
  
  if len(regions_ref) != len(regions_tst):
    print('ref has ',len(regions_ref), '//header/domain/regions/region, but test only',len(regions_tst))
    return -99
    
  for i in range(len(regions_ref)):
    print(f" - Region: {regions_ref[i].attrib['name']} ... {regions_ref[i].attrib == regions_tst[i].attrib}")
    if regions_ref[i].attrib != regions_tst[i].attrib:
      print(f"Region {i+1} of ref does not match test")
      return -99

  print(' > all Regions match')
  return 1
  
# # compare grid if written to info xml with 'cfs -G'
def compare_grid(ref,tst,verbose = False):
  grid_ref = ref.xpath('//header/domain/grid') 
  grid_tst = tst.xpath('//header/domain/grid') 
  
  print('test for grid ...'  + ('found' if grid_tst else 'none'))
  
  # do we have //header/domain/regions at all?
  if not grid_ref: 
    return 0 

  assert len(grid_ref) > 0

  node_ref = grid_ref[0].xpath('//nodeList/node')
  region_ref = grid_ref[0].xpath('//regionList/region')
  node_tst = grid_tst[0].xpath('//nodeList/node')
  region_tst = grid_tst[0].xpath('//regionList/region')

  assert len(node_ref) > 0
  assert len(region_ref) > 0


  if len(node_ref) != len(node_tst):
    print('ref: ',len(node_ref), 'nodes, test: ',len(node_tst))
    return -99

  print(f" - Checking {len(node_ref)} nodes ...")

  # Check each node
  for i in range(len(node_ref)):
    # print(node_tst[i].attrib)
    if node_ref[i].attrib != node_tst[i].attrib:
      print(f" - Node ID: {node_ref[i.attrib['id']]} of ref does not match test")
      return -99

  if len(region_ref) != len(region_ref):
    print('ref: ',len(region_ref), 'regions, test: ',len(region_tst))
    return -99

  # Check each region
  for j in range(len(region_ref)):
    element_ref = region_ref[j].xpath('//element')
    element_tst = region_tst[j].xpath('//element')
    
    if len(element_ref) != len(element_tst):
      print(f"region {region_ref[j].attrib['name']} of ref has {len(element_ref)} elements, test has {len(element_tst)}")
      return -99

    print(f" - Checking region {region_ref[j].attrib['name']}: {len(node_ref)} elements ...")
    # Check each element in region
    for i in range(len(element_ref)):
      # print(element_tst[i].attrib)
      if element_ref[i].attrib != element_tst[i].attrib:
        print(f" - Element ID: {element_ref[i].attrib['id']} of ref does not match test")
        return -99

  print(' > all Nodes and Elements in all regions match')
  return 1

## by default compare all attributes but number and duraton
# @param attrib if given compare only this one. Return silently with 0 when *ref* does not have attribute 
# @return 0 if nothing to check, 1 if tested and -100 if failed
def compare_opt_iterations(ref, tst, eps, skip_noise, last, attrib = None,verbose=False):
  opt_ref = ref.xpath('//optimization/process/iteration') 
  opt_tst = tst.xpath('//optimization/process/iteration')

  if attrib:
    print("test for attribute '" + attrib + "' in optimization iterations ... ")
  else:
    print('test for optimization iterations ... ' + ('found' if opt_ref else 'none')) 

  if not opt_ref:
    return 0

  if not last:
    if len(opt_ref) != len(opt_tst):
      print('error: number of written iterations differ: ref=' + str(len(opt_ref)) + ' test=' + str(len(opt_tst)))
      return -99
    else:
      print('compare for ' + str(len(opt_ref)) + ' written iterations') 

  for i in range(len(opt_ref)) if not last else [-1]:
     iter_ref = opt_ref[i]
     iter_tst = opt_tst[i]
   
     if attrib != None:
       if attrib not in iter_ref.attrib:
         return 0 # return silently
       else:
         ret = compare_opt_iterations_kernel(iter_ref, iter_tst, attrib, eps, skip_noise, verbose)
         if ret != 1:
           return ret
     else:    
       for prop in iter_ref.attrib:
         if prop not in ["number", "duration", "sub_prb_itr"] and not prop.startswith("infeas"):    
           ret = compare_opt_iterations_kernel(iter_ref, iter_tst, prop, eps, skip_noise, verbose)
           if ret != 1:
             return ret
  
  print(' - no significant difference found') 
  return 1

## helper for compare_opt_iterations
def compare_opt_iterations_kernel(iter_ref, iter_tst, prop, eps, skip_noise, verbose):
  if not prop in iter_tst.attrib:
    print('error: test hast no iteration attribute ' + prop)
    return -99

  val_ref = -1
  val_tst = -1
  
  try:
    val_ref = float(iter_ref.get(prop))
  except ValueError:
    pass # stay -1 and if ref and tst both fail, it is fine
  
  try:            
    val_tst = float(iter_tst.get(prop))
  except ValueError:
    pass 


  if has_rel_error(val_ref, val_tst, eps, prop, skip_noise, verbose):
    return -99
  
  return 1 # ok

# quick basecell test
# compares volumes, input parameters, radii
def compare_basecell_voxelized_mesh(ref, test, eps):
  
  ref_vol = ref.xpath('//basecell/volume/@value') 
  test_vol = test.xpath('//basecell/volume/@value') 
  
  print('test for basecell data ... ' + ('found' if ref_vol else 'none'))
  
  if not ref_vol: # if test is not being performed, return 0 and not -1
    return 0
  
  # test for same volume
  if has_rel_error(ref_vol, test_vol, eps):
    print('volumes do not match! vol_ref=' + str(ref_vol[0]) + ' vol_test=' + str(test_vol[0]))
    return -99
  
  ref_nx = ref.xpath('//basecell/@nx')
  test_nx = ref.xpath('//basecell/@nx')
  
  # test for same resolution 
  if has_rel_error(ref_nx, test_nx, eps, testinfo=True):
    print('resolution does not match! ref_nx=' + str(ref_nx) + ' and test_nx=' + str(test_nx))
    return -99
  
  # test for same input parameters
  ref_input = []
  test_input = []
  ref_input.append(ref.xpath('//basecell/input/@x1'))
  test_input.append(test.xpath('//basecell/input/@x1'))
  ref_input.append(ref.xpath('//basecell/input/@x2'))
  test_input.append(test.xpath('//basecell/input/@x2'))
  ref_input.append(ref.xpath('//basecell/input/@y1'))
  test_input.append(test.xpath('//basecell/input/@y1'))
  ref_input.append(ref.xpath('//basecell/input/@y2'))
  test_input.append(test.xpath('//basecell/input/@y2'))
  ref_input.append(ref.xpath('//basecell/input/@z1'))
  test_input.append(test.xpath('//basecell/input/@z1'))
  ref_input.append(ref.xpath('//basecell/input/@z2'))
  test_input.append(test.xpath('//basecell/input/@z2'))
  
  if has_rel_error(ref_input, test_input, eps, testinfo=True):
    print('input parameters {x1,x2,...} differ!')
    return -99
  
  # test cps for all three profiles
  for i in range(0,3):
    # test control polygons of profiles
    # a control polygon consists of 4 points (tuples with x and y coordinates)
    ref_cp = []
    test_cp = []
    
    # 0 degree spline
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P1/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P1/@y'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P2/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P2/@y'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P3/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P3/@y'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P4/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P4/@y'))
    
    # 90 degree spline
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P1/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P1/@y'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P2/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P2/@y'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P3/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P3/@y'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P4/@x'))
    ref_cp.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P4/@y'))
    
    # test: 0 degree spline
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P1/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P1/@y'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P2/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P2/@y'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P3/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P3/@y'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P4/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P4/@y'))
    # test:90 degree spline
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P1/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P1/@y'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P2/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P2/@y'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P3/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P3/@y'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="90.0"]/controlPolygon/P4/@x'))
    test_cp.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bspline[@degree="0.0"]/controlPolygon/P4/@y'))
  
    if has_rel_error(ref_cp, test_cp, eps, testinfo=True):
      print("control polygons of profile 0 mismatch!")
      return -99 
  
  ref_radius = []  
  ref_radius.append(ref.xpath('//basecell/radii/@rx1'))
  ref_radius.append(ref.xpath('//basecell/radii/@rx2'))
  ref_radius.append(ref.xpath('//basecell/radii/@ry1'))
  ref_radius.append(ref.xpath('//basecell/radii/@ry2'))
  ref_radius.append(ref.xpath('//basecell/radii/@rz1'))
  ref_radius.append(ref.xpath('//basecell/radii/@rz2'))
  
  test_radius = []
  test_radius.append(test.xpath('//basecell/radii/@rx1'))
  test_radius.append(test.xpath('//basecell/radii/@rx2'))
  test_radius.append(test.xpath('//basecell/radii/@ry1'))
  test_radius.append(test.xpath('//basecell/radii/@ry2'))
  test_radius.append(test.xpath('//basecell/radii/@rz1'))
  test_radius.append(test.xpath('//basecell/radii/@rz2'))
  
  if has_rel_error(ref_radius, test_radius, eps, testinfo=True):
    print("radii do not match! ref=" + str(ref_radius) + " test=" + str(test_radius)) 
    return -99
  
  # test bisectionFunction of all three profiles
  for i in range(0,3):
    ref_bisec = []
    test_bisec = []
    
    splines = ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/bicubic/bspline')
    for spline in splines:
      ref_bisec.append(spline.xpath('controlPolygon/P1/@x'))
      ref_bisec.append(spline.xpath('controlPolygon/P1/@y'))
    
    splines = test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/bicubic/bspline')
    for spline in splines:
      test_bisec.append(spline.xpath('controlPolygon/P1/@x'))
      test_bisec.append(spline.xpath('controlPolygon/P1/@y'))
      
    ref_bisec.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/polynomial/@a0'))
    ref_bisec.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/polynomial/@a1'))
    ref_bisec.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/polynomial/@a2'))
    ref_bisec.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/linear/@x1'))
    ref_bisec.append(ref.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/selection/@angle'))
    
    test_bisec.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/polynomial/@a0'))
    test_bisec.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/polynomial/@a1'))
    test_bisec.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/polynomial/@a2'))
    test_bisec.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/linear/@x1'))
    test_bisec.append(test.xpath('//basecell/profile[@dir="'+str(i)+'"]/bisectionFunction/selection/@angle'))
    
  if has_rel_error(ref_bisec, test_bisec, eps, testinfo=True):
    print("bisec functions mismatch") 
    return -99  
  
  return 1

# compare stochastic estimate outut by FEAST
def compare_stochastic_estimate(ref, tst, eps):
  # e.g. like the following, where N_estimated is just a hint but the nodes and weights are to be compared (size is 8 by default)
  #<header N_estimated="7"/>
  #<contour>
  #  <integrationNodes real="18.1621 90.4118 192.588 264.838 264.838 192.588 90.4118 18.1621" imag="51.0882 123.338 123.338 51.0882 -51.0882 -123.338 -123.338 -51.0882"/>
  #  <integrationWeights real="-15.4172 -6.38603 6.38603 15.4172 15.4172 6.38603 -6.38603 -15.4172" imag="6.38603 15.4172 15.4172 6.38603 -6.38603 -15.4172 -15.4172 -6.38603"/>
  #</contour>

  ref_est = ref.xpath('//solve_eigen/header/@N_estimated') # stochastic estimate from ref file
  tst_est = tst.xpath('//solve_eigen/header/@N_estimated') #

  ref_nds_r = ref.xpath('//solve_eigen/contour/integrationNodes/@real') 
  tst_nds_r = tst.xpath('//solve_eigen/contour/integrationNodes/@real') 
  ref_nds_i = ref.xpath('//solve_eigen/contour/integrationNodes/@imag') 
  tst_nds_i = tst.xpath('//solve_eigen/contour/integrationNodes/@imag') 

  ref_wts_r = ref.xpath('//solve_eigen/contour/integrationWeights/@real') 
  tst_wts_r = tst.xpath('//solve_eigen/contour/integrationWeights/@real') 
  ref_wts_i = ref.xpath('//solve_eigen/contour/integrationWeights/@imag') 
  tst_wts_i = tst.xpath('//solve_eigen/contour/integrationWeights/@imag') 

  print('test for stochasticEstimate ...'  + ('found' if tst_est else 'none'))

  if not ref_est: # returns 0 if stochastic estimate not found in ref file
    return 0 # test ignored

  print('ref header N_estimated=' + str(ref_est[0]) + ' vs. test header N_estimated=' + str(tst_est[0]) + ', allowerd to  may differ')

  # old binaries don't write nodes and weights
  assert ref_nds_r and ref_nds_i and ref_wts_r and ref_wts_i

  for r, t in [[ref_nds_r, tst_nds_r], [ref_nds_i, tst_nds_i], [ref_wts_r, tst_wts_r], [ref_wts_i, tst_wts_i]]:
    # xpath returns a list, therefore r[0], then we replce spaces by commas to eval as tuple 
    # and finally we get the first element of the tuple
    if has_rel_error(eval(r[0].replace(' ', ','))[0], eval(t[0].replace(' ', ','))[0], eps):
      print('contour nodes/weights do not match! ref=' + str(r) + ' test=' + str(t))
      return -99 # test failed

  return 1 # test passed

# compare the eigenvalue export to .info.xml
def compare_eigenvalues_to_infoxml(ref, tst, eps):

  evs_ref = ref.xpath('//sequenceStep/eigenValue/result/eigenvalues')
  evs_tst = tst.xpath('//sequenceStep/eigenValue/result/eigenvalues')

  print('test for eigenvalues to .info.xml ...'  + ('found' if evs_ref else 'none'))
  
  if not evs_ref:
    return 0

  for i in range(len(evs_ref)):
    re_ref = evs_ref[i].xpath('mode/@real')
    im_ref = evs_ref[i].xpath('mode/@imag')

    re_tst = evs_tst[i].xpath('mode/@real')
    im_tst = evs_tst[i].xpath('mode/@imag')

    # test for real part
    if has_rel_error(re_ref, re_tst, eps): # returns -1 if real part is not the same
      print(f'real parts of mode {i+1} do not match! ref=' + str(re_ref) + ' test=' + str(re_tst))
      return -99 # test failed
    
    # test for imag part
    if has_rel_error(im_ref, im_tst, eps): # returns -1 if real part is not the same
      print(f'imag parts of mode {i+1} do not match! ref=' + str(im_ref) + ' test=' + str(im_tst))
      return -99 # test failed
    
  return 1 # test passed

# compare the export of contur integration nodes and weights to .info.xml
def compare_contour_to_infoxml(ref, tst, eps):
  from numpy import fromstring

  # load data
  ## nodes
  nodes_ref = ref.xpath('//OLAS/acoustic-mechanic/solve_eigen/contour/integrationNodes')
  nodes_tst = tst.xpath('//OLAS/acoustic-mechanic/solve_eigen/contour/integrationNodes')
  ## weights
  weights_ref = ref.xpath('//OLAS/acoustic-mechanic/solve_eigen/contour/integrationWeights')
  weights_tst = tst.xpath('//OLAS/acoustic-mechanic/solve_eigen/contour/integrationWeights')

  # test real part
  ## nodes
  ref_re = fromstring(nodes_ref.xpath('@real') , dtype=float, sep=' ')
  tst_re = fromstring(nodes_tst.xpath('@real') , dtype=float, sep=' ')
  for i in range(len(tst_re)):
    if has_rel_error(ref_re[i], tst_re[i], eps):
      print(f'real parts of node {i+1} do not match! ref= {ref_re[i]} tst= {tst_re[i]}')
      return -99
  ## weights
  ref_re = fromstring(weights_ref.xpath('@real') , dtype=float, sep=' ')
  tst_re = fromstring(weights_tst.xpath('@real') , dtype=float, sep=' ')
  for i in range(len(tst_re)):
    if has_rel_error(ref_re[i], tst_re[i], eps):
      print(f'real parts of weight {i+1} do not match! ref= {ref_re[i]} tst= {tst_re[i]}')
      return -99
    
  # test imag part
  ## nodes
  ref_im = fromstring(nodes_ref.xpath('@imag') , dtype=float, sep=' ')
  tst_im = fromstring(nodes_tst.xpath('@imag') , dtype=float, sep=' ')
  for i in range(len(tst_im)):
    if has_rel_error(ref_im[i], tst_im[i], eps):
      print(f'imag parts of node {i+1} do not match! ref= {ref_im[i]} tst= {tst_im[i]}')
      return -99
  ## weights
  ref_im = fromstring(weights_ref.xpath('@imag') , dtype=float, sep=' ')
  tst_im = fromstring(weights_tst.xpath('@imag') , dtype=float, sep=' ')
  for i in range(len(tst_im)):
    if has_rel_error(ref_im[i], tst_im[i], eps):
      print(f'imag parts of weight {i+1} do not match! ref= {ref_im[i]} tst= {tst_im[i]}')
      return -99
    
  return 1


# compare the system matrices obtained after GetRidOfZeros function with reference onse
def compare_get_rid_of_zeros_hardcoded(arg_ref):
  from scipy.io import mmread
  
  path = os.path.dirname(arg_ref)
  # only for get rid of zeros
  if (os.path.basename(os.path.normpath(path)) == "GetRidOfZeros"):
    
    print('test for getRidOfZeros ... running')
    
    # ref matrices
    ref_auto_filename = os.path.join(path, "auto.ref.mtx")
    ref_yes_filename = os.path.join(path, "yes.ref.mtx")
    ref_no_filename = os.path.join(path, "no.ref.mtx")

    # test results matrices
    test_auto_filename = os.path.join(path, "auto_sys_0_0.mtx")
    test_yes_filename = os.path.join(path, "yes_step_1_sys_0_0.mtx")
    test_no_filename = os.path.join(path, "no_step_1_sys_0_0.mtx")

    ref_auto = mmread(ref_auto_filename)
    ref_yes = mmread(ref_yes_filename)
    ref_no = mmread(ref_no_filename)

    tst_auto = mmread(test_auto_filename)
    tst_yes = mmread(test_yes_filename)
    tst_no = mmread(test_no_filename)

    # (ref_auto != tst_auto) generates a sparse matrix of the same shape as ref_auto
    # and tst_auto where each element indicates whether the corresponding elements
    # in ref_auto and tst_auto are not equal.
    # If nnz of the generated matrix is zero, then the matrices are equal

    # test for "auto"
    if ((ref_auto != tst_auto).nnz == 0):
      print("\"auto\" option does not work")
      return -99 # test failed

    # test for "yes"
    if ((ref_yes != tst_yes).nnz == 0):
      print("\"yes\" option does not work")
      return -99 # test failed

    # test for "no"
    if ((ref_no != tst_no).nnz == 0):
      print("\"no\" option does not work")
      return -99 # test failed

    return 1 # test passed
  else:
    print('test for getRidOfZeros ... none')
    
    return 0 # not tested


### Script part ###
parser = argparse.ArgumentParser(description='this script is used to compare CFS++ info.xml files')
parser.add_argument("reference", help="the reference info.xml file")
parser.add_argument("test", help="the info.xml file to compare against the reference")
parser.add_argument("--eps", help="eps for relative comparison", type=float, default=1e-6)
parser.add_argument("--skip_noise", nargs='*', help="suppress too small values (e.g. rigid modes)", type=float)
parser.add_argument("--last", help="ignore different number of iterations and use last", action='store_true')
parser.add_argument("-v", "--verbose", help="show (on some tests) what is compared", action='store_true')
args = parser.parse_args()

# skip_noise is optional but it is easier to have it with 0 or 1 arguments
assert not args.skip_noise or len(args.skip_noise) == 0 or len(args.skip_noise) == 1
skip_noise = 1e-10 if (not args.skip_noise or len(args.skip_noise) == 0) else float(args.skip_noise[0])   

print("compare_info_xml: ref '" + args.reference + "'")
print("                  tst '" + args.test + "'")
print('options: eps:',args.eps,', skip_noise:',args.skip_noise, '->',skip_noise,', ignore mismatch in iterations (--last):',('yes' if args.last else 'no'))
print('')  

if not os.path.exists(args.reference):
  print("error: reference file '" + args.reference + "' does not exist")
  os.sys.exit(1)
if not os.path.exists(args.test):
  print("error: test file '" + args.test + "' does not exist")
  os.sys.exit(1)

ref = lxml.etree.parse(args.reference, lxml.etree.XMLParser(remove_comments=True))
tst = lxml.etree.parse(args.test, lxml.etree.XMLParser(remove_comments=True))

# for a passed test return 1, if test does not apply return 0, if it fails return large negative (-99).
count = 0

# test first for bandgaps, these are very hard to check otherwise, e.g. because of many inactive modes
count += compare_opt_iterations(ref, tst, args.eps, skip_noise, args.last, attrib = 'bandgap_3_4')
if count == 0:
  # we don't compare eigenvalue stuff only if there was no bandgap within the ev stuff
  count += compare_eigenfrequencies(ref, tst, args.eps, skip_noise)
  count += compare_bloch_waves(ref, tst, args.eps, skip_noise)
count += compare_homogenized_tensor(ref, tst, args.eps, skip_noise, args.last, verbose=args.verbose)
count += compare_basecell_voxelized_mesh(ref, tst, args.eps)
count += compare_output_hystOperator(ref, tst, skip_noise=skip_noise, eps=args.eps)
count += compare_iter_coupled_pde(ref,tst,verbose=args.verbose)
count += compare_stochastic_estimate(ref,tst,args.eps)
count += compare_eigenvalues_to_infoxml(ref, tst, args.eps)
count += compare_get_rid_of_zeros_hardcoded(args.reference)
# if not compare_basecell_voxelized_mesh(ref, tst, eps):
#   count += compare_all_profile_functions(ref,tst,args.eps)
if count == 0: # test opt iterations only if no other test possible
  count += compare_regions(ref,tst,verbose=args.verbose)
  count += compare_grid(ref,tst,verbose=args.verbose)
  count += compare_opt_iterations(ref, tst, args.eps, skip_noise, args.last,verbose=args.verbose)
else:
  print('skip comparing iterations')
  
if count > 0:
  print(f'++ no error detected. count = {count} ++')        
  os.sys.exit(0)
if count < 0:
  print(f'** error: a test failed. count = {count}  **')        
  os.sys.exit(1)
if count == 0:
  print(f'** error: no known datastructure found to compare. **')        
  os.sys.exit(2)
