#!/usr/bin/env python
from mesh_tool import *
from optimization_tools import *
import argparse
import types
import numpy

## structure the mesh file such that it is easier to parse
# note that for the elements, the lines might be split for too long long node lists
# abaqus has a really stupid file format here :(
def structure_abaqus(all_lines): 
  pos = -1
  #assume *NODE first to be ended by *something. Only once!
  node_start = -1
  node_end   = -1
  #assume *ELEMENT to end with *something (which might be the next *ELEMENT). Multiple!
  element_start = []
  element_end   = []
  #strange HWNAME COMP sections to be ended by **something not HWNAME COMP
  comp_start = []
  comp_end = []
  for line in all_lines:
    pos += 1
    if line.startswith('*NODE'):
      assert(node_start == -1)
      node_start = pos # including *NODE!
      continue
    if line.startswith('*') and node_start > 0 and node_end == -1:
      node_end = pos # including!
      continue
    # check element end before element start such that an element group can be ended by the start of the next
    if line.startswith('*') and len(element_start) == len(element_end) + 1:
      element_end.append(pos) # this is included  
      #conditionally no continue
    if line.startswith('*ELEMENT'):
      assert(len(element_start) == len(element_end))
      element_start.append(pos) # including *ELEMENT
      continue
    if line.startswith('**HWNAME COMP'):
      if len(comp_end) < len(comp_start):
        comp_end.append(pos)
      comp_start.append(pos)
      continue
    if line.startswith('**') and not line.startswith('**HWCOLOR') and len(comp_start) == len(comp_end) + 1:
      comp_end.append(pos)
      continue
    
  assert(len(element_start) == len(element_end))  
  return node_start, node_end, element_start, element_end, comp_start, comp_end

# find the maximum number (first element) within a range
def find_max_number(all_lines, start_search, end_search):
  max_number = -1
  
  for i in range(start_search, end_search):
    token = all_lines[i].split()[0].replace(',','')
    #print token
    max_number = max((max_number, int(token)))

  assert(max_number > -1)
  return max_number

def parse_element_header(line):
  split = line.replace('\n', ' ').split(',')
  assert(len(split) == 3)
  assert(split[0] == '*ELEMENT')
  
  type = -1
  assert(split[1].startswith('TYPE='))
  abaqus_type = split[1][5:]
  if abaqus_type == 'C3D6':
    type = WEDGE6
  if abaqus_type == 'C3D8':
    type = HEXA8
  assert(type <> -1)   

  assert(split[2].startswith('ELSET='))
  name =  split[2][6:]

  return type, name

def parse_comp_header(line):
  return line.split()[3]

def parse_abaqus(input):
  lines = open(input,"r").readlines()
  
  mesh = Mesh()

  # read and map nodes
  node_start, node_end, element_start, element_end, comp_start, comp_end = structure_abaqus(lines)
  
  # we need to map from abaqus node number to 0-based index
  max_node_nr = find_max_number(lines, node_start+1, node_end)
  
  # by index of original abaqus node number we map to the 0-based index
  abaqus_nr_to_idx = [-1] * (max_node_nr+1) # ignore *NODE
  # by 0-based index we map to the original abaqus node number
  idx_to_abaqus_nr = [-1] * (node_end - node_start - 1) # remove *NODE
  
  for i in range(node_end - node_start - 1):
    line = lines[i + node_start + 1].replace(',',' ').split()
    #print line
    if len(line) > 0 and not line[0].startswith('*'):
      nr = int(line[0])
      x = float(line[1])
      y = float(line[2])
      z = float(line[3])
      mesh.nodes.append((x, y, z))
      abaqus_nr_to_idx[nr] = i
      idx_to_abaqus_nr[i]  = nr
 
  mesh.nodes = np.array(mesh.nodes)    
      
  
  # read elements
  for e in range(len(element_start)):
    type, name = parse_element_header(lines[element_start[e]])
    numnodes = nodes_by_type(type)
    first = element_start[e] + 1 # header
    last  = element_end[e] # would index to the following *something
    # we cannot loop with for as we may increment j within the loop
    j = first
    while j < last:
      e = Element()
      e.region = name
      e.type = type
      # we have the really bad problem, that too long lines are
      # split in two lines with an eol :((
      line = lines[j]
      if line.endswith(',\n'):
        line += lines[j+1]
        j += 1
        
      tokens = line.replace(',', ' ').split()
      for n in range(1,len(tokens)): # skip number
        nr = int(tokens[n]) 
        idx = abaqus_nr_to_idx[nr]
        assert(idx > -1)
        e.nodes.append(idx)
      assert(len(e.nodes) == numnodes)
      mesh.elements.append(e)
      j += 1  
  
  #process named nodes
  # in the robot art, the *KINEMATIC COUPLING,REF NODE=711809 are not within elements, give them extra name 'ref'
  ref = []  
  for c in range(len(comp_start)):
    name = parse_comp_header(lines[comp_start[c]])
    nodes = []
    for n in range(comp_start[c], comp_end[c]):
      line = lines[n]
      if line.startswith('*KINEMATIC COUPLING,REF NODE='):
        nr = int(line.replace('\n','').split('=')[-1])
        ref.append(abaqus_nr_to_idx[nr])
      else:
        # don't know what in 5006,1,3 the second and third value means
        if not line.startswith('*'):
          nr = int(line.split(',')[0])
          nodes.append(abaqus_nr_to_idx[nr])         
    mesh.bc.append((name, nodes))
  mesh.bc.append(('ref', ref))  
  # print mesh.nodes
  # print idx_to_abaqus_nr
  return mesh

parser = argparse.ArgumentParser()
parser.add_argument("input", help="abaqus input file")
parser.add_argument("output", help="openCFS mesh file")


args = parser.parse_args()
if not os.path.exists(args.input):
  print 'input file not found: ' + args.input
  sys.exit() 

mesh = parse_abaqus(args.input)

write_ansys_mesh(mesh, args.output)