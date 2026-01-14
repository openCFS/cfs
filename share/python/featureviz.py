#!/usr/bin/env python
# This tool visualized feature mapping .density.xml pills like spaghetti but without inner

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from PIL import Image

import cfs_utils as ut


class Global:
  # the default value are for easy debugging via import
  def __init__(self):
    self.shapes = []         # array of pills
    self.n = [10, 10, 1]       # [nx, ny, nz]
    self.transition = -1.0   # optional transition zone width, -1 means not set

glob = Global()


# minimal and maximal are vectors.
def create_figure(res, minimal, maximal, scale):
  # calculate dpi such that the maximum dimension equals res / dpi
  dpi = res / 100.0 * scale

  fig = plt.figure(dpi=100, figsize=(round(dpi * maximal[0]), round(dpi * maximal[1])))
  ax = fig.add_subplot(111)
  ax.set_xlim(minimal[0], maximal[0])
  ax.set_ylim(minimal[1], maximal[1])
  ax.set_aspect('equal', adjustable='box')
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print(s)   

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'tab:orange', 'tab:brown', 'cornflowerblue', 'lime', 'tab:gray']
# transforms ids from 0 to 6 to color codes 'b' to 'k'. Only for matplotlib, not for vtk!
def matplotlib_color_coder(id):
  return colors[id % len(colors)]
    
class Pill: 
  # @param id starting from 0
  # @param base sum of all len(optvar()) of all previous shapes. Starts with 0
  def __init__(self,id,base,P,Q,p):
    self.id = id
    self.base = base
    self.color = matplotlib_color_coder(id)
    self.set(P,Q,p)

  def set(self, P,Q,p):
    
    self.P = np.array(P) # start coordinate
    self.Q = np.array(Q) # end coordinate
    self.p = p # profile is full structure thickness
    
    self.U = self.Q - self.P # base line p -> q
    norm_U = np.linalg.norm(self.U)
    if norm_U < 1e-20:
      self.U[:] = [1e-20,0]
      norm_U = np.linalg.norm(self.U)
    self.U0 = self.U / norm_U
    self.V = np.array((-self.U[1],self.U[0])) # normal to u
    self.V0 = self.V / np.linalg.norm(self.V)
   
    #print('Pill.set P',P,'Q',Q,'p',p,'U',self.U,'U0',self.U0,'V',self.V,'V0',self.V0)

  # @return p_x,p_y,q_x,q_y,p
  def optvar(self):
    return [self.P[0],self.P[1],self.Q[0],self.Q[1],self.p]

  # print the shape info
  def __str__(self):
    return "id=" + str(self.id) + " P=" + str(self.P) + " Q=" + str(self.Q) + " p=" + str(self.p) + " color=" + self.color  

  # shape info for a given index
  def to_string(self, idx):
    return "shape=" + str(self.id) + " color=" + str(self.color) + " idx=" + str(idx) 
   
# reads 2D and returns list of Spaghetti and domain. Also sets some values to glob!
# @param radius if given overwrites the value from the xml header
# @return list of spaghettis and either [[min_x, min_y], [max_x, max_y]] if in density.xml or None
def read_xml(filename, density = False, set = None, quiet = False):
 
  xml = ut.open_xml(filename)

  header = xml.xpath('/cfsErsatzMaterial/header')[0] 
  #ut.dump(header,'.')
  glob.n[0] = int(ut.xpath(header, './mesh/@x'))
  glob.n[1] = int(ut.xpath(header, './mesh/@y'))
  glob.n[2] = int(ut.xpath(header, './mesh/@z'))
  assert(glob.n[2] == 1)

  # modern cfs has a mesh domain, very modern optimization domain
  domain = None #
  pn = header.xpath('./optimizationDomain') # added 08.2025
  if len(pn) == 0:
    pn = header.xpath('./coordinateSystems/domain') # added 12.2022
  if len(pn) == 1: # we assume a single coordinate system
    dic = pn[0].attrib
    domain = [[float(dic['min_x']),float(dic['min_y'])],[float(dic['max_x']),float(dic['max_y'])]]
    if not quiet:
      print('read mesh discretization',glob.n,'domain bounds', domain)

  # this is the symmetric transition zone (void to solid) around the boundery
  pn = header.xpath('./featureMapping') # added 10.2025
  if len(pn) == 1 and 'transition' in pn[0].attrib:
    glob.transition = float(pn[0].attrib['transition'])

  shapes = []
  sq = 'last()' if set == None else '@id="' + str(set) + '"'

  set_data_res = xml.xpath('/cfsErsatzMaterial/set[' + sq + ']')
  if len(set_data_res) != 1:
    print('error: cannot open set', sq)
    sys.exit()
  data = set_data_res[0]  

  rho = None
  if density:
    physical = []
    for child in data:
      if(child.tag == "element"):
        physical.append(float(child.get("physical")))
    if len(physical) != np.prod(glob.n):
      print('cannot visualize density: there are',len(physical),'rho_e in the .density.xml set but mesh is',glob.n[0],'x',glob.n[1])
    else:
      # Fortran order (x changes fastest) 
      rho = np.reshape(physical,(glob.n[0],glob.n[1]), order='F')

  while True: # exit with break
    idx = len(shapes)
    base = './shapeParamElement[@shape="' + str(idx) + '"]'
    if len(data.xpath(base)) == 0:
      #print('no shapeParamElement with @shape=',idx)
      break

    # our xpath helper in optimization tools expects exactly one hit
    Px = float(ut.xpath(data, base + '[@dof="x"][@tip="start"]/@design'))
    Py = float(ut.xpath(data, base + '[@dof="y"][@tip="start"]/@design'))
    Qx = float(ut.xpath(data, base + '[@dof="x"][@tip="end"]/@design'))                   
    Qy = float(ut.xpath(data, base + '[@dof="y"][@tip="end"]/@design'))
    # width of noodle is 2*w -> don't confuse with P
    p  = float(ut.xpath(data, base + '[@type="profile"]/@design'))

    base = sum([len(s.optvar()) for s in shapes])
    pill = Pill(id=idx, base=base, P=(Px,Py), Q=(Qx,Qy), p=p)
    shapes.append(pill)
    # print('# read pill', pill)
      
  return shapes, domain, rho

def plot_rho(fig, domain, rho):
  if rho is None:
    return
  assert len(rho.shape) == 2
  # domain is [[min_x, min_y], [max_x, max_y]]
  nx, ny = rho.shape
  bx = domain[0][0] # base
  by = domain[0][1]
  dx = (domain[1][0] - bx) / nx
  dy = (domain[1][1] - by) / ny
  ax = fig.gca()
  for y in range(ny):
    for x in range(nx):
      ax.fill([bx+x*dx, bx+(x+1)*dx, bx+(x+1)*dx,bx+x*dx], [by+y*dy, by+y*dy, by+(y+1)*dy,by+(y+1)*dy], color=str(max(1-rho[x,y],0)), zorder=0)

def plot_outline(s, sub, p, fill, linestyle='-'):
  # code based on Patrick Jung's pltviz.py 

  L = np.hypot(s.U[0], s.U[1])
  if L < 1e-12 or s.p <= 0:
    return

  # s.U is P-Q, s.U0 is normed
  # s.V is normal and s.V0 normed normal

  # angle of the normal
  phi_n = np.arctan2(s.V0[1], s.V0[0])

  def arc(cx, cy, a0, a1, m=64):
      a = np.linspace(a0, a1, m)
      return np.c_[cx + p*np.cos(a), cy + p*np.sin(a)]

  # points of the contour (CCW)
  P_arc = arc(s.P[0], s.P[1], phi_n, phi_n + np.pi)            # half circle around P
  Q_minus = np.array([[s.Q[0] - s.V0[0]*p, s.Q[1] - s.V0[1]*p]]) # P -> Q
  Q_arc = arc(s.Q[0], s.Q[1], phi_n + np.pi, phi_n + 2*np.pi) 
  P_plus = np.array([[s.P[0] + s.V0[0]*p, s.P[0] + s.V0[1]*p]])  # back to origin

  pts = np.vstack([P_arc, Q_minus, Q_arc, P_plus])

  # construct path
  codes = [Path.MOVETO] + [Path.LINETO]*(len(pts)-2) + [Path.CLOSEPOLY]
  verts = np.vstack([pts[:-1], pts[0]])  # close
  path = Path(verts,codes)
  if fill:  
    sub.add_patch(PathPatch(path, facecolor=s.color, edgecolor="none", alpha=0.2, zorder=10))

  # contour
  sub.add_patch(PathPatch(path, facecolor="none", edgecolor=s.color, linewidth=2, zorder=10, capstyle="round", joinstyle="round", linestyle=linestyle))

def plot_data(res, shapes, detail, domain, ghost):
  # domain is typically [[0.0, 0.0], [1.0, 1.0]]
  LL = np.array(domain[0]) # lower left
  UR = np.array(domain[1]) # upper right
  minimal = LL - .5 * ghost * (UR-LL)
  maximal = UR + .5 * ghost * (UR-LL)

  # find the maximum dimension
  Dx = (maximal[0] - minimal[0])
  Dy = (maximal[1] - minimal[1])
  scale = 1 / max(Dx, Dy)

  lineopacity = 0.5 if args.gray else 1 # opacity value for plotting lines and points    

  fig, sub = create_figure(res, minimal, maximal, scale)

  if detail >= 1: # we plot the domain line alo for ghost=0 as we might have --noaxis
    # draw real domain
    sub.add_line(plt.Line2D((LL[0], UR[0], UR[0], LL[0], LL[0]), (LL[1], LL[1], UR[1], UR[1], LL[1]), color='black'))

  for s in shapes:
    # print('plot s',s)
    if s.p < 1e-10:  # omit zero-width shapes
      continue

    c = s.color

    if detail > 0:
      # start and endpoint
      fig.gca().add_artist(plt.Circle(s.P, 0.1 * s.p, alpha=lineopacity, color=c)) 
      fig.gca().add_artist(plt.Circle(s.Q, 0.1 * s.p, alpha=lineopacity, color=c))
      tail = '_{' + str(s.id) + '}$'
      if detail > 2:
        plt.annotate('$P' + tail, s.P, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = c)
        plt.annotate('$Q' + tail, s.Q, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = c)
        # we denote the profile 
        CC = s.P + .5 * s.U
        CL = CC - s.p * s.V0
        CR = CC + s.p * s.V0

        sub.add_line(plt.Line2D((CL[0],CR[0]),(CL[1],CR[1]), alpha=lineopacity, color= c))
        plt.annotate('$p' + tail, CC, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = s.color)

    # outline of the pill
    plot_outline(s, sub, s.p, fill=True)
    if detail >= 2:
      plot_outline(s, sub, s.p - .5*glob.transition, fill=False, linestyle=':')
      plot_outline(s, sub, s.p + .5*glob.transition, fill=False, linestyle=':')

  return fig

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml")
  parser.add_argument("--set", help="set within a .density.file", type=int)
  parser.add_argument('--interactive', help="interactively traverse sets in .density.xml", action='store_true')
  parser.add_argument('--save', help="save the image to the given name with the given format. Might be png, pdf, eps, vtp")
  parser.add_argument('--animate', help="save all sets as animated gif", action='store_true')
  parser.add_argument('--detail', help="level of technical details for spaghetti plot", choices=[0, 1, 2, 3], default=1, type=int)
  parser.add_argument('--density', help="visualize additionally the physical density", action='store_true')
  parser.add_argument('--ghost', help="add ghost domain in fraction of org size", type=float, default=0.0)
  parser.add_argument('--gray', help="plot grayscale image", action='store_true')

  parser.add_argument('--noshow', help="don't show the image", action='store_true')  
  parser.add_argument('--noaxis', help="dont plot axis and ticks", action='store_true')

  args = parser.parse_args()
  if not os.path.exists(args.input):
    print("error: cannot find '" + args.input + "'")
    os.sys.exit()
  file = args.input 
  
  colors = ['tab:gray'] if args.gray else ['b', 'g', 'r', 'c', 'm', 'y', 'tab:orange', 'tab:brown', 'cornflowerblue', 'lime', 'tab:gray']
  
  if args.interactive or args.animate:
    xml = ut.open_xml(file)
    sts = xml.xpath('/cfsErsatzMaterial/set/@id')
    if args.interactive: 
      idx = 0
      print('sets to traverse interactively',sts)
      while True:
        #plt.show(block=False)
        shapes, domain, rho = read_xml(args.input, args.density, str(sts[idx]), quiet=True)
        fig = plot_data(800,shapes,args.detail,domain, args.ghost)
        if args.density:
          plot_rho(fig, domain, rho)
        fig.show()  

        action = input('you see ' + str(set[idx]) + ' press (b)ack, (c)ancel/q/x or any other key + Enter for next: ')
        if action in ['c','q','x']:
          break
        idx += -1 if action == 'b' else +1
        idx = idx % len(sts) # wraps also negative numbers to 0 ... len(sts)-1
        plt.close(fig)
    if args.animate:
      frames = []
      name = file[:-len('.density.xml')] 
      print('create animated gif',name + '.gif','... ', end='',  flush=True)
      for s in sts:
        print(str(s) + ' ',end='', flush=True)
        shapes, domain, rho = read_xml(file, args.density, str(s), quiet=True)
        fig = plot_data(800,shapes,args.detail,domain, args.ghost)
        if args.density:
          plot_rho(fig, domain, rho)
        fig.gca().set_title(name + ' : ' + str(s))
        # Convert plot to image
        fig.canvas.draw()
        image = np.array(fig.canvas.renderer.buffer_rgba())
        frames.append(Image.fromarray(image))
        plt.close(fig)
      # Save as GIF, duration is in ms      
      print()
      print('save ...', flush=True)
      frames[0].save(name+ '.gif', save_all=True, append_images=frames[1:], duration=200, loop=0)
    sys.exit()  

  # sets some values in glob!
  shapes, domain, rho = read_xml(file, args.density, args.set)
  fig = plot_data(800,shapes,args.detail,domain, args.ghost)
  if args.density:
    plot_rho(fig, domain, rho)
  if args.noaxis:
    plt.axis('off')
  if args.save:
    print("write '" + args.save + "'")
    fig.savefig(args.save, bbox_inches='tight')
  if not args.noshow:
    fig.show()
    input("Press Enter to terminate.")

