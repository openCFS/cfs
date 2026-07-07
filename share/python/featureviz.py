#!/usr/bin/env python
# This tool visualized feature mapping .density.xml pills like spaghetti but without inner

import os
import sys
import argparse
import tempfile
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
    self.extension = 0.0     # optional outward enlargement of the transition zone (asymmetric Bezier)
    self.alpha_q = 1.0       # penalization exponent q of the geometry variable alpha (density scales with alpha^q)

glob = Global()


# minimal and maximal are vectors.
# @param fig reuse this figure (e.g. the persistent window of the interactive set traversal)
def create_figure(res, minimal, maximal, scale, fig=None):
  # calculate dpi such that the maximum dimension equals res / dpi
  dpi = res / 100.0 * scale

  if fig is None:
    fig = plt.figure(dpi=100, figsize=(round(dpi * maximal[0]), round(dpi * maximal[1])))
  else:
    fig.clf()
  ax = fig.add_subplot(111)
  ax.set_xlim(minimal[0], maximal[0])
  ax.set_ylim(minimal[1], maximal[1])
  ax.set_aspect('equal', adjustable='box')
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print(s)   

colors = ['tab:blue', 'tab:green', 'tab:red', 'c', 'm', 'y', 'tab:orange', 'tab:brown', 'cornflowerblue', 'lime', 'tab:gray']
# transforms ids from 0 to 6 to color codes 'b' to 'k'. Only for matplotlib, not for vtk!
def matplotlib_color_coder(id):
  return colors[id % len(colors)]
    
class Pill:
  # @param id starting from 0
  # @param base sum of all len(optvar()) of all previous shapes. Starts with 0
  # @param alpha optional geometry variable in [0,1] scaling the feature density as alpha^q * rho
  def __init__(self,id,base,P,Q,p,alpha=1.0):
    self.id = id
    self.base = base
    self.color = matplotlib_color_coder(id)
    self.alpha = alpha
    self.set(P,Q,p)

  def set(self, P,Q,p):

    self.P = np.array(P) # start coordinate, 2D or 3D
    self.Q = np.array(Q) # end coordinate
    self.p = p # profile is full structure thickness

    self.U = self.Q - self.P # base line p -> q
    norm_U = np.linalg.norm(self.U)
    if norm_U < 1e-20:
      self.U[:] = 0.0
      self.U[0] = 1e-20
      norm_U = np.linalg.norm(self.U)
    self.U0 = self.U / norm_U
    if len(self.P) == 2: # the normal is only defined (and only needed for plotting) in 2D
      self.V = np.array((-self.U[1],self.U[0])) # normal to u
      self.V0 = self.V / np.linalg.norm(self.V)

    #print('Pill.set P',P,'Q',Q,'p',p,'U',self.U,'U0',self.U0,'V',self.V,'V0',self.V0)

  # @return p_x,p_y,(p_z),q_x,q_y,(q_z),p
  def optvar(self):
    return list(self.P) + list(self.Q) + [self.p]

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
  d3 = glob.n[2] > 1

  # modern cfs has a mesh domain, very modern optimization domain
  domain = None #
  pn = header.xpath('./optimizationDomain') # added 08.2025
  if len(pn) == 0:
    pn = header.xpath('./coordinateSystems/domain') # added 12.2022
  if len(pn) == 1: # we assume a single coordinate system
    dic = pn[0].attrib
    domain = [[float(dic['min_x']),float(dic['min_y'])],[float(dic['max_x']),float(dic['max_y'])]]
    if d3:
      domain[0].append(float(dic['min_z']))
      domain[1].append(float(dic['max_z']))
    if not quiet:
      print('read mesh discretization',glob.n,'domain bounds', domain)

  # this is the symmetric transition zone (void to solid) around the boundery
  pn = header.xpath('./featureMapping') # added 10.2025
  if len(pn) == 1 and 'transition' in pn[0].attrib:
    glob.transition = float(pn[0].attrib['transition'])
    if 'extension' in pn[0].attrib:
      glob.extension = float(pn[0].attrib['extension'])
    if 'alpha_q' in pn[0].attrib: # only written when the geometry variable alpha is used
      glob.alpha_q = float(pn[0].attrib['alpha_q'])

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
        physical.append(float(child.get("physical" if "physical" in child.attrib else "design")))
    if len(physical) != np.prod(glob.n):
      print('cannot visualize density: there are',len(physical),'rho_e in the .density.xml set but mesh is',glob.n)
    else:
      # Fortran order (x changes fastest)
      shape = (glob.n[0],glob.n[1],glob.n[2]) if d3 else (glob.n[0],glob.n[1])
      rho = np.reshape(physical,shape, order='F')

  while True: # exit with break
    idx = len(shapes)
    base = './shapeParamElement[@shape="' + str(idx) + '"]'
    if len(data.xpath(base)) == 0:
      #print('no shapeParamElement with @shape=',idx)
      break

    # our xpath helper in optimization tools expects exactly one hit
    P = [float(ut.xpath(data, base + '[@dof="' + dof + '"][@tip="start"]/@design')) for dof in (('x','y','z') if d3 else ('x','y'))]
    Q = [float(ut.xpath(data, base + '[@dof="' + dof + '"][@tip="end"]/@design')) for dof in (('x','y','z') if d3 else ('x','y'))]
    # width of noodle is 2*w -> don't confuse with P
    p  = float(ut.xpath(data, base + '[@type="profile"]/@design'))
    # the optional geometry variable alpha (older files and runs without alpha lack it)
    al = data.xpath(base + '[@type="alpha"]/@design')
    alpha = float(al[0]) if len(al) == 1 else 1.0

    base = sum([len(s.optvar()) for s in shapes])
    pill = Pill(id=idx, base=base, P=P, Q=Q, p=p, alpha=alpha)
    shapes.append(pill)
    # print('# read pill', pill)
      
  return shapes, domain, rho

def plot_rho(fig, domain, rho, ref=None, grid=False):
  if rho is None and ref is None:
    return
  # ax.fill(color=c) sets BOTH face and edge to c, so cells get a per-cell border that stacks into a
  # mesh; force edgecolor='none' to suppress it. --grid instead draws an explicit dark mesh.
  edge = dict(edgecolor=(0.3, 0.3, 0.3), linewidth=0.4) if grid else dict(edgecolor='none')
  # the reference is always drawn when given; the main density only when --density supplied rho
  grid = rho if rho is not None else ref
  assert len(grid.shape) == 2
  if rho is not None and ref is not None and ref.shape != rho.shape:
    print('cannot overlay reference: mesh', ref.shape, 'differs from', rho.shape)
    ref = None
  # domain is [[min_x, min_y], [max_x, max_y]]
  nx, ny = grid.shape
  bx = domain[0][0] # base
  by = domain[0][1]
  dx = (domain[1][0] - bx) / nx
  dy = (domain[1][1] - by) / ny
  ax = fig.gca()

  def cell(x, y):
    return [bx+x*dx, bx+(x+1)*dx, bx+(x+1)*dx, bx+x*dx], [by+y*dy, by+y*dy, by+(y+1)*dy, by+(y+1)*dy]

  # reference first (blue, alpha by intensity), then the new density on top (black, alpha by intensity).
  # over void (ref==0) the black alpha-over-white reduces to (1-rho) grayscale, identical to the plain
  # density, so density above void is unchanged; only cells overlapping the reference blend.
  if ref is not None:
    for y in range(ny):
      for x in range(nx):
        r = max(min(ref[x,y], 1.0), 0.0)
        if r > 0:
          xs, ys = cell(x, y)
          ax.fill(xs, ys, color=(0, 0, 1), alpha=r, zorder=0, antialiased=False, **edge)
  if rho is not None:
    for y in range(ny):
      for x in range(nx):
        n = max(min(rho[x,y], 1.0), 0.0)
        if n > 0:
          xs, ys = cell(x, y)
          ax.fill(xs, ys, color=(0, 0, 0), alpha=n, zorder=1, antialiased=False, **edge)

def plot_outline(s, sub, p, fill, linestyle='-', fill_alpha=0.2, line_alpha=1.0):
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
    sub.add_patch(PathPatch(path, facecolor=s.color, edgecolor="none", alpha=fill_alpha, zorder=10))

  # contour
  sub.add_patch(PathPatch(path, facecolor="none", edgecolor=s.color, linewidth=2, zorder=10, capstyle="round", joinstyle="round", linestyle=linestyle, alpha=line_alpha))

def plot_data(res, shapes, detail, domain, ghost, fig=None):
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

  fig, sub = create_figure(res, minimal, maximal, scale, fig)

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
      fig.gca().add_artist(plt.Circle(s.P, 0.1 * s.p, alpha=lineopacity, color=c, zorder=11))
      fig.gca().add_artist(plt.Circle(s.Q, 0.1 * s.p, alpha=lineopacity, color=c, zorder=11))
      tail = '_{' + str(s.id) + '}$' if len(shapes) > 1 else '$'
      if detail > 2:
        plt.annotate('$P' + tail, s.P, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = c)
        plt.annotate('$Q' + tail, s.Q, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = c)
        # we denote the profile 
        CC = s.P + .5 * s.U
        CL = CC - s.p * s.V0
        CR = CC + s.p * s.V0

        sub.add_line(plt.Line2D((CL[0],CR[0]),(CL[1],CR[1]), alpha=lineopacity, color= c))
        plt.annotate('$2p' + tail, CC, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = s.color)

    # the geometry variable alpha fades the pill out: the filling and the transition zone lines scale
    # with the density factor alpha^q, the solid outline is only dimmed (stays clearly visible)
    aq = s.alpha ** glob.alpha_q
    # outline of the pill
    plot_outline(s, sub, s.p, fill=True, fill_alpha=0.2*aq, line_alpha=0.3 + 0.7*s.alpha)
    if detail >= 2:
      plot_outline(s, sub, s.p - .5*glob.transition, fill=False, linestyle=':', line_alpha=aq)                  # solid side, dist = -transition/2
      plot_outline(s, sub, s.p + .5*glob.transition + glob.extension, fill=False, linestyle=':', line_alpha=aq) # void side, dist = +(transition/2 + extension)

  return fig

# reads the tabular text file (e.g. the cfs .plot.dat, '#' comment lines are skipped) for --plot
# @param spec [file, optional x column, optional y column]. Columns are 1-based indices or names
#        from the '#1:iter\t2:compliance\t...' header (a prefix suffices), defaults 1 and 2
# @return dict with iter (first column), x, y and axis labels
def read_plot_dat(spec):
  path = spec[0]
  if not os.path.exists(path):
    print("error: cannot find --plot file '" + path + "'")
    sys.exit()
  data = np.loadtxt(path, comments='#', ndmin=2)
  names = []
  with open(path) as f:
    first = f.readline().strip()
  if first.startswith('#') and ':' in first: # cfs style '#1:iter\t2:compliance\t...'
    names = [t.split(':',1)[1].strip() for t in first[1:].split('\t') if ':' in t]
    if len(names) != data.shape[1]:
      names = []

  def col(key, default):
    if key is None:
      return default
    if key.isdigit():
      return int(key) - 1
    for i, n in enumerate(names):
      if n.startswith(key):
        return i
    print("error: --plot column '" + key + "' not found in", names if names else 'headerless ' + path)
    sys.exit()

  ix = col(spec[1] if len(spec) > 1 else None, 0)
  iy = col(spec[2] if len(spec) > 2 else None, 1)
  labels = [names[i] if names else 'column ' + str(i+1) for i in (ix, iy)]
  return {'iter': data[:,0], 'x': data[:,ix], 'y': data[:,iy], 'xlabel': labels[0], 'ylabel': labels[1]}

# the .plot.dat row of a set: sets are written per iteration, so match the set id against the iter column
# @return row index or None
def plot_dat_row(plot, set_id):
  row = np.nonzero(plot['iter'] == float(set_id))[0]
  return row[0] if len(row) else None

# renders 3D pills with vtk instead of matplotlib: interactive window, or --save to .png/.pdf/.eps
# (off-screen screenshot) or .vtp (polydata for ParaView, e.g. to overlay with the density from the .cfs result)
# @param file the .density.xml, needed to traverse sets
# @param sets for --interactive the list of set ids: Right/Left keys switch sets within the window,
#        the camera persists
# @param plot optional read_plot_dat() result shown as accompanying chart window with a marker
#        following the current set
def visualize_3d(shapes, domain, args, file=None, sets=None, plot=None):
  import vtk
  from matplotlib.colors import to_rgb
  from matviz_vtk import create_capsule_polydata, generate_outline_box, show_write_vtk, save_screenshot

  if args.save and args.save.endswith('.vtp'):
    # single polydata for ParaView with per-point 'Colors' and the geometry variable as 'alpha' array
    append = vtk.vtkAppendPolyData()
    for s in shapes:
      if s.p < 1e-10: # omit zero-width shapes
        continue
      append.AddInputData(create_capsule_polydata(s.P, s.Q, s.p, color=to_rgb(s.color), scalar=s.alpha))
    append.Update()
    show_write_vtk(append.GetOutput(), 800, args.save)
    return

  renderer = vtk.vtkRenderer()
  renderer.SetBackground(1, 1, 1)
  # one actor per shape, kept when traversing sets - update() only exchanges geometry and properties
  mappers = []
  actors = []
  for _ in shapes:
    mapper = vtk.vtkPolyDataMapper()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer.AddActor(actor)
    mappers.append(mapper)
    actors.append(actor)

  def update(shapes):
    for s, mapper, actor in zip(shapes, mappers, actors):
      actor.SetVisibility(s.p >= 1e-10) # omit zero-width shapes
      if s.p < 1e-10:
        continue
      mapper.SetInputData(create_capsule_polydata(s.P, s.Q, s.p))
      actor.GetProperty().SetColor(to_rgb(s.color))
      # the geometry variable alpha scales the density with alpha^q -> shown as transparency, with a
      # floor so switched-off features stay clearly visible (like the dimmed outline in 2D); --alpha
      # is a global base opacity on top (see through overlaps, 1.0 for opaque publication shots)
      actor.GetProperty().SetOpacity(args.alpha * (0.3 + 0.7 * s.alpha ** glob.alpha_q))
  update(shapes)

  if args.detail >= 1 and domain is not None:
    box = generate_outline_box(size=list(np.array(domain[1])-np.array(domain[0])), offset=list(domain[0]))
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(box)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer.AddActor(actor)

  window = vtk.vtkRenderWindow()
  window.AddRenderer(renderer)
  window.SetSize(800, 800)
  # correct blending of intersecting translucent capsules needs depth peeling
  window.SetAlphaBitPlanes(1)
  window.SetMultiSamples(0)
  renderer.SetUseDepthPeeling(1)
  renderer.SetMaximumNumberOfPeels(50)
  renderer.SetOcclusionRatio(0.05)

  if args.save:
    window.Render()
    save_screenshot(window, args.save)
    return
  if args.noshow:
    return

  window.SetWindowName(str(file) if file else 'featureviz')
  interactor = vtk.vtkRenderWindowInteractor()
  interactor.SetRenderWindow(window)

  # accompanying chart of the .plot.dat as a strip below the 3D scene - within the same window, a
  # second window gets hidden behind terminal/other windows and a matplotlib window cannot share the
  # vtk event loop. The red marker follows the current set
  if sets and plot is not None:
    table = vtk.vtkTable()
    for name, values in ((plot['xlabel'], plot['x']), (plot['ylabel'], plot['y'])):
      arr = vtk.vtkFloatArray()
      arr.SetName(name)
      for v in values:
        arr.InsertNextValue(v)
      table.AddColumn(arr)
    marker = vtk.vtkTable()
    for name in (plot['xlabel'], plot['ylabel']):
      arr = vtk.vtkFloatArray()
      arr.SetName(name)
      arr.InsertNextValue(0.0)
      marker.AddColumn(arr)

    chart = vtk.vtkChartXY()
    line = chart.AddPlot(vtk.vtkChart.LINE)
    line.SetInputData(table, 0, 1)
    line.SetWidth(2.0)
    dot = chart.AddPlot(vtk.vtkChart.POINTS)
    dot.SetInputData(marker, 0, 1)
    dot.SetColor(255, 0, 0, 255)
    dot.SetMarkerSize(12.0)
    chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle(plot['xlabel'])
    chart.GetAxis(vtk.vtkAxis.LEFT).SetTitle(plot['ylabel'])

    chart_actor = vtk.vtkContextActor()
    chart_actor.GetScene().AddItem(chart)
    strip = 0.27 # fraction of the window height for the chart
    renderer.SetViewport(0.0, strip, 1.0, 1.0)
    chart_ren = vtk.vtkRenderer()
    chart_ren.SetBackground(1, 1, 1)
    chart_ren.SetViewport(0.0, 0.0, 1.0, strip)
    chart_ren.AddActor(chart_actor)
    window.AddRenderer(chart_ren)
    window.SetSize(800, 1100) # 800x800 scene plus the chart strip

    def update_marker(set_id): # the caller renders the window
      row = plot_dat_row(plot, set_id)
      dot.SetVisible(row is not None)
      if row is not None:
        marker.SetValue(0, 0, float(plot['x'][row]))
        marker.SetValue(0, 1, float(plot['y'][row]))
        marker.Modified()
    update_marker(sets[0])

  if sets:
    state = {'idx': 0}
    def show_set(idx):
      state['idx'] = idx % len(sets) # wraps also negative numbers
      shapes, _, _ = read_xml(file, False, str(sets[state['idx']]), quiet=True)
      update(shapes)
      window.SetWindowName(file + ' : set ' + str(sets[state['idx']]) + '/' + str(sets[-1]))
      print('set ' + str(sets[state['idx']]) + '/' + str(sets[-1]))
      if plot is not None:
        update_marker(sets[state['idx']])
      window.Render()
    def on_key(obj, event):
      key = obj.GetKeySym()
      if key in ('Right', 'n', 'space'):
        show_set(state['idx'] + 1)
      elif key in ('Left', 'b'):
        show_set(state['idx'] - 1)
    interactor.AddObserver('KeyPressEvent', on_key)
    # the terminal works like the 2D traversal prompt: a repeating timer polls stdin
    def on_timer(obj, event):
      import select
      if select.select([sys.stdin], [], [], 0)[0]:
        line = sys.stdin.readline()
        if not line: # EOF, e.g. piped stdin - stop polling
          obj.RemoveObservers('TimerEvent')
          return
        action = line.strip()
        if action in ('c', 'q', 'x'):
          obj.TerminateApp()
        else:
          show_set(state['idx'] + (-1 if action == 'b' else +1))
    interactor.AddObserver('TimerEvent', on_timer)
    window.SetWindowName(file + ' : set ' + str(sets[0]) + '/' + str(sets[-1]))

  window.Render()
  if sets:
    interactor.Initialize()
    interactor.CreateRepeatingTimer(250)
  interactor.Start() # blocks until the window is closed
  cam = renderer.GetActiveCamera()
  print('camera when closing window: {:f} {:f} {:f} {:f} {:f} {:f} {:f}'.format(*cam.GetPosition(), *cam.GetFocalPoint(), cam.GetRoll()))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml")
  parser.add_argument("--set", help="set within a .density.file", type=int)
  parser.add_argument('--interactive', help="interactively traverse sets in .density.xml", action='store_true')
  parser.add_argument('--save', help="save the image to the given name with the given format. Might be png, pdf, eps, vtp")
  parser.add_argument('--animate', help="save all sets as animated gif", action='store_true')
  parser.add_argument('--detail', help="level of technical details for spaghetti plot", choices=[0, 1, 2, 3], default=1, type=int)
  parser.add_argument('--density', help="visualize additionally the physical density", action='store_true')
  parser.add_argument('--alpha', help="3D only: global base opacity of the capsules, multiplied on the per-feature alpha", type=float, default=0.7)
  parser.add_argument('--plot', nargs='+', help="accompany --interactive with a chart following the current set: <data-file> [x-column] [y-column] with columns as 1-based index or header name (defaults 1 2, i.e. iter and objective of a .plot.dat)")
  parser.add_argument('--ghost', help="add ghost domain in fraction of org size", type=float, default=0.0)
  parser.add_argument('--gray', help="plot grayscale image", action='store_true')
  parser.add_argument('--grid', help="draw the mesh grid on the rendered density/reference cells", action='store_true')
  parser.add_argument('--reference', help="overlay a reference .density.xml beneath the main density, shown in blue (implies --density)")

  parser.add_argument('--noshow', help="don't show the image", action='store_true')
  parser.add_argument('--detach', help="open the image in an external viewer and return immediately; the window persists after the app ends", action='store_true')
  parser.add_argument('--noaxis', help="dont plot axis and ticks", action='store_true')

  args = parser.parse_args()
  if not os.path.exists(args.input):
    print("error: cannot find '" + args.input + "'")
    os.sys.exit()
  file = args.input 
  
  ref_rho = None
  if args.reference:
    if not os.path.exists(args.reference):
      print("error: cannot find reference '" + args.reference + "'")
      os.sys.exit()
    _, _, ref_rho = read_xml(args.reference, density=True, quiet=True) # sets glob to ref mesh; main reads below reset it

  colors = ['tab:gray'] if args.gray else ['tab:blue', 'tab:green', 'tab:red', 'c', 'm', 'y', 'tab:orange', 'tab:brown', 'cornflowerblue', 'lime', 'tab:gray']
  
  if args.plot and not args.interactive:
    print('--plot is only used with --interactive')

  if args.interactive or args.animate:
    xml = ut.open_xml(file)
    d3 = int(xml.xpath('/cfsErsatzMaterial/header/mesh/@z')[0]) > 1
    if d3 and args.animate:
      print('--animate is not implemented for 3D')
      sys.exit()
    plot = read_plot_dat(args.plot) if args.plot and args.interactive else None
    sts = xml.xpath('/cfsErsatzMaterial/set/@id')
    if args.interactive:
      print('sets to traverse interactively',sts)
      print('in the focused graphics window press: Right/n/space = next set, Left/b = previous set, q = quit')
      print('or type here in the terminal: Enter = next, b + Enter = back, q + Enter = quit')
      if d3: # a single vtk window, sets are switched by key within the window or via the terminal
        shapes, domain, rho = read_xml(args.input, False, str(sts[0]), quiet=True)
        visualize_3d(shapes, domain, args, file=file, sets=sts, plot=plot)
        sys.exit()

      # 2D: one persistent figure (keeps window size/position), same dual control as 3D
      state = {'idx': 0}
      fig = None

      # accompanying chart of the .plot.dat, the red marker follows the current set
      pfig = None
      def update_marker():
        if pfig is None:
          return
        row = plot_dat_row(plot, sts[state['idx']])
        if row is not None:
          pmark.set_data([plot['x'][row]], [plot['y'][row]])
        else:
          pmark.set_data([], [])
        pfig.canvas.draw_idle()

      def show_set(idx):
        global fig
        state['idx'] = idx % len(sts) # wraps also negative numbers
        shapes, domain, rho = read_xml(args.input, args.density, str(sts[state['idx']]), quiet=True)
        fig = plot_data(800, shapes, args.detail, domain, args.ghost, fig=fig)
        if args.density or ref_rho is not None:
          plot_rho(fig, domain, rho, ref_rho, args.grid)
        fig.canvas.manager.set_window_title(file + ' : set ' + str(sts[state['idx']]) + '/' + str(sts[-1]))
        print('set ' + str(sts[state['idx']]) + '/' + str(sts[-1]))
        update_marker()
        fig.canvas.draw_idle()
      show_set(0)

      # created after the main figure so it opens on top of it - the macosx backend offers no window
      # positioning, drag it aside once, the position persists over the session
      if plot is not None:
        pfig = plt.figure(figsize=(8, 4))
        pax = pfig.add_subplot(111)
        pax.plot(plot['x'], plot['y'], '-')
        pmark, = pax.plot([], [], 'ro', markersize=8)
        pax.set_xlabel(plot['xlabel'])
        pax.set_ylabel(plot['ylabel'])
        pfig.canvas.manager.set_window_title(plot['ylabel'] + ' over ' + plot['xlabel'])
        pfig.tight_layout()
        update_marker()

      # the default keymap would put Left/Right on the toolbar view history
      kid = getattr(fig.canvas.manager, 'key_press_handler_id', None)
      if kid is not None:
        fig.canvas.mpl_disconnect(kid)
      def on_key(event):
        if event.key in ('right', 'n', ' '):
          show_set(state['idx'] + 1)
        elif event.key in ('left', 'b'):
          show_set(state['idx'] - 1)
        elif event.key == 'q':
          plt.close('all')
      fig.canvas.mpl_connect('key_press_event', on_key)
      # closing the main window (mouse) also ends the accompanying chart, else plt.show() keeps blocking
      fig.canvas.mpl_connect('close_event', lambda event: plt.close('all'))

      # the terminal works like a prompt: a timer within the gui event loop polls stdin
      def on_timer():
        import select
        if select.select([sys.stdin], [], [], 0)[0]:
          line = sys.stdin.readline()
          if not line: # EOF, e.g. piped stdin - stop polling
            timer.stop()
            return
          action = line.strip()
          if action in ('c', 'q', 'x'):
            plt.close('all')
          else:
            show_set(state['idx'] + (-1 if action == 'b' else +1))
      timer = fig.canvas.new_timer(interval=250)
      timer.add_callback(on_timer)
      timer.start()
      plt.show() # blocks until the figure is closed
    if args.animate:
      frames = []
      name = file[:-len('.density.xml')] 
      print('create animated gif',name + '.gif','... ', end='',  flush=True)
      for s in sts:
        print(str(s) + ' ',end='', flush=True)
        shapes, domain, rho = read_xml(file, args.density, str(s), quiet=True)
        fig = plot_data(800,shapes,args.detail,domain, args.ghost)
        if args.density or ref_rho is not None:
          plot_rho(fig, domain, rho, ref_rho, args.grid)
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

  if glob.n[2] > 1: # 3D is rendered with vtk, the matplotlib code below is 2D only
    if args.density:
      print('the density is not rendered in 3D - load the .cfs result and the pills .vtp in ParaView instead')
    visualize_3d(shapes, domain, args, file=file)
    sys.exit()

  fig = plot_data(800,shapes,args.detail,domain, args.ghost)
  if args.density or ref_rho is not None:
    plot_rho(fig, domain, rho, ref_rho, args.grid)
  if args.noaxis:
    plt.axis('off')
  if args.save:
    print("write '" + args.save + "'")
    fig.savefig(args.save, bbox_inches='tight')
  if not args.noshow:
    if args.detach:
      # save to a temp file and hand it to the OS viewer (like show_density.py); the window persists after we exit
      fd, tmp = tempfile.mkstemp(suffix='.png')
      os.close(fd)
      fig.savefig(tmp, bbox_inches='tight')
      Image.open(tmp).show()
    else:
      plt.show()

