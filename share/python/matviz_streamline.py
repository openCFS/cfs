from matviz_2d import *
from scipy import ndimage
import numpy
from  numpy.linalg import norm
import scipy.interpolate
import xml.etree.ElementTree



## Subclass for Fields
# this performs an interpolation, stores the data and allows access
class Data:
  def __init__(self, min, max, nx, input_space, input_data, method):

    self.NO_VAL = -12345.6

    self.min = min
    self.max = max
    self.nx = nx
    dx = (max[0] - min[0]) / nx
    self.dx = dx
    ny = int((max[1] - min[1] + 1e-3) / dx)
    self.ny = ny
    # dy might be slightly different from dx but is close
    dy = (max[1] - min[1]) / ny
    self.dy = dy

    # convenience for nx and ny
    self.ndim = (nx, ny)

    out = numpy.zeros((nx*ny, 2))
    for y in range(ny):
      for x in range(nx):
        out[y * nx + x][0] = min[0] + dx * x + 0.5 * dx
        out[y * nx + x][1] = min[0] + dy * y + 0.5 * dy

    self.data = scipy.interpolate.griddata(input_space, input_data, out, method, self.NO_VAL)

    if method != 'nearest':
      self.nearest = scipy.interpolate.griddata(input_space, input_data, out, 'nearest')



  ## tests if the coordinates are within min, max
  def within(self, x, y):
    if x < self.min[0] or x > self.max[0]:
      return False
    if y < self.min[1] or y > self.max[1]:
      return False
    return True

  # limit to be within
  def limit(self, x, y):
    # assume not to be within
    x_out = max((self.min[0], min((self.max[0], x))))
    y_out = max((self.min[1], min((self.max[1], y))))

    #print 'x=' + str(x) + ' y=' + str(y) + str((x_out, y_out))
    return (x_out, y_out)





  ## gives for the coordinates the value line. Does not further interpolate but assumes piecewise constant data.
  # has a fallback to nearest
  def getData(self, x, y):
    i = int((x - self.min[0]) / self.dx)
    #j = int((self.max[1] - y - self.min[1]) / self.dx)
    j = int((y - self.min[1]) / self.dy)
    #print 'getData x=' + str(x) + ' y=' + str(y) + ' -> i=' + str(i) + ' j=' + str(j) + ' idx=' + str(self.nx * j + i),
    v = self.data[self.nx * j + i]
    if v[0] == self.NO_VAL:
      v = self.nearest[self.nx * j + i]

    #print ' -> ' + str(v)
    return v

  def dump_data(self):
    for j in range(self.ny):
      for i in range(self.nx):
        x = self.min[0] + i * self.dx + 0.5 * self.dx
        y = self.min[1] + j * self.dy + 0.5 * self.dy
        o = self.getData(x, y)
        print('i=' + str(i) + ' j=' + str(j) + ' -> ' + str((x, y)) + ' -> ' + str(o))

  ## dumps the data as image for debug purpose
  def dump_image(self, data_idx):
    out = numpy.zeros((self.ny, self.nx))

    minv = min(self.data[:,data_idx])
    maxv = max(self.data[:,data_idx])

    print(minv)
    print(maxv)

    for j in range(self.ny):
      for i in range(self.nx):
        #v = self.data[self.nx*j + i][data_idx]
        o = self.getData(i * self.dx, j * self.dy)[data_idx]
        s = (o + -1.0 * minv) / (maxv - minv) * 255
        #print ' scale from ' + str(o) + ' -> ' + str(s) + ' base=' + str(o + -1.0 * minv) + ' scale=' + str((maxv - minv) * 255)
        out[j,i] = s

    img = Image.fromarray(out)
    #img.resize((600,600), Image.ANTIALIAS)
    img.show()


## This holds a trace
class Trace:
  #@param cell tupel if indices of macro cell where we start in the center
  #@param prominent fakes max such the this Trace is sorted firsts
  def __init__(self, points, data, cell, idx, prominent = False):
    self.points = points
    self.data   = data
    self.cell   = cell
    self.idx    = idx
    self.max    = 9.0 if prominent else max(self.data)
    self.avg    = sum(self.data) / len(self.data)
    self.prominent = prominent
    assert(min(self.data) >= 0.0)
    #print 'min trace ' + str(min(self.data)) + ' max trace ' + str(self.max)

  #give the macro fields we touch with the trace
  def touched(self, macro):
    data = numpy.zeros((macro.nx, macro.ny))
    for p in self.points:
      ix, iy = coord2index(macro, p)
      #print 'p=' + str(p) + ' min=' + str(macro.min) + ' d=' + str((macro.dy, macro.dy)) + ' i=' + str((ix, iy))

      data[ix,iy] = max((data[ix,iy], 1))

    return data



## this contains the physical fields stiff1 or stiff2 and angle in the macroscopic and fine discretization.
# The macroscopic is where we consider cells. This might be the original discretization
# The fine discretization is used for doing the streamlines. This might be the same as the macroscopic
class Fields:

  ## constructs the data sets macro and fine
  def __init__(self, coords, s, angle, macro_samples):

    centers, min, max, elem = coords
    # convert to 2D
    c, v = convert_single_data_interpolation_input(centers, s, angle)

    # this is the original data discretization
    nx = int((max[0] - min[0]) / elem[0])
    # our macro discretization can be given from the user
    dx = nx if macro_samples == None else macro_samples

    self.macro = Data(min, max, dx, c, v, 'linear')
    # no need to be much smaller than original, if coarser assume 0.2 step length.
    self.fine  = Data(min, max, numpy.min((3 * nx, 5 * dx)), c, v, 'linear')

  ## calculates the streamline from a given point up to both ends
  # @param idx 0 for s1 and 1 for s2 which changes the angle!
  # @param coord optional pair of coordinates. Take care not to have an invalid cell!!
  # @param prominent see Trace
  # return a list with two trace objects for both indices
  def streamline(self, steplength, minimal, idx, cell, coord = None, prominent = False):
    x, y = coord if coord != None else self.index2coord(cell)

    # print 'streamline x=' + str(x) + ' y=' + str(y) + ' idx=' + str(idx) + ' steplength=' + str(steplength)

    # we construct the traces for both indices and determine by the average value what is stiff1 and stiff2.
    # better would be to check for horizontal or vertical extend

    # we want a single trace but we cannot simply concatenate
    # the second trace needs to be reverted and set before first trace
    trace1 = self.directional_streamline(x,y, self.macro.dx * steplength, minimal, idx, 1.0)
    trace2 = self.directional_streamline(x,y, self.macro.dx * steplength, minimal, idx, -1.0)
    tmp = trace2[::-1] + trace1

    val = [i[1] for i in tmp]

    if len(val) > 0:
      return Trace([i[0] for i in tmp], val, cell, idx, prominent)
    else:
      return None



  ## helper for streamline
  #@param direction 1.0 for forward, -1.0 for backward
  def directional_streamline(self, x, y, steplength, minimal, idx, direction, method = "euler"):

    trace = []
    trace.append(((x, y), self.fine.getData(x, y)[0]))
    while True:
      here  = self.fine.getData(x, y)
      val   = here[0]
      assert(val >= 0.0)
      angle = here[1] + idx * numpy.pi/2 # for idx==0 nothing changes

      if method == "euler":
        #print 'x=' + str(x) + ' y=' + str(y) + ' data=' + str(here),
        x += steplength * numpy.cos(angle) * direction
        y += steplength * numpy.sin(angle) * direction
        #print ' next_x=' + str(x) + ' next_y=' + str(y) + ' val=' + str(val)
      if method == "midpoint":
        x_m = x + 0.5 * steplength * numpy.cos(angle) * direction
        y_m = y + 0.5 * steplength * numpy.sin(angle) * direction
        mid = self.fine.getData(x_m, y_m)
        angle_m = mid[1] + idx * numpy.pi/2 # for idx==0 nothing changes
        x += steplength * numpy.cos(angle_m) * direction
        y += steplength * numpy.sin(angle_m) * direction
      if val < minimal: # stop stream if we go into void
        break;
      if len(trace) > 2000: # prevent circular streams
        break
      if self.fine.within(x, y):
        trace.append(((x, y), val))
      else:
        #print 'limit ' + str(self.fine.limit(x, y))
        trace.append((self.fine.limit(x, y), val))
        break
    return trace if len(trace) > 2 else []

  ## from macro index give cell center coordinates
  def index2coord(macro, cell):
    x = (cell[0] + 0.5) * macro.macro.dx
    y = (cell[1] + 0.5) * macro.macro.dy
    return x, y

def coord2index(data, p):
  ix = (p[0] - 1e-6 - data.min[0]) / data.dx
  iy = (p[1] - 1e-6 - data.min[1]) / data.dy
  return ix, iy


## draws a trace
# @param scale scale stiffness to color. 1.0 for normal or 2.0 if max stiff = 0.5
def draw_trace(fig, trace):
  # we search for portions in the trace with similar value which need to be at least two segments wide.
  # Small value errors are ok. value equals color

  vertices = trace.points
  values   = trace.data

  start = 0
  pos = 1
  val = values[0]
  while pos < len(vertices):
    # up to the value changes or we reach the end
    if abs(values[pos] - val) > 0.1 or pos == len(vertices)-1:
      assert(pos > start)
      #print 's=' + str(start) + ' e=' + str(pos) + ' v=' + str(val)
      path = Path(vertices[start:pos+1])
      #print path
      gray = max((1.0 - val, 0.0))
      if gray > 1.0:
        print('invalid color v=' + str(val) + ' -> ' + str(gray))
        gray = 0.5
      patch = matplotlib.patches.PathPatch(path, edgecolor=str(gray), facecolor='none', lw=1)
      fig.add_patch(patch)

      start = pos-1 # make sure there is no gap between two segments
      val = values[pos]
    pos += 1


## draws a thick trace
# @param scale scale stiffness to color. 1.0 for normal or 2.0 if max stiff = 0.5
def draw_thick_trace(fig, dx, trace, minimal):
  # two adjacent coordinates and the intermediate value give a bar with the thickness


  # we search for portions in the trace with similar value which need to be at least two segments wide.
  # Small value errors are ok. value equals color

  vertices = trace.points
  values   = trace.data

  for i in range(0,len(vertices)-1,1):
    #print i
    if max(values[i:i+2]) > minimal:
      draw_frustum(fig, dx, vertices[i:i+2], values[i:i+2])



## draw a frustum defined by the start and end point of the center line and the thickness at the start and end
def draw_frustum(fig, dx, coord, thick):
  assert(len(coord) == len(thick))
  assert(len(coord) == 2)

  # construct orthogonal vector by solving for the scalar product
  # p1: start of center line
  # p2: end of center line
  # v1 = p2-p1
  # p3: new point on orthogonal vector, one component set, the other solved
  # v2 = p3-p1  -> v1*v2 = 0
  x1 = coord[0][0]
  x2 = coord[1][0]
  y1 = coord[0][1]
  y2 = coord[1][1]

  v1 = [x2-x1, y2-y1]
  v2 = numpy.array((y2-y1, x1-x2))
  v2n = norm(v2)

  if v2n < 1e-6:
    return

  v2 = v2 * 1./ v2n

  #print 'coord: ' + str(coord) + ' thickness: ' + str(thick) + ' v2: ' + str(v2)# + ' sp=' + str(numpy.dot(v1, v2))

  #v2 = numpy.zeros(2)

  p1 = [x1,y1]
  p2 = [x2,y2]

  t1 = 0.5*dx*thick[0]
  t2 = 0.5*dx*thick[1]

  verts = [p1+t1*v2, p2+t2*v2, p2-t2*v2, p1-t1*v2, (0,0)]
  codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]

  path = Path(verts, codes)
  patch = matplotlib.patches.PathPatch(path, edgecolor='black', facecolor='black', lw=1)
  fig.add_patch(patch)

  #o1 = p1+t1*v2
  #o2 = p2+t2*v2
  #o3 = p2-t2*v2
  #o4 = p1-t1*v2

  #print 'p1=' + str(p1) + ' p2=' + str(p2) + ' v2=' + str(v2) + ' t1=' + str(t1) + ' t2=' + str(t2)
  #print 'o1=' + str(o1) + ' o2=' + str(o2) + ' o3=' + str(o3) + ' o4=' + str(o4)

  #print 'plot ' + str(x1) + ', ' + str(y1) + ' with points, ', str(x2) + ', ' + str(y2) + ' with points, ',
  #print str(o1[0]) + ', ' + str(o1[1]) + ' with points, ',
  #print str(o2[0]) + ', ' + str(o2[1]) + ' with points, ',
  #print str(o3[0]) + ', ' + str(o3[1]) + ' with points, ',
  #print str(o4[0]) + ', ' + str(o4[1]) + ' with points'


# helper for forcing streamlines
def force_prominent_streamline(traces, fields, minimal, step, x, y):
  idx = 0 if fields[0].macro.getData(x,y)[0] > fields[1].macro.getData(x,y)[0] else 1
  cell = (coord2index(fields[idx].macro, (x, y)))
  trace = fields[idx].streamline(step, minimal[idx], idx, cell, coord = (x, y), prominent = True)
  if trace != None:
    traces.append(trace)


def show_streamline(coords, design, dir, scale, s1_minimal, style, step, s1_samples, s2_samples, max_traces_per_cell, res, do_save, info, force):

  s1 = design['s1']
  s2 = design['s2']
  angle = design['angle']

  assert(not (s1_samples == None and s2_samples != None))

  centers, min, max, elem = coords

  if scale > 0.0:
    s1 *= scale
    s2 *= scale
    s1_minimal *= scale

  # we scale minimal in such a way, that for different sampling the thinnest lines are drawn in the same thickness
  minimal = (s1_minimal, s1_minimal if s2_samples == None else (1.0 * s2_samples / s1_samples) * s1_minimal)

  # separate fields due to possibly separate sampling
  fields = (Fields(coords, s1, angle, s1_samples), Fields(coords, s2, angle, s2_samples if s2_samples != None else s1_samples))

  fig, sub = create_figure(min, max, res, do_save)


  #interpret horizontal as s1. Note that we sort must not sort, this results in jumps in the angle in the two-load case
  dirs = [1,0] # 1 is weaker, draw first and second the darker s1 lines
  if dir == 'horizontal':
    dirs = [0]
  if dir == 'vertical':
    dirs = [1]


  #generate all > minimal traces, draw (or not draw) them later
  traces = []

  if force == 'right_lower':
    for x in numpy.arange(0.9, 1.0, 0.01):
       tmp = fields[0].directional_streamline(x,0.0, fields[0].macro.dx * step, minimal[0], 0, -1.0)
       val = [i[1] for i in tmp]
       cell = (coord2index(fields[0].macro, (x, 0.0)))
       if len(val) > 0:
         traces.append(Trace([i[0] for i in tmp], val, cell, 0, prominent=True))
       else:
         print('zero forced trace at x=' + str(x))

  if force == 'rhombus': # assumes 0,0 -> 2,1 two-load case
    assert(not s1_samples == None)
    for x in numpy.arange(0.0, 0.5 * max[1], 0.2 * fields[0].macro.dx):
      force_prominent_streamline(traces, fields, minimal, step, x, 0.5 * max[1] - x)
      force_prominent_streamline(traces, fields, minimal, step, x, 0.5 * max[1] + x)
      force_prominent_streamline(traces, fields, minimal, step, max[0] - 0.5 * max[1] + x, x)
      force_prominent_streamline(traces, fields, minimal, step, max[0] - 0.5 * max[1] + x, max[1] - x - 1e-6)
      force_prominent_streamline(traces, fields, minimal, step, 0.5 * max[0] - 0.25 * max[1] + x, 0.5 * max[1])

  if False:
    for idx in dirs:
      field = fields[idx]
      for j in range(field.macro.ny):
        for i in range(field.macro.nx):
          trace = field.streamline(step, minimal[idx], idx, cell = (i,j))
          if trace != None:
            traces.append(trace)
  else:
    trace = fields[0].streamline(step, minimal[0], 0, cell = (0.2/fields[0].macro.nx, 0.2/fields[0].macro.ny))
    traces.append(trace)

  # sort by max field
  traces = sorted(traces, key=lambda trace: trace.max, reverse=True)

  # here we count the macro cells touched by the traces. For both indices two fields
  cells = (numpy.zeros(fields[0].macro.ndim), numpy.zeros(fields[1].macro.ndim))

  # conditionally draw the traces
  count = [0,0]
  for trace in traces:
    mycells = cells[trace.idx]
    field   = fields[trace.idx]

    # draw the trace only if there are not too much lines yet at the trace origin or if the drawing is enforced
    if trace.prominent or mycells[trace.cell] <= max_traces_per_cell:

      mycells += trace.touched(field.macro)
      count[trace.idx] += 1

      if style == 'line':
        draw_trace(sub, trace)
      else:
        draw_thick_trace(sub, field.macro.dx, trace, minimal[trace.idx])

  # finally some statistics

  print('drawn traces for s1: ' + str(count[0]) + ' and s2: ' + str(count[1]))
  if info != None:
    traces = xml.etree.ElementTree.SubElement(info, "drawnTraces")
    traces.set("s1", str(count[0]))
    traces.set("s2", str(count[1]))

  # see how much we macro cells are not drawn
  void_count = [0, 0]
  mat_count  = [0, 0]
  void_sum = [0.0, 0.0]
  mat_sum  = [0.0, 0.0]

  for idx in [0, 1]:
    mycells = cells[idx]
    field = fields[idx]
    for i in range(field.macro.nx):
      for j in range(field.macro.ny):
        x, y = field.index2coord((i, j))
        val = field.macro.getData(x, y)[0]
        if mycells[i,j] == 0:
          void_sum[idx] += val
          void_count[idx] += 1
        else:
          mat_sum[idx] += val
          mat_count[idx] += 1

  # print void_sum
  # print mat_sum
  print('below minimal cells (fraction) s1: ' + str(float(void_count[0])/(void_count[0] + mat_count[0])) + ' s2: ' + str(float(void_count[1])/(void_count[1] + mat_count[1])))
  print('below minimal material (fraction) s1: ' + str(void_sum[0]/(void_sum[0] + mat_sum[0])) + ' s2: ' + str(void_sum[1]/(void_sum[1] + mat_sum[1])))
  if info != None:
    stream = xml.etree.ElementTree.SubElement(info, "streamline")
    stream.set("max_traces_per_sell", str(max_traces_per_cell))
    cells = xml.etree.ElementTree.SubElement(stream, "belowMinimalCells")
    cells.set("s1", str(float(void_count[0])/(void_count[0] + mat_count[0])))
    cells.set("s2", str(float(void_count[1])/(void_count[1] + mat_count[1])))
    cells = xml.etree.ElementTree.SubElement(stream, "belowMinimalMaterial")
    cells.set("s1", str(void_sum[0]/(void_sum[0] + mat_sum[0])))
    cells.set("s2", str(void_sum[1]/(void_sum[1] + mat_sum[1])))

  return (fig, sub)
