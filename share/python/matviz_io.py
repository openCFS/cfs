from PIL import Image
from cfs_utils import open_xml, xpath
from hdf5_tools import get_element, has_element
from matviz_vtk import show_write_vtk
import numpy
import os.path
import xml.etree.ElementTree
import xml.dom.minidom

import matplotlib
## reads design_stiff*, design_shear* and design_rotAngle* for 2D and 3D. Fills other stuff by defaults

# considers density read
# @param angle array of anglex, angley, anglez
# @return dict with design variables (s1, s2, s3, sh1, angle)
def read_design(hdf_file, args, dim_2D, centers, TWO_SCALE):
  res = dict()
  # rot means, that we only show rotation according to rotAngle, e.g. for piezoelectric polarization
  sh1 = numpy.ones((len(centers),1)) * .5 # fix for no shearing
  if args.parametrization == 'hom_rect':
    if args.show in TWO_SCALE:
      try:
        res['microparams'] = [get_element(hdf_file, "design_microparam{}_{}".format(i + 1, args.access), args.h5_region, args.h5_step) for i in range(TWO_SCALE[args.show])]
      except:
        try:
          res['s1'] = get_element(hdf_file, "design_stiff1_" + args.access, args.h5_region, args.h5_step)
          res['s2'] = get_element(hdf_file, "design_stiff2_" + args.access, args.h5_region, args.h5_step)
          res['s3'] = get_element(hdf_file, "design_stiff3_" + args.access, args.h5_region, args.h5_step) if not dim_2D else numpy.ones((len(centers),1)) * .1
          res['sh1'] = get_element(hdf_file, "design_shear1_" + args.access, args.h5_region, args.h5_step) if (args.show == "hom_sheared_rot_cross" or args.show == "hom_cross_bar") else sh1
        except:
          # this is for two scale results with just one parameter
          detail = "physicalPseudoDensity" if args.access == 'smart' else "mechPseudoDensity"
          if args.h5_region == 'all':
            s1 = [[None]]
            for region in hdf_file['/Mesh/Regions']:
              s = get_element(hdf_file, detail, region, args.h5_step)
              s1 = numpy.concatenate((s1,s))
            res['s1'] = s1[1:]
          else:
            res['s1'] = get_element(hdf_file, detail, args.h5_region, args.h5_step)
          res['s2'] = res['s1']
          res['s3'] = res['s1']
    else:
      detail = "physicalPseudoDensity" if args.access == 'smart' else "mechPseudoDensity"
      if args.h5_region == 'all':
        s1 = [[None]]
        for region in hdf_file['/Mesh/Regions']:
          s = get_element(hdf_file, detail, region, args.h5_step)
          s1 = numpy.concatenate((s1,s))
        res['s1'] = s1[1:]
      else:
        res['s1'] = get_element(hdf_file, detail, args.h5_region, args.h5_step)
      res['s2'] = res['s1']
      res['s3'] = res['s1']
  elif args.parametrization == 'trans-iso':
    s1 = get_element(hdf_file, "design_emodul-iso_" + args.access, args.h5_region, args.h5_step)
    s2 = get_element(hdf_file, "design_emodul_" + args.access, args.h5_region, args.h5_step)
    try:
      theta = get_element(hdf_file, "design_poisson_" + args.access, args.h5_region, args.h5_step)
    except:
      theta = 0.0
      print('could not read theta (design_poisson_' + args.access + '), setting to ' + str(theta))
    m = 2.0 * numpy.max([numpy.max(s1), numpy.max(s2)])
    s1 *= 1 / (m * (1 - theta))
    s2 *= 1 / (m * (1 - theta))
    res['s1'] = s1
    res['s2'] = s2
    res['s3'] = numpy.ones((len(centers), 1)) * .1  # fix for 3D
  elif args.parametrization == 'ortho':
    t11 = get_element(hdf_file, "design_tensor11_" + args.access, args.h5_region, args.h5_step)
    t12 = get_element(hdf_file, "design_tensor12_" + args.access, args.h5_region, args.h5_step)
    t22 = get_element(hdf_file, "design_tensor22_" + args.access, args.h5_region, args.h5_step)
    t33 = get_element(hdf_file, "design_tensor33_" + args.access, args.h5_region, args.h5_step)
    s1 = t11*t11+t12*t12
    s2 = t12*t12+t22*t22
    m = 2.0*numpy.max([numpy.max(s1), numpy.max(s2)])
    s1 *= 1/m
    s2 *= 1/m
    res['s1'] = s1
    res['s2'] = s2
    res['s3'] = numpy.ones((len(centers),1)) * .1 # fix for 3D
  elif args.parametrization == "hom_iso":
    # isotropic homogenized basecell e.g. lufo fuller or V7 base cell
    res['s1'] = get_element(hdf_file, "design_stiff1_" + args.access, args.h5_region, args.h5_step)
    res['s2'] = res['s1']
    res['s3'] = res['s1']
  elif args.parametrization == 'simp':
    detail = "physicalPseudoDensity" if args.access == 'smart' else "mechPseudoDensity"
    if args.h5_region == 'all':
      s1 = [[None]]
      for region in hdf_file['/Mesh/Regions']:
        s = get_element(hdf_file, detail, region, args.h5_step)
        s1 = numpy.concatenate((s1,s))
      res['s1'] = s1[1:]
      shape(res['s1'])
    else:
      res['s1'] = get_element(hdf_file, detail, args.h5_region, args.h5_step)
    res['s2'] = res['s1']
    res['s3'] = res['s1']
    res['angle'] = numpy.zeros(((len(res['s1']), 3)))
    return res
  if has_element(hdf_file, "design_density_" + args.access):
    print("args.h5_step:" + str(args.h5_step))
    rho = get_element(hdf_file, "design_density_" + args.access, args.h5_region, args.h5_step)
    rho = pow(rho, float(args.penalty))
    res['s1'] *= rho
    res['s2'] *= rho
    res['s3'] *= rho
    print("scale stiffness values by design_density_" + args.access + " with average value " + str(numpy.mean(rho)) + " and penalty " + str(args.penalty))

  angle = numpy.zeros(((len(list(res.values())[0]), 3)))

  if args.show == "hom_rot_cross" or args.show == "hom_sheared_rot_cross" or args.show == "stream" or args.show == 'hom_rect_mod':
    try:
      if dim_2D:
        try:
          angle[:,0] = get_element(hdf_file, "design_rotAngle_" + args.access, args.h5_region, args.h5_step)[:,0]
        except:
          print('could not read design_rotAngle_' + args.access + ', trying design_rotAngle_plain')
          angle[:,0] = get_element(hdf_file, "design_rotAngle_plain", args.h5_region, args.h5_step)[:,0]
      else:
        angle[:, 0] = get_element(hdf_file, "design_rotAngleX_" + args.access, args.h5_region, args.h5_step)[:, 0]
        angle[:, 1] = get_element(hdf_file, "design_rotAngleY_" + args.access, args.h5_region, args.h5_step)[:, 0]
        angle[:, 2] = get_element(hdf_file, "design_rotAngleZ_" + args.access, args.h5_region, args.h5_step)[:, 0]
    except Exception as e:
      print('could not read angle, ignore it: ', e)
  res['angle'] = angle

  return res

def read_hom_tensor_from_info_xml(infoXmlName):
  assert(os.path.exists(infoXmlName))
  xml = open_xml(infoXmlName)
  dim = xpath(xml, "//grid/@dimensions")
  matrix = xpath(xml, "//iteration[last()]/homogenizedTensor/tensor/real/text()") # "text()" must be added due to lxml, otherwise matrix is just a string <Element real ...>
  res = list(map(float, matrix.split())) # convert list with string elements to list with float elements
  res = numpy.asarray(res)            # convert list to array
  if dim == '2':
    dim_2D = True
    res = res.reshape(3,3)         # reshaping array
    hom_tensor = [res[0][0],res[1][1],res[2][2],res[1][2],res[0][2],res[0][1]]
  else:
    assert(dim == '3')
    dim_2D = False
    res = res.reshape(6,6)         # reshaping array
    hom_tensor = [res[0][0],res[1][1],res[2][2],res[3][3],res[4][4],res[5][5],res[1][2],res[0][2],res[0][1],res[2][3],res[1][3],res[0][3],res[3][4],res[2][4],res[1][4],res[0][4],res[4][5],res[3][5],res[2][5],res[1][5],res[0][5]]
    # hom_tensor = [res[0][0],res[0][1],res[1][1],res[0][2],res[1][2],res[2][2],res[0][3],res[1][3],res[2][3],res[3][3],res[0][4],res[1][4],res[2][4],res[3][4],res[4][4],res[0][5],res[1][5],res[2][5],res[3][5],res[4][5],res[5][5]]

  return hom_tensor, dim_2D

# # show or write either Image or polydata
# @param viz either Image or polydata
# @save filename for output
# @return the volume fraction if determined or None
def show_or_write(viz, args):
  assert(viz is not None)
  global infoxml
  volume = None
  if isinstance(viz, Image.Image):
    # print 'array' + str(numpy.array(viz))
    volume = 1 - numpy.average(numpy.array(viz)) / 255
    print('volume fraction from image : ' + str(volume))
    if infoxml != None:
      vol = xml.etree.ElementTree.SubElement(infoxml, "volume")
      vol.set("imageMaterial", str(volume))
    if args.save:
      viz.save(args.save)
    else:
      viz.show()
  elif isinstance(viz, tuple):
    fig = viz[0]
    sub = viz[1]
    if args.save:
      extent = sub.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
      print('write file: ' + args.save)
      fig.savefig(args.save, bbox_inches=extent)
      if args.save.split('.')[-1] != 'pdf':
        # I was not able to render a memory image first, make an array out of the data and determine the grayness
        # So read again from file :(
        tmp = Image.open(args.save).convert('L')  # make gray, otherwise data has the dimension x*y*4 (rgb + alpha)
        dat = numpy.array(tmp)
        volume = len(numpy.where(dat.reshape(dat.size, 1) < 128)[0]) / float(dat.size)  # count fields below 128 which is become black with a threshold of 0.5
        print('volume fraction from image : ' + str(volume))
        if infoxml != None:
          vol = xml.etree.ElementTree.SubElement(infoxml, "volume")
          vol.set("imageMaterial", str(volume))
    else:
      matplotlib.pyplot.show()
      ax = fig.add_axes([0, 0, 1, 1])
      fig.show()  # Jannis: this is a temporary workaround as matplotlib.pyplot.show() does nothing for me

    matplotlib.pyplot.close(fig)

  else:
    show_write_vtk(viz, args.res, args.save, camera_settings=args.cam)

  return volume

def create_info_xml():
  global infoxml
  infoxml = xml.etree.ElementTree.Element("matviz")

def write_info_xml(file):
  print('write info xml file: ' + file)
  out = open(file, "w")
  out.write(xml.dom.minidom.parseString(xml.etree.ElementTree.tostring(infoxml)).toprettyxml())
  out.close()


def write_angle_data(file, angle, data):
  fid = open(file, "w")
  fid.write('#angle\tdata\n')
  if numpy.ndim(data) == 1:
    angle = [angle]
    data = [data]
  
  for angles, datas in zip(angle,data):
    for ang, dat in zip(angles,datas):
      fid.write(str(ang) + '\t' + str(dat) + '\n')
  fid.close()

# in this variable we can store meta-information to be exported as xml file
infoxml = None