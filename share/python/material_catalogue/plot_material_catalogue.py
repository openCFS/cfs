#!/usr/bin/env python
import os
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from cfs_utils import *

eps = 1e-6


def vol_cyl(x):
  return np.pi * x ** 2


def get_derivatives(coeff, samples, s1, s2, s3):
  idx_1, idx_2, idx_3 = get_interpolation_interval_idx(samples, s1, s2, s3)
  # flattened interval index
  iv_idx = np.ravel_multi_index((idx_1[0], idx_2[0], idx_3[0]), (len(samples[0]) - 1, len(samples[1]) - 1, len(samples[2]) - 1))

  # deriv has 64 entries
  # order: f(s1,s2,s3), dEda, dEdb, dEdc, dEdadb, dEdadc, dEdbdc, dEdadbdc
  deriv = coeff[iv_idx]
  # if np.isclose(s1, 0) and np.isclose(s3, 0) and np.isclose(s2, 1 / 3.0):
#   print("s1,s2,s3:", s1, s2, s3, "  iv_index:", iv_idx, deriv[0:8])
  return deriv[0:8], deriv[8:16], deriv[16:24], deriv[24:32], deriv[32:40], deriv[40:48], deriv[48:56], deriv[56:64] 
  

def tangent_line(slope, x, x1, y1):
  return slope * (x - x1) + y1


def get_interpolation_interval_idx(samples, s1, s2, s3):
  assert(0 <= s1 <= 1)
  assert(0 <= s2 <= 1)
  assert(0 <= s3 <= 1)
  # preset with last interval in samples list 
  bounds_s1 = bounds_s2 = bounds_s3 = [len(samples[0]) - 2, len(samples[0]) - 1]
  # s1
  for i in range(len(samples[0]) - 2):
    # find first larger sample value
    if samples[0][i] <= s1 + eps <= samples[0][i + 1]:
      bounds_s1 = [i, i + 1]
    if samples[1][i] <= s2 + eps <= samples[1][i + 1]:
      bounds_s2 = [i, i + 1]
    if samples[2][i] <= s3 + eps <= samples[2][i + 1]:
      bounds_s3 = [i, i + 1]  

#   print("\ns1: ", s1, " bounds:", bounds_s1)
#   print("s2: ", s2, " bounds:", bounds_s2)
#   print("s3: ", s3, " bounds:", bounds_s3)
  return bounds_s1, bounds_s2, bounds_s3


# for each parameter (s1,s2,s3), project it onto interval [0,1] and return normalized parameter value
# ls1: lower bound of original interval that contains s1, e.g. 0.1 for interval  [0.1,0.2] and s1 = 0.13
# da: length of original interval of s1
def normalize_parameter(s1, s2, s3, ls1, ls2, ls3, da, db, dc):
  s1_norm = (s1 - ls1) / da
  s2_norm = (s2 - ls2) / db
  s3_norm = (s3 - ls3) / dc
  
#   print("s1,s2,s3:", s1, s2, s3, " ls1,ls2,ls3:", ls1, ls2, ls3, " da, db, dc:", da, db, dc, " s1_norm, s2_norm, s3_norm:", s1_norm, s2_norm, s3_norm)
  
  return s1_norm, s2_norm, s3_norm


# returns interpolated tensor at point [s1,s2,s3]
def eval_interpolation(dim, samples, coeffs, da, db, dc, s1, s2, s3):
  assert(len(coeffs) == 4)
  k11 = eval_interpolated_point(dim, samples, coeffs[0], da, db, dc, s1, s2, s3)
  k22 = eval_interpolated_point(dim, samples, coeffs[1], da, db, dc, s1, s2, s3) 
  k33 = eval_interpolated_point(dim, samples, coeffs[2], da, db, dc, s1, s2, s3)
  vol = eval_interpolated_point(dim, samples, coeffs[3], da, db, dc, s1, s2, s3)
#   print("s1,s2,s3:", s1, s2, s3, "    k11,k22,k33:", k11, k22, k33, " vol:", vol)  
#   assert(k11 >= -eps and k22 >= -eps and k33 >= -eps and vol >= -eps)
  
#   print("s1,s2,s3:",s1,s2,s3," interp volume=", vol)
  return np.array([k11, k22, k33]), vol


def eval_interpolated_point(dim, samples, coeffs, da, db, dc, s1, s2, s3):
   # maps s1, s2, s3 onto interval [0,1]
  # s1_norm = (s1 - lower) / (upper-lower)
  idx_1, idx_2, idx_3 = get_interpolation_interval_idx(samples, s1, s2, s3)
  s1_norm, s2_norm, s3_norm = normalize_parameter(s1, s2, s3, samples[0][idx_1[0]], samples[1][idx_2[0]], samples[2][idx_3[0]], da, db, dc)
  
  if not(0 <= s1_norm <= 1 + eps and 0 <= s2_norm <= 1 + eps and 0 <= s3_norm <= 1 + eps):
    print(s1_norm, s2_norm, s3_norm)
  assert(0 <= s1_norm <= 1 + eps and 0 <= s2_norm <= 1 + eps and 0 <= s3_norm <= 1 + eps)
#   print("idx_1, idx_2, idx_3:", idx_1, idx_2, idx_3)
#   print("da, db, dc:", da, db, dc)
#   print("s1,s2,s3 normed:", s1_norm, s2_norm, s3_norm)  
  res = 0
#   print("s1,s2,s3:", s1, s2, s3)
#   print("idx1,idx2,idx3:", idx_1, idx_2, idx_3)
#   print("dim:", dim)
  idx = np.ravel_multi_index((idx_1[0], idx_2[0], idx_3[0]), (len(samples[0]) - 1, len(samples[1]) - 1, len(samples[2]) - 1))
#   print("flattened interval idx:", idx)

  for i in range(4):
    for j in range(4):
      for k in range(4):
        # ravel_multi_index computes index of linearized array of size dim*dim
        # coeffs contains the interpolation coefficients for each sample interval
        # idx_1[0] the interval number
        # idx is the flattened index of idx_1[0] 
        cidx = np.ravel_multi_index((i, j, k), (4, 4, 4), order="F") 
        summand = coeffs[idx][cidx] * s1_norm ** i * s2_norm ** j * s3_norm ** k
        res = res + summand
#         if j == 0 and k == 0:
#         print("i,j,k:", i, j, k, " coeff_idx:", cidx , " coeff:", coeffs[idx][cidx], " summand:", summand, " s1_norm**i:", s1_norm ** i, " s2_norm**j:", s2_norm ** j, " s3_norm**k:", s3_norm ** k, " res=", res)
        
#   sys.exit()
  return res


# fix_{s1,s2,s3}: three bools, indicating which parameter to fix
# val: value of fixed parameter
# plot_dir: where to put all the plots
def plot_surfaces(fix_s1, fix_s2, fix_s3, val, samples, coeffs, plot_dir):
  assert(os.path.exists(plot_dir))
  n_samples = samples.shape[1]
  if np.sum([fix_s1, fix_s2, fix_s3]) != 1:
    sys.exit("fix only one parameter!")
    
  # name of fixed parameter
  f_name = "s1"
  if fix_s2:
    f_name = "s2"
  elif fix_s3:
    f_name = "s3"  
  
  # parameter range  
  pr = np.arange(0, 1.01, 0.01)  
  k11 = np.zeros((len(pr), len(pr)))
  k22 = np.zeros((len(pr), len(pr)))
  k33 = np.zeros((len(pr), len(pr)))
  vol = np.zeros((len(pr), len(pr)))
#   for s2 in [0,1/3.,2/3.,1]:
#     vals = []
  for sc in np.arange(0, 1.1, 0.1):
    for i, sa in enumerate(pr):
      for j, sb in enumerate(pr):
        if fix_s1:
          tens, v = eval_interpolation(n_samples - 1, samples, coeffs, da, db, dc, sc, sa, sb)
        elif fix_s2:
          tens, v = eval_interpolation(n_samples - 1, samples, coeffs, da, db, dc, sa, sc, sb)
#           if j == 0 and np.isclose(sc, 0):
#             print("s1,s2,s3:", sa, sc, sb, " k11:", tens[0, 0])
        else:
          assert(fix_s3)
          tens, v = eval_interpolation(n_samples - 1, samples, coeffs, da, db, dc, sa, sb, sc)
          
        k11[i, j] = tens[0]
        k22[i, j] = tens[1]
        k33[i, j] = tens[2]
        vol[i, j] = v
#         if not (v >= -eps):
#           print("sa,sb,sc:", sa, sb, sc, " vol:", v)
#         assert(v >= -eps)
        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_title(f_name + " = " + str(round(sc, 3)))
    if fix_s3:
      ax.set_xlabel("s1")
      ax.set_ylabel("s2")
    elif fix_s2:
      ax.set_xlabel("s1")
      ax.set_ylabel("s3")
    else:
      assert(fix_s1)
      ax.set_xlabel("s2")
      ax.set_ylabel("s3")
        
    ax.set_zlabel("k11")
    # ax.plot_surface(s1, s1, np.array(vals),label="k11")
    x, y = np.meshgrid(pr, pr, indexing="ij")
    from matplotlib import cm
    ax.plot_surface(x, y, k11, label="k11", cmap=cm.viridis)
#     plt.show()
#     if np.isclose(s3, 0.2):
#       plt.show()
      
    plt.savefig(plot_dir + 'k11_' + f_name + '-' + str(round(sc, 3)) + '.png', dpi=300)
#     if fix_s2:
#       np.savetxt("k11.txt", k11)
#       fig2, axes2 = plt.subplots()
#       axes2.plot(pr, k11[:, 0])
#       plt.show()
    # delete stored data set only
    del ax.collections[0]     
    
    ax.plot_surface(x, y, k22, label="k22", cmap=cm.viridis)
    ax.set_zlabel("k22")
    plt.savefig(plot_dir + 'k22_' + f_name + '-' + str(round(sc, 3)) + '.png', dpi=300)
    # delete stored data set only
    del ax.collections[0]
    
    ax.plot_surface(x, y, k33, label="k33", cmap=cm.viridis)
    ax.set_zlabel("k33")
    plt.savefig(plot_dir + 'k33_' + f_name + '-' + str(round(sc, 3)) + '.png', dpi=300)
    del ax.collections[0]
    
    ax.plot_surface(x, y, vol, label="vol", cmap=cm.viridis)
    ax.set_zlabel("vol")
    plt.savefig(plot_dir + 'vol' + f_name + '-' + str(round(sc, 3)) + '.png', dpi=300)  
    plt.close()

  plt.close()


# helpfer func for plot_lines
# depending on which parameter (s1,s2,s3) we want to vary (flag v_s1/s2/s3),
# return s1,s2,s3 in correct order by assuming that s1 is the varying one  
def get_ordered_param(v_s1, v_s2, v_s3, s1, s2, s3):
  order = ()
  assert(np.sum([v_s1, v_s2, v_s3]) == 1)
  if v_s1:
    order = (s1, s2, s3)
  elif v_s2:
    order = (s3, s1, s2)
  else:
    assert(v_s3)  
    order = (s2, s3, s1)

  return order


# flags v_s1, v_s2, v_s3: which parameter to vary from 0 to 1    
def plot_lines(v_s1, v_s2, v_s3, samples, coeffs, coeff_deriv, plot_dir):
  # respective index in array with derivatives
  deriv_idx = 0
  order = None
  title = ""
    
  pr = np.arange(0, 1.01, 0.01)
  
  # assume equidistant sampling
  for sc in samples[0]:
    for sb in samples[1]:
      fig, axes = plt.subplots(4)
      
      val11 = []
      val22 = []
      val33 = []
      vol = []
      ref = []
      filename = plot_dir
      if v_s1:
        deriv_idx = 1
        title = "s1=[0,1] s2= " + str(round(sb, 3)) + " s3=" + str(round(sc, 3))
        filename += 's1-x_s2-' + str(round(sb, 3)) + '_s3-' + str(round(sc, 3)) + '.png'
      elif v_s2:
        deriv_idx = 2
        title = "s1= " + str(round(sc, 3)) + " s2=[0,1]" + " s3=" + str(round(sb, 3))
        filename += 's2-x_s1-' + str(round(sc, 3)) + '_s3-' + str(round(sb, 3)) + '.png'
      else:
        assert(v_s3)
        deriv_idx = 3
        title = "s1= " + str(round(sb, 3)) + " s2=" + str(round(sc, 3)) + " s3=[0,1]"
        filename += 's3-x_s1-' + str(round(sb, 3)) + '_s2-' + str(round(sc, 3)) + '.png'
      assert(deriv_idx > 0)  
      for sa in pr:
        order = get_ordered_param(v_s1, v_s2, v_s3, sa, sb, sc)
        tens, v = eval_interpolation(n_samples - 1, samples, coeffs, da, db, dc, *order)
          
        val11.append(tens[0])
        val22.append(tens[1])
        val33.append(tens[2])
        vol.append(v)
        ref.append(vol_cyl(sa / 2.0))
        if args.heat:
          assert(-eps <= tens[0] <= 1 + eps)
          assert(-eps <= tens[1] <= 1 + eps)
          assert(-eps <= tens[2] <= 1 + eps)
  
      values = [val11, val22, val33, vol]
      # title over all subplots
      fig.suptitle(title)
      # adjust vertical space between subplots
      plt.subplots_adjust(hspace=0.3)
      plt.xlabel("s" + str(deriv_idx))
      for axi in range(4):
        axes[axi].set_ylim((0, 1.0))
        if axi == 3:  # volume
          axes[axi].plot(pr, values[axi], label="volume")
        else:
          axes[axi].plot(pr, values[axi], label="k" + str(axi + 1) + str(axi + 1))
      
      tmp = []
#       fig2, ax2 = plt.subplots(1)
      ###### plot tangent lines at samples
      # as s3 = 0, we skip third dimension
      for i, sa in enumerate(samples[0]):  # [::2]):
        order = get_ordered_param(v_s1, v_s2, v_s3, sa, sb, sc)
        # f, dEda, dEdb, dEdc, dEdadb, dEdadc, dEdbdc, dEdadbdc = get_derivatives(coeff_11_deriv, samples, s1, s2, s3)
        out_11 = get_derivatives(coeff_deriv[0], samples, *order)
        out_22 = get_derivatives(coeff_deriv[1], samples, *order)
        out_33 = get_derivatives(coeff_deriv[2], samples, *order)
        out_vol = get_derivatives(coeff_deriv[3], samples, *order)
        
        # which first order derivative to plot depending on 'order'
        # out[1] -> dEda, out[2] -> dEdb, out[3] -> dEdc
        dEds_11 = out_11[deriv_idx]
        dEds_22 = out_22[deriv_idx]
        dEds_33 = out_33[deriv_idx]
        dEds_vol = out_vol[deriv_idx]
        
        # we don't have an interval with lower bound = 1.0
        # thus, in case s1/s2/s3 = 1.0 -> take interval with upper bound = 1.0 and respective function value
        offset_sa = 0 if not np.isclose(order[0], 1.0) else 1
        offset_sb = 0 if not np.isclose(order[1], 1.0) else 1
        offset_sc = 0 if not np.isclose(order[2], 1.0) else 1
#         print("offset:", offset_sa, offset_sb, offset_sc)
        
        # calc flattened array index (0,...,8)
        idx = np.ravel_multi_index((offset_sa, offset_sb, offset_sc), (2, 2, 2), order="F")
#         print("offset_idx:", idx)
        # function value at point s1
        function_values = (out_11[0][idx], out_22[0][idx], out_33[0][idx], out_vol[0][idx])
#         if np.isclose(sb, 0) and np.isclose(sc, 0) and np.isclose(sa, 1 / 3.0):
#           print("s1,s2,s3:" + str(sb) + "," + str(sc) + "," + str(sa) + " deriv_idx:" + str(deriv_idx) + "  k_11=" + str(k_11) + "  k_22=" + str(k_22) + "  k_33=" + str(k_33))
#           sys.exit()
        if args.heat:
          for f in function_values:
            assert(-eps <= f <= 1 + eps)
          
        slopes = (dEds_11[idx], dEds_22[idx], dEds_33[idx], dEds_vol[idx])
        tmp.append(slopes[0])
        
        # Define x data range for tangent line
        xrange = np.linspace(sa - 0.01, sa + 0.01, 10)  
        tangents = []
        # slope of dEds at point sa for coeff_11
        for id in range(len(slopes)):
          tangents.append(tangent_line(slopes[id], xrange, sa, function_values[id]))
        
        for aidx, ax in enumerate(axes):
          label_s = "samples" if i == 0 and aidx == 0 else ""
          label_t = "tangents" if i == 0 and aidx == 0 else ""
          ax.scatter(sa, function_values[aidx], marker="x", color="green", label=label_s)
          ax.plot(xrange, tangents[aidx], linewidth=2, color="red", label=label_t)
          ax.legend(loc="upper left")
        
#         ax2.plot(sa, function_values[0], marker="x", color="green", label=label_s)
#         ax2.plot(xrange, tangents[0], linewidth=2, color="red", label=label_t)  
#         print("sa,sb,sc:", sa, sb, sc, "f:", function_values[0], " dEda:", slopes[0])
#         ax2.plot(xrange, tangents[0], linewidth=2, color="red", label=label_t)  
  #     plt.show()
#       if v_s1 :  # and np.isclose(sb, 0) and np.isclose(sc, 0):
#         ax2.plot(pr, values[0], label="k11")
#         ax2.plot(pr, ref, linewidth=2, color="red", label="ref")  
#         print("slopes:", tmp)
#         plt.show()    
#         fig2.show()
#         fig.show()
      fig.savefig(filename, dpi=300)
      plt.clf()
      plt.close()
  
    
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("xml", help="name of xml file with interpolation coefficients", default="")
  parser.add_argument("--dim", help="geometrical dimension (default:3)", default=3, required=False, type=int)
  parser.add_argument("--eval", help="return interpolated tensor at given point, e.g. 0.1, 0.2, 0.3")
  parser.add_argument("--heat", help="flag for heat application (default=false)", action="store_true", default=False)
  args = parser.parse_args()
  
  if not os.path.exists(args.xml):
    sys.exit("Cannot find file " + args.xml)
    
  if args.dim != 3:
    sys.exit("Currently only implemented for 3d")
    
  if args.eval:
    point = args.eval.split(',')
    point = [float(p) for p in point]
    assert(len(point) == 3)  
  
  # equidistant design parameter samples in x,y and z direction
  xml = open_xml(args.xml)
  n_samples = int(xpath(xml, "//a/matrix/@dim1"))
  samples = np.zeros((args.dim, n_samples))
  
  # read parameter values for s1
  s_str = xpath(xml, "//a/matrix/real/text()")
  s = np.asarray(list(map(float, s_str.split())))
  assert(len(s) > 1)    
  samples[0, :] = s
  da = s[1] - s[0]
  
  # read parameter values for s2
  s_str = xpath(xml, "//b/matrix/real/text()")
  s = np.asarray(list(map(float, s_str.split())))    
  samples[1, :] = s
  db = s[1] - s[0]
  
  # read parameter values for s3
  s_str = xpath(xml, "//c/matrix/real/text()")
  s = np.asarray(list(map(float, s_str.split())))    
  samples[2, :] = s
  dc = s[1] - s[0]
  
  # for each interval between to sample values, we need one interpolation function
  # tricubic interpolation: for each sample we need 64 coefficients
  n_coeffs = 64
  # number of intervals between samples
  n_intvals = (n_samples - 1) ** 3
  dim = args.dim
  
  # k11, k22, k33, vol
  coeffs = []
  coeffs_deriv = []
  
  # we have 3 material coefficients per tensor
  matrix = xpath(xml, "//coeff11/matrix/real/text()")    
  coeff_11 = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(coeff_11.shape[0] == n_intvals * n_coeffs)
  coeffs.append(coeff_11.reshape((n_intvals, n_coeffs)))
  # all the derivatives used for calculating coeff_11
  # contains 8 vectors where each has 8 entries - order: f(x), dEda, dEdb, dEdc, dEdadb, dEdadc, dEdbdc, dEdadbdc
  matrix = xpath(xml, "//coeff11_deriv/matrix/real/text()")    
  coeff_11_deriv = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(coeff_11_deriv.shape[0] == n_intvals * n_coeffs)
  coeffs_deriv.append(coeff_11_deriv.reshape((n_intvals, n_coeffs)))
  
  # we have 3 material coefficients per tensor
  matrix = xpath(xml, "//coeff22/matrix/real/text()")    
  coeff_22 = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(coeff_22.shape[0] == n_intvals * n_coeffs)
  coeffs.append(coeff_22.reshape((n_intvals, n_coeffs)))
  # all the derivatives used for calculating coeff_22
  # contains 8 vectors where each has 8 entries - order: f(x), dEda, dEdb, dEdc, dEdadb, dEdadc, dEdbdc, dEdadbdc
  matrix = xpath(xml, "//coeff22_deriv/matrix/real/text()")    
  coeff_22_deriv = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(coeff_22_deriv.shape[0] == n_intvals * n_coeffs)
  coeffs_deriv.append(coeff_22_deriv.reshape((n_intvals, n_coeffs)))
  
  # we have 3 material coefficients per tensor
  matrix = xpath(xml, "//coeff33/matrix/real/text()")    
  coeff_33 = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(coeff_33.shape[0] == n_intvals * n_coeffs)
  coeffs.append(coeff_33.reshape((n_intvals, n_coeffs)))
  # all the derivatives used for calculating coeff_33
  # contains 8 vectors where each has 8 entries - order: f(x), dEda, dEdb, dEdc, dEdadb, dEdadc, dEdbdc, dEdadbdc
  matrix = xpath(xml, "//coeff33_deriv/matrix/real/text()")    
  coeff_33_deriv = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(coeff_33_deriv.shape[0] == n_intvals * n_coeffs)
  coeffs_deriv.append(coeff_33_deriv.reshape((n_intvals, n_coeffs)))
  
  matrix = xpath(xml, "//volcoeff/matrix/real/text()")
  volcoeff = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(volcoeff.shape[0] == n_intvals * n_coeffs)
  coeffs.append(volcoeff.reshape((n_intvals, n_coeffs)))
#   all the derivatives used for calculating volcoeff
#   contains 8 vectors where each has 8 entries - order: f(x), dEda, dEdb, dEdc, dEdadb, dEdadc, dEdbdc, dEdadbdc
  matrix = xpath(xml, "//volcoeff_deriv/matrix/real/text()")    
  volcoeff_deriv = np.asarray(list(map(float, matrix.split())))  # convert list with string elements to list with float elements
  assert(volcoeff_deriv.shape[0] == n_intvals * n_coeffs)
  coeffs_deriv.append(volcoeff_deriv.reshape((n_intvals, n_coeffs)))
  
#   assert(len(coeffs) == 4 and len(coeffs_deriv) == 4)
  
  cwd = os.getcwd()
  plot_dir = cwd + "/plots/"
  if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)
    
  if args.eval:
    print("tensor at point ", point, ":\n", eval_interpolated_point(args.dim, samples, coeffs[0], da, db, dc, point[0], point[1], point[2]))
    sys.exit()

  # plot surface for varying s1 and s2 - fix s3=0
  plot_surfaces(1, 0, 0, 0.0, samples, coeffs, plot_dir)
  plot_surfaces(0, 1, 0, 0.0, samples, coeffs, plot_dir)
  plot_surfaces(0, 0, 1, 0.0, samples, coeffs, plot_dir)
   
  plot_lines(1, 0, 0, samples, coeffs, coeffs_deriv, plot_dir)
  plot_lines(0, 1, 0, samples, coeffs, coeffs_deriv, plot_dir)
  plot_lines(0, 0, 1, samples, coeffs, coeffs_deriv, plot_dir)

