#!/usr/bin/python
from __future__ import print_function

import numpy as np
import tempfile
import argparse
import os
from lxml import etree

import micromodels
from mesh_tool import write_gid_mesh


def write_density(densityfile, density):
# Generates a density file out of a given density matrix.
#
# @param:
#       densityfile     name of generated .dens file
#       density         density matrix (entry has to be 0 for no material)
#

    # Set density for weak material
    density[density < 1e-10] = 1e-7;

    if np.ndim(density) == 2:

        # Rotate matrix to fit element mesh numbering (lower left to upper right, line by line)
        # to matlab linear indexing (upper left to lower right, column by column)
        density = np.transpose(np.flipud(density));

        [m,n] = np.shape(density)
        density = np.valueshape(density, (m*n,))

        # Write density file
        with open(densityfile,'wb') as f:

            print('<?xml version="1.0"?>', file=f)
            print('<cfsErsatzMaterial>', file=f)
            print('  <header>', file=f)
            print('    <mesh x="{:d}" y="{:d}" z="1"/>'.format(m, n), file=f)
            print('    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>', file=f)
            print('    <transferFunction application="mech" design="density" param="1" type="simp"/>', file=f)
            print('  </header>', file=f)
            print('  <set id="4">', file=f)

            for ii in range(m*n):
                str = '    <element nr="{:d}" type="density" design="{:e}"/>'
                print(str.format(ii+1, density[ii]), file=f)

            print('  </set>', file=f)
            print('</cfsErsatzMaterial>', file=f)

    else:

        [m,n,o] = np.shape(density)
        density = np.valueshape(density, (m*n*o,))
    
        # Write density file
        fid = fopen(densityfile,'wt', file=f)
    
        print('<?xml version="1.0"?>', file=f)
        print('<cfsErsatzMaterial>', file=f)
        print('  <header>', file=f)
        print('    <mesh x="{:d}" y="{:d}" z="{:d}"/>'.format(m, n, o), file=f)
        print('    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>', file=f)
        print('    <transferFunction application="mech" design="density" param="1" type="simp"/>', file=f)
        print('  </header>', file=f)
        print('  <set id="4">', file=f)
    
        for ii in range(m*n*o):
            str = '    <element nr="{:d}" type="density" design="{:e}"/>'
            print(str.format(ii+1, density[ii]), file=f)

        print('  </set>', file=f)
        print('</cfsErsatzMaterial>', file=f)

def read_tensor(infoxmlfile):
    xml = etree.parse(infoxmlfile)
    data = xml.xpath('//homogenizedTensor/tensor/real/text()')
    values = list(map(float, data[0].split())) # convert list with string elements to list with float elements
    if len(values) == 9:
        tensor = (values[0], values[1], values[2], values[4], values[5], values[8])
    elif len(values) == 36:
        tensor = (values[0], values[1], values[2], values[7], values[8], values[14], values[21], values[28], values[35])
    
    os.remove(temporaryfile + '.info.xml')
    
    return tensor

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('params', metavar='params', type=float, nargs='+',
                        help='parameter values for microstructure')
    parser.add_argument('--res', type=int, required=True,
                        help='resolution of catalog')
    parser.add_argument('-x', type=int, default=20,
                        help='resolution of microcell in x-direction')
    parser.add_argument('-y', type=int, default=20,
                        help='resolution of microcell in y-direction')
    parser.add_argument('-z', type=int, default=20,
                        help='resolution of microcell in y-direction')
    parser.add_argument('-c', '--catalogfile',
                        help='catalogfile to write')
    parser.add_argument('-f', '--generatorfun', dest='function_name',
                        help='function to generate micro cell')
    
    args = parser.parse_args()
    
    paramfile = '/home/daniel/code/cfs/share/matlab/material_catalogue/+Homogenization/CFS_Working_Directory/inv_tensor.xml'
    
    # write mesh
    temporaryfile = tempfile.NamedTemporaryFile().name[5:]
    meshfile = os.path.abspath('.') + '/' + temporaryfile + '.mesh'

    micro_generator_function = getattr(micromodels, args.function_name)
    
    params = [float(param)/args.res for param in args.params]
    mesh, volume = micro_generator_function(params, nx=args.x, ny=args.y)
    write_ansys_mesh(mesh, meshfile)

    os.system('cfsso.rel -m ' + meshfile + ' -p ' + paramfile + ' ' + temporaryfile)

    # delete mesh
    #os.remove(meshfile)
    os.remove(temporaryfile + '.plot.dat')

    # copy tensor to catalog
    tensor = read_tensor(temporaryfile + '.info.xml')
    with open(args.catalogfile, 'wb') as f:
        str = '{:5f} ' * len(args.params)
        f.write(str.format(*args.params))
        str = '{:e} ' * len(tensor) + '\n'
        f.write(str.format(*tensor))
