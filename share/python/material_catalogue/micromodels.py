#!/usr/bin/python
from __future__ import print_function

import numpy as np
import math
import os.path

from mesh_tool import Mesh, Element, QUAD4

def generate_cross(param, nx=128, ny=128, sparse=True):
    
    s1 = param[0]
    s2 = param[1]
    
    volume = s1 + s2 - s1 * s2
    
    # Nodes are separated into three parts:
    # 0 to remainder / remainder to remainder + s / remainder + s to 1
    # Thus the second part has exactly the width s.
    
    yremainder = (1-s1)/2.0;
    xremainder = (1-s2)/2.0;
    
    # points in each part
    npx1 = int(math.ceil(xremainder*nx)+1)
    npx2 = int(math.ceil(s2*nx)+1)
    npx3 = npx1
    npy1 = int(math.ceil(yremainder*ny)+1)
    npy2 = int(math.ceil(s1*ny)+1)
    npy3 = npy1
    
    # coordinates of points in each part
    x1 = np.linspace(0,xremainder,npx1)
    x2 = np.linspace(x1[-1],x1[-1]+s2,npx2)
    x3 = np.linspace(x2[-1],1,npx3)
    x = np.concatenate((x1,x2[1:],x3[1:]))

    y1 = np.linspace(0,yremainder,npy1)
    y2 = np.linspace(y1[-1],y1[-1]+s1,npy2)
    y3 = np.linspace(y2[-1],1,npy3)
    y = np.concatenate((y1,y2[1:],y3[1:]))
    
    # new number of elements
    npx = len(x)
    npy = len(y)
    nx = npx - 1
    ny = npy - 1
    numNodes = npx*npy
    
    mesh = Mesh(nx,ny)

    # node coordinates
    xcoords = np.tile(x,(1,npy))[0]
    ycoords = np.reshape(np.tile(y,(npx,1)),(len(y)*npx,),'F')
    nodes = np.array((xcoords, ycoords))
    for i, row in enumerate(np.transpose(nodes)):
        mesh.nodes[i] = row
    
    # boundary nodes
    if not sparse or s1 == 1 or s2 == 1:
        left = range(0, numNodes-1, npx)
        right = range(npx-1, numNodes, npx)
        bottom = range(npx)
        up = range(numNodes-npx, numNodes)
    else:
        left = range((npy1-1)*npx, (npy1+npy2-2)*npx+1, npx)
        right = range(npy1*npx-1, (npy1+npy2-1)*npx, npx)
        bottom = range(npx1-1, npx1+npx2-1)
        up = range(npx*(npy-1)+npx1-1, npx*(npy-1)+npx1+npx2-1)
    mesh.bc.append(("left", left))
    mesh.bc.append(("right", right))
    mesh.bc.append(("bottom", bottom))
    mesh.bc.append(("up", up))
    
    # elements
    num2DElements = nx*ny
    for yy in range(ny):
        for xx in range(nx):
            e = Element()
            e.type = QUAD4

            ll = npx*yy + xx # lower-left of current element
            e.nodes = ((ll, ll+1, ll+npx+1, ll+npx))
            
            if yy < len(y1)-1 or yy > len(y1)+len(y2)-3:
                if xx < len(x1)-1 or xx > len(x1)+len(x2)-3:
                    e.region = 'void'
                else:
                    e.region = 'mech'
            else:
                e.region = 'mech'
            
            if sparse and e.region == 'void':
                continue

            mesh.elements.append(e)

    return mesh, volume

def generate_sheared_cross(param, nx=128, ny=128, sparse=True):
    
    s1 = param[0]
    s2 = param[1]
    s3 = (.5 - param[2]) * np.pi
    
    volume = s1 + s2 - s1 * s2
    
    # Nodes are separated into three parts:
    # 0 to remainder / remainder to remainder + s / remainder + s to 1
    # Thus the second part has exactly the width s.
    
    yremainder = (1-s1)/2.0;
    xremainder = (1-s2)/2.0;
    
    # points in each part
    npx1 = int(math.ceil(xremainder*nx)+1)
    npx2 = int(math.ceil(s2*nx)+1)
    npx3 = npx1
    npy1 = int(math.ceil(yremainder*ny)+1)
    npy2 = int(math.ceil(s1*ny)+1)
    npy3 = npy1
    
    # coordinates of points in each part
    x1 = np.linspace(0,xremainder,npx1)
    x2 = np.linspace(x1[-1],x1[-1]+s2,npx2)
    x3 = np.linspace(x2[-1],1,npx3)
    x = np.concatenate((x1,x2[1:],x3[1:]))

    y1 = np.linspace(0,yremainder,npy1)
    y2 = np.linspace(y1[-1],y1[-1]+s1,npy2)
    y3 = np.linspace(y2[-1],1,npy3)
    y = np.concatenate((y1,y2[1:],y3[1:]))
    
    # new number of elements
    npx = len(x)
    npy = len(y)
    nx = npx - 1
    ny = npy - 1
    numNodes = npx*npy
    
    mesh = Mesh(nx,ny)

    # node coordinates
    if s3 == 0.0:
        # no shearing -> just repeat
        xcoords = np.tile(x,(1,npy))[0]
    else:
        # shearing -> add offset
        xcoords = np.zeros((numNodes,))
        for ii in range(npy):
            for jj in range(npx):
                idx = jj + ii*npx
                xcoords[idx] = x[jj] + y[ii] * np.tan(s3)
    ycoords = np.reshape(np.tile(y,(npx,1)),(len(y)*npx,),'F')
    nodes = np.array((xcoords, ycoords))
    for i, row in enumerate(np.transpose(nodes)):
        mesh.nodes[i] = row
    
    # boundary nodes
    if not sparse or s1 == 1 or s2 == 1:
        left = range(0, numNodes-1, npx)
        right = range(npx-1, numNodes, npx)
        bottom = range(npx)
        up = range(numNodes-npx, numNodes)
    else:
        left = range((npy1-1)*npx, (npy1+npy2-2)*npx+1, npx)
        right = range(npy1*npx-1, (npy1+npy2-1)*npx, npx)
        bottom = range(npx1-1, npx1+npx2-1)
        up = range(npx*(npy-1)+npx1-1, npx*(npy-1)+npx1+npx2-1)
    mesh.bc.append(("left", left))
    mesh.bc.append(("right", right))
    mesh.bc.append(("bottom", bottom))
    mesh.bc.append(("up", up))
    
    # elements
    num2DElements = nx*ny
    for yy in range(ny):
        for xx in range(nx):
            e = Element()
            e.type = QUAD4

            ll = npx*yy + xx # lower-left of current element
            e.nodes = ((ll, ll+1, ll+npx+1, ll+npx))
            
            if yy < len(y1)-1 or yy > len(y1)+len(y2)-3:
                if xx < len(x1)-1 or xx > len(x1)+len(x2)-3:
                    e.region = 'void'
                else:
                    e.region = 'mech'
            else:
                e.region = 'mech'
            
            if sparse and e.region == 'void':
                continue

            mesh.elements.append(e)

    return mesh, volume
