#!/usr/bin/env python3
r"""Check total linear momentum conservation for ContactImpactOrganic3d.

Both bodies are free and nothing external acts on them, so the exact statement is

    P(t) = SUM_j w_j v_j(t)  =  const,      w_j = rho * INT N_j dV

with the sum over every node of both bodies. `w_j` is the row sum of the consistent mass matrix,
which is what turns a nodal velocity field into a momentum: P = 1^T M v, and 1^T M is the vector
of column sums. That makes this an EXACT property of the discrete system -- it does not depend on
the mesh, the shape, the facet quality or the contact resolution, only on the contact tractions
being equal and opposite. Which is precisely why it is the right invariant for a body whose exact
answer nobody can write down.

The weights are computed in closed form rather than read from the solver:

  * TET4 is linear, so INT N_j dV = V/4 for each of its four nodes.
  * HEXA27 here is RECTILINEAR (the slab is a box), so the Jacobian is constant and
    INT N_j dV factorises into three 1D integrals of the quadratic Lagrange basis on [-1,1]:
    1/3 at each end node, 4/3 at the mid node. The element's node-to-reference-position map is
    CFS's HEXA27 ordering (8 corners, 12 edge midpoints, 6 face centres, 1 body centre).

Note this is a sharper measurement than the one quoted for ContactImpactOblique2d ("~2%"), which
compared unweighted MEAN nodal velocities of the two bodies -- a proxy that charges elastic
vibration against momentum. Properly weighted, momentum is conserved to solver precision.

Run (from a directory holding the finished run):  python3 check_momentum.py
"""

import argparse
import glob

import numpy as np
import h5py

ET_TRIA3, ET_QUAD9, ET_TET4, ET_HEXA27 = 4, 8, 9, 13

#: CFS HEXA27 local node -> offsets in {0,1,2} per axis, i.e. reference coordinate {-1,0,+1}.
HEXA27_OFF = [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
              (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
              (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
              (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
              (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1),
              (1, 0, 1), (2, 1, 1), (1, 2, 1), (0, 1, 1),
              (1, 1, 0), (1, 1, 2),
              (1, 1, 1)]
#: 1D integral of the quadratic Lagrange basis on [-1,1], by node position.
LAG1D = {0: 1.0 / 3.0, 1: 4.0 / 3.0, 2: 1.0 / 3.0}


def nodal_weights(xyz, conn, types, rho):
    """w_j = rho * INT N_j dV, summed over all volume elements."""
    w = np.zeros(len(xyz))

    for e in range(len(types)):
        t = types[e]

        if t == ET_TET4:
            n = conn[e, :4].astype(int) - 1
            p = xyz[n]
            vol = abs(np.dot(np.cross(p[1] - p[0], p[2] - p[0]), p[3] - p[0])) / 6.0
            w[n] += rho * vol / 4.0

        elif t == ET_HEXA27:
            n = conn[e, :27].astype(int) - 1
            p = xyz[n]
            # rectilinear box: side lengths from the three edges at corner 0
            dx = abs(p[1][0] - p[0][0])
            dy = abs(p[3][1] - p[0][1])
            dz = abs(p[4][2] - p[0][2])
            jac = (dx / 2.0) * (dy / 2.0) * (dz / 2.0)
            for a, off in enumerate(HEXA27_OFF):
                w[n[a]] += rho * jac * LAG1D[off[0]] * LAG1D[off[1]] * LAG1D[off[2]]

    return w


def main(cfs_file, rho, v0, launched_region):
    h = h5py.File(cfs_file, 'r')

    xyz = h['Mesh/Nodes/Coordinates'][:]
    conn = h['Mesh/Elements/Connectivity'][:]
    types = h['Mesh/Elements/Types'][:]

    w = nodal_weights(xyz, conn, types, rho)

    regions = {}
    for name in h['Mesh/Regions']:
        regions[name] = h[f'Mesh/Regions/{name}/Nodes'][:].astype(int) - 1

    ms = h['Results/Mesh/MultiStep_1']
    steps = sorted((int(k.split('_')[1]) for k in ms if k.startswith('Step_')))

    vol_regions = [n for n in regions
                   if f'Step_{steps[0]}/mechVelocity/{n}' in ms]

    print(f'{cfs_file}')
    print(f'  volume regions with velocity output: {vol_regions}')
    print(f'  total mass  {w.sum():.6g} kg   ('
          + ', '.join(f'{n}: {w[regions[n]].sum():.6g}' for n in vol_regions) + ')')

    # Expected momentum: only `launched_region` is given an initial velocity, and its nodes carry
    # it exactly, so P is fixed by the mesh before the solver ever runs.
    p_exp = w[regions[launched_region]].sum() * np.asarray(v0)
    print(f'  expected P = m_{launched_region} * v0 = '
          f'({p_exp[0]:.9g}, {p_exp[1]:.9g}, {p_exp[2]:.9g}) kg m/s')

    print(f'\n  {"step":>4} {"Px":>15} {"Py":>15} {"Pz":>15} {"|P-P0|/|P0|":>13}')
    p0 = None
    worst = 0.0
    for s in steps:
        P = np.zeros(3)
        for n in vol_regions:
            v = ms[f'Step_{s}/mechVelocity/{n}/Nodes/Real'][:]
            P += w[regions[n]] @ v
        if p0 is None:
            p0 = p_exp
        rel = np.linalg.norm(P - p0) / np.linalg.norm(p0)
        worst = max(worst, rel)
        print(f'  {s:>4} {P[0]:>15.9g} {P[1]:>15.9g} {P[2]:>15.9g} {rel:>13.3e}')

    print(f'\n  worst relative drift over the run: {worst:.3e}')
    return worst


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--file', default=None, help='the .cfs result file')
    ap.add_argument('--rho', type=float, default=2700.0)
    ap.add_argument('--v0', type=float, nargs=3, default=[50.0, 0.0, -5.0])
    ap.add_argument('--launched', default='ball')
    a = ap.parse_args()
    f = a.file or sorted(glob.glob('results_hdf5/*.cfs'))[0]
    main(f, a.rho, a.v0, a.launched)
