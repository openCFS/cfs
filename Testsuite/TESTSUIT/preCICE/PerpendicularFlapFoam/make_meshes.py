#!/usr/bin/env python3
"""Generate the two meshes of the PerpendicularFlapFoam test - the REAL OpenFOAM-coupled
FSI test (pimpleFoam + openfoam-adapter as "Fluid" vs. the cfs mechanic "Solid" flap):

  * PerpendicularFlapFoam_mesh.h5        : openCFS mesh - the "Solid" flap only
    (clamped base, three wetted interface surfaces, plane strain in z).
  * fluid-openfoam/system/blockMeshDict  : the OpenFOAM fluid mesh around the flap.

Both are derived from the SAME block spec below so the wetted interfaces coincide
geometrically (nearest-neighbour mapping of Stress/Displacement is unambiguous).

Geometry (tiny 2D channel, one cell in z / plane strain):

  y=YT +---------+--+---------+
       |    C    |D |    E    |     flap: x in [XF0,XF1], y in [Y0,YF]
  y=YF +---------+==+---------+     ==== flap top (SolidSide_Interface3)
       |    A    |##|    B    |     ## = Solid (NOT part of the fluid mesh)
  y=Y0 +---------+##+---------+
      x=XL     x=XF0/XF1     x=XR

Element budget (<400): fluid 16x12-2x8 = 176 hexes (OpenFOAM) + solid 2x8 = 16 hexes.

Outputs are committed; re-run this script only when changing the discretisation
(then also regenerate the .h5ref).
"""
import os
import numpy as np
import h5py

# --- shared block spec (keep h5 and blockMeshDict consistent!) ---
XL, XF0, XF1, XR = -1.5, -0.05, 0.05, 1.5   # channel / flap x-coordinates
Y0, YF, YT = 0.0, 1.0, 2.0                   # bottom / flap top / channel top
Z0, Z1 = 0.0, 0.1                            # one cell deep (2D, plane strain)
NXL, NXF, NXR = 7, 2, 7                      # fluid x-cells: left / flap column / right
NYL, NYU = 8, 4                              # fluid y-cells: below / above flap top
NSX, NSY = 2, 8                              # solid flap cells in x / y

ET_HEXA8, ET_QUAD4 = 11, 6
HERE = os.path.dirname(os.path.abspath(__file__))


def lin(a, b, n):
    return [a + (b - a) * i / n for i in range(n + 1)]


def build_h5():
    # solid flap lattice
    sxs, sys, szs = lin(XF0, XF1, NSX), lin(Y0, YF, NSY), [Z0, Z1]
    SNX, SNY, SNZ = len(sxs), len(sys), len(szs)

    def snid(ix, iy, iz):
        return iz * SNX * SNY + iy * SNX + ix + 1

    coords = np.array([(x, y, z) for z in szs for y in sys for x in sxs], dtype=np.float64)

    conn_rows, types = [], []
    solid_e, if1_e, if2_e, if3_e, fixed_e, fb_e = [], [], [], [], [], []

    def add(ns, et, store):
        conn_rows.append(ns + [0] * (8 - len(ns)))
        types.append(et)
        store.append(len(conn_rows))  # 1-based element id

    # solid hexes
    for cz in range(SNZ - 1):
        for cy in range(SNY - 1):
            for cx in range(SNX - 1):
                add([snid(cx, cy, cz), snid(cx + 1, cy, cz),
                     snid(cx + 1, cy + 1, cz), snid(cx, cy + 1, cz),
                     snid(cx, cy, cz + 1), snid(cx + 1, cy, cz + 1),
                     snid(cx + 1, cy + 1, cz + 1), snid(cx, cy + 1, cz + 1)],
                    ET_HEXA8, solid_e)

    # surface regions (loop over CELLS: NSX/NSY, not the node counts!)
    for cy in range(NSY):          # Interface1: -x side ; Interface2: +x side
        add([snid(0, cy, 0), snid(0, cy + 1, 0),
             snid(0, cy + 1, 1), snid(0, cy, 1)], ET_QUAD4, if1_e)
        add([snid(NSX, cy, 0), snid(NSX, cy + 1, 0),
             snid(NSX, cy + 1, 1), snid(NSX, cy, 1)], ET_QUAD4, if2_e)
    for cx in range(NSX):          # Interface3: flap top (y=YF, node row iy=NSY)
        add([snid(cx, NSY, 0), snid(cx + 1, NSY, 0),
             snid(cx + 1, NSY, 1), snid(cx, NSY, 1)], ET_QUAD4, if3_e)
    for cx in range(NSX):          # SolidFixed: clamped base (y=Y0)
        add([snid(cx, 0, 0), snid(cx + 1, 0, 0),
             snid(cx + 1, 0, 1), snid(cx, 0, 1)], ET_QUAD4, fixed_e)
    for cy in range(NSY):          # frontAndBack_solid: z=Z0 and z=Z1 (plane strain)
        for cx in range(NSX):
            add([snid(cx, cy, 0), snid(cx + 1, cy, 0),
                 snid(cx + 1, cy + 1, 0), snid(cx, cy + 1, 0)], ET_QUAD4, fb_e)
            add([snid(cx, cy, 1), snid(cx + 1, cy, 1),
                 snid(cx + 1, cy + 1, 1), snid(cx, cy + 1, 1)], ET_QUAD4, fb_e)

    conn = np.array(conn_rows, dtype=np.uint32)
    types = np.array(types, dtype=np.int32)
    assert conn.max() == len(coords), \
        "connectivity references node %d but grid has %d nodes" % (conn.max(), len(coords))
    n_elem = conn.shape[0]
    n_hex = len(solid_e)
    n_quad = n_elem - n_hex

    def rnodes(elem_ids):
        s = set()
        for e in elem_ids:
            for n in conn[e - 1]:
                if n != 0:
                    s.add(int(n))
        return np.array(sorted(s), dtype=np.uint32)

    regions = {
        "Solid":                (3, np.array(solid_e, dtype=np.int32), rnodes(solid_e)),
        "SolidSide_Interface1": (2, np.array(if1_e, dtype=np.int32),   rnodes(if1_e)),
        "SolidSide_Interface2": (2, np.array(if2_e, dtype=np.int32),   rnodes(if2_e)),
        "SolidSide_Interface3": (2, np.array(if3_e, dtype=np.int32),   rnodes(if3_e)),
        "SolidFixed":           (2, np.array(fixed_e, dtype=np.int32), rnodes(fixed_e)),
        "frontAndBack_solid":   (2, np.array(fb_e, dtype=np.int32),    rnodes(fb_e)),
    }

    mesh_path = os.path.join(HERE, "PerpendicularFlapFoam_mesh.h5")
    with h5py.File(mesh_path, "w") as f:
        fi = f.create_group("FileInfo")
        fi.create_dataset("Content", data=np.array([1, 2, 5], dtype=np.int32))
        sdt = h5py.string_dtype(encoding="utf-8")
        fi.create_dataset("Creator", data=np.array(["PerpendicularFlapFoam make_meshes.py"], dtype=object), dtype=sdt)
        fi.create_dataset("Date", data=np.array(["generated"], dtype=object), dtype=sdt)
        fi.create_dataset("Version", data=np.array(["0.9"], dtype=object), dtype=sdt)

        mesh = f.create_group("Mesh")
        mesh.attrs["Dimension"] = np.int32(3)
        nodes = mesh.create_group("Nodes")
        nodes.attrs["NumNodes"] = np.int32(coords.shape[0])
        nodes.create_dataset("Coordinates", data=coords)

        elems = mesh.create_group("Elements")
        elems.attrs["NumElems"] = np.int32(n_elem)
        elems.attrs["Num1DElems"] = np.int32(0)
        elems.attrs["Num2DElems"] = np.int32(n_quad)
        elems.attrs["Num3DElems"] = np.int32(n_hex)
        elems.attrs["Num_HEXA8"] = np.int32(n_hex)
        elems.attrs["Num_QUAD4"] = np.int32(n_quad)
        elems.attrs["QuadraticElems"] = np.int32(0)
        elems.create_dataset("Connectivity", data=conn)
        elems.create_dataset("Types", data=types)

        mesh.create_group("Groups")
        regroup = mesh.create_group("Regions")
        for name, (dim, el, nd) in regions.items():
            g = regroup.create_group(name)
            g.attrs["Dimension"] = np.int32(dim)
            g.create_dataset("Elements", data=el)
            g.create_dataset("Nodes", data=nd)

    print("wrote %s: %d nodes, %d elems (%d hex solid, %d quad)"
          % (mesh_path, coords.shape[0], n_elem, n_hex, n_quad))


def build_blockmeshdict():
    # corner vertex grid: 4 x-lines x 3 y-lines x 2 z-lines
    vx, vy, vz = [XL, XF0, XF1, XR], [Y0, YF, YT], [Z0, Z1]

    def v(ix, iy, iz):
        return iz * 12 + iy * 4 + ix

    verts = [(x, y, z) for z in vz for y in vy for x in vx]

    def hexblock(ix0, iy0, ix1, iy1, nx, ny):
        c = [v(ix0, iy0, 0), v(ix1, iy0, 0), v(ix1, iy1, 0), v(ix0, iy1, 0),
             v(ix0, iy0, 1), v(ix1, iy0, 1), v(ix1, iy1, 1), v(ix0, iy1, 1)]
        return "    hex (%s) (%d %d 1) simpleGrading (1 1 1)" % (" ".join(map(str, c)), nx, ny)

    blocks = [
        hexblock(0, 0, 1, 1, NXL, NYL),  # A: left of flap
        hexblock(2, 0, 3, 1, NXR, NYL),  # B: right of flap
        hexblock(0, 1, 1, 2, NXL, NYU),  # C: upper-left
        hexblock(1, 1, 2, 2, NXF, NYU),  # D: above flap
        hexblock(2, 1, 3, 2, NXR, NYU),  # E: upper-right
    ]

    def face(a, b, c, d):
        return "            (%d %d %d %d)" % (a, b, c, d)

    patches = """    s_inlet
    {
        type patch;
        faces
        (
%s
%s
        );
    }
    s_outlet
    {
        type patch;
        faces
        (
%s
%s
        );
    }
    s_walls
    {
        type wall;
        faces
        (
%s
%s
%s
%s
%s
        );
    }
    s_FluidSide_Interface1
    {
        type wall;
        faces
        (
%s
        );
    }
    s_FluidSide_Interface2
    {
        type wall;
        faces
        (
%s
        );
    }
    s_FluidSide_Interface3
    {
        type wall;
        faces
        (
%s
        );
    }
    s_frontAndBack_fluid
    {
        type empty;
        faces
        (
%s
%s
%s
%s
%s
%s
%s
%s
%s
%s
        );
    }""" % (
        face(v(0, 0, 0), v(0, 0, 1), v(0, 1, 1), v(0, 1, 0)),   # inlet A
        face(v(0, 1, 0), v(0, 1, 1), v(0, 2, 1), v(0, 2, 0)),   # inlet C
        face(v(3, 0, 0), v(3, 1, 0), v(3, 1, 1), v(3, 0, 1)),   # outlet B
        face(v(3, 1, 0), v(3, 2, 0), v(3, 2, 1), v(3, 1, 1)),   # outlet E
        face(v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)),   # bottom A
        face(v(2, 0, 0), v(3, 0, 0), v(3, 0, 1), v(2, 0, 1)),   # bottom B
        face(v(0, 2, 0), v(0, 2, 1), v(1, 2, 1), v(1, 2, 0)),   # top C
        face(v(1, 2, 0), v(1, 2, 1), v(2, 2, 1), v(2, 2, 0)),   # top D
        face(v(2, 2, 0), v(2, 2, 1), v(3, 2, 1), v(3, 2, 0)),   # top E
        face(v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)),   # flap -x side (A east)
        face(v(2, 0, 0), v(2, 0, 1), v(2, 1, 1), v(2, 1, 0)),   # flap +x side (B west)
        face(v(1, 1, 0), v(1, 1, 1), v(2, 1, 1), v(2, 1, 0)),   # flap top (D south)
        face(v(0, 0, 0), v(0, 1, 0), v(1, 1, 0), v(1, 0, 0)),   # z0 A
        face(v(2, 0, 0), v(2, 1, 0), v(3, 1, 0), v(3, 0, 0)),   # z0 B
        face(v(0, 1, 0), v(0, 2, 0), v(1, 2, 0), v(1, 1, 0)),   # z0 C
        face(v(1, 1, 0), v(1, 2, 0), v(2, 2, 0), v(2, 1, 0)),   # z0 D
        face(v(2, 1, 0), v(2, 2, 0), v(3, 2, 0), v(3, 1, 0)),   # z0 E
        face(v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)),   # z1 A
        face(v(2, 0, 1), v(3, 0, 1), v(3, 1, 1), v(2, 1, 1)),   # z1 B
        face(v(0, 1, 1), v(1, 1, 1), v(1, 2, 1), v(0, 2, 1)),   # z1 C
        face(v(1, 1, 1), v(2, 1, 1), v(2, 2, 1), v(1, 2, 1)),   # z1 D
        face(v(2, 1, 1), v(3, 1, 1), v(3, 2, 1), v(2, 2, 1)),   # z1 E
    )

    out = os.path.join(HERE, "fluid-openfoam", "system", "blockMeshDict")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        f.write("""/* generated by make_meshes.py - edit the block spec there, not here */
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

scale 1;

vertices
(
%s
);

blocks
(
%s
);

boundary
(
%s
);
""" % ("\n".join("    (%g %g %g)" % c for c in verts),
       "\n".join(blocks),
       patches))
    print("wrote %s (fluid cells: %d)" % (out, NXL * NYL + NXR * NYL + (NXL + NXF + NXR) * NYU))


if __name__ == "__main__":
    build_h5()
    build_blockmeshdict()
