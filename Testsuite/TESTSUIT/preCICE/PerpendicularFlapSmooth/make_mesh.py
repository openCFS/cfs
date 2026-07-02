#!/usr/bin/env python3
"""Generate a tiny FSI+smooth mesh (openCFS HDF5) and the mock partner spec for the
PerpendicularFlapSmooth preCICE test - the full coupling combination, OpenFOAM replaced
by the mock_fluid partner:

  * Solid (mechanic FSI)  : a clamped flap block. Reads the fluid traction on its wetted
                            interface (mechNormalStressPrecice) and writes its displacement.
  * Fluid (prescribed Smooth): a block whose whole-volume nodal deformation (DisplacementVolume)
                            is written by the partner and injected into the SmoothPDE (no solve).
  * Coupling              : parallel-implicit with IQN-ILS (iterative) -> exercises the
                            checkpoint/re-iteration path and two PDEs per participant.

Geometry (small + plane-strain in z): a flap of width ws clamped at y=0, with a Fluid block
on its +x side; they share the interface plane x=ws.

  y=H  +--Solid--+----Fluid----+
       |    |    |   |    |     |
  y=0  +====+====+===+====+=====+   (==== = SolidFixed, clamped)
       x=0  x=ws            x=ws+wf
            ^ SolidSide_Interface (Solid +x face, wetted)

Outputs (next to this script):
  * PerpendicularFlapSmooth_mesh.h5  : openCFS HDF5 mesh (regions Solid, Fluid + surfaces)
  * mock_fluid.spec                  : partner meshes Fluid-Faces1 (write Stress),
                                       Fluid-Nodes1 (read Displacement),
                                       Fluid-Volume (write DisplacementVolume)
"""
import os
import numpy as np
import h5py

# --- geometry / discretisation (tiny on purpose) ---
ws, wf = 0.01, 0.02      # solid width, fluid width [m]
H, d   = 0.04, 0.01      # flap height, depth [m]
NY, NZ = 4, 1            # cells in y, z ; x: 1 solid cell + 2 fluid cells

xs = [0.0, ws, ws + wf / 2.0, ws + wf]            # nx = 4 (x=ws shared at index 1)
ys = [i * (H / NY) for i in range(NY + 1)]         # nyc = 5
zs = [i * (d / NZ) for i in range(NZ + 1)]         # nzc = 2
NX, NYC, NZC = len(xs), len(ys), len(zs)
NXS = 1                                             # solid occupies cell column cx = 0

ET_HEXA8, ET_QUAD4 = 11, 6
HERE = os.path.dirname(os.path.abspath(__file__))


def nid(ix, iy, iz):
    return iz * NX * NYC + iy * NX + ix + 1


def build():
    coords = []
    for iz in range(NZC):
        for iy in range(NYC):
            for ix in range(NX):
                coords.append((xs[ix], ys[iy], zs[iz]))
    coords = np.array(coords, dtype=np.float64)

    conn_rows, types, regions = [], [], {}
    solid_e, fluid_e, iface_e, fixed_e, fb_e = [], [], [], [], []

    def add_hex(ns, store):
        conn_rows.append(ns + [0] * (8 - len(ns)))
        types.append(ET_HEXA8)
        store.append(len(conn_rows))           # 1-based element id

    def add_quad(ns, store):
        conn_rows.append(ns + [0] * (8 - len(ns)))
        types.append(ET_QUAD4)
        store.append(len(conn_rows))

    # volume hexes (Solid: cx<NXS, Fluid otherwise)
    for cz in range(NZC - 1):
        for cy in range(NYC - 1):
            for cx in range(NX - 1):
                h = [nid(cx, cy, cz),     nid(cx + 1, cy, cz),
                     nid(cx + 1, cy + 1, cz), nid(cx, cy + 1, cz),
                     nid(cx, cy, cz + 1), nid(cx + 1, cy, cz + 1),
                     nid(cx + 1, cy + 1, cz + 1), nid(cx, cy + 1, cz + 1)]
                add_hex(h, solid_e if cx < NXS else fluid_e)

    # SolidSide_Interface: +x face of the solid cell (ix = NXS = 1)
    for cz in range(NZC - 1):
        for cy in range(NYC - 1):
            add_quad([nid(NXS, cy, cz), nid(NXS, cy + 1, cz),
                      nid(NXS, cy + 1, cz + 1), nid(NXS, cy, cz + 1)], iface_e)
    # SolidFixed: -y face (iy = 0) of the solid column
    for cz in range(NZC - 1):
        for cx in range(NXS):
            add_quad([nid(cx, 0, cz), nid(cx + 1, 0, cz),
                      nid(cx + 1, 0, cz + 1), nid(cx, 0, cz + 1)], fixed_e)
    # frontAndBack_solid: -z (iz=0) and +z (iz=NZC-1) faces of the solid column
    for cy in range(NYC - 1):
        for cx in range(NXS):
            add_quad([nid(cx, cy, 0), nid(cx + 1, cy, 0),
                      nid(cx + 1, cy + 1, 0), nid(cx, cy + 1, 0)], fb_e)
            add_quad([nid(cx, cy, NZC - 1), nid(cx + 1, cy, NZC - 1),
                      nid(cx + 1, cy + 1, NZC - 1), nid(cx, cy + 1, NZC - 1)], fb_e)

    conn = np.array(conn_rows, dtype=np.uint32)
    types = np.array(types, dtype=np.int32)
    n_elem = conn.shape[0]
    n_hex = len(solid_e) + len(fluid_e)
    n_quad = n_elem - n_hex

    def rnodes(elem_ids):
        s = set()
        for e in elem_ids:
            for n in conn[e - 1]:
                if n != 0:
                    s.add(int(n))
        return np.array(sorted(s), dtype=np.uint32)

    regions = {
        "Solid":               (3, np.array(solid_e, dtype=np.int32),  rnodes(solid_e)),
        "Fluid":               (3, np.array(fluid_e, dtype=np.int32),  rnodes(fluid_e)),
        "SolidSide_Interface": (2, np.array(iface_e, dtype=np.int32),  rnodes(iface_e)),
        "SolidFixed":          (2, np.array(fixed_e, dtype=np.int32),  rnodes(fixed_e)),
        "frontAndBack_solid":  (2, np.array(fb_e,    dtype=np.int32),  rnodes(fb_e)),
    }

    mesh_path = os.path.join(HERE, "PerpendicularFlapSmooth_mesh.h5")
    with h5py.File(mesh_path, "w") as f:
        fi = f.create_group("FileInfo")
        fi.create_dataset("Content", data=np.array([1, 2, 5], dtype=np.int32))
        sdt = h5py.string_dtype(encoding="utf-8")
        fi.create_dataset("Creator", data=np.array(["PerpendicularFlapSmooth make_mesh.py"], dtype=object), dtype=sdt)
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

    # --- mock partner spec ---
    iface_nodes = regions["SolidSide_Interface"][2]
    node_xyz = [tuple(coords[n - 1]) for n in iface_nodes]
    face_xyz = []
    for e in iface_e:
        ns = [n for n in conn[e - 1] if n != 0]
        c = np.mean([coords[n - 1] for n in ns], axis=0)
        face_xyz.append(tuple(c))
    fvol_nodes = regions["Fluid"][2]
    fvol_xyz = [tuple(coords[n - 1]) for n in fvol_nodes]

    spec_path = os.path.join(HERE, "mock_fluid.spec")
    with open(spec_path, "w") as s:
        s.write("# mock Fluid partner for PerpendicularFlapSmooth (generated by make_mesh.py)\n")
        s.write("participant Fluid\n")
        s.write("dimensions 3\n")
        s.write("rampTime 0.05\n")
        # surface traction [Pa] pushing the flap in +x, ramped (kept small -> small-strain deflection)
        s.write("mesh Fluid-Faces1 write Stress 1.0e3 0.0 0.0\n%d\n" % len(face_xyz))
        for x, y, z in face_xyz:
            s.write("%.10g %.10g %.10g\n" % (x, y, z))
        # reads the solid displacement back (ignored - analytic partner)
        s.write("mesh Fluid-Nodes1 read Displacement\n%d\n" % len(node_xyz))
        for x, y, z in node_xyz:
            s.write("%.10g %.10g %.10g\n" % (x, y, z))
        # whole-volume mesh deformation [m] injected into the SmoothPDE, ramped
        s.write("mesh Fluid-Volume write DisplacementVolume 1.0e-3 5.0e-4 0.0\n%d\n" % len(fvol_xyz))
        for x, y, z in fvol_xyz:
            s.write("%.10g %.10g %.10g\n" % (x, y, z))
        s.write("end\n")

    print("wrote %s: %d nodes, %d elems (%d hex + %d quad)" % (mesh_path, coords.shape[0], n_elem, n_hex, n_quad))
    print("  Solid hexes=%d Fluid hexes=%d | interface: %d nodes %d faces | Fluid-Volume: %d nodes"
          % (len(solid_e), len(fluid_e), len(node_xyz), len(face_xyz), len(fvol_xyz)))


if __name__ == "__main__":
    build()
