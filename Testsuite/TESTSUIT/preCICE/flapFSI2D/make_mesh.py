#!/usr/bin/env python3
"""Generate a tiny 2D cantilever-beam mesh (openCFS HDF5) for the flapFSI2D
preCICE coupling test, plus the interface geometry the mock "fluid" partner
needs to provide matching preCICE meshes.

Physics / geometry (kept deliberately minimal and traceable):

    y=H  +---+---+---+---+---+      <- "interface": fluid pushes here (top edge)
         |   |   |   |   |   |         (nx LINE2 surface elements)
    y=0  +---+---+---+---+---+
        x=0                 x=L
        ^
        |  "fixed": clamped left edge (ny LINE2 surface elements)

  * Volume region "beam"      : nx*ny QUAD4 elements (mechanic PDE lives here)
  * Surface region "interface": top edge, where the mock writes the fluid
                                traction (mechNormalStress) and reads back the
                                solid displacement (mechDisplacement)
  * Surface region "fixed"    : left edge, clamped (Dirichlet)

The mesh is intentionally coarse (a handful of elements) so the reference file
stays well below 1 MB and the simulation runs in well under a second.

Outputs (written next to this script):
  * flapFSI2D_mesh.h5  : openCFS HDF5 mesh (input for the solver)
  * mock_fluid.spec    : interface node + face-centre coordinates for the mock
"""
import os
import numpy as np
import h5py

# ---------------------------------------------------------------------------
# Geometry / discretisation
# ---------------------------------------------------------------------------
L = 0.01      # beam length            [m]
H = 0.001     # beam height            [m]
NX = 5        # elements along length
NY = 2        # elements along height

# openCFS element type codes (source/Domain/ElemMapping/Elem.hh)
ET_LINE2 = 2
ET_QUAD4 = 6

HERE = os.path.dirname(os.path.abspath(__file__))


def node_id(i, j):
    """1-based node id of grid point (i in 0..NX, j in 0..NY)."""
    return j * (NX + 1) + i + 1


def build():
    # --- nodes (coordinates are always stored 3D, z = 0) ---
    coords = []
    for j in range(NY + 1):
        for i in range(NX + 1):
            coords.append((i * (L / NX), j * (H / NY), 0.0))
    coords = np.array(coords, dtype=np.float64)

    # --- elements: volume QUAD4 first, then surface LINE2 (interface, fixed) ---
    quads, interface_lines, fixed_lines = [], [], []
    for j in range(NY):
        for i in range(NX):
            quads.append((node_id(i, j), node_id(i + 1, j),
                          node_id(i + 1, j + 1), node_id(i, j + 1)))
    for i in range(NX):                       # top edge (j = NY)
        interface_lines.append((node_id(i, NY), node_id(i + 1, NY)))
    for j in range(NY):                       # left edge (i = 0)
        fixed_lines.append((node_id(0, j), node_id(0, j + 1)))

    n_quad = len(quads)
    n_if = len(interface_lines)
    n_fix = len(fixed_lines)
    n_elem = n_quad + n_if + n_fix

    # connectivity padded to 4 columns (LINE2 rows -> [n1, n2, 0, 0])
    conn = np.zeros((n_elem, 4), dtype=np.uint32)
    types = np.zeros((n_elem,), dtype=np.int32)
    eid = 0
    beam_elems, interface_elems, fixed_elems = [], [], []
    for q in quads:
        conn[eid, :] = q
        types[eid] = ET_QUAD4
        beam_elems.append(eid + 1)            # 1-based element id
        eid += 1
    for l in interface_lines:
        conn[eid, :2] = l
        types[eid] = ET_LINE2
        interface_elems.append(eid + 1)
        eid += 1
    for l in fixed_lines:
        conn[eid, :2] = l
        types[eid] = ET_LINE2
        fixed_elems.append(eid + 1)
        eid += 1

    def region_nodes(elem_ids):
        s = set()
        for e in elem_ids:
            for n in conn[e - 1]:
                if n != 0:
                    s.add(int(n))
        return np.array(sorted(s), dtype=np.uint32)

    regions = {
        "beam":      (2, np.array(beam_elems, dtype=np.int32),      region_nodes(beam_elems)),
        "interface": (1, np.array(interface_elems, dtype=np.int32), region_nodes(interface_elems)),
        "fixed":     (1, np.array(fixed_elems, dtype=np.int32),     region_nodes(fixed_elems)),
    }

    # --- write openCFS HDF5 mesh ---
    mesh_path = os.path.join(HERE, "flapFSI2D_mesh.h5")
    with h5py.File(mesh_path, "w") as f:
        fi = f.create_group("FileInfo")
        fi.create_dataset("Content", data=np.array([1, 2, 5], dtype=np.int32))
        sdt = h5py.string_dtype(encoding="utf-8")
        fi.create_dataset("Creator", data=np.array(["flapFSI2D make_mesh.py"], dtype=object), dtype=sdt)
        fi.create_dataset("Date", data=np.array(["generated"], dtype=object), dtype=sdt)
        fi.create_dataset("Version", data=np.array(["0.9"], dtype=object), dtype=sdt)

        mesh = f.create_group("Mesh")
        mesh.attrs["Dimension"] = np.int32(2)
        nodes = mesh.create_group("Nodes")
        nodes.attrs["NumNodes"] = np.int32(coords.shape[0])
        nodes.create_dataset("Coordinates", data=coords)

        elems = mesh.create_group("Elements")
        elems.attrs["NumElems"] = np.int32(n_elem)
        elems.attrs["Num1DElems"] = np.int32(n_if + n_fix)
        elems.attrs["Num2DElems"] = np.int32(n_quad)
        elems.attrs["Num3DElems"] = np.int32(0)
        elems.attrs["Num_QUAD4"] = np.int32(n_quad)
        elems.attrs["Num_LINE2"] = np.int32(n_if + n_fix)
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

    # --- interface geometry for the mock (Fluid) partner ---
    # Fluid-Nodes  : the interface nodes (top edge)            -> read Displacement
    # Fluid-Faces  : centres of the interface LINE2 elements    -> write Stress
    if_nodes = regions["interface"][2]
    node_xy = [(coords[n - 1, 0], coords[n - 1, 1]) for n in if_nodes]
    face_xy = []
    for e in interface_elems:
        n1, n2 = conn[e - 1, 0], conn[e - 1, 1]
        face_xy.append((0.5 * (coords[n1 - 1, 0] + coords[n2 - 1, 0]),
                        0.5 * (coords[n1 - 1, 1] + coords[n2 - 1, 1])))

    spec_path = os.path.join(HERE, "mock_fluid.spec")
    with open(spec_path, "w") as s:
        s.write("# mock fluid partner spec for flapFSI2D (generated by make_mesh.py)\n")
        s.write("participant Fluid\n")
        s.write("dimensions 2\n")
        # downward fluid traction [Pa], ramped linearly over rampTime [s].
        # ~2e4 Pa on this short steel beam gives a clean ~1e-6 m tip deflection
        # (small-strain / linear regime), growing linearly over the 5 windows.
        s.write("loadVector 0.0 -2.0e4\n")
        s.write("rampTime 5.0e-3\n")
        s.write(f"mesh Fluid-Faces write Stress\n{len(face_xy)}\n")
        for x, y in face_xy:
            s.write(f"{x:.10g} {y:.10g}\n")
        s.write(f"mesh Fluid-Nodes read Displacement\n{len(node_xy)}\n")
        for x, y in node_xy:
            s.write(f"{x:.10g} {y:.10g}\n")
        s.write("end\n")

    print(f"wrote {mesh_path}  ({coords.shape[0]} nodes, {n_elem} elems: "
          f"{n_quad} QUAD4 + {n_if + n_fix} LINE2)")
    print(f"wrote {spec_path}  ({len(node_xy)} interface nodes, {len(face_xy)} faces)")


if __name__ == "__main__":
    build()
