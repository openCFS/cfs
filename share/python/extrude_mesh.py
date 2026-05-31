#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""
extrude_mesh.py
=======================

Extrude a 2D openCFS HDF5 mesh (.cfs) into a 3D mesh by stacking copies
of the 2D nodes at user-supplied z-coordinates and connecting consecutive
layers with prismatic 3D elements. (It also supports extrude 1D line to 2D surface)

Element-type extrusion (linear order):
    LINE2  ->  QUAD4
    TRIA3  ->  WEDGE6
    QUAD4  ->  HEXA8

Element-type extrusion (quadratic order):
    LINE3  ->  QUAD8
    TRIA6  ->  WEDGE15
    QUAD8  ->  HEXA20

The mode is auto-detected: if the 2D mesh contains any quadratic surface
element (LINE3, TRIA6, QUAD8), the extruder switches to quadratic mode.
In quadratic mode, an extra "mid-layer" of nodes is inserted between
every pair of user z-layers — these supply the vertical mid-edge nodes
that HEXA20 / WEDGE15 / QUAD8 elements need.

The HEX8/WEDGE6/QUAD4 connectivity follows the VTK ordering
(bottom face first, then top face) which is what the openCFS CFSReader
ParaView plugin expects (see vtkCFSReader.cxx::AddElements / GetCellIdType).
The HEXA20/WEDGE15/QUAD8 connectivity follows the analogous VTK
quadratic ordering: corners first (bot then top), then bot mid-edges,
top mid-edges, and finally the vertical mid-edges.

Output schema is identical to a native openCFS mesh file, so the resulting
file is read directly by the CFSReader plugin.

Node-ID preservation
--------------------
The original 2D mesh sits at z=0 and the script is careful to keep its node
ids unchanged: 2D node id ``i`` becomes 3D node id ``i`` (same x, y, z=0).
This is true regardless of whether the extrusion is in +z, -z, or both.
Other layers occupy ids beyond ``n_nodes_2D``.  This makes it straightforward
to look up 3D nodes from 2D node ids when projecting 2.5D results.

z-list conventions
------------------
The original 2D mesh sits at z=0, so z=0 is always part of the extruded
layer set: the user only specifies the *additional* z-layers.  The list
is sanitized as follows:

  - empty list / no z given / only [0.0]   ->  no-op, warning, no output written
  - 0.0 missing                            ->  inserted automatically
  - duplicates / FP-noise around 0         ->  collapsed
  - any order                              ->  sorted ascending
  - negative values                        ->  fine, extrusion goes in -z

Examples
    [0.5, 1.0]            ->  layers [0, 0.5, 1.0]   (2 element layers, +z)
    [-1.0]                ->  layers [-1, 0]         (1 element layer, -z)
    [-1, -0.5, 0.5, 1.0]  ->  layers [-1,-0.5,0,0.5,1] (4 layers, both directions)

Usage
-----
    $ python3 extrude_mesh.py Mesh2DLinElem.cfs Mesh3DLinElem.cfs -1.0 -0.5 0.5 1.0
    $ python3 extrude_mesh.py Mesh2DQuadElem.cfs Mesh3DQuadElem.cfs -1.0 -0.5 0.5 1.0

or, programmatically:

    from extrude_mesh import extrude_mesh
    extrude_mesh("mesh2D.cfs", "mesh3D.cfs", [0.5, 1.0, 1.5, 2.0])

The script currently extrudes mesh geometry only.  Nodal/element results
under /Results are not touched yet (a warning is printed if they exist).

Author
------
Likun Luo (likun.luo@tuwien.ac.at)

Created
-------
2026-05-12
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

# -----------------------------------------------------------------------------
# openCFS element type IDs (mirrors hdf5Common.h::ElemType)
# -----------------------------------------------------------------------------
ET_UNDEF, ET_POINT = 0, 1
ET_LINE2, ET_LINE3 = 2, 3
ET_TRIA3, ET_TRIA6 = 4, 5
ET_QUAD4, ET_QUAD8, ET_QUAD9 = 6, 7, 8
ET_TET4, ET_TET10 = 9, 10
ET_HEXA8, ET_HEXA20, ET_HEXA27 = 11, 12, 13
ET_PYRA5, ET_PYRA13, ET_PYRA14 = 14, 15, 19
ET_WEDGE6, ET_WEDGE15, ET_WEDGE18 = 16, 17, 18

NUM_NODES = {
    ET_UNDEF: 0, ET_POINT: 1,
    ET_LINE2: 2, ET_LINE3: 3,
    ET_TRIA3: 3, ET_TRIA6: 6,
    ET_QUAD4: 4, ET_QUAD8: 8, ET_QUAD9: 9,
    ET_TET4: 4, ET_TET10: 10,
    ET_HEXA8: 8, ET_HEXA20: 20, ET_HEXA27: 27,
    ET_PYRA5: 5, ET_PYRA13: 13, ET_PYRA14: 14,
    ET_WEDGE6: 6, ET_WEDGE15: 15, ET_WEDGE18: 18,
}

# Spatial dimension of each element type (used for Num1D/2D/3DElems counters).
DIM_OF_TYPE = {
    ET_UNDEF: 0, ET_POINT: 0,
    ET_LINE2: 1, ET_LINE3: 1,
    ET_TRIA3: 2, ET_TRIA6: 2,
    ET_QUAD4: 2, ET_QUAD8: 2, ET_QUAD9: 2,
    ET_TET4: 3, ET_TET10: 3,
    ET_HEXA8: 3, ET_HEXA20: 3, ET_HEXA27: 3,
    ET_PYRA5: 3, ET_PYRA13: 3, ET_PYRA14: 3,
    ET_WEDGE6: 3, ET_WEDGE15: 3, ET_WEDGE18: 3,
}

# Linear extrusion mapping. (Applies when no quadratic 2D type is present.)
EXTRUSION_MAP_LINEAR = {
    ET_LINE2: ET_QUAD4,
    ET_TRIA3: ET_WEDGE6,
    ET_QUAD4: ET_HEXA8,
}

# Quadratic extrusion mapping. Activated when the 2D mesh contains any of
# {LINE3, TRIA6, QUAD8}.
EXTRUSION_MAP_QUADRATIC = {
    ET_LINE3: ET_QUAD8,
    ET_TRIA6: ET_WEDGE15,
    ET_QUAD8: ET_HEXA20,
}

# 2D types that are quadratic (presence of any of these triggers quadratic mode).
QUADRATIC_2D_TYPES = {ET_LINE3, ET_TRIA6, ET_QUAD8}

# Number of *corner* nodes for each 2D element type (the rest are mid-edges).
# Used both to identify which 2D nodes become "vertical mid-edge" suppliers
# in quadratic mode and to build the corner-only Nodes list at mid-layers.
NUM_CORNERS_2D = {
    ET_LINE2: 2, ET_LINE3: 2,
    ET_TRIA3: 3, ET_TRIA6: 3,
    ET_QUAD4: 4, ET_QUAD8: 4, ET_QUAD9: 4,
}

# Full list of Num_<TYPE> attributes that openCFS writes on /Mesh/Elements.
# Order does not matter (HDF5 attributes are unordered) but every name must
# be present so the file looks identical to one produced by openCFS itself.
ALL_TYPE_ATTRS = {
    ET_UNDEF:    "Num_UNDEF",
    ET_POINT:    "Num_POINT",
    ET_LINE2:    "Num_LINE2",
    ET_LINE3:    "Num_LINE3",
    ET_TRIA3:    "Num_TRIA3",
    ET_TRIA6:    "Num_TRIA6",
    ET_QUAD4:    "Num_QUAD4",
    ET_QUAD8:    "Num_QUAD8",
    ET_QUAD9:    "Num_QUAD9",
    ET_TET4:     "Num_TET4",
    ET_TET10:    "Num_TET10",
    ET_HEXA8:    "Num_HEXA8",
    ET_HEXA20:   "Num_HEXA20",
    ET_HEXA27:   "Num_HEXA27",
    ET_PYRA5:    "Num_PYRA5",
    ET_PYRA13:   "Num_PYRA13",
    ET_PYRA14:   "Num_PYRA14",
    ET_WEDGE6:   "Num_WEDGE6",
    ET_WEDGE15:  "Num_WEDGE15",
    ET_WEDGE18:  "Num_WEDGE18",
}

# Attribute names openCFS writes for which there is no ET_* enum.  openCFS
# writes the counter with value 0 anyway, so we mirror that to stay
# byte-faithful with native openCFS output.
NO_ENUM_TYPE_ATTRS = ("Num_POLYGON", "Num_POLYHEDRON")

# -----------------------------------------------------------------------------
# Connectivity rules (bottom face / top face split, matching VTK ordering)
# -----------------------------------------------------------------------------

def _build_3d_connectivity_row(t3: int, bot: np.ndarray,
                               top: np.ndarray) -> np.ndarray:
    """Return the connectivity row for one extruded 3D element (1-based node IDs).

    `bot` and `top` are the 1-based 3D node ids that come from the bottom
    and top z-layer respectively, in the *same order* as the original 2D
    element connectivity.
    """
    row = np.empty(NUM_NODES[t3], dtype=np.int32)
    if t3 == ET_HEXA8:
        # VTK_HEXAHEDRON: bottom 0..3, top 4..7 (top in same order as bottom)
        row[0:4] = bot
        row[4:8] = top
    elif t3 == ET_WEDGE6:
        # VTK_WEDGE: bottom 0..2, top 3..5
        row[0:3] = bot
        row[3:6] = top
    elif t3 == ET_QUAD4:
        # Extruded LINE2 -> QUAD4. The 2D edge has nodes (n0,n1).
        # CCW from outside: n0_bot, n1_bot, n1_top, n0_top
        row[0] = bot[0]
        row[1] = bot[1]
        row[2] = top[1]
        row[3] = top[0]
    else:
        raise NotImplementedError(f"3D element type {t3} not implemented")
    return row

def _build_3d_connectivity_row_quad(t3: int,
                                    bot: np.ndarray,
                                    top: np.ndarray,
                                    mid: np.ndarray) -> np.ndarray:
    """Return the connectivity row for one extruded quadratic 3D element (1-based node IDs).

    Parameters
    ----------
    t3   : 3D element type ID (HEXA20, WEDGE15, QUAD8).
    bot  : 1-based 3D node IDs at the bottom *full* layer, in the same
           order as the original 2D element connectivity (corners first,
           then 2D mid-edges).
    top  : same, for the top *full* layer.
    mid  : 1-based 3D node IDs at the *mid* layer between bot and top,
           one per 2D *corner* (so len(mid) == NUM_CORNERS_2D[t2]).
    Returns
    -------
    np.ndarray of length ``NUM_NODES[t3]``, dtype int32, 1-based node IDs.

    The orderings below match the VTK quadratic conventions
    (VTK_QUADRATIC_HEXAHEDRON, VTK_QUADRATIC_WEDGE, VTK_QUADRATIC_QUAD)
    which the openCFS CFSReader plugin consumes directly without
    remapping (see vtkCFSReader.cxx::GetCellIdType).
    """
    row = np.empty(NUM_NODES[t3], dtype=np.int32)
    if t3 == ET_HEXA20:
        # 2D QUAD8 conn order: [c0 c1 c2 c3 m01 m12 m23 m30].
        # HEXA20 ordering:
        #   0..3   bot corners
        #   4..7   top corners
        #   8..11  bot mid-edges (edges 0-1, 1-2, 2-3, 3-0)
        #   12..15 top mid-edges (edges 4-5, 5-6, 6-7, 7-4)
        #   16..19 vertical mid-edges (edges 0-4, 1-5, 2-6, 3-7)
        row[0:4]   = bot[0:4]
        row[4:8]   = top[0:4]
        row[8:12]  = bot[4:8]
        row[12:16] = top[4:8]
        row[16:20] = mid[0:4]
    elif t3 == ET_WEDGE15:
        # 2D TRIA6 conn order: [c0 c1 c2 m01 m12 m20].
        # WEDGE15 ordering:
        #   0..2   bot triangle corners
        #   3..5   top triangle corners
        #   6..8   bot mid-edges (edges 0-1, 1-2, 2-0)
        #   9..11  top mid-edges (edges 3-4, 4-5, 5-3)
        #   12..14 vertical mid-edges (edges 0-3, 1-4, 2-5)
        row[0:3]   = bot[0:3]
        row[3:6]   = top[0:3]
        row[6:9]   = bot[3:6]
        row[9:12]  = top[3:6]
        row[12:15] = mid[0:3]
    elif t3 == ET_QUAD8:
        # 2D LINE3 conn order: [c0 c1 m01].  Extruded to a vertical QUAD8.
        # CCW corners (from outside) match the linear LINE2->QUAD4 pattern:
        #   bot[0], bot[1], top[1], top[0]
        # Then VTK_QUADRATIC_QUAD mid-edges 4..7:
        #   4: edge 0-1 (bot)        = bot mid (bot[2])
        #   5: edge 1-2 (right vert) = mid[1]
        #   6: edge 2-3 (top)        = top mid (top[2])
        #   7: edge 3-0 (left vert)  = mid[0]
        row[0] = bot[0]
        row[1] = bot[1]
        row[2] = top[1]
        row[3] = top[0]
        row[4] = bot[2]
        row[5] = mid[1]
        row[6] = top[2]
        row[7] = mid[0]
    else:
        raise NotImplementedError(
            f"3D quadratic element type {t3} not implemented"
        )
    return row

# -----------------------------------------------------------------------------
# Helpers for writing openCFS-style attributes / groups
# -----------------------------------------------------------------------------

def _set_uint(obj, name: str, value: int) -> None:
    obj.attrs.create(name, np.uint32(value))

def _set_int(obj, name: str, value: int) -> None:
    obj.attrs.create(name, np.int32(value))

def _normalize_z_list(z_in, atol: float = 1e-12) -> np.ndarray:
    """Sanitize the user-supplied z list.

    Rules (z=0 is the original 2D plane and is therefore always included):
        - None / empty list   -> [0.0]                  (caller will treat as no-op)
        - any near-zero value -> snapped to exactly 0.0
        - 0.0 always inserted, then duplicates collapsed (with `atol`)
        - sorted ascending; negative values are fine
    """
    if z_in is None:
        z = np.array([], dtype=np.float64)
    else:
        z = np.asarray(z_in, dtype=np.float64).ravel()
    if z.size > 0:
        z = np.where(np.isclose(z, 0.0, atol=atol), 0.0, z)
    z = np.concatenate([z, [0.0]])
    z = np.sort(z)
    if z.size > 1:
        keep = np.concatenate([[True], ~np.isclose(np.diff(z), 0.0, atol=atol)])
        z = z[keep]
    return z

# -----------------------------------------------------------------------------
# Quadratic-mode helpers (public — also used by the 2.5D postprocessor)
# -----------------------------------------------------------------------------

def is_quadratic_2d_mesh(types_array) -> bool:
    """True if the 2D mesh contains any quadratic surface element type."""
    s = {int(t) for t in np.unique(np.asarray(types_array))}
    return bool(s & QUADRATIC_2D_TYPES)


def compute_corner_mask_2d(types_array, conn_array,
                           n_nodes_2d: int) -> np.ndarray:
    """Per-2D-node boolean mask of "is used as a corner by some element".

    A node is flagged True if any element references it in one of its
    first ``NUM_CORNERS_2D[type]`` connectivity slots.  Mid-edge-only
    nodes (used only in higher slots) get False.

    In a well-formed quadratic mesh, every 2D node is exclusively a
    corner or exclusively a mid-edge.  In quadratic mode, only the
    corner subset is replicated on the inserted *mid-layers* of the 3D
    mesh; mid-edge-only 2D nodes have no vertical-edge counterpart.
    """
    types_array = np.asarray(types_array, dtype=np.int64)
    conn_array  = np.asarray(conn_array,  dtype=np.int64)
    is_corner = np.zeros(n_nodes_2d, dtype=bool)
    for e in range(types_array.size):
        t = int(types_array[e])
        nc = NUM_CORNERS_2D.get(t, NUM_NODES.get(t, 0))
        for j in range(nc):
            nid = int(conn_array[e, j])             # 1-based; 0 = padding
            if nid > 0:
                is_corner[nid - 1] = True
    return is_corner


def compute_layer_storage(z_user_normalized: np.ndarray,
                          quadratic: bool) -> dict:
    """Plan how the 3D node layers are stored.

    Parameters
    ----------
    z_user_normalized
        Output of ``_normalize_z_list``: sorted ascending, includes 0.0,
        no duplicates.  Length L >= 2.
    quadratic
        If True, an extra mid-layer is inserted between every pair of
        consecutive user z-layers (these mid-layers carry the vertical
        mid-edge nodes of HEXA20 / WEDGE15 / QUAD8 / LINE3 elements).

    Returns
    -------
    dict with keys:
        z_storage      (np.ndarray) z values in storage order; storage
                       position 0 always corresponds to z=0.
        layer_order    (list[int])  layer_order[storage_pos] = index into
                                    z_all (the full sorted-ascending z
                                    array, including any mid-layers).
        layer_pos      (dict)       inverse of ``layer_order``.
        layer_is_full  (list[bool]) True at every storage position whose
                                    z corresponds to a user z-layer
                                    (in linear mode, all True).
        z_all          (np.ndarray) all z values in ascending order
                                    (= z_user_normalized in linear mode,
                                    or z's interleaved with midpoints).
        n_full         (int)        number of full layers (= L).
        n_total        (int)        total number of layers (= L in
                                    linear mode, 2L-1 in quadratic).
    """
    z = np.asarray(z_user_normalized, dtype=np.float64)
    L = z.size
    if L < 2:
        raise ValueError("z_user_normalized must have at least 2 entries "
                         "(including the implicit z=0 layer).")
    if quadratic:
        z_all = np.empty(2 * L - 1, dtype=np.float64)
        z_all[0::2] = z
        z_all[1::2] = 0.5 * (z[:-1] + z[1:])
        is_full_in_all = [(k % 2 == 0) for k in range(z_all.size)]
    else:
        z_all = z.copy()
        is_full_in_all = [True] * L

    # Find z=0 in z_all (must be present and unique by construction).
    k_zero_arr = np.where(np.isclose(z_all, 0.0, atol=1e-12))[0]
    assert k_zero_arr.size == 1, "internal error: z=0 missing from z_all"
    k_zero = int(k_zero_arr[0])

    layer_order = [k_zero] + [k for k in range(z_all.size) if k != k_zero]
    layer_pos   = {k: pos for pos, k in enumerate(layer_order)}
    z_storage     = z_all[layer_order]
    layer_is_full = [is_full_in_all[k] for k in layer_order]

    return {
        "z_storage":     z_storage,
        "layer_order":   layer_order,
        "layer_pos":     layer_pos,
        "layer_is_full": layer_is_full,
        "z_all":         z_all,
        "n_full":        L,
        "n_total":       z_all.size,
    }


# -----------------------------------------------------------------------------
# Main extrusion routine
# -----------------------------------------------------------------------------


def extrude_mesh(input_path: str | Path,
                     output_path: str | Path,
                     z=None):
    """Extrude a 2D openCFS mesh into a 3D mesh.

    The mode (linear vs quadratic) is auto-detected from the 2D element
    types: if any of LINE3 / TRIA6 / QUAD8 is present, quadratic mode
    is used and an extra "mid-layer" of nodes is inserted between every
    pair of user z-layers.  Otherwise the existing linear mapping is
    used and only the user z-layers become node layers.

    Parameters
    ----------
    input_path  : path to the 2D ``.cfs`` file
    output_path : path of the 3D ``.cfs`` file to write (overwritten)
    z           : list/array of additional z-coordinates for the node
                  layers.  z=0 (the original 2D plane) is always included
                  automatically, so passing e.g. ``[0.5, 1.0]`` produces
                  three layers at z = 0, 0.5, 1.0.  Empty / ``None`` /
                  ``[0.0]`` is treated as a no-op and prints a warning.
                  Negative and mixed-sign values are supported.

    Returns
    -------
    False on no-op (z-list empty / only contains 0.0), or a metadata
    dict on success with keys: ``quadratic``, ``z_user``, ``z_storage``,
    ``z_all``, ``layer_order``, ``layer_pos``, ``layer_is_full``,
    ``n_nodes2D``, ``n_nodes3D``, ``n_elem_layers``.  The dict is truthy,
    so callers using ``if extrude_mesh(...)`` continue to work.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)

    z = _normalize_z_list(z)
    if z.size < 2:
        print("WARNING: nothing to extrude — z-list is empty or contains only "
              "0.0 (the original 2D plane). No output file written.\n"
              "         Provide at least one non-zero z value, e.g. "
              "  ...  mesh3D.cfs 0.5 1.0 1.5")
        return False

    with h5py.File(input_path, "r") as fin, h5py.File(output_path, "w") as fout:
        # -------------------------------------------------------------------
        # Read 2D mesh
        # -------------------------------------------------------------------
        mesh_in = fin["/Mesh"]
        if int(mesh_in.attrs["Dimension"]) != 2:
            raise ValueError(
                f"Input mesh has Dimension={int(mesh_in.attrs['Dimension'])}, "
                f"this script only extrudes 2D meshes."
            )

        coords2 = mesh_in["Nodes/Coordinates"][:]            # (N, 3), z=0
        conn2   = mesh_in["Elements/Connectivity"][:]        # (E, M) 1-based
        types2  = mesh_in["Elements/Types"][:].astype(np.int32)

        n_nodes2 = coords2.shape[0]
        n_elems2 = conn2.shape[0]

        # -------------------------------------------------------------------
        # Auto-detect mode and pick the matching extrusion table
        # -------------------------------------------------------------------
        quadratic   = is_quadratic_2d_mesh(types2)
        extrude_map = EXTRUSION_MAP_QUADRATIC if quadratic else EXTRUSION_MAP_LINEAR
        valid_types = set(extrude_map)
        for t in np.unique(types2):
            if int(t) not in valid_types:
                mode = "quadratic" if quadratic else "linear"
                supported = sorted(valid_types)
                raise NotImplementedError(
                    f"2D element type ID {int(t)} is not extrudable in "
                    f"{mode} mode (supported in {mode}: {supported}).\n"
                    f"  linear  : LINE2={ET_LINE2}, "
                    f"TRIA3={ET_TRIA3}, QUAD4={ET_QUAD4}\n"
                    f"  quadratic: LINE3={ET_LINE3}, "
                    f"TRIA6={ET_TRIA6}, QUAD8={ET_QUAD8}\n"
                    f"  (mode is auto-detected from the 2D element types; "
                    f"mixing linear surface and quadratic surface elements "
                    f"is not supported)."
                )

        # -------------------------------------------------------------------
        # Layer storage planning.
        #
        # In linear mode:  z_all == z, n_total == L,  layer_is_full == [True]*L
        # In quadratic mode:
        #   z_all is z with midpoint layers interleaved (length 2L-1).
        #   Storage position 0 always corresponds to z=0 so 2D node id i
        #   stays at 3D node id i.
        # -------------------------------------------------------------------
        layout         = compute_layer_storage(z, quadratic)
        z_storage      = layout["z_storage"]
        layer_order    = layout["layer_order"]
        layer_pos      = layout["layer_pos"]
        layer_is_full  = layout["layer_is_full"]
        z_all          = layout["z_all"]
        n_total        = layout["n_total"]   # total node layers (full + mid)
        n_full         = layout["n_full"]    # = L (number of user z-layers)
        n_elem_layers  = n_full - 1          # element layers (same in both modes)

        if quadratic:
            is_corner_2d = compute_corner_mask_2d(types2, conn2, n_nodes2)

        # -------------------------------------------------------------------
        # Build 3D nodes.  Every storage layer holds a full copy of the 2D
        # (x, y) coordinates.  In quadratic mode, the mid-layers contain
        # extra mid-edge-only 2D nodes that are *not* referenced by any
        # element nor by any region's Nodes list — they sit in the file as
        # harmless extra coordinate entries.  This keeps the node-id math
        # uniform: 3D node id = pos * n_nodes2 + (2D node id).
        # -------------------------------------------------------------------
        n_nodes3 = n_nodes2 * n_total
        coords3 = np.empty((n_nodes3, 3), dtype=np.float64)
        for pos in range(n_total):
            sl = slice(pos * n_nodes2, (pos + 1) * n_nodes2)
            coords3[sl, 0] = coords2[:, 0]
            coords3[sl, 1] = coords2[:, 1]
            coords3[sl, 2] = z_storage[pos]

        # -------------------------------------------------------------------
        # Build 3D connectivity.
        #
        # An "element layer" connects two consecutive *full* layers.  In
        # quadratic mode, the mid-layer between them supplies the four
        # vertical mid-edge nodes that HEXA20 / WEDGE15 / QUAD8 / LINE3
        # elements need.
        # -------------------------------------------------------------------
        types3_per_2delem = np.array([extrude_map[int(t)] for t in types2],
                                     dtype=np.int32)
        max_nodes3 = int(max(NUM_NODES[int(t)] for t in types3_per_2delem))

        n_elems3 = n_elems2 * n_elem_layers
        conn3  = np.zeros((n_elems3, max_nodes3), dtype=np.int32)
        types3 = np.empty(n_elems3, dtype=np.int32)

        for k in range(n_elem_layers):
            if quadratic:
                # In z_all, full layers sit at even indices, mids at odd.
                bot_in_all = 2 * k
                top_in_all = 2 * (k + 1)
                mid_in_all = 2 * k + 1
                offset_bot = layer_pos[bot_in_all] * n_nodes2
                offset_top = layer_pos[top_in_all] * n_nodes2
                offset_mid = layer_pos[mid_in_all] * n_nodes2
            else:
                offset_bot = layer_pos[k]     * n_nodes2
                offset_top = layer_pos[k + 1] * n_nodes2

            for e in range(n_elems2):
                t2 = int(types2[e])
                t3 = extrude_map[t2]
                n2_count = NUM_NODES[t2]
                n2_ids = conn2[e, :n2_count]                  # 1-based
                bot = (offset_bot + n2_ids).astype(np.int32)
                top = (offset_top + n2_ids).astype(np.int32)

                e3 = k * n_elems2 + e                          # 0-based row

                if quadratic:
                    nc = NUM_CORNERS_2D[t2]
                    # Vertical mid-edge nodes only at the *corner* slots.
                    mid = (offset_mid + n2_ids[:nc]).astype(np.int32)
                    row = _build_3d_connectivity_row_quad(t3, bot, top, mid)
                else:
                    row = _build_3d_connectivity_row(t3, bot, top)
                conn3[e3, :NUM_NODES[t3]] = row
                types3[e3] = t3

        # -------------------------------------------------------------------
        # Element-type counts (for /Mesh/Elements attributes)
        # -------------------------------------------------------------------
        type_count = {tid: 0 for tid in NUM_NODES}
        for t in types3:
            type_count[int(t)] += 1

        n1d = sum(c for tid, c in type_count.items() if DIM_OF_TYPE[tid] == 1)
        n2d = sum(c for tid, c in type_count.items() if DIM_OF_TYPE[tid] == 2)
        n3d = sum(c for tid, c in type_count.items() if DIM_OF_TYPE[tid] == 3)

        # -------------------------------------------------------------------
        # Build 3D regions
        #
        # In linear mode, every storage layer contributes the full rn2.
        # In quadratic mode, *full* layers contribute the full rn2 but
        # *mid* layers contribute only the corner subset of rn2.  Either
        # way, storage position 0 (z=0) comes first, so the original 2D
        # region nodes are still the prefix of the 3D region's Nodes list.
        # -------------------------------------------------------------------
        def _expand_2d_node_list(n2_arr: np.ndarray) -> np.ndarray:
            """Replicate a 2D node-id array across all storage layers."""
            chunks: list[np.ndarray] = []
            if quadratic:
                corner_mask = is_corner_2d[n2_arr - 1]
                n2_corners = n2_arr[corner_mask]
                for pos in range(n_total):
                    offset = pos * n_nodes2
                    if layer_is_full[pos]:
                        chunks.append(offset + n2_arr)
                    else:
                        chunks.append(offset + n2_corners)
            else:
                # Iterate in storage order so z=0 chunk is first.
                for k in layer_order:
                    chunks.append(layer_pos[k] * n_nodes2 + n2_arr)
            return np.concatenate(chunks).astype(np.int32)

        regions_in = mesh_in["Regions"]
        regions_out_data = {}
        for region_name in regions_in:
            r_in = regions_in[region_name]
            re2 = r_in["Elements"][:].astype(np.int32)        # 1-based
            rn2 = r_in["Nodes"][:].astype(np.int32)           # 1-based
            r_dim_in = int(r_in.attrs.get("Dimension", 0))
            # element layers stack the 2D element index space (unchanged)
            re3 = np.concatenate(
                [(k * n_elems2) + re2 for k in range(n_elem_layers)]
            ).astype(np.int32)
            rn3 = _expand_2d_node_list(rn2)
            regions_out_data[region_name] = {
                "elements": re3,
                "nodes":    rn3,
                "dim":      r_dim_in + 1,                     # 2D area -> 3D volume
            }

        # -------------------------------------------------------------------
        # Build 3D groups (mirroring 2D groups; may be empty)
        # -------------------------------------------------------------------
        groups_out_data = {}
        if "Groups" in mesh_in:
            for grp_name in mesh_in["Groups"]:
                g_in = mesh_in[f"Groups/{grp_name}"]
                gd = {}
                if "Nodes" in g_in:
                    n2 = g_in["Nodes"][:].astype(np.int32)
                    gd["nodes"] = _expand_2d_node_list(n2)
                if "Elements" in g_in:
                    e2 = g_in["Elements"][:].astype(np.int32)
                    gd["elements"] = np.concatenate(
                        [(k * n_elems2) + e2 for k in range(n_elem_layers)]
                    ).astype(np.int32)
                groups_out_data[grp_name] = gd

        # -------------------------------------------------------------------
        # Warn about unhandled results
        # -------------------------------------------------------------------
        if "Results" in fin and len(fin["Results"]) > 0:
            print("WARNING: input file contains /Results data which is "
                  "NOT extruded by this script (geometry only).")

        # -------------------------------------------------------------------
        # ----------------------- Write the 3D file ------------------------
        # -------------------------------------------------------------------

        # FileInfo: copy verbatim if it exists, otherwise create a minimal one
        if "FileInfo" in fin:
            fin.copy("FileInfo", fout)
        else:
            fi = fout.create_group("FileInfo")
            for name, value in [
                ("Content",  "Grid"),
                ("Creator",  "extrude_mesh.py"),
                ("Date",     datetime.now().isoformat()),
                ("Version",  "1.0"),
            ]:
                fi.create_dataset(name, data=np.bytes_(value))

        # /Mesh
        mesh_out = fout.create_group("Mesh")
        _set_uint(mesh_out, "Dimension", 3)

        # /Mesh/Nodes
        nodes_out = mesh_out.create_group("Nodes")
        _set_uint(nodes_out, "NumNodes", n_nodes3)
        nodes_out.create_dataset("Coordinates", data=coords3)

        # /Mesh/Elements
        elems_out = mesh_out.create_group("Elements")
        elems_out.create_dataset("Connectivity", data=conn3)
        elems_out.create_dataset("Types",        data=types3)

        _set_uint(elems_out, "NumElems",     n_elems3)
        _set_uint(elems_out, "Num1DElems",   n1d)
        _set_uint(elems_out, "Num2DElems",   n2d)
        _set_uint(elems_out, "Num3DElems",   n3d)
        _set_int (elems_out, "QuadraticElems", 1 if quadratic else 0)

        for type_id, attr_name in ALL_TYPE_ATTRS.items():
            _set_uint(elems_out, attr_name, type_count.get(type_id, 0))
        for attr_name in NO_ENUM_TYPE_ATTRS:
            _set_uint(elems_out, attr_name, 0)

        # /Mesh/Groups (preserve, even if empty)
        groups_out = mesh_out.create_group("Groups")
        for grp_name, gd in groups_out_data.items():
            g = groups_out.create_group(grp_name)
            if "nodes" in gd:
                g.create_dataset("Nodes", data=gd["nodes"])
            if "elements" in gd:
                g.create_dataset("Elements", data=gd["elements"])

        # /Mesh/Regions
        regions_out = mesh_out.create_group("Regions")
        for region_name, rd in regions_out_data.items():
            r = regions_out.create_group(region_name)
            _set_uint(r, "Dimension", rd["dim"])
            r.create_dataset("Elements", data=rd["elements"])
            r.create_dataset("Nodes",    data=rd["nodes"])

        # /Results: preserve as empty placeholder, matching the source schema
        fout.create_group("Results")

        # /UserData: copy if present, otherwise create empty
        if "UserData" in fin:
            fin.copy("UserData", fout)
        else:
            fout.create_group("UserData")

    # ----------------------------------------------------------------------
    # Summary
    # ----------------------------------------------------------------------
    mode_label = "quadratic" if quadratic else "linear"
    print(f"Wrote {output_path}")
    print(f"  mode              : {mode_label}")
    print(f"  user z layers     : {n_full}  ({z.tolist()})")
    if quadratic:
        print(f"  total z layers    : {n_total}  "
              f"(={n_full} full + {n_total - n_full} mid)")
    print(f"  element layers    : {n_elem_layers}")
    print(f"  nodes    2D -> 3D : {n_nodes2} -> {n_nodes3}")
    print(f"  elements 2D -> 3D : {n_elems2} -> {n_elems3}")
    print(f"  Num1D / Num2D / Num3D : {n1d} / {n2d} / {n3d}")

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Extrude a 2D openCFS .cfs mesh into a 3D mesh.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Notes:\n"
            "  z=0 is the original 2D plane and is included automatically;\n"
            "  give only the *additional* z-layers.  Order does not matter,\n"
            "  negative and mixed-sign values are supported.\n"
            "\n"
            "Examples:\n"
            "  python3 %(prog)s mesh2D.cfs mesh3D.cfs 0.5 1.0 1.5 2.0\n"
            "  python3 %(prog)s mesh2D.cfs mesh3D.cfs -1 -0.5 0.5 1\n"
            "  python3 %(prog)s mesh2D.cfs mesh3D.cfs            # no-op + warning\n"
        ),
    )
    p.add_argument("input",  help="input 2D .cfs file")
    p.add_argument("output", help="output 3D .cfs file")
    p.add_argument("z", nargs="*", type=float,
                   help="additional z-coordinates for the node layers "
                        "(z=0 is included automatically; any order; "
                        "negatives ok). With no values given, the script "
                        "warns and exits without writing output.")
    args = p.parse_args()
    extrude_mesh(args.input, args.output, args.z)