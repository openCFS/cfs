#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
"""
postprocess_2d5_to_3d.py
========================

Postprocess one or more 2.5D openCFS result files into a single 3D
openCFS result file whose schema is indistinguishable from a normal 3D
harmonic analysis.

Each invocation performs ONE inverse Fourier transform at ONE 3D
excitation frequency and writes (or appends) ONE Step into the output
file.  The input files together hold the wavenumber spectrum for that
single excitation frequency:

  * 1 input file  → the file holds the entire spectrum.
  * N input files → each file holds a *partial* spectrum (a slice of
    the k_z sweep).  This supports distributed-storage workflows where
    multiple solver subjobs each compute part of the spectrum.  The
    script merges the slices (sort by k_z, drop near-coincident
    duplicate samples) and runs the IFT on the merged spectrum.

Old workflow (postprocess_2d5.py):
    one 2.5D file -> append MultiStep_2 with the IFT solution at user z's,
    on the *2D* mesh (the z's are step values).

New workflow (this script):
    list of 2.5D files (slices of one excitation freq's spectrum)
    + user-supplied z-list (only on first call)
    + a single excitation_freq (Hz)
    -> one 3D extruded mesh,
       /Results/Mesh/MultiStep_1 with the new Step inserted in
       monotonic-ascending StepValue order,
       Step value = the 3D excitation frequency (Hz),
       results are nodal pressures on the 3D mesh.

The result is a "regular-looking" 3D harmonic .cfs file.  All the 2.5D
machinery (wavenumber spectrum merging, inverse Fourier integration) is
hidden in the postprocessing step, so any downstream tool that consumes
3D harmonic openCFS output (ParaView CFSReader, etc.) just works.

How nodes line up
-----------------
The 3D mesh is built by ``extrude_mesh.extrude_mesh``, which
preserves 2D node ids:

    2D node id i  ->  3D node id i        (z = 0 layer, storage position 0)
    2D node id i  ->  3D node id (layer_pos[k] * n_nodes2D + i) for layer k

Inside a 3D *region's* ``Nodes`` dataset, the layers are stored back-to-back
in the same order as ``layer_order`` from the extrusion script:

    [   layer_order[0] chunk   |   layer_order[1] chunk   |  ...  ]
       (n_region_nodes_2D)        (n_region_nodes_2D)

Therefore, when we evaluate the inverse Fourier transform at a list of
z-values, the only mapping work needed is to:

  1. compute the IFT at z = z_normalized[layer_order[i]] for each storage
     position i,
  2. concatenate those 2D solution chunks in storage order,

which produces exactly the 1D vector that openCFS expects to live next to
the 3D region's ``Nodes`` dataset.

Inverse Fourier transform
-------------------------
For a 2.5D acoustic simulation with single-sided positive-k_z spectrum
(k_z >= 0) and an even spectrum in k_z, the physical 3D field at z is

    p(x, y, z) = (1 / pi) * integral_{0}^{k_max} P(x, y, k_z) cos(k_z z) dk_z

We approximate the integral with the trapezoidal rule on the discrete
wavenumber samples.  The result is identical to the formula already used
in ``postprocess_2d5.py`` when the wavenumber grid is uniform; we use
generic trapezoidal weights internally so non-uniform wavenumber grids
(e.g. produced by merging differently-sampled subjob outputs) are also
handled correctly.

Spectrum merging across files
-----------------------------
When more than one input file is given, each file's MultiStep_1 sweep is
read (its StepValues, in Hz, give k_z = 2π·f_step/c0 for each sample).
All samples are concatenated, sorted ascending by k_z, and any pair
within ``|Δk| ≤ 1e-9`` is treated as a duplicate — the entry from the
earliest-listed file wins, the rest are dropped (with a warning).  The
single merged ``(k_z, P_re, P_im)`` table per region is then fed to the
trapezoidal IFT.  All input files must share the same 2D mesh and the
same set of result regions.

Create vs. append
-----------------
If ``output_file`` does not exist, the script extrudes the 2D mesh of
the first input file into a 3D mesh using ``z_list`` (so ``z_list`` is
required) and writes a fresh ``MultiStep_1`` with the single new Step.

If ``output_file`` already exists and contains a 3D mesh built by
``extrude_mesh.extrude_mesh``, the script reuses that mesh
(``z_list`` is optional in this case).  The new step is inserted into
``MultiStep_1`` at the position that keeps StepValues monotonic-
ascending — existing Steps may be renumbered.  If a step at the same
``excitation_freq`` already exists, the call is a no-op (warning).

Usage
-----

    # CLI: create from a single full-spectrum input file
    python3 postprocess_2d5_to_3d.py \\
        -o results3D.cfs \\
        --z '0.1,0.2,0.3,0.4,0.5' \\
        --frequency 45.52 \\
        DANI2d5_45p52Hz_full.cfs

    # CLI: create from three subjob outputs covering parts of the spectrum
    python3 postprocess_2d5_to_3d.py \\
        -o results3D.cfs \\
        --z '0.1,0.2,0.3,0.4,0.5' \\
        --frequency 45.52 \\
        DANI2d5_45p52Hz_part1.cfs \\
        DANI2d5_45p52Hz_part2.cfs \\
        DANI2d5_45p52Hz_part3.cfs

    # CLI: append a different excitation frequency
    python3 postprocess_2d5_to_3d.py \\
        -o results3D.cfs \\
        --frequency 112.0 \\
        DANI2d5_112Hz_part1.cfs DANI2d5_112Hz_part2.cfs

    # Python:
    from postprocess_2d5_to_3d import postprocess_2d5_to_3d
    postprocess_2d5_to_3d(
        input_files=['DANI2d5_45p52Hz_part1.cfs',
                     'DANI2d5_45p52Hz_part2.cfs'],
        output_file='DANI3D.cfs',
        excitation_freq=45.52,                 # single scalar
        z_list=[0.1, 0.2, 0.3, 0.4, 0.5],      # required only on first call
    )

Author
------
Likun Luo (likun.luo@tuwien.ac.at)

Created
-------
2026-05-12

License
-------
BSD-3-Clause
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np

# Re-use the extrusion routine and z-list normalization that we developed
# in extrude_mesh.py.  The 3D mesh of the output file is produced
# by extrude_mesh() applied to the 2D mesh of the first input file.
# We also import the helpers that describe the layer-storage layout (which
# now knows about mid-layers in quadratic mode) and that classify 2D nodes
# as corners vs. mid-edges.
from extrude_mesh import (
    extrude_mesh,
    _normalize_z_list,
    is_quadratic_2d_mesh,
    compute_corner_mask_2d,
    compute_layer_storage,
)


# -----------------------------------------------------------------------------
# Defaults: air at ~20 °C, c = 343 m/s
# -----------------------------------------------------------------------------
DEFAULT_RHO = 1.205        # kg/m^3
DEFAULT_K   = 1.41767e5    # Pa  (bulk modulus)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def _inspect_existing_3d_mesh(path: Path) -> dict:
    """Read the layer layout from an existing 3D openCFS mesh file.

    Used when ``output_file`` already exists and contains a 3D mesh: we
    reuse that mesh and append result steps to it rather than rebuilding.

    Returns
    -------
    dict with keys:
        quadratic        bool
        n_nodes_2d       int
        z_user           ndarray, the user z-list this mesh was built from
                         (z=0 included; sorted ascending; no mid-layers).
        z_storage        ndarray, layer z's in storage order (z=0 first,
                         then ascending; matches what extrude_mesh
                         produced).
        layer_is_full    list[bool]
        layer_order, layer_pos, z_all, n_full, n_total
                         same as ``compute_layer_storage`` returns.
        region_names     set[str], regions present in /Mesh/Regions
        coords_2d        ndarray (n_nodes_2d, 3), the z=0 layer of the
                         mesh — used to sanity-check the input file's
                         2D mesh matches.

    Raises
    ------
    ValueError if the file does not contain a 3D mesh, or its layout
    does not look like something built by ``extrude_mesh`` (e.g.
    node count not divisible by the unique-z count, or coordinates
    don't match the storage convention used by this script).
    """
    with h5py.File(path, 'r') as f:
        if '/Mesh' not in f:
            raise ValueError(
                f"Existing output file '{path}' has no /Mesh group; "
                f"cannot reuse it as a 3D mesh."
            )
        dim = int(f['/Mesh'].attrs.get('Dimension', 0))
        if dim != 3:
            raise ValueError(
                f"Existing output file '{path}' has /Mesh/Dimension={dim}; "
                f"need a 3D mesh (Dimension=3) to append result steps."
            )
        coords = f['/Mesh/Nodes/Coordinates'][:]
        quadratic = bool(int(
            f['/Mesh/Elements'].attrs.get('QuadraticElems', 0)
        ))
        region_names = set(f['/Mesh/Regions'].keys()) \
            if '/Mesh/Regions' in f else set()

    n_total_nodes_3d = coords.shape[0]
    if n_total_nodes_3d == 0:
        raise ValueError(f"Existing output file '{path}' has no nodes.")

    # Unique z layers, ascending. ``extrude_mesh`` writes exact
    # float64 z's, so np.unique returns the layer set without FP slop.
    z_all_sorted = np.unique(coords[:, 2])
    n_total = z_all_sorted.size
    if n_total < 2:
        raise ValueError(
            f"Existing output file '{path}' has only {n_total} unique z "
            f"value(s); need at least 2 layers."
        )
    if n_total_nodes_3d % n_total != 0:
        raise ValueError(
            f"Existing output file '{path}': total node count "
            f"{n_total_nodes_3d} is not divisible by the {n_total} unique "
            f"z layers found in /Mesh/Nodes/Coordinates.  This file does "
            f"not look like it was built by extrude_mesh."
        )
    n_nodes_2d = n_total_nodes_3d // n_total

    if quadratic:
        if n_total % 2 != 1:
            raise ValueError(
                f"Existing output file '{path}' has QuadraticElems=1 but "
                f"{n_total} unique z layers (expected an odd number: "
                f"L user layers + L-1 mid layers)."
            )
        z_user = z_all_sorted[::2]
    else:
        z_user = z_all_sorted

    # Reconstruct the storage layout we *would* have written, then verify
    # it matches the actual file (storage position 0 → z=0, etc.).
    layout = compute_layer_storage(z_user, quadratic=quadratic)

    actual_z_storage = coords[::n_nodes_2d, 2]   # one z per layer (first node)
    if not np.allclose(actual_z_storage, layout['z_storage'], atol=1e-9):
        raise ValueError(
            f"Existing output file '{path}': node z-ordering does not match "
            f"the convention used by extrude_mesh (z=0 storage layer "
            f"first, then ascending).  Got per-layer z's = "
            f"{actual_z_storage.tolist()}; expected = "
            f"{layout['z_storage'].tolist()}.  This file may have been "
            f"produced by a different tool; please remove it or use a "
            f"different output filename."
        )

    coords_2d = coords[:n_nodes_2d, :].copy()    # the z=0 layer (storage pos 0)

    return {
        'quadratic':     quadratic,
        'n_nodes_2d':    n_nodes_2d,
        'z_user':        z_user,
        'z_storage':     layout['z_storage'],
        'layer_is_full': layout['layer_is_full'],
        'layer_order':   layout['layer_order'],
        'layer_pos':     layout['layer_pos'],
        'z_all':         layout['z_all'],
        'n_full':        layout['n_full'],
        'n_total':       layout['n_total'],
        'region_names':  region_names,
        'coords_2d':     coords_2d,
    }


def _trap_weights(k: np.ndarray) -> np.ndarray:
    """Trapezoidal weights ``w_i`` such that ``sum_i w_i f_i ~= integral f(k)dk``.

    Works for arbitrary (sorted) ``k`` grids.  For the uniform case this
    reproduces ``[dk/2, dk, dk, ..., dk, dk/2]``.
    """
    if k.size < 2:
        raise ValueError("Need at least two wavenumber samples for the trapezoidal rule.")
    w = np.empty_like(k, dtype=np.float64)
    w[0]    = (k[1]    - k[0])    * 0.5
    w[-1]   = (k[-1]   - k[-2])   * 0.5
    w[1:-1] = (k[2:]   - k[:-2])  * 0.5
    return w


def _read_partial_spectrum_one_file(h5_in: h5py.File,
                                    c0: float,
                                    region: str
                                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a single 2.5D file's partial wavenumber spectrum for one region.

    The 2.5D solver stores its discrete sweep in ``MultiStep_1``.  Each
    ``Step_i`` carries one axial-frequency value (the StepValue, in Hz)
    and the complex nodal pressure ``P(k_z = 2π·f_step / c0)`` on the 2D
    region.

    Returns
    -------
    k     : (n_steps,) ndarray of axial wavenumbers k_z in this file
    p_re  : (n_steps, n_region) ndarray, real part of P at each k
    p_im  : (n_steps, n_region) ndarray, imag part of P at each k
    """
    ms = h5_in['Results/Mesh/MultiStep_1']
    n_steps = int(ms.attrs['LastStepNum'])
    if n_steps < 1:
        raise ValueError(
            f"Input file has LastStepNum={n_steps}; need at least 1 step."
        )

    k    = np.empty(n_steps, dtype=np.float64)
    p_re = None
    p_im = None
    for i in range(1, n_steps + 1):
        step = h5_in[f'Results/Mesh/MultiStep_1/Step_{i}']
        f_step = float(step.attrs['StepValue'])
        k[i - 1] = (2.0 * np.pi * f_step) / c0

        re = step[f'acouPressure/{region}/Nodes/Real'][:].ravel()
        im = step[f'acouPressure/{region}/Nodes/Imag'][:].ravel()
        if p_re is None:
            n_region = re.size
            p_re = np.empty((n_steps, n_region), dtype=np.float64)
            p_im = np.empty((n_steps, n_region), dtype=np.float64)
        p_re[i - 1, :] = re
        p_im[i - 1, :] = im
    return k, p_re, p_im


def _load_combined_spectrum(input_files: list[Path],
                            c0: float,
                            regions: list[str],
                            duplicate_k_atol: float = 1e-9
                            ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Read the wavenumber spectrum across ALL input files and merge it.

    Each input file may hold a *partial* spectrum (a subset of the
    wavenumber sweep — e.g. one batch job covers k_z=0..1, another
    1.025..2, etc.).  This function concatenates all of them, sorts by
    k_z ascending, and removes any near-coincident duplicate k_z samples
    (with the corresponding pressure values) so the trapezoidal IFT sees
    a single clean monotonic grid.

    Parameters
    ----------
    input_files
        List of 2.5D .cfs file paths to merge.
    c0
        Speed of sound (used to convert StepValue [Hz] back to k_z).
    regions
        Region names whose nodal pressure spectra should be loaded.
    duplicate_k_atol
        Two wavenumber samples are considered the same if
        |k_a - k_b| <= duplicate_k_atol.  When that happens, the entry
        from the *first* file (in the order given) is kept and the
        others are dropped.  A warning is printed.

    Returns
    -------
    dict mapping region -> (k, P_re, P_im) where k has shape (N,) and
    P_re/P_im have shape (N, n_region_2d).  N is the union of all
    distinct wavenumbers across the input files.

    Raises
    ------
    ValueError
        if the files have inconsistent region sizes, or if the merged
        spectrum has fewer than 2 distinct wavenumbers (the trapezoidal
        rule needs at least two samples).
    """
    # Per-region accumulator: list of (origin_file_index, k, P_re, P_im).
    per_region: dict[str, list[tuple[int, np.ndarray, np.ndarray, np.ndarray]]] = {
        r: [] for r in regions
    }

    for file_idx, path in enumerate(input_files):
        with h5py.File(path, 'r') as fin:
            for region in regions:
                k, p_re, p_im = _read_partial_spectrum_one_file(fin, c0, region)
                per_region[region].append((file_idx, k, p_re, p_im))

    # Validate consistent region sizes across files.
    for region, parts in per_region.items():
        n_region_first = parts[0][2].shape[1]
        for file_idx, _, p_re, _ in parts[1:]:
            if p_re.shape[1] != n_region_first:
                raise ValueError(
                    f"Region '{region}' has {n_region_first} 2D nodes in "
                    f"'{input_files[0].name}' but {p_re.shape[1]} in "
                    f"'{input_files[file_idx].name}'.  All input files "
                    f"must share the same 2D mesh and result regions."
                )

    # Concatenate, sort, dedup.  The same wavenumber grid is shared by all
    # regions in a 2.5D file (it's just MultiStep_1 step values), so we
    # detect duplicate-k once on the first region and apply the same
    # selection mask to every region — that way the n_k axes stay aligned.
    merged: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    # Build the master k-array (concatenation of all files' k-arrays, in
    # the order they appeared).  The dedup mask we compute on the first
    # region is reused for the other regions to keep them in lock-step.
    parts_first = per_region[regions[0]]
    file_idx_per_sample = np.concatenate([
        np.full(k.size, fi, dtype=np.int64) for fi, k, _, _ in parts_first
    ])
    k_all = np.concatenate([k for _, k, _, _ in parts_first])

    # Sort by k.
    sort_idx = np.argsort(k_all, kind='stable')
    k_sorted              = k_all[sort_idx]
    file_idx_sorted       = file_idx_per_sample[sort_idx]

    # Mark duplicates: keep the first occurrence (smallest file_idx
    # at this k), drop subsequent ones within tolerance.
    keep_mask = np.ones(k_sorted.size, dtype=bool)
    duplicate_pairs: list[tuple[float, str, str]] = []  # (k, kept_file, dropped_file)
    i = 0
    n = k_sorted.size
    while i < n:
        j = i + 1
        while j < n and (k_sorted[j] - k_sorted[i]) <= duplicate_k_atol:
            # Drop j unless its file came earlier than the kept one.
            if file_idx_sorted[j] < file_idx_sorted[i]:
                # j wins; drop i (and let i become j for the next round)
                keep_mask[i] = False
                duplicate_pairs.append((
                    float(k_sorted[i]),
                    input_files[file_idx_sorted[j]].name,
                    input_files[file_idx_sorted[i]].name,
                ))
                # advance i to j; j+1 takes over
                i = j
            else:
                keep_mask[j] = False
                duplicate_pairs.append((
                    float(k_sorted[j]),
                    input_files[file_idx_sorted[i]].name,
                    input_files[file_idx_sorted[j]].name,
                ))
            j += 1
        i = j

    if duplicate_pairs:
        print(f"  note: dropped {len(duplicate_pairs)} duplicate wavenumber "
              f"sample{'s' if len(duplicate_pairs) != 1 else ''} from the "
              f"merged spectrum (|Δk| ≤ {duplicate_k_atol:g}):")
        # Show up to 5 examples; group the rest by count.
        for k_val, kept_name, dropped_name in duplicate_pairs[:5]:
            print(f"    k = {k_val:.6g}: kept '{kept_name}', "
                  f"dropped '{dropped_name}'")
        if len(duplicate_pairs) > 5:
            print(f"    ... and {len(duplicate_pairs) - 5} more.")

    final_idx = sort_idx[keep_mask]
    k_final   = k_all[final_idx]
    if k_final.size < 2:
        raise ValueError(
            f"After merging the input files' spectra, only "
            f"{k_final.size} distinct wavenumber sample"
            f"{' is' if k_final.size == 1 else 's are'} left; the "
            f"trapezoidal IFT needs at least 2."
        )

    # Now build the per-region P arrays using the same final_idx.
    for region, parts in per_region.items():
        p_re_all = np.concatenate([p_re for _, _, p_re, _   in parts], axis=0)
        p_im_all = np.concatenate([p_im for _, _, _,   p_im in parts], axis=0)
        merged[region] = (k_final,
                          p_re_all[final_idx, :],
                          p_im_all[final_idx, :])

    return merged


# -----------------------------------------------------------------------------
# Inverse Fourier transform
# -----------------------------------------------------------------------------

def _trapezoidal_ift(wavenumbers: np.ndarray,
                     p_re: np.ndarray,
                     p_im: np.ndarray,
                     z_eval: np.ndarray):
    """Trapezoidal inverse Fourier transform.

    Computes, for each z in ``z_eval`` and each 2D node in the region,

        p(z) = (1/π) Σ_i  w_i · cos(k_i · z) · P(k_i)

    where ``w_i`` are the trapezoidal weights of the (possibly non-uniform)
    wavenumber grid and ``P(k_i)`` is the complex nodal pressure at the
    i-th wavenumber sample.

    Parameters
    ----------
    wavenumbers  : (n_k,) wavenumbers k_z, sorted ascending
    p_re, p_im   : (n_k, n_region) real and imag parts of P at each k
    z_eval       : (n_z,) z-coordinates at which to evaluate the IFT

    Returns
    -------
    sol_re, sol_im : ndarrays of shape (n_z, n_region)
    """
    n_k = wavenumbers.size
    if n_k < 2:
        raise ValueError(
            f"Need ≥ 2 wavenumber samples for trapezoidal IFT, got {n_k}."
        )
    if p_re.shape[0] != n_k or p_im.shape[0] != n_k:
        raise ValueError(
            f"Pressure spectrum shape {p_re.shape} doesn't match "
            f"wavenumber length {n_k}."
        )

    weights = _trap_weights(wavenumbers)            # (n_k,)
    inv_pi  = 1.0 / np.pi
    coeff   = (inv_pi * weights)[None, :] * \
              np.cos(z_eval[:, None] * wavenumbers[None, :])   # (n_z, n_k)

    # Single matrix multiply per real/imag.  Memory: O(n_z * n_region).
    sol_re = coeff @ p_re                          # (n_z, n_region)
    sol_im = coeff @ p_im
    return sol_re, sol_im


# -----------------------------------------------------------------------------
# Output schema helpers
# -----------------------------------------------------------------------------

def _set_uint(obj, name: str, value: int) -> None:
    obj.attrs.create(name, np.uint32(value))


def _ensure_results_branch(fout: h5py.File) -> h5py.Group:
    """Make sure ``/Results/Mesh`` exists in the output file.

    ``extrude_mesh`` already creates an empty ``/Results`` group; we
    add ``/Results/Mesh`` here (with the ``ExternalFiles`` attribute that
    openCFS always writes on it).
    """
    results = fout.require_group('Results')
    mesh = results.require_group('Mesh')
    if 'ExternalFiles' not in mesh.attrs:
        _set_uint(mesh, 'ExternalFiles', 0)
    return mesh


def _copy_result_description(rd_src: h5py.Group,
                             rd_dst_parent: h5py.Group,
                             excitation_freqs: np.ndarray) -> None:
    """Copy the per-result schema (DOFNames, EntityNames, ...) from source
    to destination, replacing the StepNumbers/StepValues datasets with the
    new excitation-frequency list.

    rd_src         : .../MultiStep_1/ResultDescription/<resultName>   (input)
    rd_dst_parent  : .../MultiStep_1/ResultDescription                (output)
    excitation_freqs : the new step values (3D excitation frequencies in Hz)
    """
    name = rd_src.name.rsplit('/', 1)[-1]
    rd_out = rd_dst_parent.create_group(name)

    n_steps = excitation_freqs.size

    for key in rd_src:
        ds_in = rd_src[key]
        if key == 'StepNumbers':
            # Match the source's rank (real openCFS files typically use rank-1
            # here; the ProjectionOfSolutionOn3DMesh.py recipe uses rank-2).
            new_shape = (n_steps,) + tuple(ds_in.shape[1:]) if ds_in.ndim >= 1 \
                        else (n_steps,)
            data = np.arange(1, n_steps + 1,
                             dtype=ds_in.dtype if np.issubdtype(ds_in.dtype, np.integer)
                             else np.uint32).reshape(new_shape)
            rd_out.create_dataset('StepNumbers', data=data)
        elif key == 'StepValues':
            new_shape = (n_steps,) + tuple(ds_in.shape[1:]) if ds_in.ndim >= 1 \
                        else (n_steps,)
            data = excitation_freqs.astype(np.float64).reshape(new_shape)
            rd_out.create_dataset('StepValues', data=data)
        else:
            # bit-identical copy of every other dataset
            rd_src.copy(key, rd_out, name=key)


def _rewrite_result_description_steps(rd_group: h5py.Group,
                                      step_numbers: np.ndarray,
                                      step_values: np.ndarray) -> None:
    """Replace ResultDescription/StepNumbers and StepValues in place.

    rd_group       : .../MultiStep_1/ResultDescription/<resultName>
    step_numbers   : (N,) integer array, the *complete* final step list
    step_values    : (N,) float array,   the *complete* final step list

    Both datasets are replaced from scratch (not extended) so that a
    single call brings them into sync with whatever final ordering of
    Step_i groups we wrote.  The dataset rank/dtype is preserved (some
    files use rank-1, the legacy recipe uses rank-2).
    """
    for ds_name, new_data in [('StepNumbers', step_numbers),
                              ('StepValues',  step_values)]:
        ds_in = rd_group[ds_name]
        old_shape = ds_in.shape
        old_dtype = ds_in.dtype
        n = new_data.size
        new_shape = (n,) + tuple(old_shape[1:]) if len(old_shape) > 0 else (n,)
        del rd_group[ds_name]
        rd_group.create_dataset(ds_name,
                                data=new_data.reshape(new_shape).astype(old_dtype))


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------

def postprocess_2d5_to_3d(
    input_files: Sequence[str | Path],
    output_file: str | Path,
    excitation_freq: float,
    z_list: Iterable[float] | None = None,
    rho: float = DEFAULT_RHO,
    K: float = DEFAULT_K,
) -> None:
    """Build (or extend) a 3D harmonic .cfs file from 2.5D wavenumber-spectrum files.

    This call performs *one* inverse Fourier transform at *one* 3D
    excitation frequency and writes (or appends) *one* result step to
    the output file.  The input files together hold the wavenumber
    spectrum for that single excitation frequency:

      * 1 input file  → the file holds the entire spectrum.
      * N input files → each file holds a *partial* spectrum (a subset
        of the k_z sweep — typical of distributed-storage workflows
        where multiple solver subjobs each compute a slice of the
        spectrum).  The script merges the slices into one global
        spectrum (sorted by k_z, with near-coincident duplicate
        samples deduplicated) before doing the IFT.

    The script has two modes, picked automatically from whether
    ``output_file`` already exists:

    1. **Create mode** (output file does not exist).  The 2D mesh of
       the *first* input file is extruded into a 3D mesh using
       ``z_list``.  A fresh ``/Results/Mesh/MultiStep_1`` is created
       holding one Step at the chosen excitation frequency.

    2. **Append mode** (output file already exists with a 3D mesh
       built by ``extrude_mesh``).  The mesh and its z layers are
       reused (``z_list`` is optional).  The new step is inserted into
       ``MultiStep_1`` at the position that keeps StepValues
       monotonic-ascending (existing Steps may be renumbered).  If a
       step at the same excitation frequency already exists, the call
       is a no-op (warning printed).

    Parameters
    ----------
    input_files       list of paths to 2.5D .cfs files; each must
                      contain /Mesh (2D) and /Results/Mesh/MultiStep_1
                      with (a possibly partial) wavenumber spectrum
                      for the *same* excitation frequency.  All files
                      must share the same 2D mesh and result regions.
    output_file       path of the 3D .cfs file to create or extend.
    excitation_freq   the single 3D excitation frequency (Hz) of this
                      call.  Required.  Must be a scalar float.
    z_list            z-coordinates at which the inverse Fourier
                      transform is evaluated.  Required in *create
                      mode*; optional in *append mode* (the mesh's
                      layers are reused).  z=0 is included implicitly.
                      In *append mode*, if ``z_list`` is given but
                      does not match the mesh's z layers, a warning
                      is printed and the mesh's layers are used.
    rho, K            air density and bulk modulus; speed of sound is
                      c0 = sqrt(K/rho).
    """
    # ----- validate / resolve inputs --------------------------------------------
    if not input_files:
        raise ValueError("No input files given.")
    input_files = [Path(p) for p in input_files]
    for p in input_files:
        if not p.is_file():
            raise FileNotFoundError(f"Input file not found: {p}")

    output_file = Path(output_file)

    # ----- excitation frequency (single scalar) ---------------------------------
    if excitation_freq is None:
        raise ValueError(
            "excitation_freq is required: provide a single 3D excitation "
            "frequency (Hz) — the input files together hold its wavenumber "
            "spectrum."
        )
    try:
        excitation_freq = float(excitation_freq)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"excitation_freq must be a scalar float, got {excitation_freq!r}."
        ) from exc
    if not np.isfinite(excitation_freq):
        raise ValueError(
            f"excitation_freq must be a finite number, got {excitation_freq}."
        )

    # ----- decide create vs append mode -----------------------------------------
    output_exists = output_file.exists()
    existing_mesh: dict | None = None
    if output_exists:
        existing_mesh = _inspect_existing_3d_mesh(output_file)

    # ----- detect quadratic mode from the first input file ---------------------
    # All input files must share the same 2D mesh; we read the geometry from
    # the first one and validate the others below when we load their spectra.
    with h5py.File(input_files[0], 'r') as fin0:
        types2_pre  = fin0['/Mesh/Elements/Types'][:].astype(np.int32)
        conn2_pre   = fin0['/Mesh/Elements/Connectivity'][:].astype(np.int64)
        n_nodes2_in = int(fin0['/Mesh/Nodes/Coordinates'].shape[0])
        coords2_in  = fin0['/Mesh/Nodes/Coordinates'][:].astype(np.float64)
    quadratic_mode = is_quadratic_2d_mesh(types2_pre)
    if quadratic_mode:
        is_corner_2d = compute_corner_mask_2d(types2_pre, conn2_pre, n_nodes2_in)
    else:
        is_corner_2d = None

    # Validate that every other input file shares the same 2D mesh
    # (node count + (x,y) coordinates) — if they don't, merging their
    # spectra into one IFT is meaningless.
    for path in input_files[1:]:
        with h5py.File(path, 'r') as fin:
            n2 = int(fin['/Mesh/Nodes/Coordinates'].shape[0])
            xy = fin['/Mesh/Nodes/Coordinates'][:].astype(np.float64)
        if n2 != n_nodes2_in:
            raise ValueError(
                f"2D node count mismatch between '{input_files[0].name}' "
                f"({n_nodes2_in}) and '{path.name}' ({n2}).  All input "
                f"files must share the same 2D mesh."
            )
        if not np.allclose(xy[:, :2], coords2_in[:, :2], atol=1e-9):
            raise ValueError(
                f"2D coordinate mismatch between '{input_files[0].name}' "
                f"and '{path.name}'.  All input files must share the same "
                f"2D mesh."
            )

    # ----- z-layer plan: from existing mesh if appending, else from z_list -----
    if existing_mesh is not None:
        if existing_mesh['quadratic'] != quadratic_mode:
            raise ValueError(
                f"Mode mismatch: existing 3D mesh has QuadraticElems="
                f"{int(existing_mesh['quadratic'])}, but input file "
                f"'{input_files[0].name}' is "
                f"{'quadratic' if quadratic_mode else 'linear'}.  Cannot append."
            )
        if existing_mesh['n_nodes_2d'] != n_nodes2_in:
            raise ValueError(
                f"2D node-count mismatch: existing 3D mesh implies "
                f"{existing_mesh['n_nodes_2d']} 2D nodes, but input file "
                f"'{input_files[0].name}' has {n_nodes2_in}.  The input file "
                f"does not appear to share a 2D mesh with the existing 3D mesh."
            )
        if not np.allclose(coords2_in[:, :2],
                           existing_mesh['coords_2d'][:, :2], atol=1e-9):
            raise ValueError(
                f"2D coordinate mismatch: the (x,y) of input file "
                f"'{input_files[0].name}' does not match the z=0 layer of the "
                f"existing 3D mesh in '{output_file.name}'.  Refusing to append."
            )

        # Optional warning: user-provided z_list inconsistent with mesh's z's.
        if z_list is not None:
            z_norm_user = _normalize_z_list(list(z_list))
            if (z_norm_user.size != existing_mesh['z_user'].size
                    or not np.allclose(z_norm_user,
                                       existing_mesh['z_user'], atol=1e-9)):
                user_set = {round(float(v), 9) for v in z_norm_user}
                mesh_set = {round(float(v), 9) for v in existing_mesh['z_user']}
                in_mesh_only = sorted(mesh_set - user_set)
                in_user_only = sorted(user_set - mesh_set)
                msg_lines = [
                    "WARNING: --z (or z_list=) does not match the z layers of "
                    f"the existing 3D mesh in '{output_file.name}'.",
                    f"  user-provided z (sorted, 0 implied) : "
                    f"{[round(float(v), 9) for v in z_norm_user]}",
                    f"  mesh's z (sorted, 0 included)       : "
                    f"{[round(float(v), 9) for v in existing_mesh['z_user']]}",
                ]
                if in_mesh_only:
                    msg_lines.append(f"  in mesh but not in user list        : "
                                     f"{in_mesh_only}")
                if in_user_only:
                    msg_lines.append(f"  in user list but not in mesh        : "
                                     f"{in_user_only}")
                msg_lines.append(
                    "  → using the mesh's z layers (the IFT must be evaluated "
                    "exactly at the existing layer set so the new step's data "
                    "lines up with the existing 3D node ordering)."
                )
                print('\n'.join(msg_lines))

        layer_is_full      = existing_mesh['layer_is_full']
        z_in_storage_order = existing_mesh['z_storage']
        n_total_layers     = existing_mesh['n_total']
    else:
        if z_list is None:
            raise ValueError(
                f"Output file '{output_file}' does not exist yet, so a new "
                f"3D mesh has to be built — please provide z_list (the z "
                f"positions of the 3D node layers; z=0 is implicit)."
            )
        z_norm = _normalize_z_list(list(z_list))
        if z_norm.size < 2:
            raise ValueError(
                "z_list is empty or contains only 0.0; need at least one "
                "non-zero z to build a 3D mesh."
            )
        layout = compute_layer_storage(z_norm, quadratic=quadratic_mode)
        layer_is_full      = layout['layer_is_full']
        z_in_storage_order = layout['z_storage']
        n_total_layers     = layout['n_total']

    # ----- create the 3D mesh if needed -----------------------------------------
    if existing_mesh is None:
        print(f"[1/2] Extruding 2D mesh to 3D from {input_files[0].name} ...")
        extrude_mesh(input_files[0], output_file, list(z_list))
    else:
        print(f"[info] Reusing existing 3D mesh in '{output_file.name}' "
              f"({existing_mesh['n_total']} z layers, "
              f"{'quadratic' if quadratic_mode else 'linear'}).")

    # ----- speed of sound --------------------------------------------------------
    c0 = float(np.sqrt(K / rho))

    # ----- discover regions with results & their 2D node counts -----------------
    with h5py.File(input_files[0], 'r') as fin0:
        mesh_in = fin0['/Mesh']
        ms1 = fin0['Results/Mesh/MultiStep_1']
        step1 = ms1['Step_1']
        if 'acouPressure' not in step1:
            raise ValueError(
                f"Input '{input_files[0].name}' has no acouPressure result "
                f"in MultiStep_1/Step_1; nothing to inverse-transform."
            )
        result_names = ['acouPressure']
        result_regions: dict[str, list[str]] = {}
        for rn in result_names:
            result_regions[rn] = list(step1[rn].keys())
        region_nnodes_2d:     dict[tuple[str, str], int] = {}
        region_corner_in_rn2: dict[tuple[str, str], np.ndarray] = {}
        for rn, regs in result_regions.items():
            for region in regs:
                rn2 = mesh_in[f'Regions/{region}/Nodes'][:].astype(np.int64)
                region_nnodes_2d[(rn, region)] = int(rn2.size)
                if quadratic_mode:
                    region_corner_in_rn2[(rn, region)] = is_corner_2d[rn2 - 1]
                else:
                    region_corner_in_rn2[(rn, region)] = np.ones(
                        rn2.size, dtype=bool
                    )

    # All other input files must declare the same set of result regions.
    for path in input_files[1:]:
        with h5py.File(path, 'r') as fin:
            step1_in = fin['Results/Mesh/MultiStep_1/Step_1']
            for rn in result_names:
                if rn not in step1_in:
                    raise ValueError(
                        f"'{path.name}': missing /Results/.../Step_1/{rn}."
                    )
                regs_here = set(step1_in[rn].keys())
                if regs_here != set(result_regions[rn]):
                    raise ValueError(
                        f"'{path.name}': result regions for '{rn}' differ "
                        f"from the first file. First: {result_regions[rn]}. "
                        f"This: {sorted(regs_here)}."
                    )

    # In append mode, every result region must already exist in the existing
    # 3D mesh (otherwise the new step has nowhere to write results into).
    if existing_mesh is not None:
        missing = [r for r in result_regions['acouPressure']
                   if r not in existing_mesh['region_names']]
        if missing:
            raise ValueError(
                f"Cannot append: result region(s) {missing} are present in "
                f"input file '{input_files[0].name}' but not in the existing "
                f"3D mesh '{output_file.name}'."
            )

    # ----- prepare /Results/Mesh/MultiStep_1 in the output ----------------------
    with h5py.File(output_file, 'r+') as fout:
        mesh_out = _ensure_results_branch(fout)

        # ---------------- read existing steps (if any) -------------------------
        existing_steps: list[tuple[int, float]] = []  # (step_num, step_value)
        if 'MultiStep_1' in mesh_out:
            ms_out = mesh_out['MultiStep_1']
            n_existing_steps = int(ms_out.attrs.get('LastStepNum', 0))
            existing_analysis = ms_out.attrs.get('AnalysisType', None)
            if isinstance(existing_analysis, (bytes, np.bytes_)):
                existing_analysis = existing_analysis.decode('ascii', errors='ignore')
            if existing_analysis is not None and existing_analysis != 'harmonic':
                raise ValueError(
                    f"Existing MultiStep_1 has AnalysisType="
                    f"{existing_analysis!r}; can only append to a 'harmonic' "
                    f"analysis."
                )
            for i in range(1, n_existing_steps + 1):
                step_grp_name = f'Step_{i}'
                if step_grp_name not in ms_out:
                    raise ValueError(
                        f"Existing MultiStep_1 has LastStepNum={n_existing_steps} "
                        f"but '{step_grp_name}' is missing."
                    )
                sv = float(ms_out[step_grp_name].attrs['StepValue'])
                existing_steps.append((i, sv))
        else:
            ms_out = mesh_out.create_group('MultiStep_1')
            # IMPORTANT: write 'AnalysisType' as a variable-length, null-
            # terminated, ASCII string.  Two pitfalls to avoid here:
            #   * ``attrs.create(name, np.bytes_(s))`` -> fixed-length NULLPAD,
            #     which segfaults openCFS / ParaView's CFS reader
            #     (H5LTget_attribute_string -> std::string(char*) walks past
            #     the end of an unterminated buffer).
            #   * ``attrs[name] = str`` would work, but h5py picks UTF-8 as
            #     the character set; openCFS metadata uses ASCII, so we set
            #     it explicitly via ``string_dtype(encoding='ascii')``.
            # The dtype below yields VLEN + NULLTERM + CSET_ASCII.
            ms_out.attrs.create(
                'AnalysisType', 'harmonic',
                dtype=h5py.string_dtype(encoding='ascii'),
            )
            n_existing_steps = 0
            # Initialise ResultDescription from the first input file's schema.
            # (StepNumbers/StepValues will be overwritten below from the
            # canonical combined-sorted list.)
            rd_parent = ms_out.create_group('ResultDescription')
            for rn in result_names:
                with h5py.File(input_files[0], 'r') as fin0_re:
                    _copy_result_description(
                        fin0_re[f'Results/Mesh/MultiStep_1/ResultDescription/{rn}'],
                        rd_parent,
                        np.array([excitation_freq], dtype=np.float64),
                    )

        # ---------------- check for collision with existing freqs --------------
        existing_freqs_arr = np.array([sv for _, sv in existing_steps],
                                      dtype=np.float64)
        already_in_existing = (
            existing_freqs_arr.size > 0
            and bool(np.any(np.isclose(existing_freqs_arr, excitation_freq,
                                       atol=1e-9, rtol=0)))
        )

        # ---------------- build the final combined-sorted Step list ------------
        plan: list[dict] = []
        for old_num, freq in existing_steps:
            plan.append({'freq': freq,
                         'kind': 'existing',
                         'old_step_num': old_num})
        if not already_in_existing:
            plan.append({'freq': excitation_freq, 'kind': 'new'})
        plan.sort(key=lambda x: x['freq'])
        n_total_final = len(plan)
        for new_num, item in enumerate(plan, start=1):
            item['new_step_num'] = new_num

        # ---------------- early-exit messaging ---------------------------------
        if already_in_existing:
            print(f"Skipping: excitation frequency {excitation_freq:g} Hz is "
                  f"already present in '{output_file.name}'.  "
                  f"Nothing was written.")
            return

        # ---------------- plan + execute existing-step renames -----------------
        moves = [
            (item['old_step_num'], item['new_step_num'])
            for item in plan
            if item['kind'] == 'existing'
            and item['old_step_num'] != item['new_step_num']
        ]
        if moves:
            print(f"Reordering: renaming {len(moves)} existing step"
                  f"{'s' if len(moves) != 1 else ''} so the final Step list "
                  f"is monotonic in StepValue:")
            for old, new in moves:
                old_freq = next(v for k, v in existing_steps if k == old)
                print(f"  Step_{old} ({old_freq:g} Hz)  →  Step_{new}")
            for old, _ in moves:
                ms_out.move(f'Step_{old}', f'__tmp_Step_{old}')
            for old, new in moves:
                ms_out.move(f'__tmp_Step_{old}', f'Step_{new}')

        # ---------------- load + merge the spectrum across all input files ----
        new_item = next(item for item in plan if item['kind'] == 'new')
        new_step_idx = new_item['new_step_num']

        print(f"[2/2] Loading wavenumber spectrum from {len(input_files)} "
              f"file{'s' if len(input_files) != 1 else ''} for "
              f"f_exc = {excitation_freq:g} Hz …")
        all_regions = sorted({r for rn in result_names for r in result_regions[rn]})
        spectra = _load_combined_spectrum(input_files, c0, all_regions)

        total_k_per_region = {r: spectra[r][0].size for r in all_regions}
        print(f"  merged spectrum has "
              f"{', '.join(f'{n} k-samples for region {r}' for r, n in total_k_per_region.items())}")

        # Note about non-uniform sampling (still supported but unusual).
        for region, (k, _, _) in spectra.items():
            dk = np.diff(k)
            if dk.size > 0 and np.ptp(dk) > 1e-9 * np.abs(dk).max():
                print(f"  note: non-uniform wavenumber grid in merged "
                      f"spectrum for region '{region}' "
                      f"(Δk varies from {dk.min():.3e} to {dk.max():.3e}); "
                      f"using non-uniform trapezoidal weights.")
                break  # one note is enough

        # ---------------- write the new Step at its final slot -----------------
        step_name = f'Step_{new_step_idx}'
        if step_name in ms_out:
            raise RuntimeError(
                f"{step_name} already exists in '{output_file.name}' after "
                f"the rename phase — this is a bug.  Please report with a "
                f"copy of the file."
            )
        print(f"  writing {step_name} (f_exc = {excitation_freq:g} Hz) ...")
        step_out = ms_out.create_group(step_name)
        step_out.attrs.create('StepValue', np.float64(excitation_freq))

        for rn in result_names:
            for region in result_regions[rn]:
                k_merged, p_re, p_im = spectra[region]
                sol_re, sol_im = _trapezoidal_ift(
                    k_merged, p_re, p_im, z_in_storage_order
                )

                n_region_2d = region_nnodes_2d[(rn, region)]
                if sol_re.shape[1] != n_region_2d:
                    raise RuntimeError(
                        f"Region '{region}' has {sol_re.shape[1]} 2D "
                        f"result entries but the 2D mesh region has "
                        f"{n_region_2d} nodes."
                    )

                # Build the 3D pressure vector chunk-by-chunk in the
                # exact same order the extruder built the region's
                # /Mesh/Regions/<name>/Nodes list:
                #   full storage layer  -> all n_region_2d entries
                #   mid  storage layer  -> only the corner subset
                corner_mask = region_corner_in_rn2[(rn, region)]
                chunks_re: list[np.ndarray] = []
                chunks_im: list[np.ndarray] = []
                for pos in range(n_total_layers):
                    if layer_is_full[pos]:
                        chunks_re.append(sol_re[pos, :])
                        chunks_im.append(sol_im[pos, :])
                    else:
                        chunks_re.append(sol_re[pos, corner_mask])
                        chunks_im.append(sol_im[pos, corner_mask])
                full_re = np.concatenate(chunks_re)
                full_im = np.concatenate(chunks_im)

                step_out.create_dataset(
                    f'{rn}/{region}/Nodes/Real',
                    data=full_re.reshape((-1, 1)),
                    dtype=np.float64,
                )
                step_out.create_dataset(
                    f'{rn}/{region}/Nodes/Imag',
                    data=full_im.reshape((-1, 1)),
                    dtype=np.float64,
                )

        # ---------------- update MultiStep_1 metadata -------------------------
        # Always rebuild from the canonical combined-sorted list — this keeps
        # LastStepNum / LastStepValue / ResultDescription perfectly in sync
        # with whatever Step_i are now in the file.
        if n_total_final > 0:
            final_step_numbers = np.arange(1, n_total_final + 1, dtype=np.uint32)
            final_step_values  = np.array([item['freq'] for item in plan],
                                          dtype=np.float64)

            if 'LastStepNum' in ms_out.attrs:
                del ms_out.attrs['LastStepNum']
            _set_uint(ms_out, 'LastStepNum', n_total_final)
            if 'LastStepValue' in ms_out.attrs:
                del ms_out.attrs['LastStepValue']
            ms_out.attrs.create('LastStepValue',
                                np.float64(final_step_values[-1]))

            if 'ResultDescription' in ms_out:
                for rn in result_names:
                    rd = ms_out[f'ResultDescription/{rn}']
                    _rewrite_result_description_steps(
                        rd,
                        step_numbers=final_step_numbers,
                        step_values=final_step_values,
                    )

        # Stamp FileInfo.
        if 'FileInfo' in fout:
            fi = fout['FileInfo']
            for ds_name, value in [
                ('Date',    datetime.now().isoformat()),
                ('Creator', 'postprocess_2d5_to_3d.py'),
                ('Content', 'Grid+Results'),
            ]:
                if ds_name in fi:
                    del fi[ds_name]
                fi.create_dataset(ds_name, data=np.bytes_(value))

    # ----- summary ---------------------------------------------------------------
    print(f"Wrote {output_file}")
    if n_total_final > 0:
        all_freqs_str = (f"{plan[0]['freq']:g} … {plan[-1]['freq']:g}"
                         if n_total_final > 1
                         else f"{plan[0]['freq']:g}")
        print(f"  3D excitation frequencies (Hz) : "
              f"{all_freqs_str}  (added 1, total = {n_total_final})")
    print(f"  z layers (storage order, m)    : "
          f"{', '.join(f'{v:g}' for v in z_in_storage_order)}")
    print(f"  result regions                  : "
          f"{', '.join(result_regions['acouPressure'])}")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def _floats_csv(s: str) -> list[float]:
    """Parse a single string like '0.1,0.2,0.3' or '0.1 0.2 0.3' into floats.

    Used as the ``type=`` for argparse flags that take a list of numbers, so
    the parser cannot greedily consume the positional input file names.
    """
    parts = s.replace(',', ' ').split()
    try:
        return [float(x) for x in parts]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"could not parse '{s}' as a list of floats: {exc}"
        )


def _glue_negative_value_flags(argv: list[str],
                               flags: tuple[str, ...] = ('--z', '--frequency')
                               ) -> list[str]:
    """Collapse ``--z VAL`` into ``--z=VAL`` for the listed flags.

    Argparse refuses values that start with ``-`` and cannot be parsed as
    a single number (e.g. ``-0.5,-0.8,-1.2``), because they look like
    another flag.  We sidestep that by gluing the value onto the flag
    with ``=`` before argparse sees it.  The ``=`` form is always
    accepted, no matter what the value looks like.
    """
    out: list[str] = []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in flags and i + 1 < len(argv):
            out.append(f'{tok}={argv[i + 1]}')
            i += 2
        else:
            out.append(tok)
            i += 1
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser(
        description=("Build (or extend) a 3D harmonic .cfs file from one or "
                     "more 2.5D wavenumber-spectrum .cfs files, all "
                     "corresponding to a single 3D excitation frequency."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "What this script does\n"
            "  Each invocation performs ONE inverse Fourier transform at\n"
            "  ONE 3D excitation frequency and writes (or appends) ONE\n"
            "  Step in the output file.  The input files together hold\n"
            "  the wavenumber spectrum for that single excitation freq:\n"
            "    * 1 input file  → file holds the entire spectrum.\n"
            "    * N input files → each file holds a *partial* spectrum\n"
            "      (a slice of the k_z sweep).  This supports distributed\n"
            "      solver workflows where multiple subjobs each compute\n"
            "      part of the spectrum.  The script merges the slices\n"
            "      (sort by k_z, drop near-coincident duplicate samples)\n"
            "      and runs the IFT on the merged spectrum.\n"
            "\n"
            "Usage modes (auto-detected from whether OUTPUT exists):\n"
            "\n"
            "  CREATE mode (output file does not exist):\n"
            "    --z is REQUIRED.  A new 3D mesh is built by extruding the\n"
            "    2D mesh of the first input file along the user-supplied z\n"
            "    layers, then a fresh MultiStep_1 with one Step is written.\n"
            "\n"
            "  APPEND mode (output file already exists with a 3D mesh):\n"
            "    --z is OPTIONAL.  The mesh's z layers are reused.  The\n"
            "    new step is inserted into MultiStep_1 at the position\n"
            "    that keeps StepValues monotonic-ascending (existing\n"
            "    Steps may be renumbered).  If a step at the same\n"
            "    --frequency already exists, the call is a no-op.  If\n"
            "    --z is given but does not match the mesh's z layers,\n"
            "    a warning is printed and the mesh's layers are used.\n"
            "\n"
            "Argparse note:\n"
            "  --z and --frequency take a single value each (--z is a\n"
            "  comma- or space-separated list, --frequency is one float).\n"
            "  Either spelling works, even with a leading minus:\n"
            "      --z '-0.5,-0.8,-1.2'\n"
            "      --z=-0.5,-0.8,-1.2\n"
            "      --frequency=-12.5     (negative freqs in particular)\n"
            "\n"
            "Examples:\n"
            "  # Single input holds the full spectrum:\n"
            "  python3 %(prog)s -o DANI3D.cfs \\\n"
            "      --z '-0.5,-0.8,-1.2,-4,-6' --frequency 45.52 \\\n"
            "      DANI2d5_45p52Hz_full.cfs\n"
            "\n"
            "  # Three subjob outputs covering 0..10, 10.25..20, 20.25..30:\n"
            "  python3 %(prog)s -o DANI3D.cfs \\\n"
            "      --z '-0.5,-0.8,-1.2' --frequency 45.52 \\\n"
            "      DANI2d5_45p52Hz_part1.cfs \\\n"
            "      DANI2d5_45p52Hz_part2.cfs \\\n"
            "      DANI2d5_45p52Hz_part3.cfs\n"
            "\n"
            "  # Append another excitation frequency to an existing 3D file:\n"
            "  python3 %(prog)s -o DANI3D.cfs --frequency 112.0 \\\n"
            "      DANI2d5_112Hz_part1.cfs DANI2d5_112Hz_part2.cfs\n"
        ),
    )
    p.add_argument('-o', '--output', required=True,
                   help='output 3D .cfs file (created if missing, appended to '
                        'if it already exists with a 3D mesh).')
    p.add_argument('--z', type=_floats_csv, default=None,
                   metavar='"Z1,Z2,..."',
                   help='comma- or space-separated additional z-layers (m). '
                        'z=0 is included implicitly.  REQUIRED when OUTPUT '
                        'does not exist; optional (used as a consistency '
                        'check) when OUTPUT already has a 3D mesh.')
    p.add_argument('--frequency', type=float, required=True,
                   metavar='HZ',
                   help='the single 3D excitation frequency (Hz) of this '
                        'call (e.g. 45.52).  All input files must hold '
                        'wavenumber-spectrum slices for THIS frequency.  '
                        'REQUIRED.')
    p.add_argument('--rho', type=float, default=DEFAULT_RHO,
                   help=f'air density [kg/m^3] (default: {DEFAULT_RHO}).')
    p.add_argument('--K', type=float, default=DEFAULT_K,
                   help=f'air bulk modulus [Pa] (default: {DEFAULT_K}).')
    p.add_argument('inputs', nargs='+',
                   help='one or more 2.5D .cfs files holding the wavenumber '
                        'spectrum for --frequency (single full spectrum, or '
                        'multiple partial slices to be merged).')
    args = p.parse_args(_glue_negative_value_flags(sys.argv[1:]))

    postprocess_2d5_to_3d(
        input_files=args.inputs,
        output_file=args.output,
        excitation_freq=args.frequency,
        z_list=args.z,
        rho=args.rho,
        K=args.K,
    )