#!/usr/bin/env python3
"""h5info.py – inspect a .cfs HDF5 file produced by SimOutputHDF5.

   created by pure vibe coding (Claude Sonnet 4.6) as of 05.2026 by Fabian

Usage:
    python h5info.py result.cfs
    python h5info.py result.cfs --depth 3            # limit recursion depth

    # Extract a single string dataset as raw content (e.g. redirect to file):
    python h5info.py result.cfs --extract UserData/MaterialFile > mat.xml

    # Extract a group or numeric dataset as diff-friendly text:
    python h5info.py result.cfs --extract              # everything
    python h5info.py result.cfs --extract UserData
    python h5info.py result.cfs --extract Mesh
    python h5info.py result.cfs --extract Results/Mesh/MultiStep_1/Step_0/mechDisplacement/mech/Nodes/Real

    # Diff two files (Mesh + Results only):
    python h5info.py --diff a.cfs b.cfs
    python h5info.py --diff a.cfs b.cfs --n_diff 10
"""

import argparse
import sys

try:
    import h5py
    import numpy as np
except ImportError:
    sys.exit("Error: h5py and numpy are required.  Install with:  pip install h5py numpy")


def human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            if unit == "B":
                return f"{n:.0f} {unit}"
            prec = 0 if n >= 100 else (1 if n >= 10 else 2)
            return f"{n:.{prec}f} {unit}"
        n /= 1024
    return f"{n:.2f} PB"


def dtype_str(ds: "h5py.Dataset") -> str:
    return str(ds.dtype)


def chunk_info(ds: "h5py.Dataset") -> str:
    if ds.chunks is None:
        return "contiguous"
    # number of chunks per dimension
    num_chunks = tuple(
        int(np.ceil(s / c)) if s > 0 else 0
        for s, c in zip(ds.shape, ds.chunks)
    )
    chunk_bytes = int(np.prod(ds.chunks)) * ds.dtype.itemsize
    total_chunks = int(np.prod(num_chunks)) if ds.shape else 0
    return (
        f"chunked  chunk={ds.chunks}  {human_bytes(chunk_bytes)}/chunk"
        f"  n_chunks={num_chunks}  total={total_chunks}"
    )


def _single_chain_dataset(grp: "h5py.Group"):
    """If grp is a chain of single-item groups ending in a dataset, return
    (depth, dataset) where depth is the number of intermediate groups consumed.
    Otherwise return (0, None)."""
    depth = 0
    obj = grp
    while isinstance(obj, h5py.Group) and len(obj) == 1:
        key = next(iter(obj))
        obj = obj[key]
        depth += 1
    if isinstance(obj, h5py.Dataset):
        return depth, obj
    return 0, None


def visit(name: str, obj, args, prefix: str = "", current_depth: int = 0):
    """Recursively print groups and datasets."""
    if args.depth is not None and current_depth > args.depth:
        return

    indent = "  " * current_depth

    if isinstance(obj, h5py.Group):
        n_items = len(obj)

        # Collapse chain of single-item groups leading to one dataset
        chain_depth, ds = _single_chain_dataset(obj)
        if chain_depth >= 1 and ds is not None:
            shape = ds.shape
            nbytes = int(np.prod(shape)) * ds.dtype.itemsize if shape else ds.dtype.itemsize
            layout = chunk_info(ds)
            local_ds = ds.name.split("/")[-1]
            print(
                f"{indent}{chain_depth} * [GROUP] -> "
                f"[DATASET] {ds.name.lstrip('/')}  "
                f"shape={shape}  dtype={dtype_str(ds)}  "
                f"size={human_bytes(nbytes)}  {layout}"
            )
            return

        item_word = "item" if n_items == 1 else "items"
        print(f"{indent}[GROUP]  {name or '/'}  ({n_items} {item_word})")
        for key in obj:
            child = obj[key]
            child_name = f"{name}/{key}" if name else key
            visit(child_name, child, args, prefix, current_depth + 1)

    elif isinstance(obj, h5py.Dataset):
        shape = obj.shape
        # Variable-length string: itemsize is the pointer (8 B), not the content
        if h5py.check_string_dtype(obj.dtype):
            raw = obj[()]
            if isinstance(raw, np.ndarray):
                values = [v.decode(errors="replace") if isinstance(v, bytes) else str(v) for v in raw.flat]
            elif isinstance(raw, (bytes, np.bytes_)):
                values = [raw.decode(errors="replace")]
            else:
                values = [str(raw)]
            content_bytes = sum(len(v.encode()) for v in values)
            vl_str = h5py.check_string_dtype(obj.dtype)
            kind = "vl-string" if vl_str.length is None else f"fixed-string({vl_str.length})"
            print(
                f"{indent}[DATASET] {name}  "
                f"shape={shape}  dtype={kind}  "
                f"size={human_bytes(content_bytes)}"
            )
        else:
            nbytes = int(np.prod(shape)) * obj.dtype.itemsize if shape else obj.dtype.itemsize
            layout = chunk_info(obj)
            print(
                f"{indent}[DATASET] {name}  "
                f"shape={shape}  dtype={dtype_str(obj)}  "
                f"size={human_bytes(nbytes)}  {layout}"
            )


def _deflate_level(ds: "h5py.Dataset"):
    """Return deflate level for ds, or None if not compressed."""
    try:
        dcpl = ds.id.get_create_plist()
        for i in range(dcpl.get_nfilters()):
            code, _flags, data, _name = dcpl.get_filter(i)
            if code == h5py.h5z.FILTER_DEFLATE:
                return int(data[0]) if data else None
    except Exception:
        pass
    return None


def collect_stats(grp: "h5py.Group"):
    """Return stats tuple for all objects under grp."""
    total_bytes = 0
    ondisk_bytes = 0
    meta_bytes = 0
    total_chunks = 0
    total_chunk_bytes = 0
    n_datasets = 0
    n_groups = 0
    has_compression = False
    deflate_levels = set()

    def _walk(obj):
        nonlocal total_bytes, ondisk_bytes, meta_bytes, total_chunks, total_chunk_bytes
        nonlocal n_datasets, n_groups, has_compression
        if isinstance(obj, h5py.Dataset):
            n_datasets += 1
            shape = obj.shape
            if h5py.check_string_dtype(obj.dtype):
                raw = obj[()]
                if isinstance(raw, np.ndarray):
                    nbytes = sum(len((v.decode(errors="replace") if isinstance(v, bytes) else str(v)).encode()) for v in raw.flat)
                elif isinstance(raw, (bytes, np.bytes_)):
                    nbytes = len(raw)
                else:
                    nbytes = len(str(raw).encode())
            else:
                nbytes = int(np.prod(shape)) * obj.dtype.itemsize if shape else obj.dtype.itemsize
            total_bytes += nbytes
            # VL strings are never compressed; heap content is not in get_storage_size()
            # so use logical size as on-disk approximation for them
            if h5py.check_string_dtype(obj.dtype):
                ondisk_bytes += nbytes
            else:
                ondisk_bytes += obj.id.get_storage_size()
            try:
                info = h5py.h5o.get_info(obj.id)
                meta_bytes += info.hdr.space.total
            except Exception:
                pass
            level = _deflate_level(obj)
            if level is not None:
                has_compression = True
                deflate_levels.add(level)
            else:
                deflate_levels.add(None)
            if obj.chunks is not None:
                num_chunks = int(np.prod([
                    int(np.ceil(s / c)) if s > 0 else 0
                    for s, c in zip(shape, obj.chunks)
                ]))
                chunk_bytes = int(np.prod(obj.chunks)) * obj.dtype.itemsize
                total_chunks += num_chunks
                total_chunk_bytes += num_chunks * chunk_bytes
        elif isinstance(obj, h5py.Group):
            n_groups += 1
            try:
                info = h5py.h5o.get_info(obj.id)
                meta_bytes += info.hdr.space.total
            except Exception:
                pass
            for key in obj:
                _walk(obj[key])

    _walk(grp)
    # Summarize deflate level: single value if uniform, 'mixed' if not
    compressed_levels = deflate_levels - {None}
    if not compressed_levels:
        deflate_summary = None
    elif len(compressed_levels) == 1:
        deflate_summary = next(iter(compressed_levels))
    else:
        deflate_summary = "mixed"
    return total_bytes, ondisk_bytes, meta_bytes, total_chunks, total_chunk_bytes, n_datasets, n_groups, has_compression, deflate_summary


def print_stats(f: "h5py.File"):
    TOP_LEVEL = ["FileInfo", "Mesh", "Results", "UserData", "DataBase"]

    rows = []
    grand_bytes = 0
    grand_ondisk = 0
    grand_meta = 0
    grand_chunks = 0
    grand_chunk_bytes = 0
    grand_datasets = 0
    grand_groups = 0
    any_compression = False

    for section in TOP_LEVEL:
        if section not in f:
            continue
        sb, so, sm, sc, scb, nd, ng, hc, dlevel = collect_stats(f[section])
        grand_bytes += sb
        grand_ondisk += so
        grand_meta += sm
        grand_chunks += sc
        grand_chunk_bytes += scb
        grand_datasets += nd
        grand_groups += ng
        any_compression = any_compression or hc
        avg = (scb / sc) if sc > 0 else 0
        rows.append((section, sb, so, sm, sc, avg, nd, ng, hc, dlevel))

    # HDF5 structural overhead: B-tree nodes, symbol tables, superblock, name heaps
    import os
    file_size = os.path.getsize(f.filename)
    struct_overhead = file_size - grand_ondisk

    col_w = max(len(r[0]) for r in rows) + 2
    col_w = max(col_w, len("<Meta-Data>") + 2)
    show_ratio = any_compression

    print()
    print("=" * 105)
    print("STATISTICS")
    print("=" * 105)
    hdr = f"  {'Section':<{col_w}}  {'Logical':>10}  {'On-disk':>10}  {'Metadata':>10}  {'Datasets':>9}  {'Groups':>7}  {'Chunks':>8}  {'Avg chunk':>10}  {'Deflate':>8}"
    if show_ratio:
        hdr += f"  {'Ratio':>7}"
    print(hdr)
    print("-" * 105)
    for section, sb, so, sm, sc, avg, nd, ng, hc, dlevel in rows:
        avg_str = human_bytes(int(avg)) if sc > 0 else "n/a"
        dl_str = str(dlevel) if dlevel is not None else "-"
        line = f"  {section:<{col_w}}  {human_bytes(sb):>10}  {human_bytes(so):>10}  {human_bytes(sm):>10}  {nd:>9}  {ng:>7}  {sc:>8}  {avg_str:>10}  {dl_str:>8}"
        if show_ratio:
            ratio = so / sb if sb > 0 else 0
            line += f"  {ratio:>6.2f}x"
        print(line)
    # structural overhead row
    meta_line = f"  {'<Meta-Data>':<{col_w}}  {'':>10}  {human_bytes(struct_overhead):>10}  {'':>10}  {'':>9}  {'':>7}  {'':>8}  {'':>10}  {'':>8}"
    if show_ratio:
        meta_line += f"  {'':>7}"
    print(meta_line)
    print("-" * 105)
    grand_avg = (grand_chunk_bytes / grand_chunks) if grand_chunks > 0 else 0
    line = (f"  {'TOTAL':<{col_w}}  {human_bytes(grand_bytes):>10}  {human_bytes(file_size):>10}"
            f"  {human_bytes(grand_meta):>10}  {grand_datasets:>9}  {grand_groups:>7}"
            f"  {grand_chunks:>8}  {human_bytes(int(grand_avg)):>10}  {'':>8}")
    if show_ratio:
        ratio = file_size / grand_bytes if grand_bytes > 0 else 0
        line += f"  {ratio:>6.2f}x"
    print(line)
    print("=" * 105)


def _fmt_scalar(v) -> str:
    """Format a single scalar value for diff-friendly text output."""
    if isinstance(v, (float, np.floating)):
        return f"{v:.17g}"
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode(errors="replace")
    return str(v)


def _print_dataset(path: str, ds: "h5py.Dataset", force_header: bool = False):
    """Print a dataset in diff-friendly text form.

    If force_header is True (used during group traversal) a '### path' header
    is always printed, even for scalar string datasets, so the dataset name is
    visible.  Without force_header a scalar string dataset targeted directly is
    written as raw content only, which keeps  --extract UserData/MaterialFile
    suitable for redirection to a file.
    """
    raw = ds[()]

    # Scalar string: raw content only when directly targeted (force_header=False)
    if ds.shape == () and ds.dtype.kind in ("S", "O") and not force_header:
        content = raw.decode() if isinstance(raw, (bytes, np.bytes_)) else str(raw)
        sys.stdout.write(content)
        if not content.endswith("\n"):
            sys.stdout.write("\n")
        return

    arr = np.asarray(raw)
    print(f"### {path}")
    print(f"# shape={ds.shape}  dtype={ds.dtype}")

    # For scalar string with header: just print the content, not the raw bytes repr
    if ds.shape == () and ds.dtype.kind in ("S", "O"):
        content = raw.decode() if isinstance(raw, (bytes, np.bytes_)) else str(raw)
        print(content)
        print()
        return

    if arr.ndim == 0:
        print(f"  {_fmt_scalar(arr.item())}")
    elif arr.ndim == 1:
        for i, v in enumerate(arr.flat):
            print(f"  [{i}]  {_fmt_scalar(v)}")
    else:
        # Iterate over all leading dimensions; last dimension is on one line
        for idx in np.ndindex(arr.shape[:-1]):
            row = arr[idx]
            vals = "  ".join(_fmt_scalar(v) for v in row)
            idx_str = ",".join(str(i) for i in idx)
            print(f"  [{idx_str}]  {vals}")
    print()


def _count_leaves(obj) -> int:
    """Count the number of datasets reachable from obj."""
    if isinstance(obj, h5py.Dataset):
        return 1
    total = 0
    for key in obj:
        total += _count_leaves(obj[key])
    return total


def _visit_for_extract(obj, path: str):
    """Recursively visit an HDF5 object and print datasets in diff-friendly form."""
    if isinstance(obj, h5py.Dataset):
        _print_dataset(path, obj, force_header=True)
    elif isinstance(obj, h5py.Group):
        for key in sorted(obj.keys()):
            child_path = f"{path}/{key}" if path else key
            _visit_for_extract(obj[key], child_path)


def _print_extract_header(path: str, obj):
    """Print a compact statistics header for a multi-dataset extract."""
    sb, so, _sm, sc, _scb, nd, ng, hc, dlevel = collect_stats(obj)
    dl_str = f"  deflate={dlevel}" if dlevel is not None else ""
    print(f"# extract: {path or '/'}")
    print(f"# datasets={nd}  groups={ng}  logical={human_bytes(sb)}"
          f"  on-disk={human_bytes(so)}"
          f"  chunks={sc}{dl_str}")
    print()


def extract_path(f: "h5py.File", path: str):
    """Extract an arbitrary HDF5 path and print it to stdout.

    - Group → compact stats header, then all datasets (sorted, diff-friendly)
    - Scalar string dataset → raw content only (suitable for redirect to file)
    - Any other single dataset → compact stats header, then diff-friendly text
    """
    path = path.strip("/")
    try:
        obj = f["/" + path] if path else f["/"]
    except KeyError:
        sys.exit(f"Error: path '{path}' not found in this file.")

    # Scalar string scalar: raw only, no header (redirect-to-file use case)
    if isinstance(obj, h5py.Dataset) and obj.shape == () and obj.dtype.kind in ("S", "O"):
        _print_dataset(path, obj)
        return

    # Multiple leaves (group) or a single non-string dataset: show header first
    n_leaves = _count_leaves(obj)
    if n_leaves > 1 or isinstance(obj, h5py.Group):
        _print_extract_header(path, obj)

    _visit_for_extract(obj, path)


# ---------------------------------------------------------------------------
# Diff logic
# ---------------------------------------------------------------------------

DIFF_SECTIONS = ["Mesh", "Results"]


def _collect_datasets(grp: "h5py.Group", prefix: str = "") -> dict:
    """Return {relative_path: h5py.Dataset} for all datasets under grp."""
    result = {}
    for key in grp:
        child = grp[key]
        path = f"{prefix}/{key}" if prefix else key
        if isinstance(child, h5py.Dataset):
            result[path] = child
        elif isinstance(child, h5py.Group):
            result.update(_collect_datasets(child, path))
    return result


def _ds_meta(ds: "h5py.Dataset") -> str:
    nbytes = int(np.prod(ds.shape)) * ds.dtype.itemsize if ds.shape else ds.dtype.itemsize
    return f"shape={ds.shape}  dtype={ds.dtype}  size={human_bytes(nbytes)}"


def _diff_dataset(path: str, ds_a: "h5py.Dataset", ds_b: "h5py.Dataset",
                  n_diff: int, eps: float) -> list:
    """Compare two datasets; return list of human-readable issue strings."""
    issues = []

    if ds_a.shape != ds_b.shape:
        issues.append(f"  shape mismatch: {ds_a.shape} vs {ds_b.shape}")
        return issues
    if ds_a.dtype != ds_b.dtype:
        issues.append(f"  dtype mismatch: {ds_a.dtype} vs {ds_b.dtype}")

    # String datasets: compare raw bytes
    if ds_a.dtype.kind in ("S", "O"):
        if not np.array_equal(ds_a[()], ds_b[()]):
            issues.append("  content differs (string dataset)")
        return issues

    # Numeric datasets
    try:
        arr_a = np.asarray(ds_a[()])
        arr_b = np.asarray(ds_b[()])
    except Exception as e:
        issues.append(f"  could not read data: {e}")
        return issues

    # Use tolerance for floating-point, exact match for integer/bool
    if eps > 0 and np.issubdtype(arr_a.dtype, np.floating):
        diff_mask = np.abs(arr_a - arr_b) > eps
    else:
        diff_mask = arr_a != arr_b
    n_total = diff_mask.size
    n_bad = int(np.sum(diff_mask))
    if n_bad == 0:
        return issues

    tol_str = f"  (|a-b| > {eps})" if eps > 0 and np.issubdtype(arr_a.dtype, np.floating) else ""
    issues.append(f"  {n_bad}/{n_total} values differ{tol_str}")
    if n_diff > 0:
        shown = 0
        for flat_idx in range(n_total):
            if shown >= n_diff:
                break
            idx = np.unravel_index(flat_idx, arr_a.shape)
            va = arr_a[idx]
            vb = arr_b[idx]
            if diff_mask[idx]:
                idx_str = ",".join(str(i) for i in idx)
                issues.append(f"    [{idx_str}]  {_fmt_scalar(va)} -> {_fmt_scalar(vb)}")
                shown += 1
        remaining = n_bad - shown
        if remaining > 0:
            issues.append(f"    ... {remaining} more differing values")
    return issues


def run_diff(file_a: str, file_b: str, n_diff: int, eps: float):
    """Diff two .cfs files over Mesh and Results sections."""
    try:
        fa = h5py.File(file_a, "r")
    except OSError as e:
        sys.exit(f"Cannot open '{file_a}': {e}")
    try:
        fb = h5py.File(file_b, "r")
    except OSError as e:
        fa.close()
        sys.exit(f"Cannot open '{file_b}': {e}")

    print(f"A: {file_a}")
    print(f"B: {file_b}")

    any_section_diff = False

    with fa, fb:
        for section in DIFF_SECTIONS:
            in_a = section in fa
            in_b = section in fb

            if not in_a and not in_b:
                continue
            if not in_a:
                print(f"\n{section}: MISSING in A")
                any_section_diff = True
                continue
            if not in_b:
                print(f"\n{section}: MISSING in B")
                any_section_diff = True
                continue

            grp_a = fa[section]
            grp_b = fb[section]

            ds_a = _collect_datasets(grp_a)
            ds_b = _collect_datasets(grp_b)

            keys_a = set(ds_a)
            keys_b = set(ds_b)
            only_a = sorted(keys_a - keys_b)
            only_b = sorted(keys_b - keys_a)
            common = sorted(keys_a & keys_b)

            section_issues = []

            for k in only_a:
                section_issues.append((k, [f"only in A  ({_ds_meta(ds_a[k])})"] ))
            for k in only_b:
                section_issues.append((k, [f"only in B  ({_ds_meta(ds_b[k])})"] ))

            for k in common:
                issues = _diff_dataset(k, ds_a[k], ds_b[k], n_diff, eps)
                if issues:
                    section_issues.append((k, issues))

            if not section_issues:
                sb_a, *_ = collect_stats(grp_a)
                print(f"{section}: identical  ({len(common)} datasets, {human_bytes(sb_a)})")
            else:
                any_section_diff = True
                print(f"\n{section}: {len(section_issues)} dataset(s) differ")
                for k, msgs in sorted(section_issues):
                    # First message goes on the same line as the path
                    print(f"  {k}: {msgs[0].strip()}")
                    for m in msgs[1:]:
                        print(f"    {m.strip()}")

    print()
    print("identical" if not any_section_diff else "DIFFER")
    if any_section_diff:
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect or diff .cfs HDF5 files produced by SimOutputHDF5."
    )
    parser.add_argument("file", nargs="?", help="Path to the .cfs HDF5 file")
    parser.add_argument("--depth", type=int, default=None, metavar="N",
        help="Maximum depth to recurse into groups (default: unlimited)."
    )
    parser.add_argument("--stats", action="store_true", default=False,
        help="Print only the statistics table, skipping the tree view."
    )
    parser.add_argument(
        "--extract", metavar="PATH", nargs="?", const="",
        help="Extract an HDF5 path and print to stdout; suppresses all other "
             "output.  Omit PATH (or use '/') to extract everything.  "
             "A scalar string dataset (e.g. UserData/MaterialFile) is "
             "written as raw content so the result can be redirected to a file. "
             "A group or numeric dataset is printed in diff-friendly text form "
             "(one value/row per line, full precision, sorted paths)."
    )
    parser.add_argument(
        "--diff", metavar="FILE_B",
        help="Diff this file against FILE_B (Mesh and Results sections only)."
    )
    parser.add_argument(
        "--n_diff", type=int, default=5, metavar="N",
        help="Maximum differing values to show per dataset when using --diff "
             "(default: 5, 0 = counts only)."
    )
    parser.add_argument(
        "--eps", type=float, default=1e-12, metavar="EPS",
        help="Tolerance for floating-point comparisons when using --diff "
             "(default: 1e-12, 0 = exact)."
    )
    args = parser.parse_args()

    if args.diff:
        if not args.file:
            parser.error("the following argument is required: file (use: h5info.py a.cfs --diff b.cfs)")
        run_diff(args.file, args.diff, args.n_diff, args.eps)
        return

    if not args.file:
        parser.error("the following argument is required: file")

    try:
        f = h5py.File(args.file, "r")
    except OSError as e:
        sys.exit(f"Cannot open '{args.file}': {e}")

    with f:
        if args.extract is not None:
            extract_path(f, args.extract)
            return

        if args.stats:
            print_stats(f)
            return

        print(f"File: {args.file}")
        print(f"HDF5 library version: {h5py.version.hdf5_version}")
        print(f"h5py version:         {h5py.version.version}")
        print()
        visit("", f, args, current_depth=0)
        print_stats(f)


if __name__ == "__main__":
    main()
