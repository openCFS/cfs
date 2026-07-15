#!/usr/bin/env python
import argparse
import density
import cfs_utils
import subprocess
import pathlib
import os
import h5py
import numpy as np

def evaluate(inpath: pathlib.PurePath,
             outpath: pathlib.PurePath,
             fstart: int,
             fstop: int,
             fnum: int) -> None:
    xml = cfs_utils.open_xml(inpath)
    # if it does not specify multi frequency add it
    if fnum > 1:
        cfs_utils.replace(xml, '//cfs:costFunction/@multiple_excitation', "true")
    else:
        cfs_utils.replace(xml, '//cfs:costFunction/@multiple_excitation', "false")
    cfs_utils.replace(xml, '//cfs:optimizer/@type', "evaluate")
    cfs_utils.replace(xml, '//cfs:commit/@mode', "each_forward")
    cfs_utils.replace(xml, '//cfs:numFreq', fnum)
    cfs_utils.replace(xml, '//cfs:startFreq', fstart)
    cfs_utils.replace(xml, '//cfs:stopFreq', fstop)
    cfs_utils.replace(xml, '//cfs:sampling', "linear")
    xml.write(outpath)


def enable_anim(cfs_file: str | os.PathLike,
                fstart: int,
                fstop: int,):
    """This is a workaround, to enable animation of the "evaluate" results.
    We manually overwrite the attributes, to mimimic the output of HarmonicDriver, default export is TimeStep..."""
    with h5py.File(cfs_file, 'a') as f:
        ms = "Results/Mesh/MultiStep_1"
        dms = f[ms]
        dms.attrs["AnalysisType"] = "harmonic"
        dms.attrs["LastStepValue"] = fstop

        # change step values to frequencies
        freqs = None
        rd = ms + "/" + "ResultDescription"
        for r in f[rd]:
            sv = rd + "/" + r + "/" + "StepValues"
            dsv = f[sv]
            if freqs is None:
                freqs = np.linspace(fstart, fstop, len(dsv), endpoint=True)
                print(f"frequencies: {freqs}")
            # set frequencies as StepValues
            dsv[()] = freqs

        # extract maximum values
        dms = f[ms]
        fidx = 0
        for s in dms:
            if not s.startswith("Step_"):
                continue
            ds = f[ms + "/" + s]
            ds.attrs["StepValue"] = freqs[fidx]
            fidx += 1


def threshold(th=0.5, **kwargs):
    """Threshold function, with option to invert.

    Args:
        th: Threshold value, if < 0 the densities are inverted after thresholding with abs(th). Defaults to 0.5.

    Returns:
        Tuple[float, float]: Thresholded densities.
    """
    rho = kwargs["rho"]
    prho = kwargs["prho"]
    rho_out = 1
    if prho < abs(th):
        rho_out = 0
    if th < 0:
        rho_out = 1 - rho_out
    return rho_out, rho_out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="evaluate",
        description="Evaluate the optimized problem "
        "using the evaluate function.")

    parser.add_argument("name", help="Path to the cfs file.")
    parser.add_argument("-fl", "--flower", type=float,
                        help="Start frequency.")
    parser.add_argument("-fu", "--fupper", type=float,
                        help="Start frequency.")
    parser.add_argument("-fn", "--fnum", type=int,
                        help="Number of frequencies.")
    parser.add_argument("-t", default=1, type=int,
                        help="Threads to use for cfs simulations.")
    parser.add_argument("-th", nargs="?", type=float, const=0.5,
                        help="Enable threshold with the given value (default 0.5).")
    parser.add_argument("-m", "--mesh", type=str,
                        help="Name of the mesh file.")
    parser.add_argument("-p", "--param", type=str,
                        help="Name of XML parameter file, defaults to {name}.xml")
    parser.add_argument("-x", "--ersatz", type=str,
                        help="Name of ersatz material density file.")
    parser.add_argument("-e", "--executable", type=str, default="cfs",
                        help="Path to cfs executable, defaults to cfs")

    args = parser.parse_args()
    path = pathlib.Path(args.name)

    cmd = [args.executable, "-t", str(args.t)]

    # parse density file
    density_path = None
    if args.ersatz is not None:
        density_path = pathlib.Path(args.ersatz)
        print(f"Density path {density_path}")

    # threshold density file
    if args.th is not None:
        if density_path is None:
            raise ValueError("No density file given.")
        new_density_path = pathlib.Path(f"{density_path.parent}/{path.name}_eval.density.xml")
        print(new_density_path)
        density.mod_density(density_path,
                            new_density_path,
                            lambda **kwargs: threshold(th=args.th, **kwargs))
        density_path = new_density_path
        print(f"Updated density path after thresholding {density_path}")

    if density_path is not None:
        cmd += ["-x", density_path.absolute()]

    # copy/create eval simulation input file
    param_path = path.with_name(f"{path.name}_eval").with_suffix(".xml")
    param_path.parent.mkdir(parents=True, exist_ok=True)  # create folder if it doesnt exist
    if args.param is not None:
        # copy and edit
        in_path = pathlib.Path(args.param)
        evaluate(in_path.with_suffix(".xml"),
                 param_path,
                 args.flower, args.fupper, args.fnum)
    else:
        # if not given create form optimization file "{name}_eval.xml"
        evaluate(path.with_suffix(".xml"),
                 param_path,
                 args.flower, args.fupper, args.fnum)
    print(f"Parameter path {param_path}")
    cmd += ["-p", param_path.absolute(), f"{path.name}_eval"]

    # parse mesh
    if args.mesh:
        mesh_path = pathlib.Path(args.mesh)
        cmd += ["-m", mesh_path.absolute()]

    # run simulation
    print(f"Starting simulation run in cwd {path.parent}")
    path.parent.mkdir(parents=True, exist_ok=True)  # create folder if it doesnt exist
    subprocess.run(cmd,
                   cwd=path.parent,
                   check=True)

    # change hdf5 so it can be animated
    enable_anim(f"{path.parent}/{path.name}_eval.cfs", args.flower, args.fupper)
