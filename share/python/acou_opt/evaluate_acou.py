#!/usr/bin/env python
import argparse
import cfs_utils
import subprocess
import pathlib
from typing import Optional
import os

from optools_density import mod_density
from optools_param import apply_evaluate, apply_step, apply_em
from optools_cfs import enable_anim


def evaluate(inpath: str | os.PathLike,
             outpath: str | os.PathLike,
             fstart: float,
             fstop: float,
             fnum: int,
             step: Optional[float]) -> None:
    xml = cfs_utils.open_xml(inpath)
    apply_evaluate(xml, fstart, fstop, fnum)
    if step is not None:
        if apply_step(xml, step):
            # handle no beta found
            raise ImportError("You specified beta, but we could not find it in the xml.")
    xml.write(outpath)


def add_em(param_path: str | os.PathLike,
           density_path: str | os.PathLike):
    xml = cfs_utils.open_xml(param_path)
    apply_em(xml, str(density_path), access="physical")
    xml.write(param_path)


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
    parser.add_argument("-beta", type=float,
                        help="Set projection beta to given value.")
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
        mod_density(density_path,
                    new_density_path,
                    lambda **kwargs: threshold(th=args.th, **kwargs))
        density_path = new_density_path
        print(f"Updated density path after thresholding {density_path}")

    # copy/create eval simulation input file
    param_path = path.with_name(f"{path.name}_eval").with_suffix(".xml")
    param_path.parent.mkdir(parents=True, exist_ok=True)  # create folder if it doesnt exist
    if args.param is not None:
        # copy and edit
        in_path = pathlib.Path(args.param)
        evaluate(in_path.with_suffix(".xml"),
                 param_path,
                 args.flower, args.fupper, args.fnum,
                 args.beta)
    else:
        # if not given create form optimization file "{name}_eval.xml"
        evaluate(path.with_suffix(".xml"),
                 param_path,
                 args.flower, args.fupper, args.fnum,
                 args.beta)
    # add density as loadErsatzMaterial so we can load the physical design
    if density_path is not None:
        add_em(param_path, density_path.absolute())
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
