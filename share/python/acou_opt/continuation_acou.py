#!/usr/bin/env python
import argparse
import os
import math
import pathlib
import subprocess
from cfs_utils import open_xml, replace, xpath
from typing import Optional, List

DETA = 0.3  # distance from eta = 0.5


def get_name(name: str | os.PathLike, value: float, final: bool = False):
    if final:
        # final has only one z to make it top
        return f"z{name}"
    else:
        if value >= 1:
            value = round(value)
        else:
            exp = math.floor(math.log10(value))
            value = round(value, abs(exp))
        return f"zz{name}_s{value}".replace('.', '')


def relative_paths(xml):
    mesh_path = xpath(xml, "//cfs:mesh/@fileName")
    # if path is given relative from home, there is probably a reason, else replace
    if not mesh_path.startswith("~"):
        replace(xml, "//cfs:mesh/@fileName", f"../{mesh_path}")  
    mat_path = xpath(xml, "//cfs:materialData/@file")
    # if path is given relative from home, there is probably a reason, else replace
    if not mat_path.startswith("~"):
        replace(xml, "//cfs:materialData/@file", f"../{mat_path}")


def cont_robust(name: str | os.PathLike,
                folder: str,
                cfs_threads: str,
                ersatz: Optional[str],
                steps: List[float]):
    print(f"STEPS = {steps}")
    detas = exp_steps(1e-3, DETA, len(steps))
    print(f"ETAS = {detas}")

    densities = ""
    for istep, (step, deta) in enumerate(zip(steps, detas)):
        xml = open_xml(f"{name}.xml")
        # switch to absolute paths
        relative_paths(xml)
        search_mode = True
        # test for robust
        if search_mode:
            try:
                replace(xml, "//cfs:filter[@robust_excitation=0]//cfs:density/@beta", step)
                # replace(xml, "//cfs:filter[@robust_excitation=0]//cfs:density/@eta", 0.5 - deta)

                replace(xml, "//cfs:filter[@robust_excitation=1]//cfs:density/@beta", step)
                # replace(xml, "//cfs:filter[@robust_excitation=1]//cfs:density/@eta", 0.5 + deta)

                replace(xml, "//cfs:filter[@robust_excitation=2]//cfs:density/@beta", step)
                # replace(xml, "//cfs:filter[@robust_excitation=2]//cfs:density/@eta", 0.5)
                print("Detected mode ROBUST.")
                search_mode = False
            except RuntimeError:
                pass
        # test for shapemap
        if search_mode:
            try:
                replace(xml, "//cfs:shapeMap/@beta", step)
                print("Detected mode SHAPEMAP.")
                search_mode = False
            except RuntimeError:
                pass
        # test for shapemap
        if search_mode:
            try:
                replace(xml, "//cfs:featureMapping/@transition", step)
                print("Detected mode FEATUREMAP.")
                search_mode = False
            except RuntimeError:
                pass
        # else: normal simp mode
        if search_mode:
            try:
                replace(xml, "//cfs:filter//cfs:density/@beta", step)
                print("Detected mode SIMP.")
                search_mode = False
            except RuntimeError:
                pass
        # handle no mode found
        if search_mode:
            raise ImportError("Could not determine the continuation mode, "
                              "please verify if your input xml file is supported.")

        file = get_name(name, step)
        if istep >= len(steps) - 1:
            file = get_name(name, step, True)

        xml.write(f"{folder}/{file}.xml")
        # store name of density file
        densities += f"{file}.density.xml "

        cmd = ["cfs", "-t", cfs_threads, file]
        if istep != 0:
            density = f"{get_name(name, steps[istep - 1])}.density.xml"
            cmd += ["-x", density]
        elif ersatz:
            cmd += ["-x", ersatz]

        subprocess.run(cmd, cwd=folder, check=True)


def exp_steps(start: float, stop: float, num: int) -> List[float]:
    if num <= 1:
        return [start]
    e = (stop / start) ** (1 / (num - 1))
    return [start * (e ** i) for i in range(num)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="continuation_acou.py",
        description="Similar to continuation.py, "
        "but a folder for the subproblems can be specified. "
        "Furthermore, robust (beta + eta) and shape-/feature mapping are supported."
    )
    parser.add_argument("name",
                        help="Path to the cfs file (without file ending)")
    parser.add_argument('-o', '--folder', nargs="?", default=False, const=True,
                        help="Run the continuation in subfolder.")
    parser.add_argument("-t", default=1, type=int,
                        help="Threads to use for cfs simulations.")
    parser.add_argument("-x", "--ersatz", type=str,
                        help="Name of ersatz material density file.")
    parser.add_argument('--factor', nargs=3, type=float, required=True,
                        help="Give START STOP NUM for exponential series. "
                        "For robust the values are mapped to eta."
                        "If NUM = 1 the START value is used for a single simulation.")

    args = parser.parse_args()
    path = pathlib.Path(args.name)

    # set folder
    folder = "."
    if args.folder and isinstance(args.folder, bool):
        folder = f"./z{path.with_suffix("")}"
    if isinstance(args.folder, str):
        folder = args.folder
    os.makedirs(folder, exist_ok=True)

    # generate steps
    start, stop, num = args.factor
    steps = exp_steps(start, stop, int(num))

    cont_robust(path.with_suffix(""),
                folder,
                str(args.t),
                args.ersatz,
                steps)
