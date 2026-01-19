#!/usr/bin/env python
import cfs_utils
import argparse
import pathlib
import os


def threshold(density: str | os.PathLike,
              old_mesh: str | os.PathLike,
              new_mesh: str | os.PathLike,
              th: float):
    xml = cfs_utils.open_xml(density)
    set = xml.find("set")
    S_opt = {}
    for elem in set.findall("element"):
        S_opt[int(elem.get("nr"))] = float(elem.get("physical"))

    with open(old_mesh, "r", encoding="utf-8") as mesh, open(new_mesh, 'w', encoding='utf-8') as thmesh:
        for line in mesh:
            if "S_opt" in line:
                e_nr = int(line.split()[0])

                if S_opt[e_nr] < th:
                    thmesh.write(line.replace("S_opt", "S_void"))
                    continue

            thmesh.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="threshold.py",
        description="Remove elements from given mesh file "
        "that have a density below a threshold in the given density file.")

    parser.add_argument("density", help="path to the .density.xml file")
    parser.add_argument("mesh", help="path to the original mesh file")
    parser.add_argument("newmesh", help="path to the new mesh file")
    parser.add_argument("-th", "--threshold", default=0.5, type=float,
                        help="threshold value", dest="th")

    args = parser.parse_args()
    density_path = pathlib.Path(args.density)
    mesh_path = pathlib.Path(args.mesh)

    # check input mesh file
    with open(mesh_path, 'r') as file:
        if file.readline().strip() != "[Info]":
            raise ValueError("Only ANSYS mesh type files supported!")

    # create threshholded mesh
    threshold(f"{str(density_path).split(".")[0]}.density.xml",
              mesh_path,
              args.newmesh,
              args.th)
