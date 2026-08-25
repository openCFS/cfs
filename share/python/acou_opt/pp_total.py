#!/usr/bin/env python
from typing import List, Literal, Dict
import pathlib
import argparse

IDLE = 0
START_PARSE = 1
FINISHED_PARSE = 2


def compute_total(path: pathlib.Path,
                  cols: List[str]) -> dict | Literal[False]:
    """ Computest total for given column names."""
    # dict used to store the results
    res: Dict[str, float | str] = {"name": path.name}
    for cn in cols:
        res[cn] = float(0)

    cdxs = [0]*len(cols)

    start_parse = False
    with open(path, "r", encoding="utf8") as f:
        for line in f:
            # reset sum on new headline
            if line.startswith("# (1)"):
                start_parse = True
                continue

            if not start_parse:
                continue

            # end of postproc section, return sum
            if line.strip() == "":
                return res

            # extract column number of wallskip comments
            if line.startswith("#"):
                # only parse cols if all cols are present
                for cn in cols:
                    if cn not in line:
                        continue

                # fill cdxs, the mapping array that
                # maps names to column numbers
                for idx, name in enumerate(line[1:].split()):
                    try:
                        cdxs[cols.index(name)] = idx
                    except ValueError:
                        pass
                continue

            for cn, cidx in zip(cols, cdxs):
                tot = res[cn]
                assert isinstance(tot, float)
                tot += float(line.split()[cidx])
                res[cn] = tot
    return False


def print_res(res: dict, header: bool = False) -> None:
    if header:
        print("#", end="")
        for key in res.keys():
            print(f"{key:<20}", end="")
        print()
    for val in res.values():
        if isinstance(val, float):
            print(f"{val:20.1f}", end="")
        else:
            print(f"{val:<20}", end="")
    print("")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="postproc_total",
        description="Compute the total of given postproc files.")
    parser.add_argument("input", nargs='*',
                        help="Selection of files (containing postproc.py output) to process (with wildcards)")
    parser.add_argument("--colnames", nargs='*',
                        help="Names of columns to sum.")
    parser.add_argument("--wallperiter", action="store_true",
                        help="Computes the wall time per iteration (requires wall and iter colnames).")
    args = parser.parse_args()

    # no colnames provided
    if not args.colnames:
        parser.print_usage()
        exit(1)

    header = True
    for file in args.input:
        res = compute_total(pathlib.Path(file), args.colnames)
        if args.wallperiter:
            try:
                res["wall/iter"] = res["wall"] / res["iter"]
            except KeyError:
                raise IOError("--walliter requires the --colnames \"wall\" and \"iter\".")
        if res:
            print_res(res, header)
            header = False
        else:
            print(f"Failed {file}")
