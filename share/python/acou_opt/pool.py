#!/usr/bin/env python
import argparse
import os
import sys
import signal
import json
import subprocess
from typing import List, TextIO, Dict
import shlex
import time
import fnmatch
from dataclasses import dataclass
import pathlib
import shutil


@dataclass
class Process():
    p: subprocess.Popen
    log: TextIO
    name: str
    max_stage: int
    stage: int = 0

    def __repr__(self) -> str:
        return f"{self.name}[{self.stage}]"


def str2cmd(s: str) -> List[str]:
    return shlex.split(s)


def run_pool(jobs: Dict[str, List[str]],
             names: List[str],
             max_processes: int,
             skip_stages: int = 0,
             error_cont: bool = False) -> None:
    idx = 0
    processes: List[Process] = []

    try:
        while True:
            # check all running processes, advance them or remove finished ones
            for proc in processes:
                if proc.p.poll() is not None:
                    es = proc.p.returncode
                    # advance stage
                    if proc.stage < proc.max_stage and (es == 0 or error_cont):
                        print(f"{names.index(proc.name) + 1}:\tadvancing:\t-> ", end="")
                        proc.stage += 1
                        p = subprocess.Popen(str2cmd(jobs[proc.name][proc.stage]),
                                             stdout=proc.log,
                                             stderr=proc.log)
                        proc.p = p
                        print(f"{proc} (EXIT = {es})", flush=True)
                    else:
                        if (proc.p.returncode == 0 or error_cont):
                            print(f"{names.index(proc.name) + 1}:\tfinished:\t{proc} (EXIT = {es})", flush=True)
                        else:
                            print(f"{names.index(proc.name) + 1}:\texiterror:\t{proc} (EXIT = {es})", flush=True)
                        proc.log.close()
                        processes.remove(proc)
                        continue

            # start processes if queue has place
            while len(processes) < max_processes and idx < len(names):
                logfile = open(f"{names[idx]}.log", "w", encoding="utf-8")
                p = subprocess.Popen(str2cmd(jobs[names[idx]][skip_stages]),
                                     stdout=logfile,
                                     stderr=logfile)

                proc = Process(p, logfile, names[idx], len(jobs[names[idx]]) - 1, stage=skip_stages)
                processes.append(proc)

                print(f"{idx + 1}/{len(names)}:\tstarting:\t{proc}", flush=True)
                idx += 1

            # if all are finished return function
            if idx + 1 >= len(names) and len(processes) == 0:
                print("DONE")
                return
            time.sleep(1)

    except KeyboardInterrupt:
        while len(processes) > 0:
            proc = processes[0]
            print(f"{len(processes)}: SIGINT -> {proc}", flush=True)
            # close subprocess gently with SIGINT
            os.kill(proc.p.pid, signal.SIGINT)
            proc.p.wait()
            proc.log.close()
            processes.remove(proc)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="pool.py",
        description="Execute mutliple jobs in parallel, "
        "while respecting the maximum number of resources. "
        "Jobs are loaded from a json file, also consecutive jobs "
        "are supported!",
        epilog="""
Only exit pool.py early with KeyboardInterrupt (Ctrl+C).
In nuhup mode use kill -2 PID

Example input JSON file:
{
"twoecho": ["echo 'Hello World'", "sleep 5", "echo 'Nap successful'"],
"erroron2nd": ["echo Hi", "false", "echo 'Dont reach!'"]
}
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input",
                        help="Input json file defining the job pool.")
    parser.add_argument("-e", "--error_cont", action="store_true",
                        help="dont stop on EXIT > 0")
    parser.add_argument("-m", "--maxprocesses", nargs=1, default=[2],
                        type=int, help="maximum number of processes")
    parser.add_argument("-n", "--name", nargs='*',
                        type=str, help="only run the jobs with given name(s)")
    parser.add_argument("-s", "--skip", default=0,
                        type=int, help="skip first n stages")
    parser.add_argument("--clean", action="store_true",
                        help="clean log files produced by pool.py")

    args = parser.parse_args()

    input_path = pathlib.Path(args.input)
    with open(input_path, "r") as f:
        data = json.load(f)

    names = list(data.keys())
    if args.name is not None:
        names = [name for name in names
                 if any(fnmatch.fnmatch(name, f) for f in args.name)]

    if not names:
        raise ValueError("No tasks found.")

    if args.clean:
        for name in names:
            # delete log file
            pathlib.Path(f"{input_path.parent}/{name}.log").unlink()
            # delete results folder
            shutil.rmtree(f"{input_path.parent}/{name}")
    else:
        run_pool(data, names, args.maxprocesses[0], args.skip, args.error_cont)
