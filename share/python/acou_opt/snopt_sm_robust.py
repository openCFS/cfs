#!/usr/bin/env python
import pyoptsparse
from typing import Optional, List, get_args
import numpy as np
import lxml
import lxml.etree
import subprocess
import os
import signal
import argparse

from dataclasses import dataclass, fields, field


def main():
    parser = argparse.ArgumentParser(prog='cfs_sm_robust.py',
                                     description='CFS-like executable to run robust shape-mapping.')

    parser.add_argument('name')
    parser.add_argument('-t', '--numThreads', type=int, default=1)
    parser.add_argument('-m', '--meshFile', type=str, default=None)
    parser.add_argument('-p', '--paramFile', type=str, default=None)
    parser.add_argument('-x', '--densityFile', type=str, default=None)
    parser.add_argument('-s', '--shift', type=float, default=0.001,
                        help="Value of symmetrically applied robust shift.")

    args = parser.parse_args()

    solve_robust(args.name,
                 t=args.numThreads,
                 m=args.meshFile,
                 p=args.paramFile,
                 x=args.densityFile,
                 shift=args.shift)


def solve_robust(name,
                 t: int = 1,
                 m: Optional[str] = None,
                 p: Optional[str] = None,
                 x: Optional[str] = None,
                 shift: float = 0.001):

    if p is None:
        p = f"{name}.xml"
    elif not p.endswith(".xml"):
        raise IOError("Specify p with .xml ending!")
    # init cfs (runs once to get .gradplot file)
    cfse = CFS(f"{name}_erode", t=t, m=m, p=p, x=x)
    cfs = CFS(name, t=t, m=m, p=p, x=x)
    cfsd = CFS(f"{name}_dilate", t=t, m=m, p=p, x=x)

    # we assume slack is the last variable
    assert cfs.variables.type[-1] == "slack"

    # create offset array for profiles used for robust
    offset = [shift if vt == "profile" else 0 for vt in cfs.variables.type]

    c = cfs.metadata.n_glob_constraints
    n = cfs.metadata.n

    # pyoptsparse function definitions ---------------------------------------
    def objfunc(xdict):
        x = xdict["xvars"]

        # compute
        cfse.compute(x - offset)
        cfs.compute(x)
        cfsd.compute(x + offset)

        funcs = {}
        funcs["obj"] = x[-1]  # slack
        # nonlinear (global) contraints
        nlconval = (cfse.constraints.value[:c] +
                    cfsd.constraints.value[:c] +
                    cfs.constraints.value[:c])
        funcs["nlcon"] = nlconval
        # linear (local) constraints
        lconval = cfs.constraints.value[c:]
        funcs["lcon"] = lconval
        fail = False

        return funcs, fail

    def sens(xdict, funcs):
        # objective gradient
        grad = np.zeros(len(cfs.variables.vid))
        grad[-1] = 1  # slack

        # Return dictionary
        return {"obj": {"xvars": grad},
                "nlcon": {"xvars": get_jacsparse("nlcon")},
                "lcon": {"xvars": get_jacsparse("lcon")}}

    def get_jacsparse(congroup=""):
        # get data
        if congroup == "nlcon":
            cid = (cfse.jacobian.cid[:c*n] +
                   [cid + c for cid in cfsd.jacobian.cid[:c*n]] +
                   [cid + 2*c for cid in cfs.jacobian.cid[:c*n]])
            vid = (cfse.jacobian.vid[:c*n] +
                   cfsd.jacobian.vid[:c*n] +
                   cfs.jacobian.vid[:c*n])
            jac = (cfse.jacobian.gradient[:c*n] +
                   cfsd.jacobian.gradient[:c*n] +
                   cfs.jacobian.gradient[:c*n])
            shape = (3 * c, n)
        elif congroup == "lcon":
            # get data
            cid = [cid - c for cid in cfs.jacobian.cid[c*n:]]
            vid = cfs.jacobian.vid[c*n:]
            jac = cfs.jacobian.gradient[c*n:]
            shape = (cfs.metadata.m - c, n)
        else:
            raise NotImplementedError(f"Value congroup={congroup} not implemented!")

        # print(f"{congroup}: [{min(cid), max(cid)}]:{len(cid)}, [{min(vid), max(vid)}]:{len(vid)}, {len(jac)}")
        # pack into sparse coo format
        return {"coo": np.array([cid, vid, jac]),
                "shape": shape}

    # init pyoptsparse
    optProb = pyoptsparse.Optimization("Robust", objfunc)
    optProb.addVarGroup("xvars",
                        cfs.metadata.n,
                        lower=[lb + o for lb, o in zip(cfs.variables.lb, offset)],
                        upper=[ub - o for ub, o in zip(cfs.variables.ub, offset)])
    # here we pass the full gradient, in sens we overwrite only the non linear
    optProb.addConGroup("nlcon",
                        3 * c,
                        lower=cfse.constraints.cl[:c] + cfsd.constraints.cl[:c] + cfs.constraints.cl[:c],
                        upper=cfse.constraints.cu[:c] + cfsd.constraints.cu[:c] + cfs.constraints.cu[:c],
                        wrt=["xvars"],
                        linear=False,
                        jac={"xvars": get_jacsparse("nlcon")})
    optProb.addConGroup("lcon",
                        cfs.metadata.m - c,
                        lower=cfs.constraints.cl[c:],
                        upper=cfs.constraints.cu[c:],
                        wrt=["xvars"],
                        linear=True,
                        jac={"xvars": get_jacsparse("lcon")})
    optProb.addObj("obj")

    # snopt settings
    options = {
        "Major iterations limit": 500,
        "Iterations limit": 1000000,
        "Minor iterations limit": 100000,
        "Verify level": -1,
        "iSumm": 6,
        "Print file": name + ".snopt"
    }
    opt = pyoptsparse.SNOPT(options=options)
    sol = opt(optProb, sens=sens)


# ====================== DATACLASS ======================
@dataclass
class Metadata:
    iteration: int = -1
    n_vars: int = -1
    n_constraints: int = -1
    n_glob_constraints: int = -1
    n_jac_nonzero: int = -1

    _name = "metadata"
    _active = False

    def parse_line(self, line) -> bool:
        # for empty line return false
        if not line:
            self._active = False
            return False

        # check heading (we enter this classes section)
        if self._name in line.lower():
            self._active = True
            return True
        if not self._active:
            return False

        # parse data
        parts = line.split()
        for f in fields(self):
            # skip private fields
            if f.name.startswith("_"):
                continue
            if f.name in line:
                setattr(self, f.name, int(parts[-1]))
                return True
        return False

    @property
    def n(self) -> int:
        if self.n_vars is not None:
            return self.n_vars
        else:
            raise ValueError("Metadata not initialized")

    @property
    def m(self) -> int:
        if self.n_constraints is not None:
            return self.n_constraints
        else:
            raise ValueError("Metadata not initialized")


@dataclass
class TypeConv:

    _name = "Override"
    _active = False
    _header_names = None

    def parse_line(self, line) -> bool:
        # for empty line return false
        if not line:
            # if it was active, we now exit this dataclass
            if self._active:
                # self.convert_to_numpy()
                self._active = False
            return False

        # check heading (we enter this classes section)
        if self._name in line.lower():
            self._active = True
            return True
        if not self._active:
            return False

        # parse header
        if line.startswith("#"):
            self._header_names = line[1:].split()
            return True
        assert self._header_names is not None

        # parse data
        parts = line.split()
        for h, v in zip(self._header_names, parts):
            self.append_by_name(h, v)
        return True

    def append_by_name(self, name: str, value: str):
        # get the list to append to
        lst = getattr(self, name)
        # get the inner type
        t = None
        for f in fields(self):
            if f.name == name:
                t = get_args(f.type)[0]
                break
        assert t is not None
        # append to list
        lst.append(t(value))


@dataclass
class Variables(TypeConv):
    vid: List[int] = field(default_factory=list)
    nr: List[int] = field(default_factory=list)
    type: List[str] = field(default_factory=list)
    design: List[float] = field(default_factory=list)
    lower_bound: List[float] = field(default_factory=list)
    upper_bound: List[float] = field(default_factory=list)
    # only for shape variables
    dof: List[int] = field(default_factory=list)
    shape: List[int] = field(default_factory=list)
    ref: List[int] = field(default_factory=list)

    _name = "variables"

    @property
    def lb(self) -> List[float]:
        return self.lower_bound

    @property
    def ub(self) -> List[float]:
        return self.upper_bound


@dataclass
class Constraints(TypeConv):
    cid: List[int] = field(default_factory=list)
    type: List[str] = field(default_factory=list)
    value: List[float] = field(default_factory=list)
    bound_type: List[str] = field(default_factory=list)
    bound: List[float] = field(default_factory=list)
    local: List[int] = field(default_factory=list)

    _name = "constraint values"

    @property
    def cl(self) -> List[float]:
        cls = [-2e19] * len(self.cid)
        for idx, (bt, b) in enumerate(zip(self.bound_type, self.bound)):
            if bt == "lowerBound":
                cls[idx] = b
        return cls

    @property
    def cu(self) -> List[float]:
        cus = [2e19] * len(self.cid)
        for idx, (bt, b) in enumerate(zip(self.bound_type, self.bound)):
            if bt == "upperBound":
                cus[idx] = b
        return cus


@dataclass
class Jacobian(TypeConv):
    cid: List[int] = field(default_factory=list)
    vid: List[int] = field(default_factory=list)
    gradient: List[float] = field(default_factory=list)

    _name = "jacobian sparse"


# ====================== CFS ======================
class CFS():
    def __init__(self, name,
                 t: int = 1,
                 m: Optional[str] = None,
                 p: Optional[str] = None,
                 x: Optional[str] = None):
        self.name = name  # name
        self._m = m  # cfs -m meshFile
        self._p = p  # cfs -p paramFile
        self._t = t  # cfs -t numThreads
        self._x = x  # cfs -x densFile (only used for first iteration to init)
        self.compute()  # init dataclasses, compute and read initial values
        # we don't need to init super because we call it in compute()

    def _param_file(self):
        if self._p is None:
            return self.name + ".xml"
        else:
            return self._p

    def _dens_file(self):
        return self.name + ".density.xml"

    def _grad_file(self):
        return self.name + ".grad.dat"

    def compute(self, x: Optional[np.ndarray] = None):
        cmd = ["cfs", "-p", self._param_file()]
        if x is None:
            # init
            if self._x is not None:
                cmd += ["-x", self._x]
        else:
            # write density file
            assert len(self.variables.design) == len(x)
            self.variables.design = list(x).copy()
            self._write_density()
            cmd += ["-x", self._dens_file()]
        if self._m is not None:
            cmd += ["-m", self._m]
        if self._t > 1:
            cmd += ["-t", str(self._t)]
        # run cfs, catch SIGINT and close cfs softly
        cmd += [self.name]
        with open(f"{self.name}.log", "w") as logfile:
            p = None
            try:
                p = subprocess.Popen(cmd, stderr=logfile, stdout=logfile)
                p.wait()
                if p.returncode:
                    raise subprocess.CalledProcessError(p.returncode, cmd)
            except Exception as e:
                if p:
                    os.kill(p.pid, signal.SIGINT)
                    p.wait()
                raise e
        # read new values
        self._read_gradplot()

    def _init_dataclasses(self):
        self.metadata = Metadata()
        self.variables = Variables()
        self.constraints = Constraints()
        self.jacobian = Jacobian()

    def _write_density(self):
        parser = lxml.etree.XMLParser(remove_comments=True)
        tree = lxml.etree.parse(self._dens_file(), parser)
        root = tree.getroot()

        # parse header
        header = root.find("header")
        if header is None:
            raise IOError("Missing '<header>' element, are you sure you loaded a .density.xml file?")

        # parse sets
        sets = root.findall("set")
        if header is None:
            raise IOError("Missing '<set>' element, there is no density data in your file.")
        assert len(sets) == 1
        set = sets[0]

        # get list of shapeParamElements
        vars = self.variables
        spes = set.findall("shapeParamElement")
        # print(f"Found {len(spes)} shapeParamElements")
        # change shapeParamElement data
        for spe in spes:
            # we set by the nr
            nr = int(spe.get("nr"))
            try:
                for vnr, vtype, vdesign in zip(vars.nr, vars.type, vars.design):
                    if vnr == nr and vtype != "slack":
                        spe.set("design", str(vdesign))
                        break
            except ValueError:
                # these elements are set by cfs (e.g. symmetry)
                pass

        # set slack
        slacks = set.findall("slack")
        assert len(slacks) == 1
        slack = slacks[0]
        sidx = vars.type.index("slack")
        slack.set("design", str(vars.design[sidx]))

        # export new file
        tree.write(self._dens_file(), xml_declaration=True)

    def _read_gradplot(self):
        # reset dataclasses
        self._init_dataclasses()
        # map to list
        data = [self.metadata, self.variables, self.constraints, self.jacobian]
        current_section: Optional[int] = None
        with open(self._grad_file(), "r") as gp:
            for lidx, line in enumerate(gp):

                line = line.strip()

                # skip empty lines at the start
                if not line and current_section is None:
                    continue

                # check if we have a section heading
                if current_section is None:
                    for k, dc in enumerate(data):
                        if dc.parse_line(line):
                            current_section = k
                            # print(f"{lidx}: Entering section {k}")
                # if no section found continue
                if current_section is None:
                    continue

                # get active data object (can be class or list)
                active = data[current_section]

                # parse lines until false
                if not active.parse_line(line):
                    # print(f"{lidx}: \tDONE")
                    current_section = None


if __name__ == "__main__":
    main()
