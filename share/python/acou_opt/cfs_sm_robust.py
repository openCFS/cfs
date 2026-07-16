#!/usr/bin/env python
import cyipopt
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

    args = parser.parse_args()

    solve_robust(args.name,
                 t=args.numThreads,
                 m=args.meshFile,
                 p=args.paramFile,
                 x=args.densityFile)


def solve_robust(name,
                 t: int = 1,
                 m: Optional[str] = None,
                 p: Optional[str] = None,
                 x: Optional[str] = None,
                 shift: float = 0.001):
    # set this so hopefully ipopt pardiso runs on multiple treads
    # cfs (hopefully) does not reset this, because it runs in a different process
    os.environ['OMP_NUM_THREADS'] = str(t)
    os.environ['MKL_NUM_THREADS'] = str(t)

    if p is None:
        p = f"{name}.xml"
    elif not p.endswith(".xml"):
        raise IOError("Specify p with .xml ending!")
    cfse = CFS(f"{name}_erode", t=t, m=m, p=p, x=x)
    cfs = CFS(name, t=t, m=m, p=p, x=x)
    cfsd = CFS(f"{name}_dilate", t=t, m=m, p=p, x=x)

    # we assume slack is the last variable
    assert cfs.variables.type[-1] == "slack"

    # create offset array for profiles
    offset = [shift if vt == "profile" else 0 for vt in cfs.variables.type]

    c = cfs.metadata.n_glob_constraints
    nlp = cyipopt.Problem(
        n=cfs.metadata.n,
        m=2 * c + cfse.metadata.m,  # +2 dynamicOutput each
        problem_obj=Robust(cfs, cfse, cfsd, offset),
        lb=[lb + o for lb, o in zip(cfs.variables.lb, offset)],
        ub=[ub - o for ub, o in zip(cfs.variables.ub, offset)],
        cl=cfse.constraints.cl[:c] + cfsd.constraints.cl[:c] + cfs.constraints.cl,
        cu=cfse.constraints.cu[:c] + cfsd.constraints.cu[:c] + cfs.constraints.cu
    )
    nlp.add_option("max_iter", 500)

    nlp.add_option("linear_solver", "pardisomkl")  # we have it so why not use it
    # using multiple threads os only for the solver
    # our NLP functions, where we call compute from, dont need to be threadsafe!

    # suggested by for non convex problems ( I don't know if we need it here)
    nlp.add_option("mu_strategy", "adaptive")

    # for acceptable stop all below conditions have to be fullfilled
    nlp.add_option("acceptable_obj_change_tol", 1e-4)
    nlp.add_option("acceptable_iter", 5)
    nlp.add_option("acceptable_tol", 1e-4)

    # for optimal stop below conditoin must be met
    nlp.add_option("tol", 1e-6)

    # gradient_check(nlp)
    nlp.add_option('output_file', f"{name}.ipopt")
    nlp.solve(cfs.variables.design)


def solve_normal(name,
                 t: int = 1,
                 m: Optional[str] = None,
                 p: Optional[str] = None,
                 x: Optional[str] = None):
    cfs = CFS(name, t=t, m=m, p=p, x=x)
    nlp = cyipopt.Problem(
        n=cfs.metadata.n,
        m=cfs.metadata.m,
        problem_obj=Normal(cfs),
        lb=cfs.variables.lb,
        ub=cfs.variables.ub,
        cl=cfs.constraints.cl,
        cu=cfs.constraints.cu,
    )
    gradient_check(nlp)
    nlp.add_option('output_file', f"{name}.ipopt")
    nlp.solve(cfs.variables.design)


def gradient_check(nlp):
    # run gradient check
    nlp.add_option("max_iter", 1)
    nlp.add_option("derivative_test", "first-order")
    # we need a larger pertubation and tol, because of rounding?
    nlp.add_option("derivative_test_tol", 5e-3)
    nlp.add_option("derivative_test_perturbation", 1e-5)


# ====================== NLP Classes ======================
class Robust():
    def __init__(self, cfs: "CFS", cfs_erode: "CFS", cfs_dilate: "CFS", offset: NDArray[np.floating]):
        self._cfs = cfs
        self._cfs_erode = cfs_erode
        self._cfs_dilate = cfs_dilate

        self._offset = offset

    def objective(self, x):
        slackidx = self._cfs.variables.type.index("slack")
        return x[slackidx]

    def gradient(self, x):
        g = np.zeros(len(self._cfs.variables.vid))
        slackidx = self._cfs.variables.type.index("slack")
        g[slackidx] = 1
        return g

    def constraints(self, x):
        """Returns the constraints."""
        # recompute for given set of variables
        c = self._cfs.metadata.n_glob_constraints
        self._cfs_erode.compute(x - self._offset)
        self._cfs.compute(x)
        self._cfs_dilate.compute(x + self._offset)
        return (self._cfs_erode.constraints.value[:c] +
                self._cfs_dilate.constraints.value[:c] +
                self._cfs.constraints.value)

    def jacobian(self, x):
        # number of additional constraints per cfs
        c = self._cfs.metadata.n_glob_constraints
        n = self._cfs.metadata.n
        grad = (self._cfs_erode.jacobian.gradient[:c*n] +
                self._cfs_dilate.jacobian.gradient[:c*n] +
                self._cfs.jacobian.gradient)
        return grad

    def jacobianstructure(self):
        # number of additional constraints per cfs
        c = self._cfs.metadata.n_glob_constraints
        n = self._cfs.metadata.n
        cid = (self._cfs_erode.jacobian.cid[:c*n] +
               [cid + c for cid in self._cfs_dilate.jacobian.cid[:c*n]] +
               [cid + 2*c for cid in self._cfs.jacobian.cid])
        vid = (self._cfs_erode.jacobian.vid[:c*n] +
               self._cfs_dilate.jacobian.vid[:c*n] +
               self._cfs.jacobian.vid)
        return (cid, vid)


class Normal():
    def __init__(self, cfs: "CFS"):
        self._cfs = cfs

    def objective(self, x):
        """Returns the scalar value of the objective given x."""
        # for slack we dont recompute
        # we have one slack objective
        slackidx = self._cfs.variables.type.index("slack")
        return x[slackidx]

    def gradient(self, x):
        """Returns the gradient of the objective with respect to x."""
        # we dont have to recompute, because the slack gradient is constant
        # gradient is zero except for dJ/dslack its one
        g = np.zeros(len(self._cfs.variables.vid))
        slackidx = self._cfs.variables.type.index("slack")
        g[slackidx] = 1
        return g

    def constraints(self, x):
        """Returns the constraints."""
        # recompute for given set of variables
        self._cfs.compute(x)
        return self._cfs.constraints.value

    def jacobian(self, x):
        """Returns the Jacobian of the constraints with respect to x."""
        return self._cfs.jacobian.gradient

    def jacobianstructure(self):
        return (self._cfs.jacobian.cid, self._cfs.jacobian.vid)


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

    def compute(self, x: Optional[np.ndarray] = None) -> bool:
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
        return True

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
