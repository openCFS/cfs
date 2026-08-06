# embedded python scripts

Only scripts which are loaded by cfs itself at runtime, via the special path key `cfs:share:python`:

```xml
<optimizer type="python">
  <python file="hessian_scipy.py" path="cfs:share:python"/>
</optimizer>
```

They are no standalone tools: they import the built-in `cfs` module (`cfs.get_gradient()`,
`cfs.feature_mapping_layout()`, ...), which exists only inside the embedded interpreter, and their
functions are called back from cfs.

A `__main__` section does not disqualify a script as long as it only replays what cfs produced.
`mma.py` for instance debugs a single subproblem offline from a pickle written by cfs
(`<option key="save_pickle" value="mma_initial.pickle"/>`), substituting a `FakeCFS` for the
built-in module - that is embedded debugging, not a tool of its own.

`PythonKernel::LoadPythonModule()` resolves `cfs:share:python` by searching this directory **first**
and `share/python` second. Both end up in `sys.path` (this one in front), hence

* a script may be moved in or out of here without touching a single xml file,
* the libraries one level up (`cfs_utils`, `optimization_tools`, `matviz_*`, ...) are imported
  plainly (`import optimization_tools as ot`),
* the manual `PYTHONPATH` stays a single entry `<cfs>/share/python`.

Keep the module names unique with respect to `share/python` - as this directory comes first in
`sys.path`, a duplicate name would shadow the file up there.

## What does not belong here

Dual use scripts, that is scripts which cfs loads but which are also a commandline tool for own data
or are imported by a manual script. They stay in `share/python`, otherwise a manual
`import spaghetti` breaks. Currently `spaghetti.py` and `spaghetti3d.py` (postprocessing/plotting
from the commandline, `combine_angles` is imported by `rotviz.py`) and `mesh_tool.py`
(imported by `basecell.py`, `benchmark.py`, ...).
