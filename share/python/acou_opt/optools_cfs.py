""" Functions for result.cfs file manipulation """
import os
import h5py
import numpy as np


def enable_anim(cfs_file: str | os.PathLike,
                fstart: int,
                fstop: int,):
    """This is a workaround, to enable animation of the "evaluate" results.
        We manually overwrite the attributes, to mimimic the output of HarmonicDriver, default export is TimeStep...

    Args:
        cfs_file: Path to the cfs output file.
        fstart: Starting frequency.
        fstop: Stopping frequency.
    """

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