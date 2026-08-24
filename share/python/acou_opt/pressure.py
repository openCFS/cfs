#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import hdf5_tools
from typing import Optional, List
import os
import warnings
import lxml
import lxml.etree
import math


FUNC = {"do": "dynamicOutput",
        "doT": "dynamicOutputTracking",
        "rw": "reflectedWave"}


def vn2p(vn: float = 1) -> float:
    """Calculate the incident wae pressure based on the normal velocity b.c. value."""
    RHO_0 = 1.204
    K = 0.142e+6
    Z0 = RHO_0 * np.sqrt(K / RHO_0)
    return -vn * Z0 / 2


def pscale() -> float:
    """Used to scale down to unity tau for waveguide width changes.
    Apply to p or tau, use scale**2 for wpoer scaled tau."""
    # TODO: this should not be hardcoded
    W_IN = 0.025
    W_OUT = 0.005
    return math.sqrt(W_OUT / W_IN)


def get_pressure(hdf5: str | os.PathLike,
                 eval: bool,
                 name: str,
                 freqsteps: Optional[np.ndarray] = None):
    """Get pressure from hdf5 file."""
    if freqsteps is None:
        # this can fail if the evaluate file is not properly written (not from harmonic driver)
        # hopefully fixed someday...
        freqsteps = np.array(hdf5_tools.get_step_values(hdf5)[0])

    p_out = []
    for idx, f in enumerate(freqsteps):
        step = idx
        if not eval:
            step = idx + 1  # evaluate results use 0-indexing, while normal simulation use 1-indexing
        # result = hdf5_tools.get_result(hdf5, result="acouPressure", region="L_out", step=step)
        result = hdf5_tools.get_result(hdf5, result="acouPressure", region="S_acou", step=step)
        nids = hdf5_tools.get_nodes_in_region(hdf5, "S_acou")
        try:
            nlout = hdf5_tools.get_nodes(hdf5, name)
        except KeyError:
            raise ValueError(f"Given nodeset name {name} not found!")
        ps_out = result[np.isin(nids, nlout)]
        if not np.allclose(np.real(ps_out), np.real(ps_out[0]), rtol=0.1, atol=0.1) or \
           not np.allclose(np.imag(ps_out), np.imag(ps_out[0]), rtol=0.1, atol=0.1):
            warnings.warn("Deviations on L_out to large, no plane wave!")
            # raise ValueError("Deviations on L_out to large, no plane wave!")
        p_out.append(np.sqrt(np.mean(ps_out * np.conjugate(ps_out))))

    return freqsteps, np.abs(np.array(p_out))


def get_trans(hdf5: str | os.PathLike,
              eval: bool,
              name: str,
              freqsteps: Optional[np.ndarray] = None):
    """Get power scaled transmission coefficient from hdf5 file.
    For vn=1 and widths hardcoded."""
    fs, ps = get_pressure(hdf5, eval, name, freqsteps)
    return fs, (ps / vn2p() * pscale())**2

def info2taus(ax,
              info_file: os.PathLike | str,
              type: str = "doT",
              normalize: bool = False):
    if type not in FUNC.keys():
        raise NotImplementedError()

    scale = 1
    if normalize:
        scale = pscale() / vn2p()

    parser = lxml.etree.XMLParser(remove_comments=True)
    tree = lxml.etree.parse(info_file, parser)
    root = tree.getroot()

    opt = root.find("optimization")
    if opt is None:
        raise IOError("Missing '<optimization>' element.")
    
    # extract frequencies, we use a set, because they can appear multiple times (multiload)
    freqs = set()
    for hexc in opt.find("header").find("excitations").findall("excitation"):
        freqs.add(float(hexc.get("frequency")))
    freqs = list(freqs)
    freqs.sort()
    print(f"Found freqs {freqs}")

    # find last iteration
    iters = opt.find("process").findall("iteration")
    idxs = tuple(int(i.get("number")) for i in iters)
    maxidx = max(idxs)
    iter = iters[idxs.index(maxidx)]

    # if iteration holds the type we are the objective value
    if iter.get(FUNC[type]) is not None or iter.get(f"{FUNC[type]}_f_2"):
        print("We are not slack -> parsing from objective value")
        slack = False
    else:
        print(f"We are slack -> parsing from {FUNC[type]} constraint value")
        slack = True

    # parse tau
    excitations = iter.findall("excitation")
    ps2 = []
    for exc in excitations:
        if slack:
            val = None
            try:
                val = float(exc.get(FUNC[type]))
            except TypeError:
                pass
            try:
                # for robust
                val = float(exc.get(f"{FUNC[type]}_f_2"))
            except TypeError:
                pass
            if val is None:
                warnings.warn(f"Could not find {FUNC[type]}")
                continue
        else:
            val = float(exc.get("objective"))

        match type:
            case "doT":
                ps2.append((1 - math.sqrt(val))*scale**2)
            case "do":
                ps2.append(val*scale**2)
            case "rw":
                if normalize:
                    scale = 1 / vn2p()
                ps2.append((vn2p()**2 - val)*scale**2)
            case _:
                raise NotImplementedError("Type not implemented!")
        # only parse the first unique set of constraints
        if len(ps2) == len(freqs):
            break

    # plot
    ax.plot(freqs, ps2, "x-", color="k", label=type)


def pressure(ax,
             hdf5s: List[str],
             labels: List[str],
             names: List[str],
             freqsteps: Optional[np.ndarray] = None,
             normalize: bool = False):
    scale = 1
    if normalize:
        scale = pscale() / vn2p()
    
    for label, hdf5_file in zip(labels, hdf5s):
        for name in names:
            freqsteps, ps = get_pressure(hdf5_file, "eval" in hdf5_file, name, freqsteps)
            if len(labels) > 1:
                lab = f"{label}_{name}"
            else:
                lab = name
            # ax_r.plot(freqsteps, np.abs((p_exc-p0)/p0)**2, "x-", label=label, markersize=2)
            ax.plot(freqsteps, np.abs(ps * scale)**2, "o-", label=lab, markersize=3)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="pressure",
        description="Plot the squared output pressure "
        "over frequency for the given simulations.")

    parser.add_argument("hdf5s", nargs="+",
                        help="One or more hdf5 simulation files to compare.")
    parser.add_argument("-fl", "--flower", type=float,
                        help="Start frequency.")
    parser.add_argument("-fu", "--fupper", type=float,
                        help="Start frequency.")
    parser.add_argument("-fn", "--fnum", type=int,
                        help="Number of frequencies.")
    parser.add_argument("-l", "--labels", nargs="+",
                        help="One or more labels for the plot.")
    parser.add_argument("-n", "--names", default=["L_out"], nargs="+",
                        help="Name(s) of the output region(s) to plot.")
    parser.add_argument("--save", default=None,
                        help="Image output path.")
    parser.add_argument("--noshow", action='store_true',
                        help="Don't show the transmission plot.")
    parser.add_argument("--nolegend", action='store_true',
                        help="Don't show the legend.")
    parser.add_argument("--normalize", action='store_true',
                        help="Scale transmission coefficient for vn=1 and hardcoded dimensions.")
    parser.add_argument("--doT", type=str,
                        help="Info file with dynamicOutputTracking objectives or constraints (only for parameter=1).")
    parser.add_argument("--do", type=str,
                        help="Info file with dynamicOutput objectives or constraints.")
    parser.add_argument("--rw", type=str,
                            help="Info file with reflectedWave objectives or constraints (only for vn=1 and hardcoded dimensions).")

    args = parser.parse_args()

    if args.labels is not None:
        assert len(args.hdf5s) == len(args.labels)
    else:
        args.labels = args.hdf5s

    freqsteps = None
    if args.flower and args.fupper and args.fnum:
        freqsteps = np.linspace(args.flower, args.fupper, args.fnum, endpoint=True)
    # print(freqsteps)

    # init plot
    mm = 1 / 25.4
    plt.rcParams.update({'font.size': 7})
    fig, ax = plt.subplots(1, sharex=True, figsize=(80*mm, 40*mm))

    # plot pressure form hdf5 files
    pressure(ax, args.hdf5s, args.labels, args.names, freqsteps, normalize=args.normalize)

    # plot from info files
    if args.doT:
        info2taus(ax, args.doT, type="doT", normalize=args.normalize)
    if args.do:
        info2taus(ax, args.do, type="do", normalize=args.normalize)
    if args.rw:
            info2taus(ax, args.rw, type="rw", normalize=args.normalize)

    # format plot       
    ax.set(ylabel="transmission", ylim=(0, 1), xlabel="frequency in Hz")
    if not args.nolegend:
        ax.set_title("pressure.py")
        ax.legend()
    ax.grid()
    fig.tight_layout()

    # save
    if args.save is not None:
        fig.savefig(args.save, dpi=300)

    # show
    if not args.noshow:
        plt.show()
