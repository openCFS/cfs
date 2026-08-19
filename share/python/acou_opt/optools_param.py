""" Functions for param.xml file manipulation """
import cfs_utils
import lxml.etree


def apply_evaluate(xml: lxml.etree.ElementTree,
                   fstart: float,
                   fstop: float,
                   fnum: int):
    """Modify the xml file used for optimization for evaluation.

    Args:
        xml: Openend and parsed xml file.
        fstart: Starting frequency.
        fstop: Stopping frequency.
        fnum : Numer of frequency steps.
    """
    # if it does not specify multi frequency add it
    if fnum > 1:
        cfs_utils.replace(xml, '//cfs:costFunction/@multiple_excitation', "true")
    else:
        cfs_utils.replace(xml, '//cfs:costFunction/@multiple_excitation', "false")
    cfs_utils.replace(xml, '//cfs:optimizer/@type', "evaluate")
    cfs_utils.replace(xml, '//cfs:commit/@mode', "each_forward")
    cfs_utils.replace(xml, '//cfs:numFreq', fnum)
    cfs_utils.replace(xml, '//cfs:startFreq', fstart)
    cfs_utils.replace(xml, '//cfs:stopFreq', fstop)
    cfs_utils.replace(xml, '//cfs:sampling', "linear")


def apply_em(xml: lxml.etree.ElementTree,
             file: str,
             access: str = "physical"):
    cfs_utils.add(xml, "//cfs:cfsSimulation", "loadErsatzMaterial", attribs={"file": file, "name": access}, before="cfs:optimization")


def apply_step(xml: lxml.etree.ElementTree,
               step: float) -> bool:
    """Set the stepsize for projection.
    This can be beta or transmission zone, the function automatically detects the right case.

    Args:
        xml: Openend and parsed xml file.
        step: Value for projection step.

    Returns:
        True if a valid case was found and step was set, False otherwise.
    """
    search_mode = True
    # test for robust
    if search_mode:
        try:
            cfs_utils.replace(xml, "//cfs:filter[@robust_excitation=0]//cfs:density/@beta", step)
            # replace(xml, "//cfs:filter[@robust_excitation=0]//cfs:density/@eta", 0.5 - deta)

            cfs_utils.replace(xml, "//cfs:filter[@robust_excitation=1]//cfs:density/@beta", step)
            # replace(xml, "//cfs:filter[@robust_excitation=1]//cfs:density/@eta", 0.5 + deta)

            cfs_utils.replace(xml, "//cfs:filter[@robust_excitation=2]//cfs:density/@beta", step)
            # replace(xml, "//cfs:filter[@robust_excitation=2]//cfs:density/@eta", 0.5)
            print("Detected mode ROBUST.")
            search_mode = False
        except RuntimeError:
            pass
    # test for shapemap
    if search_mode:
        try:
            cfs_utils.replace(xml, "//cfs:shapeMap/@beta", step)
            print("Detected mode SHAPEMAP.")
            search_mode = False
        except RuntimeError:
            pass
    # test for shapemap
    if search_mode:
        try:
            cfs_utils.replace(xml, "//cfs:featureMapping/@transition", step)
            print("Detected mode FEATUREMAP.")
            search_mode = False
        except RuntimeError:
            pass
    # else: normal simp mode
    if search_mode:
        try:
            cfs_utils.replace(xml, "//cfs:filter//cfs:density/@beta", step)
            print("Detected mode SIMP.")
            search_mode = False
        except RuntimeError:
            pass
    return not search_mode
