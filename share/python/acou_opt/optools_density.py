from typing import Callable, Tuple
import os
import lxml
import lxml.etree


def nr2pos(nr: int,
           shape: Tuple[int, int],
           dim: Tuple[float, float]) -> Tuple[float, float]:
    """Convert the element number of a 2D regular mesh to its coordinates.

    Args:
        nr: Element number.
        shape: Mesh shape (Nx, Ny).
        dim: Mesh dimensions (x, y).
    Returns:
        Tuple of x and y coordinates for given element.
    """
    # extact values
    nx, ny = shape
    X, Y = dim
    dx = X / nx
    dy = Y / ny
    # one based index nr and first x rows from left to right and then y rows up
    # add half element width to get barycenters
    x = ((nr - 1) % nx) * dx + dx / 2
    y = ((nr - 1) // nx) * dy + dy / 2
    return (x, y)


def mod_density(ref_file: os.PathLike | str,
                output_file: os.PathLike | str,
                func: Callable):
    """Modifies a density file with a given function.

    Args:
        ref_file: Reference density file.
        output_file: Output (modified) density file.
        func: The function used for modification.

    Raises:
        IOError: For reference file format issues.

    Examples:
        The func should be of the following format:
        ```python
        def func(nr: int,
                 pos: Tuple[float, float],
                 opos: Tuple[float, float],
                 rho: float,
                 prho: float)
        ```

        With aruments:
        nr: Number of design element
        pos: Position with bottom left corner assumed at (0, 0)
        opos: Real position for offset meshes
        rho: Design value
        proh: Physical design value

        You don't have to specify all arguments, a simple example for thresholding would be
        ```python
        def threshold(prho, **kwargs):
            rho_out = 1
            if prho < th:
                rho_out = 0
            return rho_out, rho_out
        ```
    """
    parser = lxml.etree.XMLParser(remove_comments=True)
    tree = lxml.etree.parse(ref_file, parser)
    root = tree.getroot()

    # parse header
    header = root.find("header")
    if header is None:
        raise IOError("Missing '<header>' element, are you sure you loaded a .density.xml file?")

    # parse sets
    sets = root.findall("set")
    if header is None:
        raise IOError("Missing '<set>' element, there is no density data in your file.")
    ids = tuple(int(s.get("id")) for s in sets)
    print(f"Found sets with id(s): {ids}")
    # find maximum set id (last iteration)
    maxid = max(ids)
    set = sets[ids.index(maxid)]
    assert set.tag == "set"

    # extract mesh data
    mesh = header.find("mesh")
    if mesh is None:
        raise IOError("Missing '<mesh>' element, are you sure your mesh is regular?")
    nx, ny, nz = (int(d) for d in (mesh.get("x"), mesh.get("y"), mesh.get("z")))
    print(f"Found total mesh shape ({nx}, {ny}, {nz})")
    dm = header.find("coordinateSystems").find("domain")
    X = float(dm.get("max_x")) - float(dm.get("min_x"))
    Y = float(dm.get("max_y")) - float(dm.get("min_y"))
    print(f"Found total mesh dimensions ({X}, {Y})")
    od = header.find("optimizationDomain")
    ox_min = float(od.get("min_x"))
    oy_min = float(od.get("min_y"))

    # get list of elements
    elements = set.findall("element")
    print(f"Found {len(elements)} pseudo-densities")
    # change element data
    for element in elements:
        nr = int(element.get("nr"))
        pos = nr2pos(nr, (nx, ny), (X, Y))
        opos = (pos[0] - ox_min, pos[1] - oy_min)  # convert to rel position
        rho, prho = func(nr=nr,
                         pos=pos,
                         opos=opos,
                         rho=float(element.get("design")),
                         prho=float(element.get("physical")))
        # no idea whats used, so set both ...
        element.set("design", str(rho))
        element.set("physical", str(prho))

    # export new file
    tree.write(output_file, xml_declaration=True)
