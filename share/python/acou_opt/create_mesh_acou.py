#!/usr/bin/env python
import argparse
import mesh_tool

# Global tolerance for element/node selection
eps = 1e-6


def ismul(a: float, b: float):
    """
    Test if a is multiple of b.
    Test if b / a gives non decimal number.
    """
    return abs(b / a - round(b / a)) < eps


def make_mesh(name, Nx, W_in, W_out, W_ch, W_abh, L_in, L_out, L_abh, t_pml: float = 0):
    assert W_ch <= W_abh, "W_out must be equal or smaller than W_abh!"
    # assert W_abh <= W_in, "W_abh must be equal or smaller than W_in!"

    ws = (W_in, W_abh, W_out)

    Y = t_pml + L_in + L_out + L_abh + t_pml  # total length in y-direction
    X = max(ws) + t_pml  # total width

    # calculate element size and Ny
    dx = max(ws) / Nx

    

    # test meshing
    for v, vname in zip([W_out, W_ch, W_abh, L_in, L_out, L_abh, t_pml],
                        ["W_out", "W_ch", "W_abh", "L_in", "L_out", "L_abh", "t_pml"]):
        if not ismul(dx, v):
            raise ValueError(f"dx = {dx} is not multiple of {vname} = {v} -> can't mesh.")

    # create mesh
    Nx = round(X / dx)
    Ny = round(Y / dx)
    print(f"Nx = {Nx}, Ny = {Ny}, dx = {dx:.4f}, dy = {dx:.4f}")
    mesh = mesh_tool.create_2d_mesh(Nx, Ny, X, Y)

    # define regions
    for i, e in enumerate(mesh.elements):
        x, y = mesh.calc_barycenter(e)

        # out pml
        if y >= L_out + L_abh + L_in + t_pml:
            if x >= max(ws):
                e.region = "S_pml"
                continue
            if x >= W_abh and x >= W_out:
                e.region = "S_pml_tout"
                continue
            if x >= W_out and not x >= W_abh:
                e.region = "S_pml_void"
                continue
            e.region = "S_pml_out"
            continue

        # outlet
        if y >= L_abh + L_in + t_pml:
            if x >= max(ws):
                e.region = "S_pml"
                continue
            if x >= W_out and not x >= W_abh:
                e.region = "S_void"
                continue
            e.region = "S_acou"
            continue

        # abh
        if y >= L_in + t_pml:
            if x >= max(ws):
                e.region = "S_pml"
                continue
            if x >= W_abh:
                e.region = "S_acou"
                continue
            if x >= W_ch:
                e.region = 'S_opt'
                continue
            e.region = "S_acou"
            continue

        # inlet
        if y >= t_pml:
            if x >= max(ws):
                e.region = "S_pml"
                continue
            if x >= W_in:
                e.region = "S_void"
                continue
            e.region = "S_acou"
            continue

        # in pml
        if x >= max(ws):
            e.region = "S_pml"
            continue
        e.region = "S_pml_in"

    # define node sets, because we use normalVelocity, we also need to create line elements and corresponding surface regions.
    Ls_in = []
    Ls_out = []
    Ls_axi = []
    Ls_top = []
    last_in_i = None
    last_out_i = None
    last_axi_i = None
    last_top_i = None
    for i, n in enumerate(mesh.nodes):
        x, y = n

        # inlet
        if abs(t_pml - y) < eps:
            if x <= W_in + eps:
                Ls_in.append(i)
                # create surface region
                if last_in_i is not None:
                    e = mesh_tool.Element(region="L_in_surf", type=mesh_tool.Ansys.LINE)
                    e.nodes = (last_in_i, i)
                    mesh.elements.append(e)
                last_in_i = i

        # outlet
        if abs(Y - t_pml - y) < eps:
            if x <= W_out + eps:
                Ls_out.append(i)
                # create surface region
                if last_out_i is not None:
                    e = mesh_tool.Element(region="L_out_surf", type=mesh_tool.Ansys.LINE)
                    e.nodes = (last_out_i, i)
                    mesh.elements.append(e)
                last_out_i = i

        # axi (this is needed for bloch periodic)
        if x < eps:
            Ls_axi.append(i)
            if last_axi_i is None:
                last_axi_i = i
                continue
            if y - t_pml < eps:
                e = mesh_tool.Element(region="L_axi_pml_in_surf", type=mesh_tool.Ansys.LINE)
            elif y - (L_in + t_pml) < eps:
                e = mesh_tool.Element(region="L_axi_in_surf", type=mesh_tool.Ansys.LINE)
            elif y - Y + t_pml < eps:
                e = mesh_tool.Element(region="L_axi_abh_surf", type=mesh_tool.Ansys.LINE)
            else:
                e = mesh_tool.Element(region="L_axi_pml_out_surf", type=mesh_tool.Ansys.LINE)
            e.nodes = (last_axi_i, i)
            mesh.elements.append(e)
            last_axi_i = i

        # top
        if abs(x - max(ws)) < eps:
            Ls_top.append(i)
            if last_top_i is None:
                last_top_i = i
                continue
            if y - t_pml < eps:
                e = mesh_tool.Element(region="L_top_pml_in_surf", type=mesh_tool.Ansys.LINE)
            elif y - (L_in + t_pml) < eps:
                e = mesh_tool.Element(region="L_top_in_surf", type=mesh_tool.Ansys.LINE)
            elif y - Y + t_pml < eps:
                e = mesh_tool.Element(region="L_top_abh_surf", type=mesh_tool.Ansys.LINE)
            else:
                e = mesh_tool.Element(region="L_top_pml_out_surf", type=mesh_tool.Ansys.LINE)
            e.nodes = (last_top_i, i)
            mesh.elements.append(e)
            last_top_i = i

    mesh.bc.append(('L_in', Ls_in))
    mesh.bc.append(("L_out", Ls_out))
    mesh.bc.append(('L_axi', Ls_axi))
    mesh.bc.append(("L_top", Ls_top))

    mesh_tool.write_ansys_mesh(mesh, name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="make_mesh.py",
                                     description="Create a regular mesh for acoustic optimization.")
    parser.add_argument("-lin", help="length of the inlet", default=0, type=float)
    parser.add_argument("-labh", help="length of the optimization region", default=0, type=float)
    parser.add_argument("-lout", help="length of the outlet", default=0, type=float)
    parser.add_argument("-win", help="width of the input", default=0, type=float)
    parser.add_argument("-wch", help="width of the channel", default=0, type=float)
    parser.add_argument("-wabh", help="width of the optimization region", default=0, type=float)
    parser.add_argument("-wout", help="width of the output", default=0, type=float)
    parser.add_argument("-tpml", help="width of the PML", default=0, type=float)
    parser.add_argument("-o", help="path of the output mesh file", type=str)
    parser.add_argument("nx", help="width resolution (elements over maximum width without pml)", type=int)

    args = parser.parse_args()

    name = args.o
    if name is None:
        name = f"mm_{args.nx}.mesh"

    make_mesh(name,
              Nx=args.nx,
              L_abh=args.labh,
              L_in=args.lin,
              L_out=args.lout,
              W_in=args.win,
              W_out=args.wout,
              W_ch=args.wch,
              W_abh=args.wabh,
              t_pml=args.tpml)
