#!/usr/bin/env python3
"""Generate the two duct meshes (gmsh 2.2 ASCII) for the openCFS<->openCFS
characteristic (impedance-matched) acoustic coupling test.

A 1D plane-wave duct is split at the coupling plane x = Lf into two openCFS
acoustic participants that exchange the characteristic invariants via preCICE:

    Fluid    : x in [0, Lf]      source at S_in (x=0), characteristicCouplingBC at S_couple (x=Lf)
    Acoustic : x in [Lf, Lf+La]  characteristicCouplingBC at S_couple (x=Lf), absorbingBCs at S_far

Both ducts use a 1x1 (a x a) square cross-section meshed with one cell across (sound-hard
S_side walls keep the wave 1D). Sized like Singlefield/Acoustics/CharacteristicBC1D so a
1-cycle 5 kHz burst (lambda = c0/f = 68.2 mm) launched at S_in reaches the interface in
~15 steps and transmits into the Acoustic duct within a ~40-step run at dt = 2e-5 s.
Keeps each reference small and the runtime well under a second.
"""
import os

a  = 0.01        # cross-section side [m]
dx = 0.0068      # element size along x [m] (~ c0*dt, c0=341, dt=2e-5)
NX_F = 15        # Fluid duct elements   -> Lf = 0.102 m
NX_A = 20        # Acoustic duct elements -> La = 0.136 m

HERE = os.path.dirname(os.path.abspath(__file__))
corners = [(0.0, 0.0), (a, 0.0), (a, a), (0.0, a)]  # CCW in the y-z plane


def write_duct(path, x0, nx, name_left, name_right):
    def nid(i, c):
        return i * 4 + c + 1

    nodes = []
    for i in range(nx + 1):
        x = x0 + i * dx
        for (y, z) in corners:
            nodes.append((x, y, z))

    elements = []   # (type, phys, geom, [nodes])
    eid = [0]

    def add(etype, phys, ns):
        eid[0] += 1
        elements.append((eid[0], etype, phys, ns))

    # left / right cap quads (phys 1 / 2) and side walls (phys 3)
    add(3, 1, [nid(0,  0), nid(0,  1), nid(0,  2), nid(0,  3)])
    add(3, 2, [nid(nx, 0), nid(nx, 1), nid(nx, 2), nid(nx, 3)])
    for i in range(nx):
        j = i + 1
        add(3, 3, [nid(i,0), nid(i,3), nid(j,3), nid(j,0)])  # y = 0
        add(3, 3, [nid(i,1), nid(i,2), nid(j,2), nid(j,1)])  # y = a
        add(3, 3, [nid(i,0), nid(i,1), nid(j,1), nid(j,0)])  # z = 0
        add(3, 3, [nid(i,3), nid(i,2), nid(j,2), nid(j,3)])  # z = a
    # volume hexes (phys 4)
    for i in range(nx):
        j = i + 1
        add(5, 4, [nid(i,0), nid(i,1), nid(i,2), nid(i,3),
                   nid(j,0), nid(j,1), nid(j,2), nid(j,3)])

    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n4\n")
        f.write('2 1 "%s"\n2 2 "%s"\n2 3 "S_side"\n3 4 "V"\n' % (name_left, name_right))
        f.write("$EndPhysicalNames\n")
        f.write("$Nodes\n%d\n" % len(nodes))
        for k, (x, y, z) in enumerate(nodes, start=1):
            f.write("%d %g %g %g\n" % (k, x, y, z))
        f.write("$EndNodes\n")
        f.write("$Elements\n%d\n" % len(elements))
        for (e, etype, phys, ns) in elements:
            f.write("%d %d 2 %d %d %s\n" % (e, etype, phys, phys, " ".join(map(str, ns))))
        f.write("$EndElements\n")
    print("wrote %s: %d nodes, %d elements (%d hexes)" % (path, len(nodes), len(elements), nx))


if __name__ == "__main__":
    Lf = NX_F * dx
    write_duct(os.path.join(HERE, "channel_fluid.msh"), 0.0, NX_F, "S_in", "S_couple")
    write_duct(os.path.join(HERE, "channel_acou.msh"),  Lf,  NX_A, "S_couple", "S_far")
