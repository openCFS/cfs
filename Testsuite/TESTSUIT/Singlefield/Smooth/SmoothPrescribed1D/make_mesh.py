#!/usr/bin/env python3
# Generate the 1D duct mesh (gmsh 2.2 ASCII) for the SmoothPrescribed1D test: a
# straight duct with a 1x1 (a x a) square cross-section, meshed with Nx hex
# elements along x (single cell across). Same generator style as
# TESTSUIT/Singlefield/Acoustics/CharacteristicBC1D/make_mesh.py; the geometry
# (x in [x0, x0+L]) matches the original coupled channel case this test was
# derived from (the fluid part occupied [0, x0]). Only the volume V is used by
# the test; the surface groups are kept for consistency with the original mesh.
# Physical groups:
#   (2,1) S_interface  near face at x = x0
#   (2,2) S_out        far face at x = x0 + L
#   (2,3) S_side       4 lateral walls
#   (3,4) V            the duct volume
import sys

x0 = 0.03     # duct start [m]
L  = 1.5      # duct length [m]
a  = 0.01     # cross-section side [m]
Nx = 75       # hex elements along x

dx = L / Nx
# corner offsets (y, z), counter-clockwise in the y-z plane
corners = [(0.0, 0.0), (a, 0.0), (a, a), (0.0, a)]

def nid(i, c):           # node id for cross-section i (0..Nx), corner c (0..3)
    return i * 4 + c + 1

nodes = []
for i in range(Nx + 1):
    x = x0 + i * dx
    for (y, z) in corners:
        nodes.append((x, y, z))

elements = []            # (type, phys, geom, [node ids])
eid = 0
def add(etype, phys, geom, ns):
    global eid
    eid += 1
    elements.append((eid, etype, phys, geom, ns))

# S_interface (x=x0) and S_out (x=x0+L): one quad each (type 3)
add(3, 1, 1, [nid(0,  0), nid(0,  1), nid(0,  2), nid(0,  3)])
add(3, 2, 2, [nid(Nx, 0), nid(Nx, 1), nid(Nx, 2), nid(Nx, 3)])

# S_side: 4 lateral quads per element (type 3), phys tag 3
for i in range(Nx):
    j = i + 1
    add(3, 3, 3, [nid(i,0), nid(i,3), nid(j,3), nid(j,0)])  # y = 0
    add(3, 3, 3, [nid(i,1), nid(i,2), nid(j,2), nid(j,1)])  # y = a
    add(3, 3, 3, [nid(i,0), nid(i,1), nid(j,1), nid(j,0)])  # z = 0
    add(3, 3, 3, [nid(i,3), nid(i,2), nid(j,2), nid(j,3)])  # z = a

# Volume hexes (type 5), phys tag 4: face at x_i then face at x_{i+1}
for i in range(Nx):
    j = i + 1
    add(5, 4, 4, [nid(i,0), nid(i,1), nid(i,2), nid(i,3),
                  nid(j,0), nid(j,1), nid(j,2), nid(j,3)])

out = sys.argv[1] if len(sys.argv) > 1 else "channel.msh"
with open(out, "w") as f:
    f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    f.write("$PhysicalNames\n4\n")
    f.write('2 1 "S_interface"\n2 2 "S_out"\n2 3 "S_side"\n3 4 "V"\n')
    f.write("$EndPhysicalNames\n")
    f.write("$Nodes\n%d\n" % len(nodes))
    for k, (x, y, z) in enumerate(nodes, start=1):
        f.write("%d %g %g %g\n" % (k, x, y, z))
    f.write("$EndNodes\n")
    f.write("$Elements\n%d\n" % len(elements))
    for (e, etype, phys, geom, ns) in elements:
        f.write("%d %d 2 %d %d %s\n" % (e, etype, phys, geom, " ".join(map(str, ns))))
    f.write("$EndElements\n")
print("wrote %s: %d nodes, %d elements (%d hexes)" % (out, len(nodes), len(elements), Nx))
