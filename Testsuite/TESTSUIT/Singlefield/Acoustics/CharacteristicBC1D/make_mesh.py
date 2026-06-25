#!/usr/bin/env python3
# Generate a small 1D plane-wave duct (gmsh 2.2 ASCII) for the characteristicCouplingBC
# regression test. A straight duct of length L with a 1x1 (a x a) square cross-section,
# meshed with Nx hex elements along x (single cell across), so the wave stays 1D
# (sound-hard side walls). Physical groups match the original channelBC mesh:
#   (2,1) S_interface  far/near source face at x = 0
#   (2,2) S_out        far face at x = L   (device under test: characteristicCouplingBC)
#   (2,3) S_side       4 lateral walls (no BC => sound-hard => keeps the wave plane)
#   (3,4) V            the duct volume (air)
#
# Sized so a 1-cycle 5 kHz burst (lambda = c0/f = 341/5000 = 68.2 mm) traverses the
# duct (L = 0.17 m ~ 2.5 lambda, ~25 elements) in ~25 time steps at dt = 2e-5 s, and
# is absorbed within a ~40-step run. Keeps the .h5ref tiny and runtime well under 1 s.
import sys

L  = 0.17     # duct length [m]
a  = 0.01     # cross-section side [m]
Nx = 25       # hex elements along x

dx = L / Nx
# corner offsets (y, z), counter-clockwise in the y-z plane (matches the original mesh)
corners = [(0.0, 0.0), (a, 0.0), (a, a), (0.0, a)]

def nid(i, c):           # node id for cross-section i (0..Nx), corner c (0..3)
    return i * 4 + c + 1

nodes = []
for i in range(Nx + 1):
    x = i * dx
    for (y, z) in corners:
        nodes.append((x, y, z))

elements = []            # (type, phys, geom, [node ids])
eid = 0
def add(etype, phys, geom, ns):
    global eid
    eid += 1
    elements.append((eid, etype, phys, geom, ns))

# S_interface (x=0) and S_out (x=L): one quad each (type 3)
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
