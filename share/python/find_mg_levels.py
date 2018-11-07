""""
This tool will help in determining the number of levels and the number of processor that might be ideal for multigrid

 For maximum performance its always good to have more levels.(Need to verify this hyposthesis)
 Therefore the levels are calculated using how many times the number of elements can be divided into half.

 the number of processor is also crucial since sometimes the division of the problem in all x,y and z directions makes
 the solver very slow

 This script can be used to find for a given problem the ideal range of procs and levels to set for the mg
 also the script suggests if increasing the mesh size by a small fraction will result in faster solution.

"""

from optimization_tools import *
import argparse
import math
import matplotlib.pyplot as plt


def valid_level(mesh_size):
    # The mesh size in each direction
    # we assume that the mesh size is same in all three directions. and indicates the number of elements
    level = 100
    for lvl in reversed(range(level)):
        divisor = math.pow(2.0, float(lvl - 1))
        if math.floor((mesh_size) / divisor) == ((mesh_size) / (divisor)):
            level = lvl
            break
    return level


# This function takes the ideal levels and now tries to find a ideal number of processor for the mesh size.
# The default processor distribution occurs in x, y and z direction of a 3d mesh. So when processor number is
# prime then we have no issues as the

def find_ideal_mesh(mesh_size, low_bound=10, upper_bound=10):

    # look in the +10 and -10 neighborhood of the mesh size to see if we can obtain a mesh which can use more levels
    level = valid_level(mesh_size)
    print(" The maximum number of level for the mesh size " + str(mesh_size) + " is " + str(level))
    mesh_changed = False
    for mesh in range(mesh_size - low_bound, mesh_size + upper_bound):
        level_range = valid_level(mesh)

        if level_range > level:
            level = level_range
            mesh_ideal = mesh
            mesh_changed = True
            
    if mesh_changed:
        print(" The ideal mesh for the given problem is " + str(mesh_ideal) + " with levels " + str(level))
        mesh_size = mesh_ideal
    else:
        print("The provided mesh " + str(
            mesh_size) + " has the maximum number of level " + str(
            level) + " in the neighbourhood. Try increasing the neighborhood of search. Default is 10")

    return level, mesh_size


def plot_mesh_levels(mesh_range):
    levels = list()
    [levels.append(valid_level(mesh)) for mesh in range(1, mesh_range)]
    l = plt
    l.plot(levels[:])
    l.xlabel('Number of Elements in a direction')
    l.ylabel('Max MG levels')
    l.title('Variaton of MG levels with Number of elements')
    l.show()
    # print(levels)

    return levels


def plot_mesh_procs(mesh_range):
    procs = list()
    [procs.append(find_procs(mesh)) for mesh in range(1, mesh_range)]

    l = plt
    l.plot(procs[:])
    l.xlabel('Mesh Size')
    l.ylabel('Max procs')
    l.title('Variaton of MG levels with Number of elements')
    l.show();

    return procs


def find_procs(mesh_size, levels=None,max_procs=100,):
    if levels is None:
        levels = valid_level(mesh_size)
    divisor = math.pow(2.0, float(levels - 1))

    # this ensures that for the coarsest mesh atleast one node is present for the  mesh size
    coarse_mesh = (mesh_size / divisor)
    if (max_procs % 3 == 0 or max_procs % 2 == 0) and (max_procs !=3 or max_procs !=2):
        for procs in reversed(range(max_procs + 1)):
            if procs % 3 != 0 and procs % 2 != 0 and (int(coarse_mesh / procs) >= 1):
                # print("The max procs to be used are " + str(procs))
                max_procs = procs
                break
        
    return max_procs


def ideal_parameters(mesh_size,max_procs=100,env=None,flexible_mesh_size=True):

    # The ideal parameters tries to give the levels and procs that will be ideal for the given mesh size and
    # depending upon the env. If you are flexible in defining the mesh size it can also find out a mesh thats
    # closer to the original mesh size but accomodates more levels
    # The cluster mode tries to go for using more procs in that we manually try to set the max possible levels at
    # half  of what is possible.
    levels = valid_level(mesh_size)
    if flexible_mesh_size:
        if env == "cluster":
            # Since we want to run in cluster we assume more procs are required
            levels,mesh_size = find_ideal_mesh(mesh_size,10,10)
            levels = int(levels-3)

        else:
            levels,mesh_size = find_ideal_mesh(mesh_size,10,10)

    procs = find_procs(mesh_size, levels, max_procs)

    return mesh_size,levels,procs

if __name__ == '__main__':
    # levels = plot_mesh_levels(121)
    # plot_mesh_procs(121)
    


    print(ideal_parameters(64,60,"cluster",True))

