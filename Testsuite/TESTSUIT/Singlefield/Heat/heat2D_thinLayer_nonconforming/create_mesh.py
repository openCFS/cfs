import gmsh

# geometry parameters

domain_width = 0.2 # m
diameter_solid = 0.072 # m
diameter_solid_inner = 0.060 # m

mesh_size_fluid = 0.016
mesh_size_solid = 0.008


def create_outter_mesh(filename):
    gmsh.initialize()

    # Outter mesh (e.g. fluid)
    
    gmsh.model.add("outter_mesh")

    # exterior boundary 
    gmsh.model.geo.add_point(0, 0, 0, mesh_size_fluid, 1)
    gmsh.model.geo.add_point(0, domain_width / 2, 0, mesh_size_fluid, 2)
    gmsh.model.geo.add_point(0, -domain_width / 2, 0, mesh_size_fluid, 3)

    gmsh.model.geo.addCircleArc(2, 1, 3, 1)
    gmsh.model.geo.addCircleArc(3, 1, 2, 2)

    gmsh.model.geo.addCurveLoop([1, 2], 1)

    # interior boundary
    gmsh.model.geo.add_point(0, 0, 0, mesh_size_fluid, 5)
    gmsh.model.geo.add_point(0, diameter_solid / 2, 0, mesh_size_fluid, 6)
    gmsh.model.geo.add_point(0, -diameter_solid / 2, 0, mesh_size_fluid, 7)

    gmsh.model.geo.addCircleArc(6, 5, 7, 5)
    gmsh.model.geo.addCircleArc(7, 5, 6, 6)

    gmsh.model.geo.addCurveLoop([5, 6], 2)

    # create surface
    gmsh.model.geo.addPlaneSurface([1, 2], 1)

    # set physical lines
    gmsh.model.geo.addPhysicalGroup(1, [1, 2], 1)
    gmsh.model.setPhysicalName(1, 1, "l_exterior")
    gmsh.model.geo.addPhysicalGroup(1, [5, 6], 5)
    gmsh.model.setPhysicalName(1, 5, "l_fluid-solid")

    # set physical surfaces
    gmsh.model.geo.addPhysicalGroup(2, [1], 6)
    gmsh.model.setPhysicalName(2, 6, "s_fluid")

    gmsh.model.geo.synchronize()

    # options    
    options = dict()
    # for compatibility with openCFS
    options['Mesh.MshFileVersion'] = 2.1
    options['Mesh.Algorithm'] = 2


    for mo, val in iter(options.items()):
        gmsh.option.setNumber(mo, val)

    # mesh
    gmsh.model.mesh.generate(2)

    gmsh.write(filename)
    gmsh.write(filename.split(".")[0]+".geo_unrolled")
    gmsh.finalize()
    return

def create_solid_mesh(filename):
    gmsh.initialize()

    # Solid mesh
    
    gmsh.model.add("solid_mesh")

    # solid to fluid boundary 
    gmsh.model.geo.add_point(0, 0, 0, mesh_size_solid, 8)
    gmsh.model.geo.add_point(0, diameter_solid / 2, 0, mesh_size_solid, 9)
    gmsh.model.geo.add_point(0, -diameter_solid / 2, 0, mesh_size_solid, 10)
    gmsh.model.geo.addCircleArc(9, 8, 10, 7)
    gmsh.model.geo.addCircleArc(10, 8, 9, 8)

    gmsh.model.geo.addCurveLoop([7, 8], 3)

    # solid inner boundary
    gmsh.model.geo.add_point(0, diameter_solid_inner / 2, 0, mesh_size_solid, 11)
    gmsh.model.geo.add_point(0, -diameter_solid_inner / 2, 0, mesh_size_solid, 12)

    gmsh.model.geo.addCircleArc(11, 8, 12, 9)
    gmsh.model.geo.addCircleArc(12, 8, 11, 10)

    gmsh.model.geo.addCurveLoop([9, 10], 4)

    # create surface (solid)
    gmsh.model.geo.addPlaneSurface([3, 4], 2)


    # set physical lines
    gmsh.model.geo.addPhysicalGroup(1, [7, 8], 7)
    gmsh.model.setPhysicalName(1, 7, "l_solid-fluid") # outter surface solid
    gmsh.model.geo.addPhysicalGroup(1, [9, 10], 8)
    gmsh.model.setPhysicalName(1, 8, "l_solid_inner") # inner surface soldi

    # set physical surfaces
    gmsh.model.geo.addPhysicalGroup(2, [2], 9)
    gmsh.model.setPhysicalName(2, 9, "s_solid")

    gmsh.model.geo.synchronize()

    # options    
    options = dict()
    # for compatibility with openCFS
    options['Mesh.MshFileVersion'] = 2.1
    options['Mesh.Algorithm'] = 2
    #options['Mesh.RecombinationAlgorithm'] = 1  # Default 1; Simple 0
    # recombine triangular elements the quads
    #options['Mesh.RecombineAll'] = 1

    for mo, val in iter(options.items()):
        #print(mo, val)
        gmsh.option.setNumber(mo, val)

    # mesh
    gmsh.model.mesh.generate(2)

    gmsh.write(filename)
    gmsh.write(filename.split(".")[0]+".geo_unrolled")
    gmsh.finalize()
    return


def merge_mesh(filename):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.merge("mesh_fluid.msh")
    gmsh.merge("mesh_solid.msh")

    lines = gmsh.model.getEntities(dim=1)
    print("Lines:", lines)

    surfaces = gmsh.model.getEntities(dim=2)
    print("Surfaces:", surfaces)

    volumes = gmsh.model.getEntities(dim=3)
    print("Volumes:", volumes)

    options = dict()
    # for compatibility with openCFS
    options['Mesh.MshFileVersion'] = 2.1
    #options['Mesh.HighOrderOptimize'] = mesh_HighOrderOptimize # for positive jacobian

    for mo, val in iter(options.items()):
        print(mo, val)
        gmsh.option.setNumber(mo, val)

    gmsh.write(filename)
    gmsh.write(filename.split(".")[0]+".geo_unrolled")

    #gmsh.fltk.run()
    gmsh.finalize()


create_outter_mesh("mesh_fluid.msh")
create_solid_mesh("mesh_solid.msh")
merge_mesh("mesh_merged.msh")
