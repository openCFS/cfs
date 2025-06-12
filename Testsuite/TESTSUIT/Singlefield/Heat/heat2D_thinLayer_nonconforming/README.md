# Readme

## Mesh Generation

Information to create the non conforming mesh for this test case with ```gmsh```.

### Prerequisites

Installed python with pip (python package manager).

Install gmsh via pip (python package manager):

```pip install --upgrade gmsh```

### Running the Script

Run python script:

```python3 create_mesh.py```

Remark: It creates first separte meshes for the fluid and solid domain (```mesh_fluid.msh``` & ```mesh_solid.msh```). Afterwards those files get merged into ```mesh_merged.msh``` which is used during simulation.

### Cleanup

The files ```mesh_fluid.msh``` & ```mesh_solid.msh``` and the respective ```*.geo_unrolled``` can be deleted, as those files are not used during simulation.
