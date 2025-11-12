# Import the needed libraries :
import numpy as np
from pyCFS import pyCFS

def get_setup(project_name="LiftMagnet_MagEdgeAdj", generate_mesh=False):

    # Set up project :
    # project_name = sim_name#"scatter_ring"
    cfs_path = "/home/eniz/Devel/CFS_BIN/build_opt_TO/bin/cfs"

    # Set simulation parameter names :
    cfs_params = []
    trelis_params = []

    # Construct cfssimulation object :
    sim = pyCFS(
        project_name,
        cfs_path,
        cfs_params_names=cfs_params,
        trelis_params_names=trelis_params,
        trelis_version="/home/eniz/Public/Cubit/bin/coreform_cubit",
        n_threads=6,
        save_hdf_results=True,
    )

    if generate_mesh:
        sim(np.array([[]]), mesh_only=True)

    return sim
