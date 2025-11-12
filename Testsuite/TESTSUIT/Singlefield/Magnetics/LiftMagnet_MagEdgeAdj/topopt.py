from sim_setup import get_setup
import numpy as np
from pyCFS.topt import topt

MU0 = 4*np.pi*1e-7
EPS0 = 8.854e-12

iron_params = {
                "magElemReluctivity": 1/(1.0 * MU0),
                "magElemPermeability": 1.0 * MU0,
                 "magElemPermittivity": 1.0 * EPS0, 
                 "magElemConductivity": 1.0, 
                }

air_params = {
              "magElemReluctivity": 1/(1e0 * MU0), 
              "magElemPermeability": 1e0 * MU0, 
              "magElemPermittivity": 1.0e0 * EPS0, 
              "magElemConductivity": 1.0,
              }

region_params = {
    "V_air": air_params,
    "V_objective_domain": air_params,
    "V_armature": air_params,
    "V_inner_coil_part": air_params,
    "V_outer_coil_part": air_params,
    "V_design_domain": air_params,
    "V_core": iron_params,
}

param_grad_map = {
                "magElemPermeability": 0,
                "magElemPermittivity": 1,
                }

# get simulation setup : 
sim = get_setup()
# sim(np.array([[]]), mesh_only=True)

# construct the topology optimizer : 
sim_topt = topt.Topt(
    sim, 
    "3d",
    "V_design_domain",
    region_params,
    param_grad_map,
    penalization_method="simp",
    # penalization_method="sigmoid",
    r_filter=1.0e-2,
)