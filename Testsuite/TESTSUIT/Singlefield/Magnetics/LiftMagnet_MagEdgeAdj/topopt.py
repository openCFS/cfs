from sim_setup import get_setup
import numpy as np
from pyCFS.topt import topt

MU0 = 4*np.pi*1e-7
EPS0 = 8.854e-12
MU_IRON = 1000

def forward_problem(sim_topt, p, normalize_grad=False):
    x = np.array([[]])
    sim_topt.update_design_parameters(p)
    sim_topt.sim(x)

    Ecp = sim_topt.sim.results[0][1]["elecFieldIntensity"]["V_objective"].squeeze().flatten()
    obj_value = -0.5*np.abs(np.dot(np.conjugate(Ecp), Ecp))
    dg_dp = sim_topt.compute_gradient(normalize=normalize_grad)
    volfrac = sim_topt.compute_volume_constraint(p)
    dvolfrac_dx = np.array([sim_topt.compute_volume_constraint_deriv()])
    
    return (obj_value, dg_dp, volfrac, dvolfrac_dx)

iron_params = {
                "magElemReluctivity": 1/(MU_IRON * MU0),
                "magElemPermeability": MU_IRON * MU0,
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
    "V_objective_domain": iron_params,
    "V_armature": iron_params,
    "V_inner_coil_part": air_params,
    "V_outer_coil_part": air_params,
    "V_design_domain": iron_params,
    "V_core": iron_params,
}

param_grad_map = {
                "magElemPermeability": 0,
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

print("-- finished setup --")

V_max = 0.7
p0 = np.ones(sim_topt.n_design_elems)*V_max
print(f"init xval[0]={p0[0]}")

x = np.array([[]])
sim_topt.update_design_parameters(p0)
sim_topt.sim(x)

E = sim_topt.sim.hist_results[0][1]["magEnergy"]["V_objective_domain"][0,0]
dg_dp = sim_topt.compute_gradient(normalize=False)

print(" -- finished sim -- ")
