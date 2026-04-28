import numpy as np

# parameters of air
rho0 = 1.204
eta0 = 18.232e-6
kappa = 25.68410e-6
gamma = 1.4
Cp = 1006.825
p0 = 100325

# load JCAL from dataset
data = np.loadtxt("JCAL_params_Basotect.csv", delimiter=",", skiprows=1)
data = data[6]  # select D=98mm and t = 50mm
t, D, phi, k0, lmbd, lmbds, sigma, alpha = data


# create freq vector
fs = np.linspace(400, 600, 3, endpoint=True)
ws = 2 * np.pi * fs

# compute JCAL
rho = alpha * rho0 / phi
rho *= (1 + (sigma * phi) / (1j * ws * rho0 * alpha) * np.sqrt(1 + 1j * (4 * alpha**2 * eta0 * rho0 * ws) / (sigma**2 * lmbd**2 * phi**2)))
K = gamma * p0 / phi 
K /= (gamma - (gamma - 1) / (1 - 1j * phi * kappa / (k0 * Cp * rho0 * ws) * np.sqrt(1 + 1j * 4 * k0**2 * Cp * rho0 * ws / (kappa * lmbds**2 * phi**2))))

np.savetxt("densReal.txt", np.vstack((fs, np.real(rho))).T, delimiter="\t")
np.savetxt("densImag.txt", np.vstack((fs, np.imag(rho))).T, delimiter="\t")
np.savetxt("cmpMReal.txt", np.vstack((fs, np.real(K))).T, delimiter="\t")
np.savetxt("cmpMImag.txt", np.vstack((fs, np.imag(K))).T, delimiter="\t")
