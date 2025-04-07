# -*- coding: utf-8 -*-
"""
1. Calculates the absorption coefficient from the the wall impedance as it is
   specified for the simulation for 100Hz, 500Hz, 1000Hz and 1500Hz.
2. Calculates the absorption coefficient in analogy to the measurement in an
   impedance tube according to Uebertragungsfunktionsmethode (ISO 10534-2)

@author: cjunger
"""
import numpy as np

rho0 = 1.205 # density
c = 343 # speed of sound
f = np.array([100, 500, 1000, 1500]) # frequencies
k = 2*np.pi*f[:]/c # wavenumber
L = 0.8 # length of tube
x1= 0.3 # distance of the first microphone to the impedanceBC
s = 0.1 # distance between microphones, for higher frequencies this has to be
        # reduced. accordig to the norm s < 0.45 lambda

# compute the values from the input quantities
Z0 = rho0*c                               # impedance of air
Zw= np.array([1.56492e+001-1.31373e+001j, # given wall impedance from files
              1.53945e+000-4.25873e+000j,
              1.85736e+000-1.50777e+000j,
              1.81165e+000-4.40121e-001j])

r_input = np.abs((Zw-Z0)/(Zw+Z0)) # given reflection coefficient
alpha_input = 1-r_input**2        # given absorbtion coefficient
print("Alpha specified in simulation input = {0}".format(alpha_input))


# get simulation results
p1 = np.loadtxt("history/ImpedanceBC_axi-acouPressure-node-106-mic1.hist", usecols=[1,2])
p2 = np.loadtxt("history/ImpedanceBC_axi-acouPressure-node-126-mic2.hist", usecols=[1,2])
p1=p1[:,0]+p1[:,1]*1j   # pressure at microphone 1
p2=p2[:,0]+p2[:,1]*1j   # pressure at microphone 2

# compute the values from simulation results
H12 = p2/p1           # transferfunction
Hh = np.exp(-1j*k*s)  # transferfunction of the outgoing wave
Hr = np.exp(1j*k*s)   # transferfunction of the incoming wave

r = np.abs(((H12-Hh)/(Hr-H12)) * np.exp(1j*2*k*(L-x1)))
alpha = 1-r**2
print("Alpha computed from simulation = {0}".format(alpha))

diff = np.abs(alpha_input-alpha)/alpha_input*100
print("Error = {}%".format(diff))