# calculate and plot current, flux and inductance
# Dominik Perchtold, TU Wien

import pylab

E1 = pylab.loadtxt('history/Coil3DVoltEdge-magEnergy-region-air.hist')
E2 = pylab.loadtxt('history/Coil3DVoltEdge-magEnergy-region-coil.hist')
I = pylab.loadtxt('history/Coil3DVoltEdge-coilCurrent-coil-AirCoreCoil.hist')
Psi = pylab.loadtxt('history/Coil3DVoltEdge-coilLinkedFlux-coil-AirCoreCoil.hist')
Phi = pylab.loadtxt('history/Coil3DVoltEdge-magFlux-surfRegion-mainFlux.hist')

N = 50

# symmetry corrections
Phi[:,1] = abs(Phi[:,1])*4
I[:,1] = I[:,1]*2
Psi[:,1] = Psi[:,1]*4
W = (E1[:,1] + E2[:,1])*8

# inductance with corrected values
LW = 2*W/I[:,1]**2 # static inductance (only valid for linear magnetics without permanent magnets)
Ld = pylab.diff(Psi[:,1])/pylab.diff(I[:,1]) # differential inductance (always valid)

# current
pylab.subplot(311)
pylab.plot(I[:,0],I[:,1])
pylab.grid(True)
pylab.xlabel('t (s)')
pylab.ylabel('i (A)')

# flux through center, flux linkage divided by N (should be almost the same)
pylab.subplot(312)
pylab.plot(I[:,1],Phi[:,1],'-',I[:,1],Psi[:,1]/N,'-')
pylab.grid(True)
pylab.xlabel('i (A)')
pylab.ylabel('Phi (Wb)')

# inductance via energy, differential inductance (should be the same)
pylab.subplot(313)
pylab.plot(I[:,1],LW,'-',I[1:,1],Ld,'-')
pylab.grid(True)
pylab.xlabel('i (A)')
pylab.ylabel('L (H)')

pylab.show()
