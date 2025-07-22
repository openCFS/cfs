# coding: utf-8
u = linspace(-3000,+3000,100)
sig = np.zeros((u.size))
sig
def spark(u):
    if u < -1000:
        return 1e8 # ideal conductor
    if u > 1000:
        return 1e8 # ideal conductor
    else:
        return 1e-15 # air as isolator

for i in arange(u.size):
    sig[i] = spark(u[i])
    
sig
np.vstack((u, sig)).T
np.savetxt('sigAir.dat', np.vstack((u, sig)).T, delimiter=' ')
