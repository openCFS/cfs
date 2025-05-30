import numpy as np
import matplotlib.pyplot as plt

frequency1 = 50 #Hz
frequency2 = frequency1/2
frequency3 = frequency1/3

timestep = 1.0/(6*frequency1)
amplitude1 = 1.0
amplitude2 = -1.5
amplitude3 = 0.5

numSteps = 18

timeVec = np.linspace(0,numSteps*timestep,numSteps,endpoint=False)
sinSig1 = amplitude1*np.sin(2*np.pi*frequency1*timeVec)
sinSig2 = amplitude2*np.sin(2*np.pi*frequency2*timeVec)
sinSig3 = amplitude3*np.sin(2*np.pi*frequency3*timeVec)

#plt.plot(timeVec,sinSig1,label='sin1')
#plt.plot(timeVec,sinSig2,label='sin2')
#plt.plot(timeVec,sinSig3,label='sin3')
plt.plot(timeVec,sinSig1+sinSig2+sinSig3,label='sin1+sin2+sin3')
plt.legend()
plt.show()

np.savetxt('inputCurrent.txt',np.c_[timeVec,sinSig1+sinSig2+sinSig3])