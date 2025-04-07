#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:53:37 2017

@author: cjunger
"""

import numpy as np
import matplotlib.pyplot as plt

dt = 0.01
time_shift = 0*dt
t = np.arange(0.01, 6.29, dt)
ref = np.cos(t+time_shift)
data = np.loadtxt('dpdt0.00.csv',delimiter=',',skiprows=1)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(data[:,44],data[:,9]-ref)
#ax.plot(data[:,44],data[:,9])
#ax.plot(t,ref,linestyle='--')
#plt.xlim([0.1,5.9])
plt.ylim([-0.02,0.02])
ax.grid(True)
#plt.show()
