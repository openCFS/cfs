#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 06:04:17 2020

@author: kroppert
"""


import matplotlib.pylab as plt
import scipy.integrate
import numpy as np
from sys import exit
#from matplotlib2tikz import save as tikz_save


# =============================================================================
# Input Section
# =============================================================================

# Basename of the cfs result file
basename = 'capacitor_quarter'


endstep = 280



fig1, ax1 = plt.subplots()
file = 'history/{0}-magEddyCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current, phase= np.loadtxt(file, skiprows=3, unpack=True)
#ax1.semilogx(f,1/current)
ax1.semilogx(f[0:endstep],phase[0:endstep])
ax1.grid(True)
ax1.set_xlabel('Frequency in Hz')
ax1.set_ylabel('Phase')
ax1.set_title('Phase')


fig3, ax3 = plt.subplots()
file = 'history/{0}-displacementCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current, phase= np.loadtxt(file, skiprows=3, unpack=True)
#ax1.semilogx(f,1/current)
ax3.semilogx(f[0:endstep],current[0:endstep])
ax3.grid(True)
ax3.set_xlabel('Frequency in Hz')
ax3.set_ylabel('displacement current')
ax3.set_title('dD/dt')



fig5, ax5 = plt.subplots()
file_lad = 'history/{0}-magEddyCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current_lad, phase= np.loadtxt(file_lad, skiprows=3, unpack=True)
file_disp = 'history/{0}-displacementCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current_disp, phase= np.loadtxt(file_disp, skiprows=3, unpack=True)
ax5.loglog(f[0:endstep],current_disp[0:endstep]/(current_disp[0:endstep]+current_lad[0:endstep]))
ax5.grid(True)
ax5.set_xlabel('Frequency in Hz')
ax5.set_ylabel('displacement current/total current')
ax5.set_title('dD/dt / I')

fig6, ax6 = plt.subplots()
file_lad = 'history/{0}-magEddyCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current_lad, phase= np.loadtxt(file_lad, skiprows=3, unpack=True)
file_disp = 'history/{0}-displacementCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current_disp, phase= np.loadtxt(file_disp, skiprows=3, unpack=True)
ax6.loglog(f[0:endstep],current_lad[0:endstep]/(current_disp[0:endstep]+current_lad[0:endstep]), label='lad')
ax6.loglog(f[0:endstep],current_disp[0:endstep]/(current_disp[0:endstep]+current_lad[0:endstep]), label='disp')
ax6.grid(True)
ax6.legend()
ax6.set_xlabel('Frequency in Hz')
ax6.set_ylabel('I_ladungsgebunden /total current')
ax6.set_title('I_ladungsgebunden / I')
#tikz_save('data.tex')


# fig4, ax4 = plt.subplots()
# file = 'history/{0}-magEddyCurrent-surfRegion-S_electrode_top.hist'.format(basename)
# f, current, phase= np.loadtxt(file, skiprows=3, unpack=True)
# #ax1.semilogx(f,1/current)
# ax4.semilogx(f,current)
# ax4.grid(True)
# ax4.set_xlabel('Frequency in Hz')
# ax4.set_ylabel('eddy current')
# ax4.set_title('eddy current')


fig2, ax2 = plt.subplots()
file = 'history/{0}-magEddyCurrent-surfRegion-S_electrode_top.hist'.format(basename)
f, current, phase= np.loadtxt(file, skiprows=3, unpack=True)
ax2.loglog(f[0:endstep],1/(current_lad[0:endstep]+current_disp[0:endstep]))
ax2.grid(True)
ax2.set_xlabel('Frequency in Hz')
ax2.set_ylabel('Impedance in Ohm')
ax2.set_title('Impedance')




