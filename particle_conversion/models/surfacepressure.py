# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

import numpy as np
from scipy.optimize import bisect
from math import pow as pow

#============#
#============#

gas_const = 8.31446   # gas constant, J/mol.K

# 100% CO2
kin = [10000.,283000,1]   # A0 in , Ea in J/mol, n in -
C1  = 1e-13               # kg/m.s.Pa.K^0.75
eta = 1                   # -
pg  = 2e6                 # Pa


# Particle
dp    = 200.0e-6  # m
Tp    = 1673.15   # K
rhop0 = 800.0     # kg/m3
dp0   = 200.0e-6  # m
wc0   = 0.80      # -

#============#
#============#

def surface():
    return np.pi * pow(dp,2)

def volume():
    return np.pi/6 * pow(dp,3)


#============#
#============#

def diffusion_rate_coeff():
    #>> diffusion rate constant, kg/s.m2.Pa
    kdiff = C1 * (pow(Tp,0.75)/dp)
    return kdiff

#============#
#============#

def kinetic_rate_coeff():
    #>> kinetic reaction coefficient
    #>> 1/S * dm/dt , kg/(s m^2 Pa^n)
    #>> 1/s.Pa^n = 1/s.Pa^n * exp[ (J/mol) / (J/mol.K * K) ]
    kinConst = kin[0] * np.exp(-1. * kin[1]/(gas_const * Tp))

    return kinConst

#============#
#============#

def Func(PS):

    #--------------------------------------
    #-- Diffusion rate == Kinetic rate

    k_diff = diffusion_rate_coeff()
    k_rate = kinetic_rate_coeff()

    #>> Fluent approach
    #>> diffRate = 1/S * dm/dt [=] kg/m2.s = kg/m2.Pa.s * Pa
    #>> kinRate  = 1/S * dm/dt [=] kg/m2.s = 1/s.Pa^n * Pa^n * kg/m3 * m
    diffRate = k_diff * (pg - PS)
    kinRate  = k_rate * pow(PS, kin[2]) 

    #--------------------------------------
    diff = kinRate - diffRate

    return diff

#============#
#============#

def Pi_surface():
    min_PS=0
    max_PS=pg
    ps = bisect(Func, min_PS, max_PS,
                xtol=1e-12, rtol=1e-15, maxiter=100, full_output=False, disp=True)

    return ps

#================#
#== START HERE ==#
#================#

pS        = Pi_surface()
kin_rate  = eta * surface() * kinetic_rate_coeff() * pow(pS, kin[2])
diff_rate = eta * surface() * diffusion_rate_coeff() * (pg - pS)
mp0        = rhop0 * volume()
mC0        = rhop0 * volume() * wc0

dmdt_check = surface() * eta * pg * kinetic_rate_coeff()*diffusion_rate_coeff()/(kinetic_rate_coeff()+diffusion_rate_coeff())

print('Kin. rate, kg/s:',kin_rate)
print('Diff. rate, kg/s:',diff_rate)
print('PS, Pa:',pS)
print('Check-rate, kg/s:',dmdt_check)
print('dXdt Kin, 1/s:',kin_rate/mC0)
print('dXdt diff, 1/s:',diff_rate/mC0)
print('--END--')

