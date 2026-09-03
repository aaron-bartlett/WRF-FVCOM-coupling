#!/usr/bin/env python

import os
import sys

znu = 0.9975             # Mass-level eta (or sigma)
znw = 1.0                # w-level eta (or sigma)
rdnw = -200.             # Eta-thickness
mu_pert = 1473.773       # Perturbation dry air mass in column
mu_base = 79728.1        # Base-state dry air mass in column
p_top = 5000.            # Model top pressure
p_pert = 1659.141        # Perturbation hydrostatic pressure
p_base = 84528.78        # Base-state hydrostatic pressure
p_hyd = 86187.92         # Hydrostatic pressure
p_sfc = 86387.24         # Surface pressure
theta_bar = 300.         # Base-state potential temperature
theta_pert = -0.4732056  # Perturbation potential temperature
geopot_bar_1 = 84528.78
geopot_pert_1 = 1659.141
geopot_bar_2 =
geopot_pert_2 =

# For equation of state

def calc_pressure_hyd(eta,mu_base,mu_pert,p_top):
    return eta*(mu_base + mu_pert) + p_top

def calc_alpha(phi1,phi2,dznw,mu_total):
    return -1*(theta1-theta2)/dznw/mu_total

def calc_pressure(theta,alpha):
    P0 = 100000.
    RD = 287.
    CP = 7.*RD*0.5
    GAMMA = CP/(CP-RD)
    return P0*(RD*theta/(P0*alpha))**GAMMA

p_calc = calc_pressure_hyd(znw,mu_base,mu_pert,p_top)
print "znw = ",znw
print "mu_base = ",mu_base
print "mu_pert = ",mu_pert
print "mu_total = ",mu_base+mu_pert
print "p_calc = ",p_calc
print "p_pert+p_base = ",p_pert+p_base
print "p_hyd = ",p_hyd
print "p_sfc = ",p_sfc

