# Tray

Possible adaption:

- self.dXdt = 0.5 * self.dXdt
    - bader pS approach
    - fluent pS approach

- self.Deff = 0.5 * self.Deff

```
# SPDX-FileCopyrightText: 
# SPDX-License-Identifier: CC0-1.0

[reactor]
###################################################
#fallback_gas = 'CO2'
#default_temp = 1673.15 #K
reactor_pressure = 20e5 #Pa
###################################################

[particle.properties]
###################################################
#.. - initial apparent particle density, Bader: 538.1 kg/s
#.. - ash density, kg/m3
#.. - true density of char, klose1988 1850..1950 kg/m3
#.. - absolute internal particle surface (initial), m^2
particle_diameter_init = 1e-3 #m
density_particle_apparent_init = 783.6 #kg/m3
density_ash_true = 1900 #kg/m3
density_carbon_true = 1900 #kg/m3
reactive_surface_init=1 #m2
###################################################

[particle.conversion_model]
###################################################
#.. parameters for the apparent density of char fraction, 
#.. if 1.0, shrinking density is applied -> regime I
#.. if 0.0, shrinking particle/rpm -> regime II/III
#.. else see input alpha
#.. Surface development:
#.. - Shrinking Density Model, SDM  --> alpha=1
#.. - Shrinking Particle Model, SPM --> alpha=0
#.. - Random Pore Model, RPM        --> alpha=1
#.. !!! Model should corresponds to the kinetics !!!
model = 'SDM'
###################################################

[particle.proxymate_analysis]
###################################################
#.. Proximate Analysis of char, kg/kg:
#.. - fixed carbon
#.. - ash fraction
fixedcarbon_init= 0.8342
ash_init = 0.1658
###################################################

[particle.reaction]
###################################################
#.. Dummys
A=1.38e+1
Ea=1.38e+8
n=0.8
psi_R0=0.0
#.. Reaction 1: C + 0.5 O2 -> CO 
A_R1   = 4.95392582e+02
Ea_R1  = 1.38e+5
n_R1   = 0.8
psi_R1 = 0.0
#.. Reaction 2: C + CO2 -> 2 CO
A_R2   = 5.87769389e+07
Ea_R2  = 2.89e+5
n_R2   = 1.0
psi_R2 = 0.0
#.. Reaction 3: C + H2O -> CO + H2
A_R3   = 2.04054428e+06  
Ea_R3  = 2.28e+5
n_R3   = 0.4
psi_R3 = 0.0
#.. structural char pore parameter
#.. equal to tortuosity * constriction-factor; 1/7.13 1/8
tortuosity-factor = 0.1
#.. constant porosity during conversion
adaptionfactor = 1
###################################################

[numerical]
###################################################
# dpm model:
# final particle conversion until the model is solved
conversion_init = 0.0 #kg/kg
conversion_final = 0.98
#mechanism = 'gri30.cti'
calctime = 5 #seconds
steps = 200
###################################################

[environment]
###################################################
#co   = ''
#co2  = ''
#h2o  = ''
#o2   = ''
#n2   = ''
#ch4  = ''
#h2   = ''
#Re   = ''
#Tg   = ''
#Tp   = ''
#law  = ''
###################################################
```
