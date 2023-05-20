<!--
SPDX-FileCopyrightText: 2023 Thomas Förster

SPDX-License-Identifier: CC-BY-4.0
-->

# Configuration File 

## Basic structure

The configuration for the modelling is stored in the `toml` format. This means values are separated by an equal sign with in a heading or subheading. **Be aware that all inputs are in SI-Units (K,kg,m,s)**


## Collection of all available options

```
├── reactor
│   ├── fallback_gas
│   ├── fallback_temp
│   └── reactor_pressure
├── particle
│   ├── properties
│   │   ├── particle_diameter_init
│   │   ├── density_particle_apparent_init
│   │   ├── density_ash_true
│   │   └── density_char_true
│   ├── conversion_model
│   │   └── model
│   ├── proxymate_analysis
│   │   ├── fixedcarbon_init
│   │   └── ash_init
│   └── reaction
│       ├── A
│       ├── Ea
│       ├── n
│       ├── psi_R0
│       ├── A_R1
│       ├── Ea_R1
│       ├── n_R1
│       ├── psi_R1
│       ├── A_R2
│       ├── Ea_R2
│       ├── n_R2
│       ├── psi_R2
│       ├── A_R3
│       ├── Ea_R3
│       ├── n_R3
│       ├── psi_R3
│       └── tortuosity-factor
├── numerical
│   ├── conversion_init
│   ├── conversion_final
│   ├── calctime
│   └── steps
└── environment
    ├── co  
    ├── co2 
    ├── h2o 
    ├── o2  
    ├── n2 
    ├── ch4
    ├── h2 
    ├── Re 
    ├── Tg 
    ├── Tp 
    └── law
```

### `[reactor]`

- `fallback_gas = CO2/O2/H2O [default: CO2]`
    : Default gas composition if no `[environment]` gas composition is given.

- `fallback_temp = <value> [default: 1273]`
    : constant reactor temperature, $\mathrm{K}$

- `reactor_pressure = <value> [default: 1e5]`
    : pressure of the reactor where the conversion takes place, $\mathrm{Pa}$

### `[particle]`

- `[particle.properties]`

    - `particle_diameter_init = <value>`
        : particle diameter, $\mathrm{m}$

    - `density_particle_apparent_init = <value>`
        : apparent (incl. porous) density of the particle, $\mathrm{kg/m^3}$

    - `density_ash_true = <value>`
        : true density of ash

    - `density_char_true = <value>`
        : true density of the initial char (incl. carbon and ash), $\mathrm{kg/m^3}$
    
- `[particle.conversion_model]`
    
    - `model = SDM/SPM/RPM/SDMp/SPMp`
        : selection of the used conversion routines
    
- `[particle.proxymate_analysis]`
    
    - `fixedcarbon_init = <value>`
        : Fixed carbon fraction, $\mathrm{kg/kg}$
    
    - `ash_init = <value>`
        : Ash content, stays constant during the conversion, $\mathrm{kg/kg}$

- `[particle.reaction]`

    - `A, Ea, n, psi_R0`
        : fall back parameters for all Arrhenius equation

    - `A_R1, Ea_R1, n_R1, psi_R1`
        : Constants for Arrhenius equation for the reaction $\mathrm{C + O_2}$

    - `A_R2, Ea_R2, n_R2, psi_R2`
        : Constants for Arrhenius equation for the reaction $\mathrm{C + CO_2}$
    
    - `A_R3, Ea_R3, n_R3, psi_R3`
        : Constants for Arrhenius equation for the reaction $\mathrm{C + H_2O}$

    - `tortuosity-factor`
        : A model parameter to estimate the development of the effective diffusion, see [following publication][vascellari_2015].


### `[numerical]`

Either the conversion goal or the calculation time (`calctime`) is met. 

- `conversion_init = <value> [default: 0.0]`
  : particle are introduced with the initial conversion. This can be useful if you have a rig which uses repeating runs, $\mathrm{kg/kg}$

- `conversion_final = <value> [default: 0.999]`
  : particle final conversion, $\mathrm{kg/kg}$
.

- `calctime = <value> [default: 5.0]`
  : calculation time, $\mathrm{s}$

- `steps = <value> [default: 200]`
  : Time step resolution. But the solver determines dynamical the step size after the initial run
    

### `[environment]`
   
You can define transient values for the gas composition (volume/mole fraction), particle and gas temperature, as well as Reynolds number and the custom ANSYS fluent input, called particle law. The particle law is important to determine the time the particle is in the conversion state, instead for example drying. The file output format is not important, only the content will be evaluated and read. File endings like `*.out`, `*.dat`, `*.xy`, etc. can be used..

The file(s) should contain two columns. The first represents the time and the second column represents the desired value.

Values will be extrapolated for each time step. Further, if the values in the file do not cover the desired time frame, it will be kept constant based on the last known value. In addition, in each time step the gas composition will be normalized!

**Example**

```
[environment]
co2 = './data/particle_1234/gas_co2.out'
o2  = './data/particle_1234/gas_o2.out'
```

`gas_co2.out`
```
0.00 0.75
1.00 0.40
```

`gas_o2.out`
```
0.00 0.25
2.00 0.60
```

# Example

## Basic Example 

In the following is a basic example with the name, e.g `base.toml`. The same file can be found in the `./config/` folder.

```
[reactor]
fallback_gas = 'O2'
reactor_pressure = 5e5

[particle.properties]
particle_diameter_init = 1e-3
density_particle_apparent_init = 783.6
density_ash_true = 1900
density_carbon_true = 1900

[particle.conversion_model]
model = 'SDM'

[particle.proxymate_analysis]
fixedcarbon_init= 0.8342
ash_init = 0.1658

[particle.reaction]
#Fallback
A      = 1
Ea     = 1
n      = 1
psi    = 1
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
tortuosity-factor = 0.1

[numerical]
calctime = 3.0

[environment]

```

## Advanced example

In the following is a more advanced setup with the name, e.g `inj.toml`. The same file can be found in the `./config/` folder.

```
[reactor]
###################################################
fallback_temp = 1643 #K
reactor_pressure = 20e5 #Pa
###################################################

[particle.properties]
###################################################
particle_diameter_init = 76e-6 #m
density_particle_apparent_init = 727 #kg/m3
density_ash_true = 2580 #kg/m3
density_char_true = 2270 #kg/m3
###################################################

[particle.conversion_model]
###################################################
model = 'SDMp'
###################################################

[particle.proxymate_analysis]
###################################################
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
n_R2   = 0.4
psi_R2 = 0.0
#.. Reaction 3: C + H2O -> CO + H2
A_R3   = 2.04054428e+06  
Ea_R3  = 2.28e+5
n_R3   = 0.4
psi_R3 = 0.0
#.. structural char pore parameter
tortuosity-factor = 0.08
###################################################

[numerical]
###################################################
# boundary conditions for conversion
#conversion_init = 0.0 #kg/kg
#conversion_final = 0.99 #kg/kg
calctime = 6 #seconds
steps = 100
###################################################

[environment]
###################################################
co   = './docs/exampledata/inj_gas_co.xy'
co2  = './docs/exampledata/inj_gas_co2.xy'
h2o  = './docs/exampledata/inj_gas_h2o.xy'
o2   = './docs/exampledata/inj_gas_o2.xy' 
n2   = './docs/exampledata/inj_gas_n2.xy'
ch4  = './docs/exampledata/inj_gas_ch4.xy'
h2   = './docs/exampledata/inj_gas_h2.xy'
Re   = './docs/exampledata/inj_p_Re.xy' 
Tg   = './docs/exampledata/inj_gas_temp.xy' 
Tp   = './docs/exampledata/inj_p_temp.xy' 
law  = './docs/exampledata/inj_p_law.xy' 
###################################################
```


# Notes

## Particle properties

- Apparent particle density tend to be in the range between $400 .. 800 \mathrm{kg/m^3}$. 
- The true density of char was estimated by Klose (Klose1988) between $1850 .. 1950 \mathrm{kg/m^3}$.


## Reaction mechanism

Be aware that each models parameters are fit accordingly. Meaning the Arrhenius parameters should correspond with model ($A_{0,i}$, $E_{a,i}$, $n_i$)

Kinetic is defined via the following equation:

$$\begin{aligned}
r_X
 &= \frac{\mathrm{d} X_C}{\mathrm{d} t}
  = -\frac{1}{m_{C,0}}\frac{\mathrm{d} m_C}{\mathrm{d} t} \\
 &= \sum\limits_{i} \eta_{i} \cdot A_{0,i} \cdot \exp\left( -\frac{E_{a,i}}{R_u \cdot T_p} \right) \cdot p_i^{n_i} \cdot (1-X_C) \cdot \frac{s_r}{s_{r,0}}
\end{aligned}$$

Be aware that the effectiveness factor $\eta_i$ is dependent on a variety of parameters.


Model assumptions for $\tfrac{s_r}{s_{r,0}}$

| Abbr.           | Model | Name                | Definition                             |
| -----           | ----- | -----               | -----                                  |
| RPM             | RPM   | Random pore         | $\sqrt{1-\varPsi_i \cdot \ln (1-X_C)}$ |
| SDM             | SDM   | Shrinking density   | $(1-X_\mathrm{C})^{-1}$                |
| SPM             | SPM   | Shrinking particle  | $(1-X_\mathrm{C})^{-\tfrac 13}$        |
| SDM<sub>p</sub> | SDMp  | Shrinking density*  | $(1-X_p)^{-1}$                         |
| SPM<sub>p</sub> | SPMp  | Shrinking particle* | $(1-X_p)^{-\tfrac 13}$                 |



## Diffusion

Properties, like the mixed diffusion coefficient, are determined by the GRI 3.0 mechanism in the software package **cantera**. No cantera solves are utilized.

<!-- 
#.. constant porosity during conversion
adaptionfactor = 1
-->

<h1>Literature</h1>

- **Klose1988**

    Klose, E.; Kuchling, T. & Born, M. <br>
    Brennstofftechnische Arbeitsmappe  <br>
    TU Freiberg, IEC, TU Freiberg, IEC, 1988

<!---- Literature ---->
[vascellari_2015]: https://doi.org/10.1016/j.fuel.2015.01.038

<!---- Links ---->

