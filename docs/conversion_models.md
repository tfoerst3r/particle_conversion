# Configuration File 

## Basic structure

The configuration for the modelling is stored in the `toml` format. This means values are separated by an equal sign with in a heading or subheading. **Be aware that all inputs are in SI-Units (K,kg,m,s)**

**Example**

```
[reactor]
reactor_pressure = 20e5 #Pa

[particle.conversion_model]
model = SDM
```

## Collection of all available options

```
├── reactor
│   ├── fallback_gas
│   ├── default_temp
│   └── reactor_pressure
├── particle
│   ├── properties
│   │   ├── particle_diameter_init
│   │   ├── density_particle_apparent_init
│   │   ├── density_ash_true
│   │   ├── density_carbon_true
│   │   └── reactive_surface_init
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
│   ├── mechanism
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

`fallback_gas = CO2/O2/H2O [default: CO2]`
- Default gas composition if no `[environment]` gas composition is given.

`default_temp = 1500 [default: 1673]`
- constant reactor temperature, $\mathrm{K}$

`reactor_pressure`
- pressure of the reactor where the conversion takes place, $\mathrm{Pa}$

### `[particle]`

- `[particle.properties]`

    - particle_diameter_init



│   │   ├── density_particle_apparent_init
│   │   ├── density_ash_true
│   │   ├── density_carbon_true
│   │   └── reactive_surface_init


### `[environment]`
   
You can define transient behaviours of the gas composition (volume/mole fraction), particle and gas temperature, as well as Reynolds number and the custom ANSYS fluent input (particle law). The particle law is important to determine the time the particle is in the conversion state, instead of for example drying. The file output format is not important, only the content will be evaluated and read. Formats like `*.out`, `*.dat`, `*.xy`, etc. can be used..

The file(s) should contain a first row which represents the time and the second line represents the desired value.

Values will be extrapolated and kept constant, based on the last known value. Further in each time step the gas composition will be normalized!

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





Possible adaption:

- self.dXdt = 0.5 * self.dXdt
    - bader pS approach
    - fluent pS approach

- self.Deff = 0.5 * self.Deff




<!-- 
#.. constant porosity during conversion
adaptionfactor = 1
-->


