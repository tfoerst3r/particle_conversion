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
: Default gas composition. It will be overwritten by `[environment]` gas compositions.

`default_temp = 1500 [default: 1673]`
: constant reactor temperature

`reactor_pressure`
: 






Possible adaption:

- self.dXdt = 0.5 * self.dXdt
    - bader pS approach
    - fluent pS approach

- self.Deff = 0.5 * self.Deff




<!-- 
#.. constant porosity during conversion
adaptionfactor = 1
-->


