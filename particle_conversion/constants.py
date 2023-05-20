# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

# Default Values
_defaults = {
    'numerical': {
        'calctime': 5.0, #s
        'steps': 200,
        'conversion_init': 0.0,
        'conversion_final': 0.999,
        'mechanism': 'gri30.cti'
    },
    'reactor': {
        'default_temp': 1273.15, #K
        'fallback_gas': 'CO2',
        'reactor_pressure': 1e5 #Pa
    },
    'environment': {
        'co'   : '',
        'co2'  : '',
        'h2o'  : '',
        'o2'   : '',
        'n2'   : '',
        'ch4'  : '',
        'h2'   : '',
        'Re'   : '',
        'Tg'   : '',
        'Tp'   : '',
        'law'  : '',
    },
    'particle':{
        'properties': {
            'particle_init_conversion': 0.0,
            'mechanism': 'gri30.cti'
        }
    }
}
