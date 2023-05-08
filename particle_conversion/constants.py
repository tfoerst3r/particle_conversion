# Default Values
_defaults = {
    'numerical': {
        'calctime': 2.0, #s
        'steps': 200,
        'conversion_init': 0.0,
        'conversion_final': 0.999,
        'mechanism': 'gri30.cti'
    },
    'reactor': {
        'default_temp': 1673.15, #K
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
            'particle_init_conversion': 0.0
        }
    }
}
