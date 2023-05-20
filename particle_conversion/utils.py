# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

"""
SPDX-FileCopyrightText: 
SPDX-License-Identifier: MIT


This file provides supporting functions:
- loading parameters from a config file
- saving the multi_dim model output
"""

import tomli
import numpy as np

def load_config(config):
    """
    Load configuration file (in toml-format) with model-settings
    @param config: configuration-file with different model-settings (see documentation for details)
    @return: dict with settings
    """
    with open(config, 'rb') as f:
        config_output = tomli.load(f)
    return config_output


## ======== ##
## ======== ##

def func_inter(timex,data):
    return float(np.interp(timex,data[0], data[1]))

## ======== ##
## ======== ##

def normalize_dict(inquiry: dict):

    if not inquiry:
        print('>>> Warning! Dictionary was empty! <<<')
        return {}

    inquiry_new = {}

    try:
        for i,value in enumerate(dict(inquiry)):
            inquiry_new[value] = inquiry[value]
    except TypeError:
        error_string = 'Not a dictionary type class!'
        raise TypeError(error_string)


    for i,(valA,valB) in enumerate(inquiry_new.items()):
        if type(valB)!=float and type(valB)!=int:
            raise ValueError(valB,'is not a number')
        if float(valB) < 0:
            print ('Input is negative. They are ignored!')
            continue

    sum = 0
    for i,(valA,valB) in enumerate(inquiry_new.items()):
        if valB < 0:
            valB = 0
        sum += valB

    for i,(valA,valB) in enumerate(inquiry_new.items()):
        inquiry_new[valA] = valB/sum

    return inquiry_new
