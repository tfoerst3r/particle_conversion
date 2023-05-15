"""
SPDX-FileCopyrightText: 
SPDX-License-Identifier: MIT


This file provides supporting functions:
- loading parameters from a config file
- saving the multi_dim model output
"""

import tomli

def load_config(config):
    """
    Load configuration file (in toml-format) with model-settings
    @param config: configuration-file with different model-settings (see documentation for details)
    @return: dict with settings
    """
    with open(config, 'rb') as f:
        config_output = tomli.load(f)
    return config_output
