"""
SPDX-FileCopyrightText: 
SPDX-License-Identifier: MIT


This file provides supporting functions:
- loading parameters from a config file
- saving the multi_dim model output
"""

import tomli

def load_config(config_file):
    """
    Load configuration file (in toml-format) with model-settings
    @param config_file: configuration-file with different model-settings (see documentation for details)
    @return: dict with settings
    """
    with open(config_file, 'rb') as f:
        config = tomli.load(f)
    return config
