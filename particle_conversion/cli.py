"""
SPDX-FileCopyrightText: 
SPDX-License-Identifier: 

Functions in this file take information from a command line user input 
and parameters from config .toml file 
"""

import argparse
import sys
from pathlib import Path

from particle_conversion.utils import load_config
from particle_conversion.models.zeroDim_reactor import ZEROreactor

def _parse_args():
    """
    Parsing the commandline arguments

    @return: the parsed arguments (argparse.Namespace)
    """

    #-- create argparse object
    parser = argparse.ArgumentParser(
        description='Estimates conversion of a particle',
        formatter_class=argparse.RawTextHelpFormatter
    )

    #-- add arguements to the CLI
    parser.add_argument(
        '--config_file', 
        type=str,
        help='Configuration file (.toml) - contains model and settings for running the calculation'
    )

    #-- check if file exists
    config_file = parser.parse_args().config_file
    if not Path(config_file).is_file():
        raise ImportError("Given config does not exists.")

    parser.add_argument(
        '--output',
        default='./output_' + Path(config_file).stem + '.csv', 
        type=str,
        help=   'Output where calculation results are stored. '
                'Destination will be created if it does not exist, existing files will be overwritten.\n'
                'For models "base", the output is a CSV-File\n'
    ) 

    return parser.parse_args()


def main():
    """
    Entrypoint of the command line interface for interacting with particle conversion model
    @return: exit code
    """

    # Parse cli arguments
    args = _parse_args()
    settings = load_config(args.config_file)
    output = args.output


    # TODO: check for input errors, maybe already in load_config()
    # check_input(settings)

    calculation = ZEROreactor(settings,output)
    calculation.calc()
    
    # important call for analysis_zeroD.py to find the EOF
    print('===EOF===')

    #output = pathlib.Path(args.output)

    #models_output_dir = ['multi_dim']
    #models_output_file = ['base', 'one_dim', 'one_temp', 'one_weight']

    #if model in models_output_dir and output.is_file():
    #    sys.exit(f'Error: Output must be directory for models {models_output_dir}')
    #if model in models_output_file and output.is_dir():
    #    sys.exit(f'Error: Output must be a file for models {models_output_file}')
   
    sys.exit(0)

