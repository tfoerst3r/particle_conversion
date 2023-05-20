# SPDX-FileCopyrightText: 2023 Thomas Förster
#
# SPDX-License-Identifier: MIT

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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    #-- add arguements to the CLI
    parser.add_argument(
        '--config', 
        type=str,
        required=True,
        help='Configuration file (.toml) - contains model and settings for running the calculation'
    )
    
    parser.add_argument(
        '--output',
        type = str,
        help = 'Output where calculation results are stored.',
        default = './output_default.csv',
    ) 
 
    config_file = parser.parse_args().config
    #default = './output_' + Path(parser.parse_args().config).stem + '.csv',
    #output_file_default = './output_' + Path(config).stem + '.csv'
   
    #-- check if file exists
    if config_file and (not Path(config_file).is_file()):
        raise ImportError("Given config does not exists.")

    return parser.parse_args()


def main():
    """
    Entrypoint of the command line interface for interacting with particle conversion model
    @return: exit code
    """

    # Parse cli arguments
    args = _parse_args()
    settings = load_config(args.config)
    output = args.output

    print('===START===')
    calculation = ZEROreactor(settings,output)
    calculation.calc()
    
    # important call for analysis_zeroD.py to find the EOF
    print('====EOF====')

    sys.exit(0)

if __name__ == "__main__":
    main()
