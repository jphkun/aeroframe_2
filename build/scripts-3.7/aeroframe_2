#!python
# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------
# Copyright 2020 Airinnova AB and the PyTornado authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ----------------------------------------------------------------------

# Author: Aaron Dettmann

"""
Command line interface
"""

import argparse
# from pathlib import Path
# import sys

from aeroframe_2.__version__ import __version__
import aeroframe_2.stdfun.run as stdrun


def main():
    """
    Command line interface
    """

    HELP_RUN = f"Settings file (entry point for {stdrun.__prog_name__})"
    # TODO: add a make example case

    parser = argparse.ArgumentParser(prog=f'{stdrun.__prog_name__} {__version__}')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-r', '--run',
                       metavar='<Settings file>',
                       type=str,
                       help=HELP_RUN)
    # TODO: add a make example command

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-v', '--verbose', action='store_true')
    group.add_argument('-d', '--debug', action='store_true')
    group.add_argument('-q', '--quiet', action='store_true')

    group = parser.add_mutually_exclusive_group()
    args = parser.parse_args()

    if args.run:
        stdrun.standard_run(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
