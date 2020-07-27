#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:48:15 2020

@author: Jean-Philippe Kuntzer
"""

import logging
import json
import os
import aeroframe_2.fileio.settings as settings

logging.basicConfig(level=logging.DEBUG)
__prog_name__ = "Aeroframe2.0"
logger = logging.getLogger(__prog_name__+"."+__name__)



def main():
    logger.info(f"{__prog_name__} started")
    
    # dir_path = os.path.dirname(os.path.realpath(__file__))
    # logger.debug(f"Current directory is: {dir_path}")
    
    # for construction purposes
    dir_path = "/home/cfse2/Documents/aeroframe_2.0/test/static/aeroframe2.0_case1.json"
    logger.debug(f"Settings file is located: {dir_path}")
    
    with open(dir_path) as json_file:
        simu_settings = settings(json.load(json_file))
    # settings.verify()
    # CFD_solver = simu_settings["CFD_solver"]
    # logger.debug(f"simulation settings raw: \n{simu_settings}")
    # logger.debug(f"Settings data: \n{type(CFD_solver)}")
    # logger.debug(f"Settings data \n: {data}")


if __name__ is "__main__":
    main()