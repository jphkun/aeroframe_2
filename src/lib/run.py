#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:48:15 2020

@author: Jean-Philippe Kuntzer
"""

import logging
import json
import argparse
import aeroframe_2.fileio.settings as Settings
import aeroframe_2.deformation.functions as aeroDef
import pytornado.stdfun.run as cfd
import pytornado.fileio as io
import os
# import SU2_CFD

logging.basicConfig(level=logging.DEBUG)
__prog_name__ = "Aeroframe2.0"
logger = logging.getLogger(__prog_name__+"."+__name__)


def main():
    # Starts simulation
    logger.info(f"{__prog_name__} started")

    # for construction purposes
    dir_path = "/home/cfse2/Documents/aeroframe_2/test/static/"
    name = "aeroframe2.0_case1.json"
    aeroframe_settings_path = dir_path + name
    logger.debug(f"Settings file is located: {aeroframe_settings_path}")
    with open(aeroframe_settings_path) as json_file:
        dictionary = json.load(json_file)
        simulation = Settings.settings(dictionary)
    # Checks if JSON corresponds to what is implemented in this code
    simulation.verify()
    if simulation.CFD_solver == "Pytornado":
        settings_file = dir_path + "/CFD/settings/Optimale.json"
        
        with open(settings_file, "r") as fp:
            settings = json.load(fp)
        # Command line simulation
        parser = argparse.ArgumentParser(description='Process some integers.')
        args = parser.parse_args()
        args.clean = True
        args.clean_only = False
        args.cpacs2json = None
        args.debug = False
        args.list_example_from_db = False
        args.make_example = False
        args.make_example_from_db = None
        args.run = 'settings/Optimale.json'
        args.verbose = True
        os.chdir(dir_path + "/CFD")
        
        # Buids CFD mesh
        lattice, vlmdata, settings, aircraft, cur_state = cfd.meshing(args)
        
        # Deforms CFD mesh
        logger.info(settings.settings["deformation"])
        if settings.settings["deformation"]:
            # Deforms the mesh and uploads the deformed one into the code
            logger.info("===== Mesh deformation function activated =====")
            mesh_def = aeroDef.Mesh_Def(lattice)
            mesh_def.deformation(settings)
            lattice.p = mesh_def.f_p
            lattice.v = mesh_def.f_v  # turns everything down
            lattice.c = mesh_def.f_c
            lattice.bound_leg_midpoints = mesh_def.f_b  # turns everything up
            lattice.n = mesh_def.f_n
            lattice.a = mesh_def.f_a
        else:
            logger.info("===== Mesh deformation function deactivated =====")
        
        # Computes CFD solution
        cfd.solver(lattice, vlmdata, settings, aircraft, cur_state)
        
    # TODO cfd Mesh
    # TODO VLM Mesh
    # TODO structure mesh


if __name__ == "__main__":
    main()
