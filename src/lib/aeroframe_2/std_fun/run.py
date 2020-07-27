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
import pickle
# import SU2_CFD

logging.basicConfig(level=logging.DEBUG)
__prog_name__ = "Aeroframe2.0"
logger = logging.getLogger(__prog_name__+"."+__name__)

def save_to_pkl(path,act_dea,lattice,vlmdata):
    if act_dea == "activated":
        name_lattice = "/lattice_defActivated.pkl"
        name_vlmdata = "data_defActivated.pkl"
    elif act_dea == "deactivated":
        name_lattice = "/lattice_defDeactivated.pkl"
        name_vlmdata = "data_defDeactivated.pkl"
    else:
        logger.error("activation or diactivation not specified")
    with open(path + name_lattice, "wb") as la:
        var = [lattice.p, # 0
               lattice.v, # 1
               lattice.c, # 2
               lattice.n, # 3
               lattice.a, # 4
               lattice.bound_leg_midpoints] # 5
        pickle.dump(var, la)
    la.close()
    
    with open(path + name_vlmdata, "wb") as d:
        pickle.dump(vlmdata, d)
    d.close()


def deform_mesh(settings,lattice,file):
    logger.info(settings.settings["deformation"])
    if settings.settings["deformation"]:
        # Deforms the mesh and uploads the deformed one into the code
        logger.info("===== Mesh deformation function activated =====")
        mesh_def = aeroDef.Mesh_Def(lattice,file)
        mesh_def.deformation(settings)
        lattice.p = mesh_def.f_p
        lattice.v = mesh_def.f_v  # turns everything down
        lattice.c = mesh_def.f_c
        lattice.bound_leg_midpoints = mesh_def.f_b  # turns everything up
        lattice.n = mesh_def.f_n
        lattice.a = mesh_def.f_a
    else:
        logger.info("===== Mesh deformation function deactivated =====")
    return lattice


def getSettings(inputFileAndPath):
    try:
        with open(inputFileAndPath) as f:
            settings = json.load(f)
        f.close()
    except FileNotFoundError:
        logger.error("<setting file>.json not found for aeroframe")
        raise FileNotFoundError
    return settings


def standard_run(args):
    # OK TODO: make it work from the command line
    # TODO do a better job with the settings files workflow
    # TODO get json file
    # TODO get all json variables
    # TODO checks all JSON variables
    # TODO checks folder structure and prints the potential errors
    # TODO try except error management with setting file
    # TODO add test to see if file exists
    # TODO add test for filepath

    # Starts simulation
    logger.info(f"{__prog_name__} started")

    # Simulation input file
    logger.debug("Setting file path:\n" + args.cwd+"/"+args.run)
    settingsFileAndPath = args.cwd + "/" + args.run
    aeroframe_2_settings = getSettings(settingsFileAndPath)

    # simulation.verify()
    logger.debug(aeroframe_2_settings)
    if aeroframe_2_settings["CFD_solver"] == "Pytornado":
        # Command line simulation
        pytornado_settings_file = args.cwd + "/CFD/settings/" + aeroframe_2_settings["CFD_settings_file"]
        def_file_path = aeroframe_2_settings["deformation_file"]
        dir_path = args.cwd

        # Buids CFD mesh
        lattice, vlmdata, settings, aircraft, cur_state, state = cfd.meshing(args,pytornado_settings_file)
        logger.debug("Meshing done")

        # Deforms CFD mesh
        file_path = dir_path + "/" + def_file_path
        logger.debug(file_path)
        # Activates or diactivates function
        if aeroframe_2_settings["deformation_activation"]:
            deform_mesh(settings,lattice,file_path)
        
        # Computes CFD solution
        cfd.solver(lattice, vlmdata, settings, aircraft, cur_state, state)

        if aeroframe_2_settings["deformation_activation"]:
            if aeroframe_2_settings["save_pkl"]:
                save_to_pkl(dir_path, "activated", lattice, vlmdata)
        else:
            if aeroframe_2_settings["save_pkl"]:
                save_to_pkl(dir_path, "deactivated", lattice, vlmdata)

    # TODO cfd Mesh
    # TODO structure mesh


# if __name__ == "__main__":
#     main()


# -TODO get p
# -TODO get v
# -TODO get c
# -TODO get b
# -TODO get n
# -TODO get a
# -TODO get what is in between. it looks like something wrong lives here
# -TODO get RHS
# -TODO get Downwash
