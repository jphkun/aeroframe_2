#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 10:48:15 2020

@author: Jean-Philippe Kuntzer


OK TODO: make it work from the command line
OK TODO get json file
OK TODO get all json variables
OK TODO checks all JSON variables
OK TODO try except error management with setting file
OK TODO add test to see if file exists
OK TODO add test for filepath


TODO: change CFD mesh
TODO: inject the framat setting file to the wrapper
TODO: Compute the error
TODO: File management with csv deformation
TODO: compute CFD moments
TODO:
TODO make a features that makes the program run a save the results and mesh.
TODO checks folder structure and prints the potential errors
TODO do a better job with the settings files workflow
TODO: get CFD forces for csv file
"""

import logging
import json
# import argparse
# import aeroframe_2.fileio.settings as Settings
import aeroframe_2.deformation.functions as aeroDef
import aeroframe_2.deformation.framatNormalVecConverter as normalVecFramat
import aeroframe_2.csdImportGeometry.importGeomerty as importGeomerty
import aeroframe_2.informationTransfer.mapping as mapping
import aeroframe_2.wrappers.framatWrapper as framatWrapper
import pytornado.stdfun.run as cfd
import numpy as np
# import pytornado.fileio as io
# import os
import pickle
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
# import SU2_CFD

# TODO take the information from the args parameters
logging.basicConfig(level=logging.DEBUG)
__prog_name__ = "aeroframe_2"
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
        var = [lattice.p,  # 0
               lattice.v,  # 1
               lattice.c,  # 2
               lattice.n,  # 3
               lattice.a,  # 4
               lattice.bound_leg_midpoints]  # 5
        pickle.dump(var, la)
    la.close()

    with open(path + name_vlmdata, "wb") as d:
        pickle.dump(vlmdata, d)
    d.close()


def deform_mesh(settings,lattice,file,aeroframe_2_settings):
    logger.info(settings.settings["deformation"])
    if settings.settings["deformation"]:
        # Deforms the mesh and uploads the deformed one into the code
        logger.info("===== Mesh deformation function activated =====")
        mesh_def = aeroDef.Mesh_Def(lattice,file,aeroframe_2_settings)
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


def meshComputation(args,aeroframe_2_settings):
    ##########################################################################
    #   Aerodynamic part mesh computation
    ##########################################################################
    # Pytornado meshing
    # All the accepted ways of writing pytornado in the json file
    pytornado = ["pytornado","Pytornado","pyTornado","PyTornado"]
    if aeroframe_2_settings["CFD_solver"] in pytornado:
        # Command line simulation
        pytornado_settings_file = args.cwd + "/CFD/settings/" + aeroframe_2_settings["CFD_settings_file"]
        # dir_path = args.cwd

        # Buids CFD mesh
        lattice, vlmdata, settings, aircraft, cur_state, state = cfd.meshing(args,pytornado_settings_file)
        logger.debug(" VLM meshing done")

    ##########################################################################
    #   Structure part mesh computation
    ##########################################################################
    # Structure beam meshing
    # All the accepted ways of writing FramAT in the json file
    framat = ["framat","Framat","FramAT","framAT"]
    if aeroframe_2_settings["CSD_solver"] in framat:
        """
        Reads the cpacs file and gets the rough geometry. Geometry properties
        are given by the user.
        """
        logger.info("FramAT is chosen as the CSD solver")
        # Importing geometry
        aircraft_path = args.cwd + "/CFD/aircraft/" + aeroframe_2_settings["aircraft_file"]
        csdGeometry = importGeomerty.CsdGeometryImport(aircraft_path,aeroframe_2_settings)
        csdGeometry.getAllPoints()
        # csdGeometry.plotSectionsPoints()

    return lattice, vlmdata, settings, aircraft, cur_state, state, csdGeometry


def meshPreprocessing(args,aeroframe_2_settings,settings,lattice,csdGeometry,vlmdata):
    """
    Deforms the VLM mesh if asked to otherwise compute the transfer matrix.
    """
    framat = ["framat","Framat","FramAT","framAT"]
    # Deforms VLM mesh if the users asks to deform from an input file
    if aeroframe_2_settings["deformation_from_file"]:
        def_file_path = aeroframe_2_settings["deformation_file"]
        file_path = args.cwd + "/" + def_file_path
        logger.debug(file_path)
        deform_mesh(settings,lattice,file_path,aeroframe_2_settings)

    # Computes the radial basis function matrix
    elif aeroframe_2_settings["CSD_solver"] in framat:
        # Computes structural mesh
        csd = framatWrapper.framat(csdGeometry)
        csd.mesh()
        # Computes transformation matrix
        transform = mapping.mapper(lattice,vlmdata,csdGeometry,csd)
        
    return settings, lattice, csd, transform

def csdComputations(transform,csd):
    # apply loads
    # run
    transform.aeroToStructure()
    csd.run(transform)
    pass

def meshDeformation(lattice,csd,transform):
    # U and T separated for each beam
    transform.structureToAero()
    # Separation of the results for each transformation matrix
    return lattice

def cfdComputation(lattice, vlmdata, settings, aircraft, cur_state, state):
    cfd.solver(lattice, vlmdata, settings, aircraft, cur_state, state)
    return lattice, vlmdata


def saveCFDresults(args,settings,lattice,vlmdata):
    if settings["deformation_activation"]:
        if settings["save_pkl"]:
            save_to_pkl(args.cwd, "activated", lattice, vlmdata)
    else:
        if settings["save_pkl"]:
            save_to_pkl(args.cwd, "deactivated", lattice, vlmdata)


def plotMesh(lattice):
    fig = plt.figure("figure 2")
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(lattice.c[:,0],
               lattice.c[:,1],
               lattice.c[:,2],
               label='wing')
    val = 15
    ax.set_xlim(-val,val)
    ax.set_ylim(-val,val)
    ax.set_zlim(-val,val)
    ax.legend()
    plt.show()


def standard_run(args):
    # Starts simulation
    logger.info(f"{__prog_name__} started")

    # Simulation input file
    logger.debug("Setting file path:\n" + args.cwd+"/"+args.run)
    settingsFileAndPath = args.cwd + "/" + args.run
    aeroframe_2_settings = getSettings(settingsFileAndPath)

    # Computes the necessary meshes
    lattice, vlmdata, settings, aircraft, cur_state, state, csdGeometry = meshComputation(args,aeroframe_2_settings)

    # Computes the deformation parameters.
    # Deforms the VLM mesh if the deformation is imported from a file
    settings, lattice, csd, transformOriginal = meshPreprocessing(args,aeroframe_2_settings,settings,lattice,csdGeometry,vlmdata)
    
    # Computes the first CFD solution
    lattice, vlmdata = cfdComputation(lattice, vlmdata, settings, aircraft, cur_state, state)

    # saves results to pkl file if asked.
    saveCFDresults(args,aeroframe_2_settings,lattice,vlmdata)
    for i in range(3):
        # Computes the first CSD solution
        # transformCurrent = copy.deepcopy(transformOriginal)
        transformCurrent = transformOriginal
        csdCurrent = csd
        latticeCurrent = copy.deepcopy(lattice)
        
        
        csdComputations(transformCurrent,csdCurrent)
        
        # Deforms CFD mesh from the results of the structure calculation
        meshDeformation(latticeCurrent,csdCurrent,transformCurrent)
        vlmDef = normalVecFramat.Mesh_Def(latticeCurrent,transformCurrent)
        vlmDef.deformation()
        logger.debug(latticeCurrent.c.shape)
        latticeCurrent.p = np.copy(vlmDef.f_p)
        latticeCurrent.v = np.copy(vlmDef.f_v)
        latticeCurrent.c = np.copy(vlmDef.f_c)
        latticeCurrent.n = np.copy(vlmDef.f_n)
        latticeCurrent.a = np.copy(vlmDef.f_a)
        latticeCurrent.b = np.copy(vlmDef.f_b)
        logger.debug("\n"*20)
        logger.debug(np.max(csd.results.get('tensors').get('comp:U')["uz"]))
        logger.debug(csd.results.get('tensors').get('comp:U')["uz"])
        logger.debug("\n"*20)
        logger.debug(vlmDef.f_c.shape)
        # plotMesh(latticeCurrent)

        # Computes the first CFD solution
        latticeCurrent, vlmdata = cfdComputation(latticeCurrent, vlmdata, settings, aircraft, cur_state, state)
    
        # saves results to pkl file if asked.
        saveCFDresults(args,aeroframe_2_settings,latticeCurrent,vlmdata)
    sys.exit()

        ######################################################################
        # # FramAT part
        # csd = framatWrapper.framat(csdGeometry)
        # csd.csdRun()
        # # csd.plotPoints()
        # logger.debug(csd_mesh.wingsSectionsCenters)
        ######################################################################
        # TODO get path
        # TODO Upload path and aircraft to the mesher
        # csd = mesher()
        ######################################################################

    # TODO 3D cfd Mesh
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
