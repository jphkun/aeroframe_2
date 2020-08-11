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
OK TODO: change CFD mesh
OK TODO: inject the framat setting file to the wrapper
OK TODO: Compute the error
OK TODO: File management with csv deformation
OK TODO do a better job with the settings files workflow
OK TODO: get CFD forces for csv file
OK TODO: compute CFD moments
TODO make a features that makes the program run a save the results and mesh.
TODO checks folder structure and prints the potential errors

"""

import logging
import json
# import argparse
# import aeroframe_2.fileio.settings as Settings
import aeroframe_2.deformation.functions as aeroDef
import aeroframe_2.deformation.framatNormalVecConverter as normalVecFramat
import aeroframe_2.csdImportGeometry.importGeomerty as importGeomerty
import aeroframe_2.informationTransfer.mapping as mapping
import aeroframe_2.informationTransfer.mappingSU2 as mappingSU2
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
import pandas as pd

# import SU2_CFD

# TODO take the information from the args parameters
logging.basicConfig(level=logging.DEBUG)
__prog_name__ = "aeroframe_2"
logger = logging.getLogger(__prog_name__+"."+__name__)


def getSettings(inputFileAndPath):
    try:
        with open(inputFileAndPath) as f:
            settings = json.load(f)
        f.close()
    except FileNotFoundError:
        logger.error("<setting file>.json not found for aeroframe")
        raise FileNotFoundError
    return settings


def pytornadoMeshing(args, aeroframeSettings):
    """
    Calls pytornado functions to build up a lattice (which is basically the
    mesh and our point of interest for aeroelasticity) and all other
    variables/instances needed by pytornado.
    """
    # Virtual command line setting file input path. the "args" contains the
    # rest of the needed information
    pytornadoSttings = args.cwd \
                       + "/CFD/settings/" \
                       + aeroframeSettings["CFD_settings_file"]
    # Buids CFD mesh
    lattice,vlmdata,settings,aircraft,cur_state,state = cfd.meshing(args,pytornadoSttings)
    # Join all the variables for ease of use.
    pytornadoVariables = [lattice, vlmdata, settings, aircraft, cur_state, state]
    logger.info("Pytornado meshing done")

    return pytornadoSttings, pytornadoVariables


def csvDeformation(args,pytornadoSettings,aeroframeSettings,pytornadoVariables):
    """
    Deforms the VLM mesh from a csv file
    """
    # Extract the vlm mesh
    lattice = pytornadoVariables[0]
    # Feeds all the needed settings in order to read the .csv deformation file
    # and then be able to deform the mesh (basically compute the new normnal
    # vector).
    mesh_def = aeroDef.Mesh_Def(args,aeroframeSettings,lattice)

    # Calls the deformation function that computes the new points positions
    mesh_def.deformation()

    # Feeds the result back to an understandable pytornado mesh. The idea is
    # not to touch any unwanted variable, hence the action of feeding back the
    # new point to the old lattice variable.
    lattice.p = mesh_def.f_p
    lattice.v = mesh_def.f_v
    lattice.c = mesh_def.f_c
    lattice.bound_leg_midpoints = mesh_def.f_b
    lattice.n = mesh_def.f_n
    lattice.a = mesh_def.f_a

    # For ease of use and code readability
    pytornadoVariables[0] = lattice
    return pytornadoVariables


def feeder(pytornadoVariables,meshDeformation):
    """
    Transfers the deformed mesh points to pytornadoVariables.
    """
    pytornadoVariables[0].p = meshDeformation.f_p
    pytornadoVariables[0].v = meshDeformation.f_v
    pytornadoVariables[0].c = meshDeformation.f_c
    pytornadoVariables[0].a = meshDeformation.f_a
    pytornadoVariables[0].bound_leg_midpoints = meshDeformation.f_b
    pytornadoVariables[0].n = meshDeformation.f_n

    return pytornadoVariables


def forcesToCsv(args,pytornadoVariables,results):
    """
    Writes the results to a csv file.
    """
    # Computes the path
    path = args.cwd + "/CSD/results/panelwiseForces.csv"
    headers = "x;y;z;fx;fy;fz"
    # Gets simulation values
    panelCoordinates = pytornadoVariables[0].c
    panelFx = results["vlmdata"].panelwise['fx']
    panelFy = results["vlmdata"].panelwise['fy']
    panelFz = results["vlmdata"].panelwise['fz']

    results = np.array([panelCoordinates[:,0],
                        panelCoordinates[:,1],
                        panelCoordinates[:,2],
                        panelFx,
                        panelFy,
                        panelFz])

    np.savetxt(path, results.T, delimiter=';', header=headers)
    logger.info("Simulation finised")


def solverPytornadoCSV(args, aeroframeSettings, acceptedNames):
    """
    The workflow is as follows:
        1) TODO: Verifies that the input file contains all the necessary
           information.
        2) Builds CFD mesh.
        3) Deforms CFD mesh.
        4) Computes CFD problem.
        5) Extract panelwise forces and saves them into a csv file.
    """
    # Step 1) Verifying input
    # TODO: verify all parameters one by one

    # Step 2) pytornado meshing
    pytornadoSettings, pytornadoVariables = pytornadoMeshing(args, aeroframeSettings)

    # Step 3) Deforms CFD mesh
    lattice = pytornadoVariables[0]
    meshDeformation = aeroDef.Mesh_Def(args,aeroframeSettings,lattice)
    meshDeformation.deformation(acceptedNames)
    pytornadoVariables = feeder(pytornadoVariables,meshDeformation)

    # Step 4) Computes the CFD problem
    pytornadoVariables, results = cfd.solver(pytornadoVariables)

    # Step 5) Saves panelwise forces results
    forcesToCsv(args,pytornadoVariables,results)
    logger.info("End of simulation")
    sys.exit()


def solverPytornadoFramat(args, aeroframeSettings, acceptedNames):
    """
    Function called when the user desires to couple the VLM solver Pytornado
    and the structure solver FramAT.

    The workflow is as follows:
        1) TODO: Verifies that the input file contains all the necessary
           information.
        2)  Builds CFD mesh.
        3)  Reads CPACS files and computes the nodes of pseudo 1D structural
            mesh.
        4)  Builds CSD instance in FramAT.
        5)  Computes the transformation matrices
        6)  Computes CFD problem.
        7)  Enters the aeroelastic loop.
        8)  Projects the loads on CSD instance.
        9)  Computes CSD solution
        10) Deforms the CFD mesh.
        11) Computes the norm of the displacement error
        12) Computes new CFD problem.
        13) loops back to point 6) if simulation has not converged.
    """
    # Step 1: File verification
    # TODO: verify all parameters one by one

    # Step 2) pytornado meshing
    pytornadoSettings, pytornadoVariables = pytornadoMeshing(args, aeroframeSettings)

    # Step 3)  Reads CPACS files and computes the nodes of pseudo 1D structural
    #          mesh. Aeroframe function pre-meshes the aircraft to get each
    #          structure node.
    preMeshedStructre = importGeomerty.CsdGeometryImport(args,aeroframeSettings)
    preMeshedStructre.getAllPoints()

    # Step 4) feeds the computed nodes to the structure solver and builds a
    # structure mesh
    csdSolverClassVar = framatWrapper.framat(preMeshedStructre)
    csdSolverClassVar.mesh()

    # Step 5) feeds the computed nodes to a mapping function which computes the
    # tranformation matrices (based on RBF)
    transform = mapping.mapper(pytornadoVariables,preMeshedStructre,csdSolverClassVar)

    # Step 6) Computes CFD problem.
    cfd.solver(pytornadoVariables)
    pytornadoVariablesInit = copy.deepcopy(pytornadoVariables)
    # Setp 7) Aeroelastic loop.
    N = aeroframeSettings["MaxIterationsNumber"]
    i = 0
    maxDisplacement = np.array([0])
    error = []
    absoluteDisplacement = []
    tol = aeroframeSettings["ConvergeanceTolerence"]
    while (i < N):
        # basic user comminication
        logger.debug("aeroelastic loop number: "+str(i))

        # Makes a copy to avoid memory linked mistakes
        transformCurrent = transform
        csdSolverClassVarCurrent = csdSolverClassVar
        # Step 8) Projects the loads on CSD instance.
        transformCurrent.aeroToStructure()

        # Step 9) Compute structure solution
        csdSolverClassVarCurrent.run(transformCurrent)

        # Step 9) deforms the CFD mesh. Computes beam deformation
        latticeCurrent = pytornadoVariablesInit[0]
        meshDeformation = aeroDef.Mesh_Def(args,aeroframeSettings,latticeCurrent)

        # Step 10) computes new aerodynamic points
        transformCurrent.structureToAero()
        meshDeformation.deformation(acceptedNames,transformCurrent)
        pytornadoVariables = feeder(pytornadoVariables,meshDeformation)
        # Step 11) Computes the norm of the displacement error
        # TODO find a clean way to do it
        # TODO one can add a convergence graph
        maxDisplacement = np.append(maxDisplacement, np.max(transform.displacements))
        error.append(np.abs(maxDisplacement[-1] - maxDisplacement[-2]))
        absoluteDisplacement.append(np.abs(maxDisplacement[-1] - maxDisplacement[0]))
        logger.info("Max error between two iteration: "+str(error))
        # Step 12) Deforms the CFD mesh.
        pytornadoVariables = cfd.solver(pytornadoVariables)
        i += 1
        if i == N-1:
            logger.warning("Simulation has reached max number of step,")
            logger.warning("convergeance is yet to determine!")
        if error[-1] <= tol:
            logger.info("Simulation has converged")
            i = N
    # Writes a file which contains the error for each timestep
    N = len(error)
    path = args.cwd + "/results.csv"
    MyFile = open(path,"w")
    MyFile.write("Relative error; max displacement")
    MyFile.write("\n")
    for i in range(N):
        MyFile.write(str(error[i]) + ";" + str(absoluteDisplacement[i]))
        MyFile.write("\n")
    MyFile.close()
    sys.exit()

def solverSU2Framat(args, aeroframeSettings, acceptedNames):
    """
    Function called when the user desires to couple the CFD solver SU2
    and the structure solver FramAT.
    """
    # TODO Step 1) Checks the entry data.
    
    # Step 2) Reads CFD data from CSV file
    logger.debug("Inside the folder")
    path = aeroframeSettings["CFDResultFile"]
    logger.debug(path)
    SU2Data = pd.read_csv(path)

    # Step 3)  Reads CPACS files and computes the nodes of pseudo 1D structural
    #          mesh. Aeroframe function pre-meshes the aircraft to get each
    #          structure node.
    preMeshedStructre = importGeomerty.CsdGeometryImport(args,aeroframeSettings)
    preMeshedStructre.getAllPoints()

    # Step 4) feeds the computed nodes to the structure solver and builds a
    # structure mesh
    csdSolverClassVar = framatWrapper.framat(preMeshedStructre)
    csdSolverClassVar.mesh()

    # Step 5) feeds the computed nodes to a mapping function which computes the
    # tranformation matrices (based on RBF)
    logger.debug("Next step is transformation")
    transform = mappingSU2.mapper(SU2Data,preMeshedStructre,csdSolverClassVar)
    transform.computesTransformationsMatrices()
    transform.aeroToStructure()

def standard_run(args):
    """
    Master function which selects the solvers
    """
    # Starts simulation
    logger.info(f"{__prog_name__} started")

    # Simulation input file
    logger.debug("Setting file path:\n" + args.cwd+"/"+args.run)
    settingsFileAndPath = args.cwd + "/" + args.run
    aeroframeSettings = getSettings(settingsFileAndPath)

    # Calls the correct solver in function of the user input. Supported choice
    # at this point of the development are:
    # 1) CFD solver: pytornado
    #    CSD solver: .csv file -> from external solver
    #                file structure should be a surface of points for each
    #                lines: x;y;z;dx;dy;dz
    # 2) CFD solver: pytornado
    #    CSD solver: FramAT
    logger.info("CFD solver will be: "+str(aeroframeSettings["CFD_solver"]))
    logger.info("CSD solver will be: "+str(aeroframeSettings["CSD_solver"]))

    # Accepted solver notations. Add solvers after. DO NOT CHANGE THE ORDER.
    # Theses variables are used in many classes and the order matters. Change
    # them only if you know what you are doing. If you want to add a new solver
    # just add the solver accepted names after the existing lists.
    pytornado = ["pytornado","Pytornado","pyTornado","PyTornado"]
    deformationFromFile = ["dff","deformationFromFile","external"]
    framat = ["framat","Framat","FramAT","framAT"]
    su2 = ["SU2","su2","Su2"]
    # Assembles all the accepted names for ease of use
    acceptedNames = [pytornado,deformationFromFile,framat,su2]
    # Reads user input
    cfdSolver = aeroframeSettings["CFD_solver"]
    csdSolver = aeroframeSettings["CSD_solver"]

    # Chooses the solvers form what the user input asks.
    if cfdSolver in acceptedNames[0] and csdSolver in acceptedNames[1]:
        logger.info("CFD with pytornado and external CSD solver")
        solverPytornadoCSV(args, aeroframeSettings, acceptedNames)
    elif cfdSolver in acceptedNames[0] and csdSolver in acceptedNames[2]:
        logger.info("CFD with pytornado and CSD with FramAT")
        solverPytornadoFramat(args, aeroframeSettings, acceptedNames)
    elif cfdSolver in acceptedNames[3] and csdSolver in acceptedNames[2]:
        logger.info("CFD with SU2 and CSD with FramAT")
        solverSU2Framat(args, aeroframeSettings, acceptedNames)
    else:
        logger.error("CFD solver or/and CSD solver not supported")
        sys.exit()
    sys.exit()
