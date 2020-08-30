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
OK TODO: Verify that the displacements do not add up

TODO: checks folder structure and prints the potential errors
TODO: In deformation from csv mode the program should read a logfile where
      there is the results of the previous simulations
TODO: set the logging level from the command line
TODO: add gravity in Pytornado and SU2
TODO: add motor push

WARNING the table accepted names is not supported in all classes
"""

import logging
import json
import aeroframe_2.deformation.functions as aeroDef
import aeroframe_2.deformation.framatNormalVecConverter as normalVecFramat
import aeroframe_2.csdImportGeometry.importGeomerty as importGeomerty
import aeroframe_2.informationTransfer.mapping as mapping
import aeroframe_2.informationTransfer.mappingSU2 as mappingSU2
import aeroframe_2.wrappers.framatWrapper as framatWrapper
import aeroframe_2.wrappers.framatWrapperSU2 as framatWrapperSU2
import pytornado.stdfun.run as cfd
import ceasiompy.SU2Run.su2run as SU2_fsi
import numpy as np
# import pickle
import sys
import matplotlib.pyplot as plt
import copy
# import pandas as pd
import os
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd


# from ceasiompy.SU2Run.func.extractloads import extract_loads
# import SU2_CFD

# TODO take the information from the args parameters
logging.basicConfig(level=logging.DEBUG)
__prog_name__ = "aeroframe_2"
logger = logging.getLogger(__prog_name__+"."+__name__)


def getSettings(inputFileAndPath):
    """
    Reads the json input file. This files contains the user input choices. If
    there is no file, the function throws an error.

    Parameters
    ----------
    inputFileAndPath : string
        Contains the path to the user input JSON file.

    Raises
    ------
    FileNotFoundError
        If the file is not found the function throws that error

    Returns
    -------
    settings : Dictionay
        Returns a dictionaly with all the user input settings like which CFD
        solver should be used, which CSD solver, wings mechanical properties,
        materials, etc..

    """
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

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.
    aeroframeSettings : Dicitonary
        Contains all the information that there is in the input JSON file, i.e.
        the which solver to use and mechanical properties of the wing

    Returns
    -------
    pytornadoSttings : string
        Contains the equivalent of a command line input command for pytornado.
        With this variable we are tricking pytornado into thinking that he is
        called from the command line which is not the case.
    pytornadoVariables : list
        Variables containing the mesh, the results and some other information,
        necessary to run a pytornado simulation.

    """

    # Virtual command line setting file input path. the "args" contains the
    # rest of the needed information
    pytornadoSttings = args.cwd + "/CFD/settings/" + aeroframeSettings["CFD_settings_file"]
    # Buids CFD mesh
    lattice,vlmdata,settings,aircraft,cur_state,state = cfd.meshing(args,pytornadoSttings)
    # Join all the variables for ease of use.
    pytornadoVariables = [lattice, vlmdata, settings, aircraft, cur_state, state]
    logger.info("Pytornado meshing done")

    return pytornadoSttings, pytornadoVariables


def csvDeformation(args,pytornadoSettings,aeroframeSettings,pytornadoVariables):
    """
    Deforms the VLM mesh from a csv file, by calling the aeroDef class and the
    method Mesh_Def.

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.
    pytornadoSettings : TYPE
        Settings variables files needed to run pytornado
    aeroframeSettings : Dicitonary
        Contains all the information that there is in the input JSON file, i.e.
        the which solver to use and mechanical properties of the wing
    pytornadoVariables : list
        Variables containing the mesh, the results and some other information,
        necessary to run a pytornado simulation.

    Returns
    -------
    pytornadoVariables : list
        Variables containing the mesh, the results and some other information,
        necessary to run a pytornado simulation.

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

    Parameters
    ----------
    pytornadoVariables : list
        Variables containing the mesh, the results and some other information,
        necessary to run a pytornado simulation.
    meshDeformation : aeroDef
        Vairalbe of class aeroDef that contains all the necessary mesh points
        to run pytornado.

    Returns
    -------
    pytornadoVariables : list
        Variables containing the mesh, the results and some other information,
        necessary to run a pytornado simulation.
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
    Writes the results to a csv file. WARNING at the moment (19.08.20) there is
    another function doing the same in CEASIOMpy but that function extracts the
    loads at the boundleg midpoint which not what the FEM mesh needs.

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.
    pytornadoVariables : list
        Variables containing the mesh, the results and some other information,
        necessary to run a pytornado simulation.
    results : pytornado result class
        Variable containg all the panelwise results, hence airspeed, forces,
        etc.

    Returns
    -------
    None.

    """
    # assembless the path
    path = args.cwd + "/CSD/results/panelwiseForces.csv"
    headers = "x;y;z;fx;fy;fz"
    # Gets simulation values
    panelCoordinates = pytornadoVariables[0].c
    panelFx = results["vlmdata"].panelwise['fx']
    panelFy = results["vlmdata"].panelwise['fy']
    panelFz = results["vlmdata"].panelwise['fz']
    panelMx = results["vlmdata"].panelwise['mx']
    panelMy = results["vlmdata"].panelwise['my']
    panelMz = results["vlmdata"].panelwise['mz']

    results = np.array([panelCoordinates[:,0],
                        panelCoordinates[:,1],
                        panelCoordinates[:,2],
                        panelFx,
                        panelFy,
                        panelFz,
                        panelMx,
                        panelMy,
                        panelMz])

    np.savetxt(path, results.T, delimiter=';', header=headers)
    logger.info("Simulation finised")


def solverPytornadoCSV(args, aeroframeSettings, acceptedNames):
    """
    The workflow is as follows:
        1) The function builds CFD mesh.
        2) Then it deforms the CFD mesh.
        3) Solves the CFD problem. (basically calls the solver)
        4) Extract panelwise forces and saves them into a csv file.

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.
    aeroframeSettings : Dicitonary
        Contains all the information that there is in the input JSON file, i.e.
        the which solver to use and mechanical properties of the wing
    acceptedNames : list
        List containing all the implemented solvers. If the user asks a solver
        which is not implemented the program will throw an error and close
        itself.

    Returns
    -------
    None.

    """

    # Step 1) pytornado meshing
    pytornadoSettings, pytornadoVariables = pytornadoMeshing(args, aeroframeSettings)

    # Step 2) Deforms CFD mesh
    lattice = pytornadoVariables[0]
    meshDeformation = aeroDef.Mesh_Def(args,aeroframeSettings,lattice)
    meshDeformation.deformation(acceptedNames)
    pytornadoVariables = feeder(pytornadoVariables,meshDeformation)

    # Step 3) Computes the CFD problem
    pytornadoVariables, results = cfd.solver(pytornadoVariables)

    # Step 4) Saves panelwise forces results
    forcesToCsv(args,pytornadoVariables,results)
    logger.info("End of simulation")
    sys.exit()


def solverPytornadoFramat(args, aeroframeSettings, acceptedNames):
    """
    Function called when the user desires to do a simulation with the VLM
    solver Pytornado and the structure solver FramAT.

    The workflow is as follows:
        1)  Builds CFD mesh.
        2)  Reads CPACS files and computes the nodes of pseudo 1D structural
            mesh.
        3)  Builds CSD instance in FramAT.
        4)  Computes the transformation matrices
        5)  Computes CFD problem.
        6)  Enters the aeroelastic loop.
        7)  Projects the loads on CSD instance.
        8)  Computes CSD solution
        9) Deforms the CFD mesh.
        10) Computes the norm of the displacement error
        11) Computes new CFD problem.
        12) loops back to point 6) if simulation has not converged.
        13) Saves the results to a CSV file

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.
    aeroframeSettings : Dicitonary
        Contains all the information that there is in the input JSON file, i.e.
        the which solver to use and mechanical properties of the wing
    acceptedNames : list
        List containing all the implemented solvers. If the user asks a solver
        which is not implemented the program will throw an error and close
        itself.

    Returns
    -------
    None.

    """

    # Step 1) pytornado meshing
    pytornadoSettings, pytornadoVariables = pytornadoMeshing(args, aeroframeSettings)

    # Step 2)  Reads CPACS files and computes the nodes of beam model structural
    #          mesh. Aeroframe function pre-meshes the aircraft to get each
    #          structure node.
    preMeshedStructre = importGeomerty.CsdGeometryImport(args,aeroframeSettings)
    preMeshedStructre.getAllPoints()

    # Step 3) feeds the computed nodes to the structure solver and builds a
    # structure mesh
    csdSolverClassVar = framatWrapper.framat(preMeshedStructre)
    csdSolverClassVar.mesh()

    # Step 4) feeds the computed nodes to a mapping function which computes the
    # tranformation matrices (based on RBF)
    logger.debug(pytornadoVariables[0].c)
    # sys.exit()
    transform = mapping.mapper(pytornadoVariables,preMeshedStructre,csdSolverClassVar)
    transform.computesTransformationsMatrices()
    
    # Step 5) Computes CFD problem.
    cfd.solver(pytornadoVariables)
    pytornadoVariablesInit = copy.deepcopy(pytornadoVariables)
    
    # Setp 6) Aeroelastic loop.
    N = aeroframeSettings["MaxIterationsNumber"]
    i = 0
    maxDisplacement = np.array([0])
    error = []
    absoluteDisplacement = []
    aeroFx = []
    aeroFy = []
    aeroFz = []
    aeroMx = []
    aeroMy = []
    aeroMz = []
    structFx = []
    structFy = []
    structFz = []
    structMx = []
    structMy = []
    structMz = []
    pVold = pytornadoVariables
    tol = aeroframeSettings["ConvergeanceTolerence"]
    while (i < N):
        # basic user comminication
        logger.debug("aeroelastic loop number: "+str(i))
        
        # Saves the aerodynamic results
        points = pVold[0].bound_leg_midpoints
        Fx = pVold[1].panelwise['fx']
        Fy = pVold[1].panelwise['fy']
        Fz = pVold[1].panelwise['fz']
        Mx = pVold[1].panelwise['mx']
        My = pVold[1].panelwise['my']
        Mz = pVold[1].panelwise['mz']
        df = pd.DataFrame(points)
        df['Fx'] = Fx
        df['Fy'] = Fy
        df['Fz'] = Fz
        df['Mx'] = Mx
        df['My'] = My
        df['Mz'] = Mz
        df.columns = ['x', 'y', 'z', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz']
        df.to_csv(args.cwd + '/CFD/_results/forces' + str(i) + '.csv')

        # Step 7) Projects the loads on CSD instance.
        transform.aeroToStructure(args,i)
        logger.debug(transform)

        # Step 8) Compute structure solution
        csdSolverClassVar.run(transform)

        # Step 9) deforms the CFD mesh. Computes beam deformation
        latticeCurrent = pytornadoVariablesInit[0]
        meshDeformation = aeroDef.Mesh_Def(args,aeroframeSettings,latticeCurrent)

        # Step 10) computes new aerodynamic points
        transform.structureToAero(args,i)
        meshDeformation.deformation(acceptedNames,transform)
        pytornadoVariables = feeder(pytornadoVariables,meshDeformation)

        # Step 11) Computes the norm of the displacement error
        # Max structure displacement form one aeroelastic iteration to the next
        maxDisplacement = np.append(maxDisplacement, np.max(transform.suz))
        error.append(np.abs(maxDisplacement[-1] - maxDisplacement[-2]))
        # Max structure displacement from undeformed state
        absoluteDisplacement.append(np.abs(maxDisplacement[-1] - maxDisplacement[0]))
        
        aeroFx.append(transform.totalAerodynamicFx)
        aeroFy.append(transform.totalAerodynamicFy)
        aeroFz.append(transform.totalAerodynamicFz)
        aeroMx.append(transform.totalAerodynamicMx)
        aeroMy.append(transform.totalAerodynamicMy)
        aeroMz.append(transform.totalAerodynamicMz)
        structFx.append(transform.totalStructureFx)
        structFy.append(transform.totalStructureFy)
        structFz.append(transform.totalStructureFz)
        structMx.append(transform.totalStructureMx)
        structMy.append(transform.totalStructureMy)
        structMz.append(transform.totalStructureMz)
        
        # logger.info("Max error between two iteration: "+str(error))
        # logger.info("Max displacement between two iteration: "+str(absoluteDisplacement))
        # logger.info('G load: '+str(transform.G))
        
        # Step 12) Deforms the CFD mesh.
        pytornadoVariables, results = cfd.solver(pytornadoVariables)
        # logger.debug(pytornadoVariables[0].bound_leg_midpoints)
        # sys.exit()
        del(csdSolverClassVar)
        csdSolverClassVar = framatWrapper.framat(preMeshedStructre)
        csdSolverClassVar.mesh()
        
        # del(transform)
        transform = mapping.mapper(pytornadoVariables,preMeshedStructre,csdSolverClassVar)
        
        transform.computesTransformationsMatrices()
        # sys.exit()
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
    MyFile.write("Relative error; max displacement;"\
                 " aero Fx; aero Fy; aero Fz; aero Mx; aero My; aero Mz;"\
                 " struct  Fx; struct  Fy; struct  Fz; struct Mx; struct My; struct Mz"\
                     )
    MyFile.write("\n")
    for i in range(N):
        MyFile.write(str(error[i]) + ";" + \
                     str(absoluteDisplacement[i]) + ";" + \
                     str(aeroFx[i]) + ";" + \
                     str(aeroFy[i]) + ";" + \
                     str(aeroFz[i]) + ";" + \
                     str(aeroMx[i]) + ";" + \
                     str(aeroMy[i]) + ";" + \
                     str(aeroMz[i]) + ";" + \
                     str(structFx[i]) + ";" + \
                     str(structFy[i]) + ";" + \
                     str(structFz[i]) + ";" + \
                     str(structMx[i]) + ";" + \
                     str(structMy[i]) + ";" + \
                     str(structMz[i]) + ";"
                     )
        MyFile.write("\n")
    MyFile.close()
    
    # Writes the forces and points at which they apply
    
    
    sys.exit()


def solverSU2Framat(args, aeroframeSettings, acceptedNames):
    """
    Function called when the user desires to couple the CFD solver SU2
    and the structure solver FramAT.

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.
    aeroframeSettings : Dicitonary
        Contains all the information that there is in the input JSON file, i.e.
        the which solver to use and mechanical properties of the wing
    acceptedNames : list
        List containing all the implemented solvers. If the user asks a solver
        which is not implemented the program will throw an error and close
        itself.

    Returns
    -------
    None.

    """
    # TODO Step ) Checks the entry data.

    # Step 1) Initialization of the loop
    #   Creation of the current loop directory done by the su2run.py file.
    # Case s2etup
    iteration = 0
    wkdir = args.cwd
    # TODO: read it form the config file
    nb_proc = aeroframeSettings["SU2_Nproc"]
    logger.debug("Configuration path: \n"+str(aeroframeSettings))
    logger.debug("nb of proc: \n"+str(nb_proc))
    logger.debug("WKDIR: \n"+str(wkdir))

    # Step 2) Runs a single SU2 simulation
    case = '/Case00_alt0_mach0.3_aoa2.0_aos0.0/'
    ###
    # WARNING
    ###
    # SU2_fsi.run_SU2_fsi(aeroframeSettings, wkdir, case, iteration)

    # Step 3)  Reads CPACS files and computes the nodes of 3D beam structural
    #          mesh. Aeroframe_2 function 'importGeomerty' pre-meshes the
    #          aircraft to get each structure node.
    preMeshedStructre = importGeomerty.CsdGeometryImport(args,aeroframeSettings)
    preMeshedStructre.getAllPoints()
    logger.info("Structure mesh points computed")

    # Step 4) feeds the computed nodes to the structure solver and builds a
    # structure mesh
    csdSolverClassVar = framatWrapper.framat(preMeshedStructre)
    csdSolverClassVar.mesh()
    logger.info("FramAT mesh computed")

    # Step 5) feeds the computed nodes to a mapping function which computes the
    # tranformation matrices (based on RBF)
    logger.debug("Next step is transformation")
    forceFile = '/CFD/Case00_alt0_mach0.3_aoa2.0_aos0.0/' + str(iteration) + '/force.csv'
    logger.debug(wkdir)
    logger.debug(forceFile)
    forceInitFilePath = wkdir + forceFile
    transform = mappingSU2.mapper(forceInitFilePath,preMeshedStructre,csdSolverClassVar)
    transform.computesTransformationsMatrices(forceInitFilePath)

    # Setp 6) Aeroelastic loop.
    N = aeroframeSettings["MaxIterationsNumber"]
    maxDisplacement = np.array([0])
    error = []
    absoluteDisplacement = []
    Gloads = []
    tol = aeroframeSettings["ConvergeanceTolerence"]
    while (iteration < N):

        # basic user comminication
        logger.debug("aeroelastic loop number: "+str(iteration))
        forceFilePath = wkdir + '/CFD/Case00_alt0_mach0.3_aoa2.0_aos0.0/' + str(iteration) + '/force.csv'
        # Makes a copy to avoid memory linked mistakes
        transformCurrent = transform
        csdSolverClassVarCurrent = csdSolverClassVar
        # Step 7) Projects the loads on CSD instance.
        transformCurrent.aeroToStructure(forceFilePath)
        # Step 8) Compute structure solution
        csdSolverClassVarCurrent.run(transformCurrent)
        transformCurrent.structureToAero(iteration,forceInitFilePath,forceFilePath)

        # Step 9) Computes convergence
        maxDisplacement = np.append(maxDisplacement, np.max(transform.displacements))
        error.append(np.abs(maxDisplacement[-1] - maxDisplacement[-2]))
        absoluteDisplacement.append(np.abs(maxDisplacement[-1] - maxDisplacement[0]))
        Gloads.append(transform.G)
        logger.info("Max error between two iteration: "+str(error))
        logger.info("Max displacement between two iteration: "+str(absoluteDisplacement))
        logger.info('G load: '+str(Gloads[iteration]))
        # WARNING do not change it's place unless you know what you are doing
        iteration += 1
        # Step 10) computes new CFD solution

        SU2_fsi.run_SU2_fsi(aeroframeSettings, wkdir, case, iteration)
        if iteration == N-1:
            logger.warning("Simulation has reached max number of step,")
            logger.warning("convergeance is yet to determine!")
        if error[-1] <= tol:
            logger.info("Simulation has converged")
            iteration = N
    # Writes a file which contains the error for each timestep
    N = len(error)
    path = args.cwd + "/results.csv"
    MyFile = open(path,"w")
    MyFile.write("Relative error; max displacement")
    MyFile.write("\n")
    for i in range(N):
        MyFile.write(str(error[i]) + ";" + \
                     str(absoluteDisplacement[i]) + ";" + \
                     str(Gloads[i]))
        MyFile.write("\n")
    MyFile.close()
    sys.exit()


def standard_run(args):
    """
    Master function which manage selects the "path". There is 3 different
    supported combinations at the moment:
        Pytornado - FramAT
        Pytornado - External CSV file
        SU2       - FramAT
    This function will call the correct functions in order to perform the type
    of simulation desired by the user.

    Parameters
    ----------
    args : sys.argv
        Unser input from the command line.

    Returns
    -------
    None.

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
