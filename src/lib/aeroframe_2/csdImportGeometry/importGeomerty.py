#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 08:45:26 2020

@author: Jean-Philippe Kuntzer

OK TODO: See how the program behaves with a single wing and correct up to here
OK TODO: See how the program behaves with a boxwing and correct up to here
OK TODO: distribute the nodes accordingly for the fuselage
OK TODO: Generate the nodes names for the fuselage
OK TODO vérifier qu'il y a assez de points
OK TODO répartir linéairement
OK TODO: Verify that the wing section is correct 
    -> there is some dispersion
OK TODO: Verify that the fuselage section is correct
OK TODO: distribute the nodes accordingly for the wing
OK TODO: add the clamped node
OK TODO: Generate the nodes names for the wing
OK TODO: Assemble the points together in a big matrix
OK TODO: Assemble the names together in a big matrix
OK TODO: compute 'A' for each points
OK TODO: compute 'Iy' for each points
OK TODO: compute 'Ix' for each points
OK TODO: compute 'J' for each points
OK TODO: Add an savety for when there is no input origin.0
OK TODO: solve the wing node name issue
OK TODO: Verify if the areas correspond

TODO: for wings find the attached nodes pairs
TODO: find a uid for each cross sections
      cs = self.model.add_feature('cross_section', uid='dummy')

(BONUS)
TODO: Fix the points distibution with the boxwing, there is an issue with tigl
TDOD: Take into account boxwings

"""

import logging
import tigl3.configuration
from tixi3 import tixi3wrapper
from tigl3 import tigl3wrapper
from tigl3 import geometry
import numpy as np
from numpy.core.umath_tests import inner1d
import math
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull,convex_hull_plot_2d
import sys
from scipy.spatial.transform import Rotation as R


logger = logging.getLogger(__name__)


class CsdGeometryImport:

    def __init__(self,inputAircraftPath,settings):
        """
        Initializes the instance and stores the values of the CPACS file to the
        handle into the class.

        Getwingcount
        """
        # Aircraft absolute path
        self.aircraftPath = inputAircraftPath
        # JSON input file
        self.settings = settings

        # Tixi and Tigl handler. Necessary to read the CPACS
        self.tixi = tixi3wrapper.Tixi3()
        self.tigl = tigl3wrapper.Tigl3()
        self.tixi.open(self.aircraftPath)
        self.tigl.open(self.tixi,"")

        # Number of fuselage instances, generally one, more could create
        # problems.
        try:
            self.nFuselage = self.tigl.getFuselageCount()
            logger.info("CPACS file fuselage count: " + str(self.nFuselage))
            if self.nFuselage > 1:
                logger.error("Too many fuselages")
                sys.exit()
        except :
            logger.info("No fuselage found in CPACS file")
            sys.exit()
        
        # Number of wings
        try:
            self.nWings = self.tigl.getWingCount()
            logger.info("CPACS file wing count: "+str(self.nWings))
        except:
            logger.error("No wings found in CPACS file")
            sys.exit()
                    
        
        # Reads JSON
        # Searches for fuselage instance
        fuselageSettings = []
        try:
            for i in range(self.nFuselage):
                fuselageSettings.append(settings["fuselage"])
            self.fuselageJsonExists = True
            if len(fuselageSettings) != self.nFuselage:
                logger.error("JSON fuselage instances number do not match CPACS file")
                sys.exit()
            else:
                logger.info("All fuselage instance(s) found in JSON file")
        except:
            self.fuselageJsonExists = False
            logger.info("No fuselage instance found in JSON")
            sys.exit()
            
        
        # Searches for wing instance
        wingsSettings = []
        try:
            for i in range(self.nWings):
                wingsSettings.append(settings["wing" + str(i+1)])
            if len(wingsSettings) != self.nWings:
                logger.error("JSON wings instances number do not match CPACS file")
                sys.exit()
            else:
                logger.info("All wing instances found in JSON file")
        except:
            logger.error("Not all CPACS wings where found in JSON setting file")
            sys.exit()


    def getOrigin(self):
        """
        Reads where the origin of the aircraft is and stores it into the class
        """
        try:
            ogX = self.tixi.getDoubleElement("/cpacs/vehicles/aircraft/model/reference/point/x")
            ogY = self.tixi.getDoubleElement("/cpacs/vehicles/aircraft/model/reference/point/y")
            ogZ = self.tixi.getDoubleElement("/cpacs/vehicles/aircraft/model/reference/point/z")
            self.origin = np.array([ogX, ogY, ogZ])
        except:
            logger.error("No origin found in the CPACS file")
            sys.exit()

    def getAllPoints(self):
        """
        Calls all the necessary functions to build a mesh
        """
        self.getOrigin()
        self.computesWingsMeshPoints()
        self.computesFuselageMeshPoints()
        self.assembleMatrices()
        # self.areasComparison()  # for debuig
        self.computesWingConnexions()
        self.plotSectionsPoints()
        
    def areasComparison(self):
        # if self.nFuselage > 0:
        #     logger.debug("Fuselage initial area:\n "+str(self.fs_m_pointsInitArea[0]))
        #     logger.debug("Fuselage final A:\n "+str(self.fs_m_pointsA[0]))
        #     logger.debug("Fuselage final Iy:\n "+str(self.fs_m_pointsIy[0]))
        #     logger.debug("Fuselage final Iz:\n "+str(self.fs_m_pointsIz[0]))
        #     logger.debug("Fuselage final J:\n "+str(self.fs_m_pointsJ[0]))
        #     logger.debug("Fuselage nodes names"+str(self.fs_m_pointsName[0]))
        # for i in range(self.nWings):
        #     logger.debug("Wing initial area:\n "+str(self.ws_m_pointsInitArea[i]))
        #     logger.debug("Wing final A:\n "+str(self.ws_m_pointsA[i]))
        #     logger.debug("Wing final Iy:\n "+str(self.ws_m_pointsIy[i]))
        #     logger.debug("Wing final Iz:\n "+str(self.ws_m_pointsIz[i]))
        #     logger.debug("Wing final J:\n "+str(self.ws_m_pointsJ[i]))
        #     logger.debug("Wing nodes names"+str(self.ws_m_pointsName[i]))
        N = len(self.aircraftNodesPoints)
        for i in range(N):
            logger.debug("Aircraft nodes:\n"+str(self.aircraftNodesPoints[i]))
            logger.debug("Aircraft nodes names:\n"+str(self.aircraftNodesNames[i]))
            logger.debug("Aircraft A:\n"+str(self.aircraftNodesA[i]))
            logger.debug("Aircraft Iy:\n"+str(self.aircraftNodesIy[i]))
            logger.debug("Aircraft Iz:\n"+str(self.aircraftNodesIz[i]))
            logger.debug("Aircraft J:\n"+str(self.aircraftNodesJ[i]))
        sys.exit()
        
    def getWingChordLinePoint(self,wingIndex,segmentIndex,eta,xsi):
        # tigl.wingGetUpperPoint(wingIndex, segmentIndex, eta -> y, xsi->x)
        up = self.tigl.wingGetUpperPoint(wingIndex,segmentIndex,eta,xsi)
        dw = self.tigl.wingGetLowerPoint(wingIndex,segmentIndex,eta,xsi)
        upper = np.array(up)
        down = np.array(dw)
        center = 0.5*(upper + down)
        return center
    
    def computePointSectionArea(self,wingIndex,segmentIndex,eta,xsi):
        """
        Computes the wing section area. No matter the ortientation of the wing
        """
        # tigl.wingGetUpperPoint(wingIndex, segmentIndex, eta -> y, xsi->x)
        # WARNING there is a slight difference in the area computed with this
        # method ans CPACSCREATOR. At the moment it is undetermined who is more
        # accurate.
        N = 20
        xsi1 = np.linspace(0,1,N)
        upper = np.empty((N,3))
        lower = np.empty((N,3))
        for i in range(N):
            u = self.tigl.wingGetUpperPoint(wingIndex,segmentIndex,eta,xsi1[i])
            l = self.tigl.wingGetLowerPoint(wingIndex,segmentIndex,eta,xsi1[i])
            upper[i] = np.array(u)
            lower[i] = np.array(l)
        v1 = upper[0]-upper[-1]
        v2 = upper[7] - lower[7]
        v1xv2 = np.cross(v1,v2)
        upper = np.flip(upper,axis=0)
        wingSectionPoints = np.concatenate((upper, lower))
        ey_0 = np.array([0,1,0])
        e_1 = v1xv2
        # Computes the cross prodct
        cross = np.cross(ey_0,e_1)
        normCross = np.linalg.norm(cross)
        cross = cross/normCross
        if normCross < 1e-8:
            # No need to rotate
            wingSectionPoints = np.delete(wingSectionPoints,1,1)
            hull = ConvexHull(wingSectionPoints)
            area = hull.volume
        else:
            ab = inner1d(ey_0,e_1)
            a = np.linalg.norm(ey_0)
            b = np.linalg.norm(e_1)
            angle = np.arccos(ab / (a*b))
            logger.debug("angle: "+str(angle))
            quat = angle*cross
            r = R.from_rotvec(quat)
            # Deletes the y column since the Convex hull will struggle with
            # a 3d plane otherwise
            wingSectionPoints = r.apply(wingSectionPoints)
            wingSectionPoints = np.delete(wingSectionPoints,1,1)
            hull = ConvexHull(wingSectionPoints)
            area = hull.volume
        
        logger.debug("Computed section area: "+str(area))

        return area
                            
    def computesWingsMeshPoints(self):
        """
        Nomencalture:
            c_ : current
            w_ : wing  (relate to current wing)
            ws_: wings (related to all wings)
            sg : segments
            sc : section
            *N_: number
            m_ : mesh
        """
        self.userAskedNNodesWings = np.zeros(self.nWings)
        self.ws_m_points = []
        self.ws_m_pointsName = []
        self.ws_m_pointsInitArea = []
        for i in range(self.nWings):
            
            # Basic wing input check
            self.userAskedNNodesWings[i] = self.settings["wing"+str(i+1)]["FEM"]["nodesFEM"]
            w_m_N_nodes = int(self.userAskedNNodesWings[i])
            if w_m_N_nodes < 2:
                logger.error("Not enough points for wing"+str(i+1)+" (min 2)")
                sys.exit()

            logger.debug("Number of wing nodes asked: "+str(w_m_N_nodes))    
            # distance from leading edge
            xsi = self.settings["wing"+str(i+1)]["mechanicalCenter"]
            logger.debug("Wing"+str(i+1)+" mechanical center is: "+str(xsi))
            wingIndex = i+1
            
            # Searches for each wing segment the start and end relative point
            w_N_sg = self.tigl.wingGetSegmentCount(i+1)
            w_N_sc = self.tigl.wingGetSectionCount(i+1)
            logger.debug("Wing"+str(i+1)+" has "+str(w_N_sg)+" segments")
            logger.debug("Wing"+str(i+1)+" has "+str(w_N_sc)+" sections")
            if w_m_N_nodes < w_N_sc:
                logger.warning("Mesh underdetermined, less mesh points than actual CPACS sections")
            
            # Gets each segments starting and ending points
            w_sg_points = np.empty((w_N_sg+1,3))
            for j in range(w_N_sg):   
                w_sg_points[j] = self.getWingChordLinePoint(wingIndex,j+1,0,xsi)
            w_sg_points[-1] = self.getWingChordLinePoint(wingIndex,j+1,1,xsi)
            logger.debug("Wing"+str(wingIndex)+" segment points:\n"+str(w_sg_points))
            
            # Gets each segments length
            w_sg_length = np.empty(w_N_sg)
            w_sg_relativePosition = np.empty(w_N_sg+1)
            w_length = 0
            for j in range(w_N_sg):
                w_sg_relativePosition[j] = w_length
                length = np.linalg.norm(w_sg_points[j] - w_sg_points[j+1])
                w_sg_length[j] = length
                w_length += length
            w_sg_relativePosition[-1] = w_length
            logger.debug("Wing"+str(wingIndex)+" segments lengths are:\n"+str(w_sg_length))
            logger.debug("Wing"+str(wingIndex)+" segments relative positions are:\n"+str(w_sg_relativePosition))
            logger.debug("Wing"+str(wingIndex)+" length is:"+str(w_length))

            # Computes mesh relative points
            w_m_relativePoints = np.linspace(0, w_length, w_m_N_nodes)
            logger.debug("Wing"+str(wingIndex)+" relative mesh points:\n"+str(w_m_relativePoints))
            
            # If the user askes more points that there in the CPACS file 
            # definitions the program automatically changes the position to the
            # closest known point. This features ensures that the simulations
            # will be made with maximal fidelity to the definintion.
            logger.debug("+"*20)
            logger.debug("wing relative pos:\n"+str(w_sg_relativePosition))
            logger.debug("mesh relative pos:\n"+str(w_m_relativePoints))
            if w_N_sc <= w_m_N_nodes:
                for j in range(w_N_sc):
                    diff = np.abs(w_m_relativePoints - w_sg_relativePosition[j])
                    index = np.argmin(diff)
                    w_m_relativePoints[index] = w_sg_relativePosition[j]
            
            logger.debug("mesh relative pos:\n"+str(w_m_relativePoints))
            
            # Computes the eta for each segment in order to get the mesh point
            # from tigl
            w_m_points = np.empty((w_m_N_nodes,3))
            w_m_pointsName = []
            w_m_pointsInitArea = np.empty(w_m_N_nodes)
            for j in range(w_m_N_nodes):
                # finds in which segment the mesh point will be
                relativePosition = w_m_relativePoints[j]
                dist = w_sg_relativePosition - relativePosition
                segmentIndex = np.argmin(np.abs(dist))+1
                # o--x-------o situations
                if dist[segmentIndex-1] < 0:
                    case = 1
                    eta = w_m_relativePoints[j] - w_sg_relativePosition[segmentIndex-1]
                    eta = (eta/w_sg_length[segmentIndex-1])
                # o--x-------o situation
                elif dist[segmentIndex-1] > 0:
                    case = 2
                    eta = w_sg_relativePosition[segmentIndex-1] - w_m_relativePoints[j]
                    segmentIndex = segmentIndex -1
                    eta = 1 - (eta/w_sg_length[segmentIndex-1])
                elif dist[segmentIndex-1] == 0.0 and segmentIndex == 1:
                    case = 3
                    eta = 0
                elif dist[segmentIndex-1] == 0.0 and segmentIndex != 1:
                    case = 4
                    eta = 1
                    segmentIndex -= 1
                else:
                    logger.error("Something wrong with CPACS file")
                    sys.exit()
                # logger.debug()
                logger.debug("case "+str(case)+"  eta = "+str(eta))

                # Gets the wing point
                w_m_points[j] = self.getWingChordLinePoint(wingIndex,segmentIndex,eta,xsi)
                name = "w_"+str(i+1)+"_n_"+str(j)
                if self.nFuselage == 0:
                    if np.abs(w_m_points[j][1]) < 1e-2:
                        name = "w_n_clamped"
                w_m_pointsName.append(name)
                # Computes section area
                area = self.computePointSectionArea(wingIndex,segmentIndex,eta,xsi)
                w_m_pointsInitArea[j] = area
                
            # In tigl3wrapper.py the symmetry is defined as such
            # class TiglSymmetryAxis(object):
            # TIGL_NO_SYMMETRY = 0
            # TIGL_X_Y_PLANE = 1
            # TIGL_X_Z_PLANE = 2
            # TIGL_Y_Z_PLANE = 3
            symmetry = self.tigl.wingGetSymmetry(i+1)
            if symmetry > 0:
                w_m_points_copy = np.copy(w_m_points)
                w_m_pointsName_copy = w_m_pointsName.copy()
                w_m_pointsInitArea_c = np.copy(w_m_pointsInitArea)
                if symmetry == 1:
                    index = 2
                elif symmetry == 2:
                    index = 1
                elif symmetry == 3:
                    index = 0
                
                # Computes symmetry
                for k in range(w_m_N_nodes):
                    w_m_points_copy[k][index] = -w_m_points[k,index]
                    w_m_pointsName_copy[k] = w_m_pointsName_copy[k] + "sym"
                    # The -1 avoids copying two times the "same" point
                w_m_points = np.concatenate((np.flip(w_m_points_copy[1:],axis=0),w_m_points))
                rev = w_m_pointsName_copy[::-1]
                w_m_pointsName = rev[:-1] + w_m_pointsName
                # logger.debug(w_m_pointsInitArea)
                # logger.debug(np.flip(w_m_pointsInitArea_c))
                w_m_pointsInitArea = np.concatenate((np.flip(w_m_pointsInitArea_c[1:],axis=0),w_m_pointsInitArea))
            
            logger.debug("Wing mesh points:\n"+str(w_m_points))
            self.ws_m_points.append(w_m_points)
            self.ws_m_pointsInitArea.append(w_m_pointsInitArea)
            self.ws_m_pointsName.append(w_m_pointsName)

    def computesFuselageMeshPoints(self):
        """
        
        """
        # Uploads user mesh informations
        if self.nFuselage > 0:
            self.userAskedNNodesFuselage = self.settings["fuselage"]["FEM"]["nodesFEM"]
            # Initializes fuselage part class variables
            f_m_N_nodes = int(self.userAskedNNodesFuselage)
            self.fs_m_points = []
            self.fs_m_pointsName = []
            self.fs_m_pointsInitArea = []
            logger.debug("Number of fuselage nodes asked: " + str(self.userAskedNNodesFuselage))
            if self.userAskedNNodesFuselage < 3:
                logger.error("Not enough fuselage points (min 3)")
                sys.exit()

            logger.debug("Number of fuselage nodes asked: "+str(f_m_N_nodes))    
            fuselageIndex = 1
            
            # Searches for each wing segment the start and end relative point
            f_N_sg = self.tigl.fuselageGetSegmentCount(fuselageIndex)
            f_N_sc = self.tigl.fuselageGetSectionCount(fuselageIndex)
            logger.debug("Fuselage has "+str(f_N_sg)+" segments")
            logger.debug("Fuselage has "+str(f_N_sc)+" sections")
            if f_m_N_nodes < f_N_sc:
                logger.warning("Fuselage mesh underdetermined, less mesh points than actual CPACS sections")
            
            # Gets each segments starting and ending points
            f_sg_points = np.empty((f_N_sg+1,3))
            f_sg_names = []
            for j in range(f_N_sg):
                f_sg_names.append(self.tigl.fuselageGetSegmentUID(1, j+1))
                f_sg_points[j] = self.tigl.fuselageGetSectionCenter(f_sg_names[j], 0)
            f_sg_points[-1] = self.tigl.fuselageGetSectionCenter(f_sg_names[-1], 1)
            logger.debug("Fuselage segment points:\n"+str(f_sg_points))
            
            # Gets each segments length
            f_sg_length = np.empty(f_N_sg)
            f_sg_relativePosition = np.empty(f_N_sg+1)
            f_length = 0
            for j in range(f_N_sg):
                f_sg_relativePosition[j] = f_length
                length = np.linalg.norm(f_sg_points[j] - f_sg_points[j+1])
                f_sg_length[j] = length
                f_length += length
            f_sg_relativePosition[-1] = f_length
            logger.debug("Fuselage segments lengths are:\n"+str(f_sg_length))
            logger.debug("Fuselage segments relative positions are:\n"+str(f_sg_relativePosition))
            logger.debug("Fuselage length is:"+str(f_length))

            # Computes mesh relative points
            f_m_relativePoints = np.linspace(0, f_length, f_m_N_nodes)
            logger.debug("Fuselage relative mesh points:\n"+str(f_m_relativePoints))
            
            # Corrects node closet to the origin to be exactly on the same
            # x postion.
            dist = np.empty(f_N_sg)
            for j in range(f_N_sg):
                dist[j] = np.linalg.norm(f_sg_points[j]-self.origin)
            index = np.argmin(dist)
            if f_sg_points[index,0] - self.origin[0] > 0:
                index -= 1
            etaClamped = (self.origin[0] - f_sg_points[index][0])/f_sg_length[index]
            diff = np.empty(f_m_N_nodes)
            
            logger.debug(f_sg_relativePosition[index])
            logger.debug(f_sg_relativePosition[index] + etaClamped)
            logger.debug(f_m_relativePoints)
            
            for j in range(f_m_N_nodes):
                diff[j] = np.abs(f_sg_relativePosition[index]+etaClamped - f_m_relativePoints[j])
            
            closest = np.argmin(diff)
            logger.debug("Origin: "+str(self.origin))
            logger.debug("closest point:"+str(f_sg_points[index]))
            logger.debug("segment length"+str(f_sg_length[index]))
            logger.debug("eta = "+str(etaClamped))
            # closest = 100000
            # etaClamped = 0

            # Computes the eta for each segment in order to get the mesh point
            # from tigl
            f_m_points = np.empty((f_m_N_nodes,3))
            f_m_pointsName = []
            f_m_pointsInitArea = np.empty(f_m_N_nodes)
            for j in range(f_m_N_nodes):
                # finds in which segment the mesh point will be
                relativePosition = f_m_relativePoints[j]
                dist = f_sg_relativePosition - relativePosition
                segmentIndex = np.argmin(np.abs(dist))+1
                # o--x-------o situations
                if dist[segmentIndex-1] < 0:
                    case = 1
                    eta = f_m_relativePoints[j] - f_sg_relativePosition[segmentIndex-1]
                    eta = (eta/f_sg_length[segmentIndex-1])
                    name = "f_n_"+str(j+1)
                    if j == closest:
                        name = "f_n_clamped"
                        eta = etaClamped
                # o--x-------o situation
                elif dist[segmentIndex-1] > 0:
                    case = 2
                    eta = f_sg_relativePosition[segmentIndex-1] - f_m_relativePoints[j]
                    segmentIndex = segmentIndex -1
                    eta = 1 - (eta/f_sg_length[segmentIndex-1])
                    name = "f_n_"+str(j+1)
                    if j == closest:
                        name = "f_n_"+str(j+1)+"_clamped"
                        eta = etaClamped
                elif dist[segmentIndex-1] == 0.0 and segmentIndex == 1:
                    case = 3
                    eta = 0
                    name = "f_n_"+str(j+1)
                    if j == closest:
                        name = "f_n_"+str(j+1)+"_clamped"
                        eta = 0
                elif dist[segmentIndex-1] == 0.0 and segmentIndex != 1:
                    case = 4
                    eta = 1
                    segmentIndex -= 1
                    name = "f_n_"+str(j+1)
                    if j == closest:
                        name = "f_n_"+str(j+1)+"_clamped"
                        eta = 1
                else:
                    logger.error("Something wrong with CPACS file")
                # logger.debug()

                # Gets the fuselage center point
                f_m_points[j] = self.tigl.fuselageGetSectionCenter(f_sg_names[segmentIndex-1], eta)
                
                f_m_pointsName.append(name)
                # Computes section area
                area = self.tigl.fuselageGetCrossSectionArea(f_sg_names[segmentIndex-1], eta)
                f_m_pointsInitArea[j] = area
            
            logger.debug("fuselage center points:\n"+str(f_m_points))
            logger.debug("fuselage center area:\n"+str(f_m_pointsInitArea))
            logger.debug("fuselage points names:\n"+str(f_m_pointsName))
            # sys.exit()
            self.fs_m_points.append(f_m_points)
            self.fs_m_pointsName.append(f_m_pointsName)
            self.fs_m_pointsInitArea.append(f_m_pointsInitArea)

    def assembleMatrices(self):
        """
        Assembles fuselage points and wing points into a "nodes" matrix/instance
        Assembles each point CPACS area in an "CPACSarea" matrix/instance
        Assembles each point names into a "nodesNames" matrix/instance
        """
        # All nodes informations
        self.aircraftNodesPoints = []
        self.aircraftNodesNames = []
        self.aircraftInitNodesAreas = [] # Will be destroid
        self.aircraftNodesA = []
        self.aircraftNodesIy = []
        self.aircraftNodesIz = []
        self.aircraftNodesJ = []
        # More general information
        self.aircraftBeamsMaterials = []
        self.aircraftConnectedNodes = []
        
        self.computeProportionFuselage()
        self.computeProportionWings()
        
        # adds fulseage infos to matrices if there is one
        if self.nFuselage > 0:
            self.aircraftNodesPoints.append(self.fs_m_points[0])
            self.aircraftNodesNames.append(self.fs_m_pointsName[0])
            self.aircraftInitNodesAreas.append(self.fs_m_pointsInitArea[0])
            self.aircraftNodesA.append(self.fs_m_pointsA[0])
            self.aircraftNodesIy.append(self.fs_m_pointsIy[0])
            self.aircraftNodesIz.append(self.fs_m_pointsIz[0])
            self.aircraftNodesJ.append(self.fs_m_pointsJ[0])

        for i in range(self.nWings):
            self.aircraftNodesPoints.append(self.ws_m_points[i])
            self.aircraftNodesNames.append(self.ws_m_pointsName[i])
            self.aircraftInitNodesAreas.append(self.ws_m_pointsInitArea[i])
            self.aircraftNodesA.append(self.ws_m_pointsA[i])
            self.aircraftNodesIy.append(self.ws_m_pointsIy[i])
            self.aircraftNodesIz.append(self.ws_m_pointsIz[i])
            self.aircraftNodesJ.append(self.ws_m_pointsJ[i])

    def computeProportionFuselage(self):
        """
        Computes the area proportions witht the following procedure
        1) get the biggest area of the fuselage
        2) get the biggest area of the wing
        3) checks if the user asks for a constant/linear/quadratic 
           distribution.
        4) Computes the local variable
        5) Feeds it back to the fuselage mesh matrices
        """
        
        # reads settings file for distribution        
        if self.nFuselage > 0:
            self.fs_m_pointsA = []
            self.fs_m_pointsIy = []
            self.fs_m_pointsIz = []
            self.fs_m_pointsJ = []
            inputMaxA = self.settings["fuselage"]["mechanicalProperties"]["A"]
            inputMaxIy = self.settings["fuselage"]["mechanicalProperties"]["Iy"]
            inputMaxIz = self.settings["fuselage"]["mechanicalProperties"]["Iz"]
            inputMaxJ = self.settings["fuselage"]["mechanicalProperties"]["J"]
            inputInterpolationType = self.settings["fuselage"]["FEM"]["mechanicalInterpolationType"]
            logger.debug(self.fs_m_pointsInitArea[0])
            index = np.argmax(self.fs_m_pointsInitArea[0])
            logger.debug(index)
            fuselageMaxArea = self.fs_m_pointsInitArea[0][index]

            if inputInterpolationType == "constant":
                exposant = 0
            elif inputInterpolationType == "linear":
                exposant = 1
            elif inputInterpolationType == "quadratic":
                exposant = 2
            else:
                logger.error("Fuselage mechanical properties distribution is")
                logger.error("wrong. Accepted values are \"constant\",")
                logger.error("\"linear\", \quadratic\"")
                sys.exit()
            
            coef = np.empty(self.userAskedNNodesFuselage)
            self.f_m_pointsA = np.empty(self.userAskedNNodesFuselage)
            self.f_m_pointsIy = np.empty(self.userAskedNNodesFuselage)
            self.f_m_pointsIz = np.empty(self.userAskedNNodesFuselage)
            self.f_m_pointsJ = np.empty(self.userAskedNNodesFuselage)
            for i in range(self.userAskedNNodesFuselage):
                coef[i] = self.fs_m_pointsInitArea[0][i]/fuselageMaxArea
                self.f_m_pointsA[i] = inputMaxA * coef[i]
                self.f_m_pointsIy[i] = inputMaxIy * coef[i]**exposant
                self.f_m_pointsIz[i] = inputMaxIz * coef[i]**exposant
                self.f_m_pointsJ[i] = inputMaxJ * coef[i]**exposant
            self.fs_m_pointsA.append(self.f_m_pointsA)
            self.fs_m_pointsIy.append(self.f_m_pointsIy)
            self.fs_m_pointsIz.append(self.f_m_pointsIz)
            self.fs_m_pointsJ.append(self.f_m_pointsJ)
        else:
            logger.warning("No fuselage found")

    def computeProportionWings(self):
        """
        Computes the area proportions witht the following procedure
        1) get the biggest area of the fuselage
        2) get the biggest area of the wing
        3) checks if the user asks for a constant/linear/quadratic 
           distribution.
        4) Computes the local variable
        5) Feeds it back to the wing mesh matrices
        """
        self.ws_m_pointsA = []
        self.ws_m_pointsIy = []
        self.ws_m_pointsIz = []
        self.ws_m_pointsJ = []
        for i in range(self.nWings):
            wingIndex = i+1
            wingNnodes = self.settings["wing"+str(wingIndex)]["FEM"]["nodesFEM"]
            if self.tigl.wingGetSymmetry(i+1) != 0:
                wingNnodes = 2*wingNnodes - 1
            inputMaxA = self.settings["wing"+str(wingIndex)]["mechanicalProperties"]["A"]
            inputMaxIy = self.settings["wing"+str(wingIndex)]["mechanicalProperties"]["Iy"]
            inputMaxIz = self.settings["wing"+str(wingIndex)]["mechanicalProperties"]["Iz"]
            inputMaxJ = self.settings["wing"+str(wingIndex)]["mechanicalProperties"]["J"]
            inputInterpolationType = self.settings["wing"+str(wingIndex)]["FEM"]["mechanicalInterpolationType"]
            index = np.argmax(self.ws_m_pointsInitArea[i])
            wingMaxArea = self.ws_m_pointsInitArea[i][index]

            if inputInterpolationType == "constant":
                exposant = 0
            elif inputInterpolationType == "linear":
                exposant = 1
            elif inputInterpolationType == "quadratic":
                exposant = 2
            else:
                logger.error("Wing mechanical properties distribution is")
                logger.error("wrong. Accepted values are \"constant\",")
                logger.error("\"linear\", \quadratic\"")
                sys.exit()
            
            coef = np.empty(wingNnodes)
            self.w_m_pointsA = np.empty(wingNnodes)
            self.w_m_pointsIy = np.empty(wingNnodes)
            self.w_m_pointsIz = np.empty(wingNnodes)
            self.w_m_pointsJ = np.empty(wingNnodes)
            for j in range(wingNnodes):
                coef[j] = self.ws_m_pointsInitArea[i][j]/wingMaxArea
                self.w_m_pointsA[j] = inputMaxA * coef[j]
                self.w_m_pointsIy[j] = inputMaxIy * coef[j]**exposant
                self.w_m_pointsIz[j] = inputMaxIz * coef[j]**exposant
                self.w_m_pointsJ[j] = inputMaxJ * coef[j]**exposant
            self.ws_m_pointsA.append(self.w_m_pointsA)
            self.ws_m_pointsIy.append(self.w_m_pointsIy)
            self.ws_m_pointsIz.append(self.w_m_pointsIz)
            self.ws_m_pointsJ.append(self.w_m_pointsJ)

    def computesWingConnexions(self):
        """
        Computes each wing pair connexions.
        """
        if self.nWings <= 1:
            logger.info("No wing connexions needed")
        elif self.nWings > 1 and self.nFuselage < 1:
            logger.error("Multiple wings with no fusekage!")
            sys.exit()
        else:
            ws_connextionsPoints = self.nWings + self.nFuselage
            # if self.nFuselage > 0:
            #     ws_connextionsPoints = self.nWings
            # else:
            #     ws_connextionsPoints = self.nWings -1
            
            # ws: wings
            # d: distance
            # L: left, R: right, c: center
            ws_d_L = []
            ws_d_L = []
            ws_d_L = []
            for i in range(ws_connextionsPoints):
                """
                for each wing there is 3 points that need to be taken into 
                account:
                    1) Left tip position[0]
                    2) center [np.floor(nNodes/2))] (note symmetric wings always have odd number of points)
                    3) right tip position[-1]
                """
                if self.nFuselage > 0:
                    wingIndex = i+1
                else:
                    wingIndex = i

                # # N = len(self.aircraftNodesPoints)
                # c = int(np.floor(len(self.aircraftNodesPoints[wingIndex])/2))
                # logger.debug("c = "+str(c))
                # ws_identifiers = []
                # for j in range(ws_connextionsPoints):
                #     if wingIndex != j:
                #         # Computes the distances between all point
                #         # logger.debug(self.aircraftNodesPoints[wingIndex][0])
                #         # logger.debug(self.aircraftNodesPoints[j])
                #         # logger.debug(self.aircraftNodesPoints[j] - self.aircraftNodesPoints[wingIndex][0])
                #         dist_l = np.linalg.norm(self.aircraftNodesPoints[j] - self.aircraftNodesPoints[wingIndex][0],axis=0)
                #         dist_c = np.linalg.norm(self.aircraftNodesPoints[j] - self.aircraftNodesPoints[wingIndex][c],axis=0)
                #         dist_r = np.linalg.norm(self.aircraftNodesPoints[j] - self.aircraftNodesPoints[wingIndex][-1],axis=0)
                #         # Gets the index of the minimal distance between all point
                #         index_l = np.argmin(dist_l)
                #         index_c = np.argmin(dist_c)
                #         index_r = np.argmin(dist_r)
                #         indexes = np.array([index_l,index_c,index_r])
                #         # gets the minimal distance
                #         minDist_l = dist_l[index_l]
                #         minDist_c = dist_c[index_c]
                #         minDist_r = dist_r[index_r]
                #         minDist = np.array([minDist_l,minDist_c,minDist_r])
                #         k = np.argmin(minDist)
                #         minimalDistance = minDist[k]
                #         tab = np.array([0,c,1e6]) # TODO Correct the mistake
                #         k = tab[k]
                #         l = indexes[k]
                #         # for wing i the minimal distance to object j is
                #         # point k of the wing and point l of the object. 
                #         # Minimal distance between these objects is: minimalDistance
                #         # i = wing number
                #         # j = object number (wing or fuselage)
                #         # k = point on wing (l,c,r)
                #         # minimalDistance: distance between both of these points
                #         identifier = np.array([wingIndex,j,k,l,minimalDistance])
                #         ws_identifiers.append(identifier)
                #         logger.debug(ws_identifiers)
                #         sys.exit()
                # logger.debug(ws_identifiers)
                # sys.exit()
    def plotSectionsPoints(self):
        # N = len(self.wingsSectionsCenters)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(self.origin[0],self.origin[1],self.origin[2])
        if self.nFuselage > 0:
            ax.scatter(self.fs_m_points[0][:,0],
                       self.fs_m_points[0][:,1],
                       self.fs_m_points[0][:,2])
        for i in range(self.nWings):
            ax.scatter(self.ws_m_points[i][:,0],
                        self.ws_m_points[i][:,1],
                        self.ws_m_points[i][:,2])
        size = 10
        ax.set_xlim(-size+self.origin[0],size+self.origin[0])
        ax.set_ylim(-size+self.origin[1],size+self.origin[1])
        ax.set_zlim(-size+self.origin[2],size+self.origin[2])
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')        
        plt.show()