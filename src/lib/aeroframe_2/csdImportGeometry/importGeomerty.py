#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 08:45:26 2020

@author: Jean-Philippe Kuntzer
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
        # TODO: can be a source of error. Should be tested with multiple CPACS
        # files.
        ogX = self.tixi.getDoubleElement("/cpacs/vehicles/aircraft/model/reference/point/x")
        ogY = self.tixi.getDoubleElement("/cpacs/vehicles/aircraft/model/reference/point/y")
        ogZ = self.tixi.getDoubleElement("/cpacs/vehicles/aircraft/model/reference/point/z")
        self.origin = np.array([ogX, ogY, ogZ])

    def getFuselageSectionCenterPointsAndArea(self):
        """
        Gets all the fuselage sections center points. This will give us a line
        of points that will be meshed afterwards for the CSD simulation.
        """
        if self.nFuselage > 0:
            logger.info("Starts computing each section of the fuselage central point")
            self.nFuseSect = self.tigl.fuselageGetSectionCount(1)
            self.nFuseSegm = self.tigl.fuselageGetSegmentCount(1)
            self.fuseSegmNames = []

            # Gets each fuselage section UID and stores it
            for i in range(1,self.nFuseSegm + 1):
                self.fuseSegmNames.append(self.tigl.fuselageGetSegmentUID(1,i))

            # Computes each fuselage section central point
            self.fuselageSectionsPoints = np.empty((self.nFuseSegm+1,3))
            self.fuselageSectionsArea = np.empty(self.nFuseSegm+1)
            for i in range(0,self.nFuseSegm):
                name = self.fuseSegmNames[i]
                point = self.tigl.fuselageGetSectionCenter(name,0.0)
                area = self.tigl.fuselageGetCrossSectionArea(name, 0.0)
                self.fuselageSectionsPoints[i] = point
                self.fuselageSectionsArea[i] = area
            point = self.tigl.fuselageGetSectionCenter(name,1.0)
            area = self.tigl.fuselageGetCrossSectionArea(name, 0.0)
            self.fuselageSectionsPoints[-1] = point
            self.fuselageSectionsArea[-1] = area
            logger.info("Fuselage sections center points and area are computed")
            logger.debug("Fuselage center points:\n"+str(self.fuselageSectionsPoints))
            logger.debug("Fuselage areas:\n"+str(self.fuselageSectionsArea))

        elif self.nFuselage > 1:
            logger.warning("More than one fuselage found in the CPACS file,")
            logger.warning("only the first one is taken into account for this")
            logger.warning("simulation")

        else:
            logger.warning("No fuselage found in the CPACS file")

    def areaWingSection(self,i,j,eta):
        """
        Computes wing section area by assuming wing is a rectagle. This is
        done this way because sections can be tilted and complexity will
        increase tremendeuly if complex Hull are used.
        """
        # The more you know
        # self.tigl.wingGetUpperPoint(wingIndex, segmentIndex, eta -> y, xsi->x)
        
        # Computes some key points, 0.25 is arbitrary but should be a good
        # estimate of biggest wing thickness
        leadingEdge = upper = np.array(self.tigl.wingGetUpperPoint(i+1,j+1,eta,0))
        trailingEdge = upper = np.array(self.tigl.wingGetUpperPoint(i+1,j+1,eta,1))
        upper = np.array(self.tigl.wingGetUpperPoint(i+1,j+1,eta,0.25))
        lower = np.array(self.tigl.wingGetLowerPoint(i+1,j+1,eta,0.25))
        
        width = np.linalg.norm(leadingEdge - trailingEdge)
        height = np.linalg.norm(upper - lower)
        
        area = width * height
        return area

    # def getWingsSectionCenterPointsAndArea(self):
    #     """
    #     Warning, the hypothesis of having a center of rotation which is not
    #     alligned with the center of pressure of the wing is done. This 
    #     hypothesis at the moment is such that the center of rotation is at
    #     1/3 of the chord.
    #     """
    #     logger.info("Starts computing each section of the wing central point")
    #     # OK TODO: add a parameter for the chord axis of rotation.
        
        
    #     if self.nWings > 0:
    #         # Number of sections in each wing
    #         self.wingsNsegm = np.empty((self.nWings))
    #         self.wingsNsect = np.empty((self.nWings))
    #         # Distance from the leading edge
    #         # WARNING This is probably also the mass line of the wing
    #         self.wingHingeAxis = np.empty((self.nWings))
    #         for i in range(self.nWings):
    #             self.wingsNsegm[i] = self.tigl.wingGetSegmentCount(i+1)
    #             self.wingsNsect[i] = self.tigl.wingGetSectionCount(i+1)
    #             # logger.debug(self.settings["wing"+str(i+1)]["mechanicalCenter"])
    #             self.wingHingeAxis[i] = self.settings["wing"+str(i+1)]["mechanicalCenter"]
            
    #         logger.debug("Number of sections per wing"+str(self.wingsNsect))
    #         logger.debug("Number of segments per wing"+str(self.wingsNsegm))
            
    #         # Gets each section and segement name
    #         self.wingsSegmentNames = []
    #         self.wingsSectionsCenters = []
    #         self.wingsSectionsArea = []

    #         for i in range(self.nWings):
    #             # segmentUID = []
    #             N = int(self.wingsNsect[i])
    #             center = np.empty((N,3))
    #             area = np.empty(N)
    #             for j in range(N-1):
    #                 # segmentUID.append(self.tigl.wingGetSegmentUID(i+1, j+1))
    #                 upper = self.tigl.wingGetUpperPoint(i+1, j+1, 0, self.wingHingeAxis[i])  # WARNING WITH THE 0
    #                 upper = np.array(upper)
    #                 lower = self.tigl.wingGetLowerPoint(i+1, j+1, 0, self.wingHingeAxis[i])
    #                 lower = np.array(lower)
    #                 center[j] = 0.5*(upper + lower)
    #                 area[j] = self.areaWingSection(i,j,0)
    #                 # logger.debug(area[j])
    #                 # logger.debug("Wing section area:"+str(self.tigl.))
    #             # Gets last part of each segment
    #             upper = self.tigl.wingGetUpperPoint(i+1, j+1, 1, self.wingHingeAxis[i])
    #             upper = np.array(upper)
    #             lower = self.tigl.wingGetLowerPoint(i+1, j+1, 1, self.wingHingeAxis[i])
    #             lower = np.array(lower)
    #             center[j+1] = 0.5*(upper + lower)
    #             area[j+1] = self.areaWingSection(i,j,1)
    #             # logger.debug("="*50)
    #             # logger.debug(center)
    #             # logger.debug(area)

    #             # If there is a symmetry, a set of point symmetric to the plane
    #             # is added
    #             symmetry = self.tigl.wingGetSymmetry(i+1)
    #             if symmetry == 0:
    #                 self.wingsSectionsCenters.append(center)
    #                 self.wingsSectionsArea.append(area)
    #             # Necessary, otherwise python will add the same vector to
    #             # self.wingsSectionsCenters leading to an error
                

                
    #             # In tigl3wrapper.py the symmetry is defined as such
    #             # class TiglSymmetryAxis(object):
    #             # TIGL_NO_SYMMETRY = 0
    #             # TIGL_X_Y_PLANE = 1
    #             # TIGL_X_Z_PLANE = 2
    #             # TIGL_Y_Z_PLANE = 3
    #             if symmetry > 0:
    #                 center_copy = np.copy(center)
    #                 area_copy = np.copy(area)
    #                 if symmetry == 1:
    #                     index = 2
    #                 elif symmetry == 2:
    #                     index = 1
    #                 elif symmetry == 3:
    #                     index = 0
    #                 # Computes symmetry
    #                 for k in range(len(center)):
    #                     center_copy[k][index] = -center[k,index]
    #                     # The -1 avoids copying two times the "same" point
    #                 centerConcat = np.concatenate((np.flip(center_copy[1:],axis=0),center))
    #                 areaConcat = np.concatenate((np.flip(area_copy[1:],axis=0),area))
    #                 self.wingsSectionsCenters.append(centerConcat)
    #                 self.wingsSectionsArea.append(areaConcat)
    #             logger.debug("="*50)
    #             logger.debug("Wing sections areas")
    #             logger.debug(self.wingsSectionsArea)
            
            

    #     else:
    #         logger.warning("No wings found in the CPACS file")
    #         sys.exit()

    def getAllPoints(self):
        # OK TODO: See how the program behaves with a single wing and correct up to here
        # OK TODO: See how the program behaves with a boxwing and correct up to here
        # OK TODO: distribute the nodes accordingly for the fuselage
        # TODO: distribute the nodes accordingly for the wing
        # TODO: add the clamped node
        # OK TODO: Generate the nodes names for the fuselage
        # TODO: Generate the nodes names for the wing
        # TODO: Assemble the points together in a big matrix
        # TODO: Assemble the names together in a big matrix
        # OK TODO vérifier qu'il y a assez de points
        # OK TODO répartir linéairement
        """
        Calls all the necessary functions to build a mesh
        """
        self.getOrigin()
        # self.getWingsSectionCenterPointsAndArea()
        self.computesWingsMeshPoints()
        # self.getFuselageSectionCenterPointsAndArea()
        self.computesFuselageMeshPoints()
        
    def getWingChordLinePoint(self,wingIndex,segmentIndex,eta,xsi):
        # tigl.wingGetUpperPoint(wingIndex, segmentIndex, eta -> y, xsi->x)
        up = self.tigl.wingGetUpperPoint(wingIndex,segmentIndex,eta,xsi)
        dw = self.tigl.wingGetLowerPoint(wingIndex,segmentIndex,eta,xsi)
        upper = np.array(up)
        down = np.array(dw)
        center = 0.5*(upper + down)
        return center
    
    def computePointSectionArea(self,wingIndex,segmentIndex,eta,xsi):
        # tigl.wingGetUpperPoint(wingIndex, segmentIndex, eta -> y, xsi->x)
        # TODO try to get something working with convex hulls.
        up = self.tigl.wingGetUpperPoint(wingIndex,segmentIndex,eta,xsi)
        lw = self.tigl.wingGetLowerPoint(wingIndex,segmentIndex,eta,xsi)
        le = self.tigl.wingGetUpperPoint(wingIndex,segmentIndex,eta,0)
        te = self.tigl.wingGetUpperPoint(wingIndex,segmentIndex,eta,1)
        upper = np.array(up)
        lower = np.array(lw)
        leadi = np.array(up)
        trail = np.array(lw)
        hight = np.linalg.norm(upper-lower)
        width = np.linalg.norm(leadi-trail)
        area = hight*width
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
        self.ws_m_pointsArea = []
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
            
            # Computes the eta for each segment in order to get the mesh point
            # from tigl
            w_m_points = np.empty((w_m_N_nodes,3))
            w_m_pointsName = []
            w_m_pointsArea = np.empty(w_m_N_nodes)
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
                w_m_pointsName.append(name)
                # Computes section area
                area = self.computePointSectionArea(wingIndex,segmentIndex,eta,xsi)
                w_m_pointsArea[j] = area
                
            # In tigl3wrapper.py the symmetry is defined as such
            # class TiglSymmetryAxis(object):
            # TIGL_NO_SYMMETRY = 0
            # TIGL_X_Y_PLANE = 1
            # TIGL_X_Z_PLANE = 2
            # TIGL_Y_Z_PLANE = 3
            symmetry = self.tigl.wingGetSymmetry(i+1)
            if symmetry > 0:
                w_m_points_copy = np.copy(w_m_points)
                w_m_pointsName_copy = w_m_pointsName
                w_m_pointsArea = np.copy(w_m_pointsArea)
                if symmetry == 1:
                    index = 2
                elif symmetry == 2:
                    index = 1
                elif symmetry == 3:
                    index = 0
                
                # Computes symmetry
                for k in range(w_m_N_nodes):
                    w_m_points_copy[k][index] = -w_m_points[k,index]
                    w_m_pointsName_copy[k] += "sym"
                    # The -1 avoids copying two times the "same" point
                w_m_points = np.concatenate((np.flip(w_m_points_copy[1:],axis=0),w_m_points))
                rev = w_m_pointsName_copy[::-1]
                w_m_pointsName = rev + w_m_pointsName
                logger.debug(w_m_pointsArea)
                logger.debug(np.flip(w_m_pointsArea))
                w_m_pointsArea = np.concatenate((np.flip(w_m_pointsArea),w_m_pointsArea))
            
            # TODO add improvement that when the number of element is higher
            # than the number of segments point. a function pushes the closer
            # point to the exact segment end to get more accurate results.
            
            logger.debug("Wing mesh points:\n"+str(w_m_points))
            self.ws_m_points.append(w_m_points)
            self.ws_m_pointsArea.append(w_m_pointsArea)
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
            self.fs_m_pointsAera = []
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
                logger.warning("Mesh underdetermined, less mesh points than actual CPACS sections")
            
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
            
            # Computes the eta for each segment in order to get the mesh point
            # from tigl
            f_m_points = np.empty((f_m_N_nodes,3))
            f_m_pointsName = []
            f_m_pointsArea = np.empty(f_m_N_nodes)
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
                # o--x-------o situation
                elif dist[segmentIndex-1] > 0:
                    case = 2
                    eta = f_sg_relativePosition[segmentIndex-1] - f_m_relativePoints[j]
                    segmentIndex = segmentIndex -1
                    eta = 1 - (eta/f_sg_length[segmentIndex-1])
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

                # Gets the fuselage center point
                f_m_points[j] = self.tigl.fuselageGetSectionCenter(f_sg_names[segmentIndex-1], eta)
                name = "f_n_"+str(j+1)
                f_m_pointsName.append(name)
                # Computes section area
                area = self.tigl.fuselageGetCrossSectionArea(f_sg_names[segmentIndex-1], eta)
                f_m_pointsArea[j] = area
            logger.debug("fuselage center points:\n"+str(f_m_points))
            logger.debug("fuselage center area:\n"+str(f_m_pointsArea))
            logger.debug("fuselage points names:\n"+str(f_m_pointsName))
            self.fs_m_points.append(f_m_points)
            self.fs_m_pointsName.append(f_m_pointsArea)
            self.fs_m_pointsAera.append(f_m_pointsName)
# =============================================================================
            # # Initialization of the fuselage mesh points
            # self.fuselageMeshPoints = np.empty((self.userAskedNNodesFuselage,3))
            # self.fuselageMeshArea = np.empty(self.userAskedNNodesFuselage)
            # self.fuselageNodesName = []
            # # gets fuselage total length:
            # self.fuselageLength = 0
            # self.segmentsLength = np.empty(len(self.fuselageSectionsPoints)-1)
            # self.segmentsRelativePoints = np.empty(len(self.fuselageSectionsPoints)-1)
            # # Computes
            # #   relative position in terms of lenth
            # #   total length
            # #   each segement lemgth
            # for i in range(len(self.fuselageSectionsPoints)-1):

            #     self.segmentsRelativePoints[i] = self.fuselageLength

            #     p1 = self.fuselageSectionsPoints[i]
            #     p2 = self.fuselageSectionsPoints[i+1]
            #     dist = np.linalg.norm(p2-p1)

            #     self.fuselageLength += dist
            #     self.segmentsLength[i] = dist
            # self.segmentsRelativePoints[-1] = self.fuselageLength

            # # sets equidistant nodes from 0 to max fuselage length
            # self.futureMeshPoints = np.linspace(0,self.fuselageLength,self.userAskedNNodesFuselage)

            # logger.debug(self.fuselageLength)
            # logger.debug(self.segmentsLength)
            # logger.debug(self.segmentsRelativePoints)
            # logger.debug(self.futureMeshPoints)

            # for i in range(len(self.futureMeshPoints)):
            #     x = self.futureMeshPoints[i]
            #     # logger.debug(self.segmentsRelativePoints[:]-x)
            #     dist = self.segmentsRelativePoints[:]-x
            #     segmentIndex = np.argmin(np.abs(dist))
            #     # logger.debug(x)
            #     # logger.debug(self.segmentsRelativePoints[i])
            #     if segmentIndex == 0:
            #         segmentUID = self.tigl.fuselageGetSegmentUID(1, 1)
            #         eta = 0
            #     elif segmentIndex == self.nFuseSegm:
            #         segmentUID = self.tigl.fuselageGetSegmentUID(1, self.nFuseSegm)
            #         eta = 1
            #     elif x > self.segmentsRelativePoints[segmentIndex]:
            #         segmentUID = self.tigl.fuselageGetSegmentUID(1, segmentIndex+1)
            #         eta = self.segmentsRelativePoints[segmentIndex+1]-x
            #         eta = eta / self.segmentsLength[segmentIndex]
            #     elif x < self.segmentsRelativePoints[segmentIndex]:
            #         segmentUID = self.tigl.fuselageGetSegmentUID(1, segmentIndex)
            #         eta = self.segmentsRelativePoints[segmentIndex]-x
            #         eta = eta / self.segmentsLength[segmentIndex-1]
                
            #     point = self.tigl.fuselageGetSectionCenter(segmentUID,eta)
            #     self.fuselageMeshPoints[i] = np.array(point)
            #     area = self.tigl.fuselageGetCrossSectionArea(segmentUID, eta)
            #     self.fuselageMeshArea[i] = area
            #     self.fuselageNodesName.append("fuse1node"+str(i))

            # logger.debug(self.fuselageMeshPoints)
            # logger.debug(self.fuselageMeshArea)

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

    def plotMeshPoints(self):
        pass