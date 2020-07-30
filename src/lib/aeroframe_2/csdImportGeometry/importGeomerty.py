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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull,convex_hull_plot_2d
import sys


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
            self.fuselageCPACSexists = True
            logger.info("CPACS file fuselage count: " + str(self.nFuselage))
            if self.nFuselage > 1:
                logger.error("Too many fuselages")
                sys.exit()
        except:
            self.fuselageCPACSexists = False
            logger.info("No fuselage found in CPACS file")
            sys.exit()
        
        # Number of wings
        try:
            self.nWings = self.tigl.getWingCount()
            self.wingsCPACSexists = True
            logger.info("CPACS file wing count: "+str(self.nWings))
        except:
            self.wingsCPACSexists = False
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
            self.wignsJsonExists = True
            # logger.debug(wingsSettings)
            if len(wingsSettings) != self.nWings:
                logger.error("JSON wings instances number do not match CPACS file")
                sys.exit()
            else:
                logger.info("All wing instances found in JSON file")
        except:
            logger.error("Not all CPACS wings where found in JSON setting file")
            self.wignsJsonExists = False
            sys.exit()
        
        # logger.debug(self.fuselageJsonExists)
        # logger.debug(self.wignsJsonExists)
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

    def getWingsSectionCenterPointsAndArea(self):
        """
        Warning, the hypothesis of having a center of rotation which is not
        alligned with the center of pressure of the wing is done. This 
        hypothesis at the moment is such that the center of rotation is at
        1/3 of the chord.
        """
        logger.info("Starts computing each section of the wing central point")
        # OK TODO: add a parameter for the chord axis of rotation.
        
        
        if self.nWings > 0:
            # Number of sections in each wing
            self.wingsNsegm = np.empty((self.nWings))
            self.wingsNsect = np.empty((self.nWings))
            # Distance from the leading edge
            # WARNING This is probably also the mass line of the wing
            self.wingHingeAxis = np.empty((self.nWings))
            for i in range(self.nWings):
                self.wingsNsegm[i] = self.tigl.wingGetSegmentCount(i+1)
                self.wingsNsect[i] = self.tigl.wingGetSectionCount(i+1)
                # logger.debug(self.settings["wing"+str(i+1)]["mechanicalCenter"])
                self.wingHingeAxis[i] = self.settings["wing"+str(i+1)]["mechanicalCenter"]
            # Gets each section and segement name
            self.wingsSegmentNames = []
            self.wingsSectionsCenters = []
            self.wingsSectionsArea = []

            for i in range(self.nWings):
                # segmentUID = []
                N = int(self.wingsNsect[i])
                center = np.empty((N,3))
                area = np.empty(N)
                for j in range(N-1):
                    # segmentUID.append(self.tigl.wingGetSegmentUID(i+1, j+1))
                    upper = self.tigl.wingGetUpperPoint(i+1, j+1, 0, self.wingHingeAxis[i])  # WARNING WITH THE 0
                    upper = np.array(upper)
                    lower = self.tigl.wingGetLowerPoint(i+1, j+1, 0, self.wingHingeAxis[i])
                    lower = np.array(lower)
                    center[j] = 0.5*(upper + lower)
                    area[j] = self.areaWingSection(i,j,0)
                    # logger.debug(area[j])
                    # logger.debug("Wing section area:"+str(self.tigl.))
                # Gets last part of each segment
                upper = self.tigl.wingGetUpperPoint(i+1, j+1, 1, self.wingHingeAxis[i])
                upper = np.array(upper)
                lower = self.tigl.wingGetLowerPoint(i+1, j+1, 1, self.wingHingeAxis[i])
                lower = np.array(lower)
                center[j+1] = 0.5*(upper + lower)
                area[j+1] = self.areaWingSection(i,j,1)
                # logger.debug("="*50)
                # logger.debug(center)
                # logger.debug(area)

                # If there is a symmetry, a set of point symmetric to the plane
                # is added
                symmetry = self.tigl.wingGetSymmetry(i+1)
                if symmetry == 0:
                    self.wingsSectionsCenters.append(center)
                    self.wingsSectionsArea.append(area)
                # Necessary, otherwise python will add the same vector to
                # self.wingsSectionsCenters leading to an error
                

                
                # In tigl3wrapper.py the symmetry is defined as such
                # class TiglSymmetryAxis(object):
                # TIGL_NO_SYMMETRY = 0
                # TIGL_X_Y_PLANE = 1
                # TIGL_X_Z_PLANE = 2
                # TIGL_Y_Z_PLANE = 3
                if symmetry > 0:
                    center_copy = np.copy(center)
                    area_copy = np.copy(area)
                    if symmetry == 1:
                        index = 2
                    elif symmetry == 2:
                        index = 1
                    elif symmetry == 3:
                        index = 0
                    # Computes symmetry
                    for k in range(len(center)):
                        center_copy[k][index] = -center[k,index]
                        # The -1 avoids copying two times the "same" point
                    centerConcat = np.concatenate((np.flip(center_copy[1:],axis=0),center))
                    areaConcat = np.concatenate((np.flip(area_copy[1:],axis=0),area))
                    self.wingsSectionsCenters.append(centerConcat)
                    self.wingsSectionsArea.append(areaConcat)
                logger.debug("="*50)
                logger.debug("Wing sections areas")
                logger.debug(self.wingsSectionsArea)
            
            

        else:
            logger.warning("No wings found in the CPACS file")
            sys.exit()

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
        self.getWingsSectionCenterPointsAndArea()
        self.getFuselageSectionCenterPointsAndArea()
        # self.computesFuselageMeshPoints()
        self.computesWingsMeshPoints()
    
    def computesWingsMeshPoints(self):
        # TODO corriger l'emplacement des points clef si moins de 7 points dans le fuselage initial
        # TODO verify the tables lengths for the fuselage!!
        self.userAskedNNodesWings = np.zeros(self.nWings)
        for i in range(self.nWings):
            self.userAskedNNodesWings[i] = self.settings["wing"+str(i+1)]["FEM"]["nodesFEM"]
            if self.userAskedNNodesWings[i] < 2:
                logger.error("Not enough wing points")
                sys.exit()
            logger.debug("Number of wing nodes asked: "+str(self.userAskedNNodesWings))    
            
            # Initialization of the fuselage mesh points
            N = int(self.userAskedNNodesWings[i])
            self.currentWingMeshPoints = np.empty((N,3))
            self.currentWingMeshArea = np.empty(N)
            self.currentWingNodesName = []
            # gets fuselage total length:
            self.currentWingLength = 0
            logger.debug(len(self.wingsSectionsCenters[i]-2))
            self.currentWingSegmentsLength = np.empty(len(self.wingsSectionsCenters[i])-1)
            self.currentWingSegmentsRelativePoints = np.empty(len(self.wingsSectionsCenters[i])-1)
            # Computes
            #   relative position in terms of lenth
            #   total length
            #   each segement lemgth
            for j in range(len(self.wingsSectionsCenters[i])-1):

                self.currentWingSegmentsRelativePoints[j] = self.currentWingLength
                p1 = self.wingsSectionsCenters[i][j]
                p2 = self.wingsSectionsCenters[i][j+1]
                dist = np.linalg.norm(p2-p1)
                logger.debug(j)
                logger.debug(dist)
                self.currentWingLength += dist
                self.currentWingSegmentsLength[j] = dist
            
            self.currentWingSegmentsRelativePoints[-1] = self.currentWingLength
            # sets equidistant nodes from 0 to max current wing length
            self.currentWingFutureMeshPoints = np.linspace(0,self.currentWingLength,N)

            logger.debug("self.currentWingLength")
            logger.debug(self.currentWingLength)
            logger.debug("self.currentWingSegmentsLength")
            logger.debug(self.currentWingSegmentsLength)
            logger.debug("self.currentWingSegmentsRelativePoints")
            logger.debug(self.currentWingSegmentsRelativePoints)
            logger.debug("self.currentWingFutureMeshPoints")
            logger.debug(self.currentWingFutureMeshPoints)
            self.currentWingNsegments = len(self.wingsSectionsCenters[i])
            for j in range(len(self.currentWingFutureMeshPoints)):
                x = self.currentWingFutureMeshPoints[j]
                # logger.debug(self.currentWingSegmentsRelativePoints[:]-x)
                dist = self.currentWingSegmentsRelativePoints[:]-x
                segmentIndex = np.argmin(np.abs(dist))
                logger.debug("segment index"+str(segmentIndex))
                if segmentIndex == 0:
                    segmentUID = self.tigl.wingGetComponentSegmentUID(i+1, 1)
                    eta = 0
                # elif segmentIndex == self.currentWingNsegments:
                #     segmentUID = self.tigl.wingGetComponentSegmentUID(i+1, self.currentWingNsegments)
                #     eta = 1
                elif x > self.currentWingSegmentsRelativePoints[segmentIndex]:
                    logger.debug("i+1 = "+str(i+1))
                    logger.debug("currentWingSegments"+str(self.currentWingNsegments))
                    logger.debug("segmentIndex = "+str(segmentIndex))
                    logger.debug(self.tigl.wingGetComponentSegmentCount(i+1))
                    segmentUID = self.tigl.wingGetComponentSegmentUID(i+1, 9)
                    # segmentUID = self.tigl.wingGetComponentSegmentUID(i+1, segmentIndex)
                # #     eta = self.currentWingSegmentsRelativePoints[segmentIndex+1]-x
                # #     eta = eta / self.currentWingSegmentsLength[segmentIndex]
                # elif x < self.currentWingSegmentsRelativePoints[segmentIndex]:
                #     segmentUID = self.tigl.fuselageGetSegmentUID(i+1, segmentIndex)
            #         eta = self.currentWingSegmentsRelativePoints[segmentIndex]-x
            #         eta = eta / self.currentWingSegmentsLength[segmentIndex-1]
                
            #     point = self.tigl.fuselageGetSectionCenter(segmentUID,eta)
            #     self.fuselageMeshPoints[j] = np.array(point)
            #     area = self.tigl.fuselageGetCrossSectionArea(segmentUID, eta)
            #     self.fuselageMeshArea[j] = area
            #     self.fuselageNodesName = "fuse"+str(j)

        #     logger.debug(self.fuselageMeshPoints)
        #     logger.debug(self.fuselageMeshArea)
    
    def computesFuselageMeshPoints(self):
        # TODO corriger l'emplacement des points clef si moins de 7 points dans le fuselage initial
        
        # Uploads user mesh informations
        if self.nFuselage > 0:
            self.userAskedNNodesFuselage = self.settings["fuselage"]["FEM"]["nodesFEM"]
            logger.debug("Number of fuselage nodes asked: " + str(self.userAskedNNodesFuselage))
            if self.userAskedNNodesFuselage < 2:
                logger.error("Not enough fuselage points")
                sys.exit()
            
            # Initialization of the fuselage mesh points
            self.fuselageMeshPoints = np.empty((self.userAskedNNodesFuselage,3))
            self.fuselageMeshArea = np.empty(self.userAskedNNodesFuselage)
            self.fuselageNodesName = []
            # gets fuselage total length:
            self.fuselageLength = 0
            self.segmentsLength = np.empty(len(self.fuselageSectionsPoints)-1)
            self.segmentsRelativePoints = np.empty(len(self.fuselageSectionsPoints)-1)
            # Computes
            #   relative position in terms of lenth
            #   total length
            #   each segement lemgth
            for i in range(len(self.fuselageSectionsPoints)-1):

                self.segmentsRelativePoints[i] = self.fuselageLength

                p1 = self.fuselageSectionsPoints[i]
                p2 = self.fuselageSectionsPoints[i+1]
                dist = np.linalg.norm(p2-p1)

                self.fuselageLength += dist
                self.segmentsLength[i] = dist
            self.segmentsRelativePoints[-1] = self.fuselageLength

            # sets equidistant nodes from 0 to max fuselage length
            self.futureMeshPoints = np.linspace(0,self.fuselageLength,self.userAskedNNodesFuselage)

            logger.debug(self.fuselageLength)
            logger.debug(self.segmentsLength)
            logger.debug(self.segmentsRelativePoints)
            logger.debug(self.futureMeshPoints)

            for i in range(len(self.futureMeshPoints)):
                x = self.futureMeshPoints[i]
                # logger.debug(self.segmentsRelativePoints[:]-x)
                dist = self.segmentsRelativePoints[:]-x
                segmentIndex = np.argmin(np.abs(dist))
                # logger.debug(x)
                # logger.debug(self.segmentsRelativePoints[i])
                if segmentIndex == 0:
                    segmentUID = self.tigl.fuselageGetSegmentUID(1, 1)
                    eta = 0
                elif segmentIndex == self.nFuseSegm:
                    segmentUID = self.tigl.fuselageGetSegmentUID(1, self.nFuseSegm)
                    eta = 1
                elif x > self.segmentsRelativePoints[segmentIndex]:
                    segmentUID = self.tigl.fuselageGetSegmentUID(1, segmentIndex+1)
                    eta = self.segmentsRelativePoints[segmentIndex+1]-x
                    eta = eta / self.segmentsLength[segmentIndex]
                elif x < self.segmentsRelativePoints[segmentIndex]:
                    segmentUID = self.tigl.fuselageGetSegmentUID(1, segmentIndex)
                    eta = self.segmentsRelativePoints[segmentIndex]-x
                    eta = eta / self.segmentsLength[segmentIndex-1]
                
                point = self.tigl.fuselageGetSectionCenter(segmentUID,eta)
                self.fuselageMeshPoints[i] = np.array(point)
                area = self.tigl.fuselageGetCrossSectionArea(segmentUID, eta)
                self.fuselageMeshArea[i] = area
                self.fuselageNodesName.append("fuse1node"+str(i))

            logger.debug(self.fuselageMeshPoints)
            logger.debug(self.fuselageMeshArea)

    def plotSectionsPoints(self):
        N = len(self.wingsSectionsCenters)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(self.origin[0],self.origin[1],self.origin[2])
        if self.nFuselage > 0:
            ax.scatter(self.fuselageSectionsPoints[:,0],
                       self.fuselageSectionsPoints[:,1],
                       self.fuselageSectionsPoints[:,2])
            # ax.scatter(self.fuselageMeshPoints[:,0],
            #            self.fuselageMeshPoints[:,1],
            #            self.fuselageMeshPoints[:,2])
        for i in range(N):
            ax.scatter(self.wingsSectionsCenters[i][:,0],
                       self.wingsSectionsCenters[i][:,1],
                       self.wingsSectionsCenters[i][:,2])
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