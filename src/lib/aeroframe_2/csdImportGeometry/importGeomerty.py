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


logger = logging.getLogger(__name__)


class CsdGeometryImport:

    def __init__(self,inputAircraftPath):
        """
        Initializes the instance and stores the values of the CPACS file to the
        handle into the class.

        Getwingcount
        """
        # Aircraft absolute path
        self.aircraftPath = inputAircraftPath

        # Tixi and Tigl handler. Necessary to read the CPACS
        self.tixi = tixi3wrapper.Tixi3()
        self.tigl = tigl3wrapper.Tigl3()
        self.tixi.open(self.aircraftPath)
        self.tigl.open(self.tixi,"")

        # Number of fuselage instances, generally one, more could create
        # problems.
        self.nFuselage = self.tigl.getFuselageCount()
        # Number of wings
        self.nWings = self.tigl.getWingCount()
        logger.info("Number of fuselage instances found: " +
                    str(self.nFuselage))
        logger.info("Number of wing instances found:     " +
                    str(self.nWings))

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

    def getFuselageSectionCenterPoints(self):
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
            self.fuselagePoints = np.empty((self.nFuseSegm+1,3))
            for i in range(0,self.nFuseSegm):
                name = self.fuseSegmNames[i]
                point = self.tigl.fuselageGetSectionCenter(name,0.0)
                self.fuselagePoints[i] = point
            point = self.tigl.fuselageGetSectionCenter(name,1.0)
            self.fuselagePoints[-1] = point
            logger.info("Fuselage points are stored")

        elif self.nFuselage > 1:
            logger.warning("More than one fuselage found in the CPACS file,")
            logger.warning("only the first one is taken into account for this")
            logger.warning("simulation")

        else:
            logger.warning("No fuselage found in the CPACS file")

    def getWingsSectionCenterPoints(self):
        """
        Warning, the hypothesis of having a center of rotation which is not
        alligned with the center of pressure of the wing is done. This 
        hypothesis at the moment is such that the center of rotation is at
        1/3 of the chord.
        """
        logger.info("Starts computing each section of the wing central point")
        # TODO: add a parameter for the chord axis of rotation.
        
        chord_rotation = 1/3
        if self.nWings > 0:
            # Number of sections in each wing
            self.wingsNsegm = np.empty((self.nWings))
            self.wingsNsect = np.empty((self.nWings))
            for i in range(self.nWings):
                self.wingsNsegm[i] = self.tigl.wingGetSegmentCount(i+1)
                self.wingsNsect[i] = self.tigl.wingGetSectionCount(i+1)

            # Gets each section and segement name
            self.wingsSegmentNames = []  # np.empty((self.nWings,))
            self.wingsSectionsCenters = []
            # self.wingsSectionsNames = np.empty((self.nWings,))

            for i in range(self.nWings):
                # segmentUID = []
                N = int(self.wingsNsect[i])
                center = np.empty((N,3))
                for j in range(N-1):
                    # segmentUID.append(self.tigl.wingGetSegmentUID(i+1, j+1))
                    upper = self.tigl.wingGetUpperPoint(i+1, j+1, 0, chord_rotation)  # WARNING WITH THE 0
                    upper = np.array(upper)
                    lower = self.tigl.wingGetLowerPoint(i+1, j+1, 0, chord_rotation)
                    lower = np.array(lower)
                    center[j] = 0.5*(upper + lower)
                # Gets last part of each segment
                upper = self.tigl.wingGetUpperPoint(i+1, j+1, 1, chord_rotation)
                upper = np.array(upper)
                lower = self.tigl.wingGetLowerPoint(i+1, j+1, 1, chord_rotation)
                lower = np.array(lower)
                center[j+1] = 0.5*(upper + lower)

                # If there is a symmetry, a set of point symmetric to the plane
                # is added
                symmetry = self.tigl.wingGetSymmetry(i+1)
                if symmetry == 0:
                    self.wingsSectionsCenters.append(center)
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
                    logger.debug(center)
                    logger.debug(center_copy)
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
                    tab = np.concatenate((np.flip(center_copy[1:],axis=0),center))
                    logger.debug(tab)
                    self.wingsSectionsCenters.append(tab)

        else:
            logger.warning("No wings found in the CPACS file")

    def getFuselageProprieties(self):
        pass

    def getWingsPropriertes(self):
        pass

    def getHorizontalTailProprieties(self):
        pass

    def getVerticalTailProprieties(self):
        pass

    def getAllPoints(self):
        """
        Calls all the necessary functions to build a mesh
        """
        self.getOrigin()
        self.getWingsSectionCenterPoints()
        self.getFuselageSectionCenterPoints()

    def plotPoints(self):
        N = len(self.wingsSectionsCenters)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(self.origin[0],self.origin[1],self.origin[2])
        ax.scatter(self.fuselagePoints[:,0],
                   self.fuselagePoints[:,1],
                   self.fuselagePoints[:,2])
        for i in range(N):
            ax.scatter(self.wingsSectionsCenters[i][:,0],
                       self.wingsSectionsCenters[i][:,1],
                       self.wingsSectionsCenters[i][:,2])
        ax.set_xlim(-5,25)
        ax.set_ylim(-15,15)
        ax.set_zlim(-15,15)
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')        
        plt.show()
