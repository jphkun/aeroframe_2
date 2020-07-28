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

logger = logging.getLogger(__name__)

class Csd_mesh:

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
        # Number of fuselage instances (generally one)
        self.nFuselage = self.tigl.getFuselageCount()
        # Number of wings
        self.nWings = self.tigl.getWingCount()
        logger.info("Number of fuselage instances found: " +
                    str(self.nFuselage))
        logger.info("Number of wing instances found:     " +
                    str(self.nWings))
    
    class fuselage:
        def __init__(fuselage):
            pass
    
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
        if self.nFuselage > 0:
            logger.info("Starts computing each section of the fuselage central point")
            self.nFuseSect = self.tigl.fuselageGetSectionCount(1)
            self.nFuseSegm = self.tigl.fuselageGetSegmentCount(1)
            self.fuseSegmNames = []
            # Gets each fuselage section UID and stores it
            for i in range(1,self.nFuseSegm + 1):
                self.fuseSegmNames.append(self.tigl.fuselageGetSegmentUID(1,i))

            # Computes each fuselage section central point
            self.fusePoints = np.empty((self.nFuseSegm+1,3))
            for i in range(0,self.nFuseSegm):
                name = self.fuseSegmNames[i]
                point = self.tigl.fuselageGetSectionCenter(name,0.0)
                self.fusePoints[i] = point
            point = self.tigl.fuselageGetSectionCenter(name,1.0)
            self.fusePoints[-1] = point
            logger.info("Fuselage points are stored")
        
        elif self.nFuselage > 1:
            logger.warning("More than one fuselage found in the CPACS file,")
            logger.warning("only the first one is taken into account for this")
            logger.warning("simulation")
        
        else:
            logger.warning("No fuselage found in the CPACS file")
    
    def getWingsSectionCenterPoints(self):
        # TODO: make it work properly at this point there is a bug with the
        # symmetry.
        # Tomorrow:
        # 1) plot all the points
        # 2) solve the wing symmetry issue
        if self.nWings > 0:
            # Number of sections in each wing
            self.wingsNsegm = np.empty((self.nWings))
            self.wingsNsect = np.empty((self.nWings))
            for i in range(self.nWings):
                self.wingsNsegm[i] = self.tigl.wingGetSegmentCount(i+1)
                self.wingsNsect[i] = self.tigl.wingGetSectionCount(i+1)
            logger.debug(self.wingsNsegm)
            logger.debug(self.wingsNsect)

            # Gets each section and segement name
            self.wingsSegmentNames = []  # np.empty((self.nWings,))
            self.wingsSectionsCenters = []
            # self.wingsSectionsNames = np.empty((self.nWings,))
            
            for i in range(self.nWings):
                segmentUID = []
                center = []
                for j in range(int(self.wingsNsegm[i])):
                    # segmentUID.append(self.tigl.wingGetSegmentUID(i+1, j+1))
                    upper = self.tigl.wingGetUpperPoint(i+1, j+1, 1, 0.25)
                    lower = self.tigl.wingGetLowerPoint(i+1, j+1, 1, 0.25)
                    center.append([0.5*(upper[0]+lower[0]),
                                   0.5*(upper[1]+lower[1]),
                                   0.5*(upper[2]+lower[2])])
                upper = self.tigl.wingGetUpperPoint(i+1, j+1, 1, 0.25)
                lower = self.tigl.wingGetLowerPoint(i+1, j+1, 1, 0.25)
                center.append([0.5*(upper[0]+lower[0]),
                               0.5*(upper[1]+lower[1]),
                               0.5*(upper[2]+lower[2])])
                self.wingsSectionsCenters.append(center)
                
                # logger.debug("Wing Symmetry "+str(self.tigl.wingGetSymmetry(i+1)))
                
                # If there is a symmetry, a set of point symmetric to the plane
                # is added
                symmetry = self.tigl.wingGetSymmetry(i+1)
                if symmetry > 0:
                    for i in range(int(self.wingsNsegm[i])):
                        logger.debug(len(center))
                        logger.debug(symmetry)
                        logger.debug(center[i][symmetry]) #= -center[i,symmetry]
                    self.wingsSectionsCenters.append(center)
            
            logger.debug(self.wingsSectionsCenters)
            
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
        
        # TODO find a way to take into account multiple wings
        # self.getWingPoints
