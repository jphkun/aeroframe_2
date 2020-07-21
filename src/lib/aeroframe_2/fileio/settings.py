#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 11:28:11 2020

@author: Jean-Philippe Kuntzer
"""
import sys
import logging

logger = logging.getLogger(__name__)

class settings():
    
    def __init__(self,dictionary):
        logger.debug("Initialisazion of the settings file")
        self.CFD_solver = dictionary["CFD_solver"]
        self.CSD_solver = dictionary["CSD_solver"]
    
    def verify(self):
        """
        Verifies if the input file is complete enougth.
        TODO: Add a python test to see if corresponding module is installed.
        TODO: Add a verification to see if SU2 is installed.
        """
        cfd_solver = self.CFD_solver
        csd_solver = self.CSD_solver
        if cfd_solver != "Pytornado" and cfd_solver != "SU2":
            logger.error("CFD solver choice must be \"Pytornado\" or \"SU2\"")
            sys.exit()
        if csd_solver != "FramAT" and \
           csd_solver != "FeniCS" and \
           csd_solver != "External":
            logger.error("CFD solver choice must be \"FramAT\", \"FeniCS\" or \"External\"")
            sys.exit()
