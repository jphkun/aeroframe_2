#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:54:08 2020

@author: Jean-Philippe Kuntzer

# TODO: clamp the wings
# TODO: correct the sections
"""
import logging
import numpy as np
from framat import Model
import scipy as sp


logger = logging.getLogger(__name__)


class framat:

    def __init__(self,geometry):  # propertiesPath,csd_mesh):
        self.model = Model()
        self.geo = geometry

    def csdRun(self):
        # Run the beam analysis
        self.loadMaterials()
        self.loadGeometryPropertiesFromJson()
        self.computesCsdMesh()
        self.applysLoad()
        self.imposeBC()
        self.postProcessing()
        self.results = self.model.run()

    def loadMaterials(self):
        logger.debug(self.geo.aircraftBeamsMaterials)
        N = self.geo.nFuselage + self.geo.nWings
        mat = []
        
        for i in range(N):
            name = self.geo.aircraftBeamsMaterials[i][0] + "_mat"
            mat.append(self.model.add_feature('material', uid=name))
            mat[i].set('E', self.geo.aircraftBeamsMaterials[i][1])
            mat[i].set('G', self.geo.aircraftBeamsMaterials[i][2])
            mat[i].set('rho', self.geo.aircraftBeamsMaterials[i][3])
        # logger.debug(mat[0])

    def loadGeometryPropertiesFromJson(self):
        # TODO add versatility to the cross_section part
        # TODO load from json file
        N = len(self.geo.aircraftNodesPoints)
        for i in range(N):
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                name = self.geo.aircraftNodesNames[i][j] + "_cross_section"
                logger.debug(name)
                A = self.geo.aircraftNodesA[i][j]
                Iy = self.geo.aircraftNodesA[i][j]
                Iz = self.geo.aircraftNodesA[i][j]
                J = self.geo.aircraftNodesA[i][j]
                cs = self.model.add_feature('cross_section', uid=name)
                cs.set('A', A)
                cs.set('Iy', Iy)
                cs.set('Iz', Iz)
                cs.set('J', J)

    def computesCsdMesh(self):
        """
        """
        self.beams = []
        
        # Number of beams
        N = len(self.geo.aircraftNodesPoints)
        for i in range(N):
            self.beams.append(self.model.add_feature('beam'))
            # Number of nodes in the current beam
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                # adds point
                point = self.geo.aircraftNodesPoints[i][j].tolist()
                # logger.debug(point)
                name = self.geo.aircraftNodesNames[i][j]
                self.beams[i].add("node",point,uid=name)
                if j+1 < int(np.floor(M/2)) + 2:
                    
                    name1 = self.geo.aircraftNodesNames[i][j]
                    name2 = self.geo.aircraftNodesNames[i][j+1]
                    uid = name2 + "_cross_section"
                    self.beams[i].add('cross_section', {'from': name1,
                                                        'to': name2,
                                                        'uid': uid})
                else:
                    name1 = self.geo.aircraftNodesNames[i][j-1]
                    name2 = self.geo.aircraftNodesNames[i][j]
                    uid = name1 + "_cross_section"
                    self.beams[i].add('cross_section', {'from': name1,
                                                        'to': name2,
                                                        'uid': uid})
                
                # self.beams[i].add('cross_section', {'at': name,'uid': uid})
            # Sets beam number of elements
            self.beams[i].set('nelem', M)
            uid = self.geo.aircraftBeamsMaterials[i][0] + "_mat"
            a = self.geo.aircraftNodesNames[i][0]
            b = self.geo.aircraftNodesNames[i][-1]
            self.beams[i].add('material', {'from': a, 'to': b, 'uid': uid})
            self.beams[i].add('orientation',{'from': a, 'to': b, 'up': [0, 0, 1]})

    def applysLoad(self):
        if self.geo.nFuselage > 0:
            name1 = self.geo.aircraftNodesNames[1][0]
            name2 = self.geo.aircraftNodesNames[1][-1]
            index = 1
            self.beams[index].add('distr_load', {'from': name1,
                                'to': name2,
                                'load': [0, 0, -2e-3, 0, 0, 0]})
            name1 = self.geo.aircraftNodesNames[2][0]
            name2 = self.geo.aircraftNodesNames[2][-1]
            index = 2
            self.beams[index].add('distr_load', {'from': name1,
                                'to': name2,
                                'load': [0, 0, 0.5, 0, 0, 0]})
            name1 = self.geo.aircraftNodesNames[3][0]
            name2 = self.geo.aircraftNodesNames[3][-1]
            index = 3
            self.beams[index].add('distr_load', {'from': name1,
                                'to': name2,
                                'load': [0, 0, 0, 0, 0, 0]})
        else:
            name1 = self.geo.aircraftNodesNames[0][0]
            name2 = self.geo.aircraftNodesNames[0][-1]
            index = 0
            self.beams[index].add('distr_load', {'from': name1,
                                'to': name2,
                                'load': [0, 0, 2, 0, 0, 0]})

    def imposeBC(self):
        # ===== BOUNDARY CONDITIONS =====
        bc = self.model.set_feature('bc')
        # name = self.geo.aircraftConnectedNodes
        # logger.debug(name)
        if self.geo.nFuselage > 0:
            name = "f_n_clamped"
            logger.debug(self.geo.aircraftNodesNames)
            
            # Connect nodes
            Nc = len(self.geo.aircraftConnectedNodes[0])
            # logger.debug(self.geo.aircraftConnectedNodes)
            logger.debug(self.geo.aircraftConnectedNodes[0])
            for i in range(Nc):
                logger.debug(self.geo.aircraftConnectedNodes[0][i])
                beamIndex1 = int(self.geo.aircraftConnectedNodes[0][i][0])
                beamIndex1Node = int(self.geo.aircraftConnectedNodes[0][i][2])
                beamIndex2 = int(self.geo.aircraftConnectedNodes[0][i][1])
                beamIndex2Node = int(self.geo.aircraftConnectedNodes[0][i][3])
                logger.debug(beamIndex1)
                logger.debug(beamIndex1Node)
                logger.debug(beamIndex2)
                logger.debug(beamIndex2Node)
                name1 = self.geo.aircraftNodesNames[beamIndex1][beamIndex1Node]
                name2 = self.geo.aircraftNodesNames[beamIndex2][beamIndex2Node]
                bc.add('connect',{'node1': name1,
                                  'node2': name2,
                                  'fix':['all']})
                logger.debug("Connects node"+name1)
                logger.debug("connects node"+name2)
                logger.debug("="*30)
        else:
            name = "w_n_clamped"
        bc.add('fix', {'node': name, 'fix': ['all']})

    def postProcessing(self):
        # ===== POST-PROCESSING =====
        # By default the analysis is run without any GUI, but to get a visual
        # representation of the results we can create a plot
        pp = self.model.set_feature('post_proc')
        pp.set('plot_settings', {'show': True})
        pp.add('plot', ['undeformed', 'deformed', 'node_uids', 'nodes', 'forces'])
        
