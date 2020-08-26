#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:54:08 2020

@author: Jean-Philippe Kuntzer

# OK TODO: clamp the wings
# TODO: correct the sections
# TODO: correct the mesh points
# TODO: correct force
"""
import logging
import numpy as np
from framat import Model
import scipy as sp
import sys


logger = logging.getLogger(__name__)


class framat:

    def __init__(self,geometry):  # propertiesPath,csd_mesh):
        self.model = Model()
        self.geo = geometry

    def mesh(self):
        # Does all the necessary step to get a useful mesh
        logger.info("FramAT meshing process started")
        self.loadMaterials()
        self.loadGeometryPropertiesFromJson()
        self.computesCsdMesh()
        logger.info("FramAT mesh completed")

    def run(self,tranform):
        """
        Function that applies BC, apply the external loads.
        """
        self.imposeBC()
        self.applysLoad(tranform)
        # TODO add a user input if he wants or not to see the results
        # self.postProcessing()
        logger.debug("Framat solver stars computing")
        logger.debug(tranform.afx)
        logger.debug(self.geo)
        self.results = self.model.run()
        logger.debug("Framat solver finised computing")

    def checkResults(self):
        logger.debug(self.results[""])

    def loadMaterials(self):
        logger.error(self.geo.aircraftBeamsMaterials)
        N = self.geo.nFuselage + self.geo.nWings

        for i in range(N):
            name = self.geo.aircraftBeamsMaterials[i][0] + "_mat"
            mat = self.model.add_feature('material', uid=name)
            mat.set('E', self.geo.aircraftBeamsMaterials[i][1])
            mat.set('G', self.geo.aircraftBeamsMaterials[i][2])
            mat.set('rho', self.geo.aircraftBeamsMaterials[i][3])
            logger.debug("mat E = "+str(mat.get("E")))
            logger.debug("mat G = "+str(mat.get("G")))
            logger.debug("mat rho = "+str(mat.get("rho")))
            logger.error(self.geo.aircraftBeamsMaterials[i][1])
            logger.error(self.geo.aircraftBeamsMaterials[i][2])
            logger.error(self.geo.aircraftBeamsMaterials[i][3])

    def loadGeometryPropertiesFromJson(self):
        # TODO add versatility to the cross_section part
        # TODO load from json file
        N = len(self.geo.aircraftNodesPoints)
        cs = []
        for i in range(N):
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                name = self.geo.aircraftNodesNames[i][j] + "_cross_section"
                logger.debug(name)
                A = self.geo.aircraftNodesA[i][j]
                Iy = self.geo.aircraftNodesIy[i][j]
                Iz = self.geo.aircraftNodesIz[i][j]
                J = self.geo.aircraftNodesJ[i][j]
                cs = self.model.add_feature('cross_section', uid=name)
                cs.set('A', A)
                cs.set('Iy', Iy)
                cs.set('Iz', Iz)
                cs.set('J', J)
                logger.debug("uid = "+str(name))
                logger.debug("Iy = "+str(Iy))
                logger.debug("Iz = "+str(Iz))
                logger.debug("J = "+str(J))

    def computesCsdMesh(self):
        """
        """
        self.beams = []
        np.set_printoptions(precision=3)
        # Number of beams
        N = len(self.geo.aircraftNodesPoints)
        for i in range(N):
            self.beams.append(self.model.add_feature('beam'))
            # Number of nodes in the current beam
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                # adds the points (nodes) to the beam
                point = self.geo.aircraftNodesPoints[i][j].tolist()
                logger.debug("point = "+str(point))
                name = self.geo.aircraftNodesNames[i][j]
                logger.debug("name = "+str(name))
                self.beams[i].add("node",point,uid=name)

                # WARNING this bit of code should be indented one more time.
                # It should be nested inside the second for loop
                # Adds the cross_section feature to each points
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

            # Sets beam number of elements
            self.beams[i].set('nelem', 1)
            uid = self.geo.aircraftBeamsMaterials[i][0] + "_mat"
            a = self.geo.aircraftNodesNames[i][0]
            b = self.geo.aircraftNodesNames[i][-1]
            logger.debug("Material uid = "+str(uid))
            self.beams[i].add('material', {'from': a, 'to': b, 'uid': uid})
            self.beams[i].add('orientation',{'from': a, 'to': b, 'up': [0, 0, 1]})

    def applysLoad(self,tranform):
        """
        WARNING We should pass the accepted names variable into the class
        """
        # Number of beams
        N = self.geo.nFuselage + self.geo.nWings
        for i in range(N):
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                name = self.geo.aircraftNodesNames[i][j]

                # Distributes loads due to inertia. m__ for mass, _f_ for
                # force, _m_ for moment, __* direction
                if self.geo.settings['G_loads']:
                    mfx = tranform.smf[i][j,0]
                    mfy = tranform.smf[i][j,1]
                    mfz = tranform.smf[i][j,2]
                    mmx = tranform.smm[i][j,0]
                    mmy = tranform.smm[i][j,1]
                    mmz = tranform.smm[i][j,2]
                else:
                    mfx = 0
                    mfy = 0
                    mfz = 0
                    mmx = 0
                    mmy = 0
                    mmz = 0

                # Distributes loads due to aerodynamics if CFD solver is Pytornado
                if (self.geo.settings['CFD_solver'] == 'Pytornado'
                    and i >= self.geo.nFuselage
                    ):
                    fx = tranform.sfx[i - self.geo.nFuselage][j]
                    fy = tranform.sfy[i - self.geo.nFuselage][j]
                    fz = tranform.sfz[i - self.geo.nFuselage][j]
                    mx = tranform.smx[i - self.geo.nFuselage][j]
                    my = tranform.smy[i - self.geo.nFuselage][j]
                    mz = tranform.smz[i - self.geo.nFuselage][j]
                elif self.geo.settings['CFD_solver'] == 'Pytornado':
                    fx = 0
                    fy = 0
                    fz = 0
                    mx = 0
                    my = 0
                    mz = 0
                logger.debug(fx)
                logger.debug(fz)
                # Distributes loads due to aerodynamics if CFD solver is SU2
                if self.geo.settings['CFD_solver'] == 'SU2':
                    fx = tranform.sfx[i][j]
                    fy = tranform.sfy[i][j]
                    fz = tranform.sfz[i][j]
                    mx = tranform.smx[i][j]
                    my = tranform.smy[i][j]
                    mz = tranform.smz[i][j]
                load = [fx+mfx, fy+mfy, fz+mfz, mx+mmx, my+mmy, mz+mmz]
                self.beams[i].add('point_load', {'at': name, 'load': load})
                # for removing them at the end of the simulation but keeping
                # the same mesh for computation time efficiency.
                self.minusLoads = []
                for t in load:
                    self.minusLoads.append(-2*t)
                # self.minusLoads = [-fx-mfx, -fy-mfy, -fz-mfz, -mx-mmx, -my-mmy, -mz-mmz]

    def eraseLoad(self,tranform):
        N = self.geo.nFuselage + self.geo.nWings
        for i in range(N):
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                name = self.geo.aircraftNodesNames[i][j]
                self.beams[i].add('point_load', {'at': name, 'load': self.minusLoads})

    def imposeBC(self):
        # ===== BOUNDARY CONDITIONS =====
        bc = self.model.set_feature('bc')
        if self.geo.nFuselage > 0:
            name = "f_n_clamped"
            logger.debug(self.geo.aircraftNodesNames)

            # # Connect nodes
            # Nc = len(self.geo.aircraftConnectedNodes[0])
            # logger.debug(self.geo.aircraftConnectedNodes[0])
            # for i in range(Nc):
            #     logger.debug(self.geo.aircraftConnectedNodes[0][i])
            #     beamIndex1 = int(self.geo.aircraftConnectedNodes[0][i][0])
            #     beamIndex1Node = int(self.geo.aircraftConnectedNodes[0][i][2])
            #     beamIndex2 = int(self.geo.aircraftConnectedNodes[0][i][1])
            #     beamIndex2Node = int(self.geo.aircraftConnectedNodes[0][i][3])
            #     logger.debug(beamIndex1)
            #     logger.debug(beamIndex1Node)
            #     logger.debug(beamIndex2)
            #     logger.debug(beamIndex2Node)
            #     name1 = self.geo.aircraftNodesNames[beamIndex1][beamIndex1Node]
            #     name2 = self.geo.aircraftNodesNames[beamIndex2][beamIndex2Node]
            #     bc.add('connect',{'node1': name1,
            #                       'node2': name2,
            #                       'fix':['all']})
            #     logger.debug("Connects node: " + name1)
            #     logger.debug("connects node: " + name2)
            #     logger.debug("="*30)
            
            # Connect nodes
            N = len(self.geo.aircraftConnectedNodes)
            logger.debug(self.geo.aircraftConnectedNodes)
            for i in range(N):
                # logger.debug(self.geo.aircraftConnectedNodes[0][i])
                beamIndex1 = int(self.geo.aircraftConnectedNodes[i][0])
                beamIndex1Node = int(self.geo.aircraftConnectedNodes[i][2])
                beamIndex2 = int(self.geo.aircraftConnectedNodes[i][1])
                beamIndex2Node = int(self.geo.aircraftConnectedNodes[i][3])
                logger.debug(beamIndex1)
                logger.debug(beamIndex1Node)
                logger.debug(beamIndex2)
                logger.debug(beamIndex2Node)
                name1 = self.geo.aircraftNodesNames[beamIndex1][beamIndex1Node]
                name2 = self.geo.aircraftNodesNames[beamIndex2][beamIndex2Node]
                bc.add('connect',{'node1': name1,
                                  'node2': name2,
                                  'fix':['all']})
                logger.debug("Connects node: " + name1)
                logger.debug("connects node: " + name2)
                logger.debug("="*30)
                
            
        else:
            name = "w_n_clamped"
        bc.add('fix', {'node': name, 'fix': ['all']})
        N = len(self.geo.aircraftNonRotatingNodes)
        # for i in range(N):
        #     name = self.geo.aircraftNonRotatingNodes[i]
        #     # bc.add('fix', {'node': name, 'fix': ['ty']})
        #     bc.add('fix', {'node': name, 'fix': ['ux','uy','uz]})

    def postProcessing(self):
        # ===== POST-PROCESSING =====
        # By default the analysis is run without any GUI, but to get a visual
        # representation of the results we can create a plot
        pp = self.model.set_feature('post_proc')
        pp.set('plot_settings', {'show': True})
        # pp.add('plot', ['undeformed', 'deformed', 'node_uids', 'nodes', 'forces'])
        pp.add('plot', ['undeformed', 'deformed', 'nodes', 'forces'])
