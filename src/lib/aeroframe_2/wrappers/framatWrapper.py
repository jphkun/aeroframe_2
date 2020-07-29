#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:54:08 2020

@author: Jean-Philippe Kuntzer
"""
import logging
import numpy as np
from framat import Model
import scipy as sp


logger = logging.getLogger(__name__)


class framat:

    def __init__(self,geometry):  # propertiesPath,csd_mesh):
        self.model = Model()
        self.geometryInput = geometry

    def loadMaterials(self):
        # TODO add versatility to the material part
        # TODO load from json file
        mat = self.model.add_feature('material', uid='dummy')
        mat.set('E', 1)
        mat.set('G', 1)
        mat.set('rho', 1)

    def loadGeometryPropertiesFromJson(self):
        # TODO add versatility to the cross_section part
        # TODO load from json file
        cs = self.model.add_feature('cross_section', uid='dummy')
        cs.set('A', 1)
        cs.set('Iy', 1)
        cs.set('Iz', 1)
        cs.set('J', 1)

    def computesCsdMesh(self):

        # OK TODO verify if there is any fuselage
        # OK TODO verify if there is any HTW
        # OK TODO verify if there is any VTW
        # OK TODO clamp the node closer to the origin

        # Tomorrow
        # TODO verify if there is any wing, throw error if nothing is found  
        # TODO clamp wings together if wings are close enough
        # TODO clamp the biggest wing if there is no fuselage
        # TODO map loads on the wings
        # TODO add nodes of find a way to map the load on the elements
        # TODO read materials from file
        # TODO read mechanical properties from file
        # TODO add versatility
        # TODO add lisibility
        
        # BONUS
        # TODO clamp together boxwing if necessary
        # TODO handel multiple fuselages case


        
        self.nFuselage = self.geometryInput.nFuselage
        self.nWings = self.geometryInput.nWings
        fuselagePoints = self.geometryInput.fuselagePoints
        wingsPoints = self.geometryInput.wingsSectionsCenters
        self.origin = self.geometryInput.origin
        # ============================
        #   Fuselage
        # ============================
        # Tests if there is a fuselage
        if self.nFuselage == 1:
            logger.info("Number of fuselages found: "+str(self.nFuselage))

            # Generates an instance for the fuselage
            self.fuselage = self.model.add_feature('beam')
            self.fuselageNnodes = len(fuselagePoints)
            self.distancesFromOg = np.empty((self.fuselageNnodes))

            # adds nodes to the beam
            for i in range(self.fuselageNnodes):
                difference = self.origin-fuselagePoints[i]
                self.distancesFromOg[i] = np.linalg.norm(difference)
                # Computes the distance from the origin
                self.fuselage.add('node',
                                  fuselagePoints[i].tolist(),
                                  uid='fuselageNode' + str(i))
            # logger.debug(self.distancesFromOg)
            self.bcNode = np.argmin(self.distancesFromOg)
            # logger.debug(self.bcNode)
            # Set the number of elements for the fuselage beam.
            self.fuselage.set('nelem', 14)

            # Adds features to the fuselage beam
            self.fuselage.add('material', {'from': 'fuselageNode0',
                                           'to': 'fuselageNode'+str(i),
                                           'uid': 'dummy'})
            self.fuselage.add('cross_section', {'from': 'fuselageNode0',
                                                'to': 'fuselageNode'+str(i),
                                                'uid': 'dummy'})
            # TODO What does the up direction means??
            self.fuselage.add('orientation', {'from': 'fuselageNode0',
                                              'to': 'fuselageNode'+str(i),
                                              'up': [0, 0, 1]})
        # ============================
        #   Wings
        # ============================
        # TODO, find a was to see if two wings are linked one to the other.
        if self.nWings > 0:
            logger.info("Number of wings found:     "+str(self.nWings))
            self.wingsBeams = []
            self.indecesNodesToConnect = []
            self.instancesToConnectTo = []
            self.indexFuselage = []
            # iterates through all the wings
            for i in range(self.nWings):

                # Recovers the current wing points WARING BE CAREFUL WITH THE
                # PLURAL wingssssPointsss and wingPointssss
                wingPoints = wingsPoints[i]
                wingNnodes = len(wingPoints)

                # Generate a beam instances
                self.wingsBeams.append(self.model.add_feature('beam'))

                ###############################################################
                # Computes the distance between the current wing and the
                # fuselage
                ###############################################################
                distanceToFuselage = sp.spatial.distance.cdist(wingsPoints[i],
                                                      fuselagePoints,
                                                      metric='euclidean')
                # Computes the distance between the current wing and all the
                # other wings
                # TODO correct the bug with the 0

                # distanceToWings = np.empty((len(wingsPoints[i]), len(wingsPoints)))
                # it = 0
                # for k in range(self.nWings):
                #     for l in range(len(wingsPoints[i])):
                #         for m in range(len(wingsPoints[k])):
                #             # Avoids auto-conncetion
                #             if i == k:
                #                 dist = 1e6
                #             else:
                #                 dist = np.linalg.norm(wingsPoints[i][l] - wingsPoints[k][m]
                #             distanceToWings[l,m] =
                # logger.debug(distance)
                
                # Finds the index of the point of minimal distance between the
                # current wing, the fuselage and the other wings
                logger.debug(distanceToFuselage)
                self.indexFuselage.append(np.where(distanceToFuselage == np.amin(distanceToFuselage)))
                logger.debug(self.indexFuselage)
                # indexWings = np.where(distanceToWings == np.amin(distanceToWings))

                # Finds the absolute minimal distance for the connection of the
                # current wing and the "outside beams".
                # indexFuselage = np.array([indexFuselage[0][0],indexFuselage[1][0]])
                # logger.debug(indexWings)
                # indexWings = np.array([indexWings[0][0],indexWings[1][0]])
                # logger.debug(indexFuselage)
                # logger.debug(indexWings)
                # logger.debug(distanceToFuselage)
                # logger.debug(distanceToWings)
                # if distanceToFuselage(indexFuselage[1]) <= distanceToWings(indexWings):
                #     self.indecesNodesToConnect.append(indexFuselage)
                #     self.instancesToConnectTo.append("fuselage")
                # else:
                #     self.indecesNodesToConnect.append(indexWings)
                #     self.instancesToConnecTo.append("wing")

                # Inserts all the known points of the mesh by creating the
                # nodes. The nodes are placed at a certain distance of the
                # leading edge this distance is imposed ahead in the code.
                for j in range(wingNnodes):
                    self.wingsBeams[i].add('node', wingPoints[j].tolist(),
                                            uid='wing' + str(i)+'node' + str(j))
                # Inserts the elements inside the current wing beam.
                self.wingsBeams[i].set('nelem', 10)
                # Adds features to the current wing beam
                self.wingsBeams[i].add('material',
                                        {'from': 'wing' + str(i)+'node0',
                                        'to': 'wing' + str(i) + 'node' + str(wingNnodes-1),
                                        'uid': 'dummy'})
                self.wingsBeams[i].add('cross_section',
                                        {'from': 'wing' + str(i)+'node0',
                                        'to': 'wing' + str(i) + 'node' + str(wingNnodes-1),
                                        'uid': 'dummy'})
                self.wingsBeams[i].add('orientation',
                                        {'from': 'wing' + str(i)+'node0',
                                        'to': 'wing' + str(i) + 'node' + str(wingNnodes-1),
                                        'up': [0, 0, 1]})

    def applysLoad(self):
        # Add some line loads [N/m] and point loads [N]
        if self.nFuselage == 1:
            self.fuselage.add('distr_load',
                              {'from': 'fuselageNode0',
                               'to': 'fuselageNode' + str(self.fuselageNnodes-1),
                               'load': [0, 0, 5e-3, 0, 0, 0]})

    def imposeBC(self):
        # ===== BOUNDARY CONDITIONS =====
        # We also must constrain our model. Below, we fix the nodes 'a' and 'd'
        bc = self.model.set_feature('bc')
        if self.nFuselage == 1:
            bc.add('fix', {'node': 'fuselageNode' + str(self.bcNode), 'fix': ['all']})
            # bc.add('connect',{'node1': '', 'node2': '', 'fix':['all']})
            for i in range(self.nWings):
                logger.debug("fuselage node: " + str(self.indexFuselage[i][1][0]))
                logger.debug("wing node:     " + str(self.indexFuselage[i][0][0]))
                fuseNode = int(self.indexFuselage[i][1][0])
                wingNode = int(self.indexFuselage[i][0][0])
                logger.debug(fuseNode)
                logger.debug(wingNode)
                bc.add('connect',{'node1': 'fuselageNode' + str(fuseNode),
                                  'node2': 'wing' + str(i)+'node' + str(wingNode),
                                  'fix':['all']})

    def postProcessing(self):
        # ===== POST-PROCESSING =====
        # By default the analysis is run without any GUI, but to get a visual
        # representation of the results we can create a plot
        pp = self.model.set_feature('post_proc')
        pp.set('plot_settings', {'show': True})
        pp.add('plot', ['undeformed', 'deformed', 'node_uids', 'nodes', 'forces'])

    def csdRun(self):
        # Run the beam analysis
        self.loadMaterials()
        self.loadGeometryPropertiesFromJson()
        self.computesCsdMesh()
        self.applysLoad()
        self.imposeBC()
        self.postProcessing()
        self.results = self.model.run()
        
