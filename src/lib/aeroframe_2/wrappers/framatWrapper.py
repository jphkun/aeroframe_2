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
        # logger.debug("Framat solver stars computing")
        # logger.debug(tranform.afx)
        # logger.debug(self.geo)
        self.results = self.model.run()
        var = self.results.get('tensors').get('comp:U')["thx"]
        # for i in var:
        #     logger.debug(i)
        
        # var = self.results.get('tensors').get('F')[3::6]
        # for i in var:
        #     logger.debug(i)
        # var = self.model.get('beam')[0].get('cross_section')
        # for i in var:
        #     logger.debug(i)
            
        # logger.debug(self.model)
        # logger.debug(self.model.get('beam')[0].get('material'))
        # logger.debug())
        # logger.debug(self.model.get('beam')[0].get('node'))
        # sys.exit()
        logger.debug("Framat solver finised computing")

    # def checkResults(self):
        # logger.debug(self.results[""])

    def loadMaterials(self):
        logger.error(self.geo.aircraftBeamsMaterials)
        N = self.geo.nFuselage + self.geo.nWings

        for i in range(N):
            name = self.geo.aircraftBeamsMaterials[i][0] + "_mat"
            mat = self.model.add_feature('material', uid=name)
            mat.set('E', self.geo.aircraftBeamsMaterials[i][1])
            mat.set('G', self.geo.aircraftBeamsMaterials[i][2])
            mat.set('rho', self.geo.aircraftBeamsMaterials[i][3])
            # logger.debug("mat E = "+str(mat.get("E")))
            # logger.debug("mat G = "+str(mat.get("G")))
            # logger.debug("mat rho = "+str(mat.get("rho")))
            # logger.error(self.geo.aircraftBeamsMaterials[i][1])
            # logger.error(self.geo.aircraftBeamsMaterials[i][2])
            # logger.error(self.geo.aircraftBeamsMaterials[i][3])

    def loadGeometryPropertiesFromJson(self):
        # TODO add versatility to the cross_section part
        # TODO load from json file
        N = len(self.geo.aircraftNodesPoints)
        cs = []
        for i in range(N):
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                name = self.geo.aircraftNodesNames[i][j] + "_cross_section"
                # logger.debug(name)
                A = self.geo.aircraftNodesA[i][j]
                Iy = self.geo.aircraftNodesIy[i][j]
                Iz = self.geo.aircraftNodesIz[i][j]
                J = self.geo.aircraftNodesJ[i][j]
                cs = self.model.add_feature('cross_section', uid=name)
                cs.set('A', np.round_(A , decimals=4))
                cs.set('Iy',np.round_(Iy, decimals=9))
                cs.set('Iz',np.round_(Iz, decimals=9))
                cs.set('J', np.round_(J , decimals=9))
                # logger.debug("uid = "+str(name))
                # logger.debug("Iy = "+str(Iy))
                # logger.debug("Iz = "+str(Iz))
                # logger.debug("J = "+str(np.round_(J , decimals=9)))
                # logger.debug(cs.get('Iy'))
            # sys.exit()

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
                # adds the points (nodes) to the beam
                point = self.geo.aircraftNodesPoints[i][j].tolist()
                point[0] = np.round_(point[0], decimals=3)
                point[1] = np.round_(point[1], decimals=3)
                point[2] = np.round_(point[2], decimals=3)
                # logger.debug("point = "+str(point))
                name = self.geo.aircraftNodesNames[i][j]
                # logger.debug("name = "+str(name))
                self.beams[i].add("node",point,uid=name)
                # logger.debug(str(j) + ' ' + str(point))
            # Sets beam number of elements
            self.beams[i].set('nelem', 2)
            uid = self.geo.aircraftBeamsMaterials[i][0] + "_mat"
            a = self.geo.aircraftNodesNames[i][0]
            b = self.geo.aircraftNodesNames[i][-1]
            # logger.debug("Material uid = "+str(uid))
            # logger.debug()
            self.beams[i].add('material', {'from': a, 'to': b, 'uid': uid})
            self.beams[i].add('orientation',{'from': a, 'to': b, 'up': [0, 0, 1]})
        del(i)
        del(j)

        for i in range(N):
            # self.beams.append(self.model.add_feature('beam'))
            # Number of nodes in the current beam
            M = len(self.geo.aircraftNodesPoints[i])
            for j in range(M):
                # Adds the cross_section feature to each points
                if j < int(np.floor(M/2)):
                    name1 = self.geo.aircraftNodesNames[i][j]
                    name2 = self.geo.aircraftNodesNames[i][j+1]
                    uid = name2 + "_cross_section"
                    # logger.debug('Before '+name1+' '+name2+' '+uid)
                    self.beams[i].add('cross_section', {'from': name1,
                                                        'to':   name2,
                                                        'uid':  uid})
                elif j < M-1:
                    # M - 1 in order to add the last point
                    name1 = self.geo.aircraftNodesNames[i][j]
                    name2 = self.geo.aircraftNodesNames[i][j+1]
                    uid = name1 + "_cross_section"
                    # logger.debug('After '+name1+' '+name2+' '+uid)
                    self.beams[i].add('cross_section', {'from': name1,
                                                        'to':   name2,
                                                        'uid':  uid})

            # sys.exit()
            


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
                    if j < int(np.floor(M/2)):
                        coef = 1
                    else:
                        coef = 1
                    mfx = np.round(tranform.smf[i][j,0])
                    mfy = np.round(tranform.smf[i][j,1])
                    mfz = np.round(tranform.smf[i][j,2])
                    mmx = coef * np.round(tranform.smm[i][j,0])
                    mmy = coef * np.round(tranform.smm[i][j,1])
                    mmz = coef * np.round(tranform.smm[i][j,2])
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
                    if j < int(np.floor(M/2)):
                        coef = 1
                    else:
                        coef = 1
                    fx = np.round(tranform.sfx[i - self.geo.nFuselage][j])
                    fy = np.round(tranform.sfy[i - self.geo.nFuselage][j])
                    fz = np.round(tranform.sfz[i - self.geo.nFuselage][j])
                    mx = np.round(tranform.smx[i - self.geo.nFuselage][j])
                    my = coef * np.round(tranform.smy[i - self.geo.nFuselage][j])
                    mz = np.round(tranform.smz[i - self.geo.nFuselage][j])
                elif self.geo.settings['CFD_solver'] == 'Pytornado':
                    fx = 0
                    fy = 0
                    fz = 0
                    mx = 0
                    my = 0
                    mz = 0
                # logger.debug(fx)
                # logger.debug(fz)
                # Distributes loads due to aerodynamics if CFD solver is SU2
                if self.geo.settings['CFD_solver'] == 'SU2':
                    if j < int(np.floor(M/2)):
                        coef = 1
                    else:
                        coef = 1
                    fx = np.round(tranform.sfx[i][j])
                    fy = np.round(tranform.sfy[i][j])
                    fz = np.round(tranform.sfz[i][j])
                    mx = np.round(tranform.smx[i][j])
                    my = coef * np.round(tranform.smy[i][j])
                    mz = np.round(tranform.smz[i][j])
                # Warning the frame of reference of the structure mesh is
                # rotated from the one of the airplane so:
                # airplaine x -> structure y
                # airplaine y -> structure x
                # airplaine z -> structure z
                if i >= self.geo.nFuselage:
                    load = [(fx-mfx), (fy-mfy), (fz-mfz), (mx-mmx), (my-mmy), (mz-mmz)]
                    self.beams[i].add('point_load', {'at': name, 'load': load})
                    
                    # logger.debug('Forces : '+ str(fx)  + ' ' + str(fy)  + ' ' + str(fz))
                    # logger.debug('Forces : '+ str(mfx) + ' ' + str(mfy) + ' ' + str(mfz))
                    # logger.debug('load: ' + str(load))
                    # logger.debug('Moments: '+ str( mx) + ' ' + str( my) + ' ' + str( mz))
                    # logger.debug('Moments: '+ str(mmx) + ' ' + str(mmy) + ' ' + str(mmz))

                
                # else:
                #     # Fuselage
                #     load = [0*fx+mfx, 0*fy+mfy, 0*fz+mfz, 0*mx+mmx, 0*my+mmy, 0*mz+mmz]
                #     self.beams[i].add('point_load', {'at': name, 'load': load})
                    
                # load = [0*fx+mfx, 0*fy+mfy, 0*fz+mfz, mx+mmx, my+mmy, mz+mmz]
                # self.beams[i].add('point_load', {'at': name, 'load': load})
                
            logger.debug('iteration finised')
            # sys.exit()
            # sys.exit()
                # for removing them at the end of the simulation but keeping
                # the same mesh for computation time efficiency.
                # self.minusLoads = []
                # for t in load:
                #     self.minusLoads.append(-2*t)
                # self.minusLoads = [-fx-mfx, -fy-mfy, -fz-mfz, -mx-mmx, -my-mmy, -mz-mmz]

    # def eraseLoad(self,tranform):
    #     N = self.geo.nFuselage + self.geo.nWings
    #     for i in range(N):
    #         M = len(self.geo.aircraftNodesPoints[i])
    #         for j in range(M):
    #             name = self.geo.aircraftNodesNames[i][j]
    #             self.beams[i].add('point_load', {'at': name, 'load': self.minusLoads})

    def imposeBC(self):
        # ===== BOUNDARY CONDITIONS =====
        bc = self.model.set_feature('bc')
        if self.geo.nFuselage > 0:
            name = "f_n_clamped"
            # logger.debug(self.geo.aircraftNodesNames)

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
            # logger.debug(self.geo.aircraftConnectedNodes)
            for i in range(N):
                # logger.debug(self.geo.aircraftConnectedNodes[0][i])
                beamIndex1 = int(self.geo.aircraftConnectedNodes[i][0])
                beamIndex1Node = int(self.geo.aircraftConnectedNodes[i][2])
                beamIndex2 = int(self.geo.aircraftConnectedNodes[i][1])
                beamIndex2Node = int(self.geo.aircraftConnectedNodes[i][3])
                # logger.debug(beamIndex1)
                # logger.debug(beamIndex1Node)
                # logger.debug(beamIndex2)
                # logger.debug(beamIndex2Node)
                name1 = self.geo.aircraftNodesNames[beamIndex1][beamIndex1Node]
                name2 = self.geo.aircraftNodesNames[beamIndex2][beamIndex2Node]
                bc.add('connect',{'node1': name1,
                                  'node2': name2,
                                  'fix':['all']})
                # logger.debug("Connects node: " + name1)
                # logger.debug("connects node: " + name2)
                # logger.debug("="*30)
                
            
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
