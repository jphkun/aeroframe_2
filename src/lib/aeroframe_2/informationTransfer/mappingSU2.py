#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:56:01 2020

@author: Jean-Philippe Kuntzer

TODO: For the computation of the moment, distance must be computed in the x,y
      plane.
TODO: verify the moment orientation
TODO: set to 0 small values to speed up the code and limit error.
"""

import logging
import numpy as np
import numpy.linalg as LA
import scipy as sp
import skspatial.objects as sk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import sys

logger = logging.getLogger(__name__)


class mapper:
    def __init__(self,SU2Data,preMeshedStructre,csdSolverClassVar):
        """
        Initialises the class and separes the wing points for enhanced quality
        results.
        """
        # For debug purposes
        self.plotting = True
        np.set_printoptions(precision=3)

        # Assembles matrices
        self.geo = preMeshedStructre
        self.SU2 = SU2Data
        logger.debug(self.SU2)
        self.geoP = preMeshedStructre.aircraftNodesPoints
        self.csd = csdSolverClassVar

        # Separates input points into each aircraft instance (fuselage, wings)
        self.aircraftPointsAndForcesCFD = []
        self.aircraftPartsNames = self.SU2["marker"].unique()
        logger.debug(self.aircraftPartsNames)
        logger.debug(self.geo.aircraftPartsUIDs)
        self.order = []
        # Constructs a dictionary for ease of use when regrouping displacements
        for i in self.aircraftPartsNames:
            index = self.geo.aircraftPartsUIDs.index(i)
            logger.debug(index)
            self.order.append(index)
            # logger.debug(self.aircraftPartsNames.index(i))
            # self.order.append(index)
            # logger.debug(self.order) 
        # # Constructs a dictionary for ease of use when regrouping displacements
        # for i in self.geo.aircraftPartsUIDs:
        #     index = np.where(self.aircraftPartsNames == i)
        #     self.order.append(index[0][0])
        #     # logger.debug(self.aircraftPartsNames.index(i))
        #     # self.order.append(index)
        #     # logger.debug(self.order)
        logger.debug(self.order)
        # sys.exit()
        
        N = len(self.geo.aircraftPartsUIDs)
        logger.debug(N)
        for i in range(N):
            logger.debug(i)
            self.aircraftPointsAndForcesCFD.append(self.SU2[self.SU2['marker'].str.contains(self.geo.aircraftPartsUIDs[i])].to_numpy())

        # if self.plotting:
        #     fig = plt.figure("figure 1")
        #     ax = fig.add_subplot(111, projection='3d')
        #     # for i in range(len(self.wingsPoints)):
        #     for i in range(len(self.aircraftPartsNames)):
        #         ax.scatter(self.aircraftPointsAndForcesCFD[i][:,1],
        #                    self.aircraftPointsAndForcesCFD[i][:,2],
        #                    self.aircraftPointsAndForcesCFD[i][:,3],
        #                    label=self.geo.aircraftPartsUIDs[i])
        #     val = 15
        #     ax.set_xlim(-val,val)
        #     ax.set_ylim(-val,val)
        #     ax.set_zlim(-val,val)
        #     ax.legend()
        #     plt.show()

    def computesTransformationsMatrices(self):
        """
        Computes the transformation matrix in order for the mesh tranformation
        to follow the principle of virtual work.
        """

        # Computes transformations matrices
        self.iM= []
        self.A = []
        self.H = []

        # For debugging
        self.dzsGlob = []
        self.dzaGlob = []

        # Computes the matrix M and then invert it for each aircraft part.
        for i in range(len(self.aircraftPartsNames)):
            # permitted choices are: G,L,TPS,HMQ,HIMQ,C0,C2,C4,C6,EH. See the 
            # function below for the definition of the shape functions.
            fun = "L"
            n = self.geoP[i].shape
            n = n[0]

            # Computes the distance matrix for the radial basis function
            # computation
            X = self.geoP[i]
            r = sp.spatial.distance.cdist(X,X,"euclidean")
            # Computes the radial basis function matrix "M" in the report
            Mbeam = self.phi(r,fun)
            self.iM.append(np.linalg.inv(Mbeam))

            # Computes the matrix Q which is also the matrix A transposed in
            # this specific case.
            m = len(self.aircraftPointsAndForcesCFD[i])
            X = self.geoP[i]
            Y = self.aircraftPointsAndForcesCFD[i][:,1:4]
            r = sp.spatial.distance.cdist(X,Y,"euclidean")
            # Computes the radial basis function matrix "Q" in the report
            Q = self.phi(r,fun)
            self.A.append(Q.T)
            self.H.append(np.matmul(self.A[i],self.iM[i]))

            logger.debug("A "+str(self.A[i].shape))
            logger.debug("iM"+str(self.iM[i].shape))
            logger.debug("H "+str(self.H[i].shape))
        
        # Computes the distances from the leading edge for each point. This
        # computation result is necessary for the moments computation.
        self.computeDistanceForMorments()
        logger.info("Finised computing distance matrix")
        #     # tests the mapping:
        #     n = self.geoP[i].shape
        #     n = n[0]
        #     dzs = np.zeros(n)
        #     for k in range(n):
        #         dzs[k] = 0.01 * self.geoP[i][k,1]**2
        #     self.dzsGlob.append(dzs)
        #     dza = np.matmul(self.H[i],self.dzsGlob[i])
        #     self.dzaGlob.append(dza)

        # logger.debug(self.dzsGlob)
        # logger.debug(self.dzaGlob)

        # # sys.exit()
        # # Plots line
        # if self.plotting:
        #     fig = plt.figure("figure 2")
        #     ax = fig.add_subplot(111, projection='3d')
        #     for p in range(len(self.aircraftPartsNames)):
        #         logger.debug(self.aircraftPointsAndForcesCFD[p].shape)
        #         # ax.scatter(self.aircraftPointsAndForcesCFD[p]["x"],
        #         #            self.aircraftPointsAndForcesCFD[p]["y"],
        #         #            self.aircraftPointsAndForcesCFD[p]["z"],
        #         #            label=self.aircraftPartsNames[i])
        #         ax.scatter(self.aircraftPointsAndForcesCFD[p][:,1],
        #                    self.aircraftPointsAndForcesCFD[p][:,2],
        #                    self.aircraftPointsAndForcesCFD[p][:,3] + self.dzaGlob[p],
        #                    label=self.geo.aircraftPartsUIDs[p])
        #     val = 15
        #     ax.set_xlim(-val,val)
        #     ax.set_ylim(-val,val)
        #     ax.set_zlim(-val,val)
        #     ax.legend()
        #     plt.show()

    def phi(self,r,fun):
        """
        Set of radial basis functions that the user can choose of. After some
        test "Wendland C2" seems to be the better choice, but this is really
        up to user preference.
        """
        eps = 1
        # r = np.linalg.norm(x1-x2)
        if fun == "G":
            # Gaussian
            phi_x = np.exp(-eps*r**2)
        elif fun == "L":
            # Linear
            phi_x = r
        elif fun == "TPS":
            # Thin plate spline
            phi_x = r**2 * np.log(r)
        elif fun == "HMQ":
            # Hardy's multiquadratic
            phi_x = (eps**2 + r**2)**0.5
        elif fun == "HIMQ":
            # Hardy's inverse multiquadratic
            phi_x = 1/(eps**2 + r**2)**0.5
        elif fun == "C0":
            # Wendland's C0
            phi_x = (1-r)**2
        elif fun == "C2":
            # Wendland's C2
            phi_x = (1-r)**4 * (4*r + 1)
        elif fun == "C4":
            # Wendland's C4
            phi_x = (1-r)**6 * (35*r**2 + 18*r + 3)
        elif fun == "C6":
            # Wendland's C6
            phi_x = (1-r)**8 * (32*r**3 + 25*r**2 + 8*r + 1)
        elif fun == "EH":
            # Euclid's hat
            phi_x = np.pi*((1/12*r**3) - r*eps**2 + 4/3*eps**3)
        return phi_x

    def aeroToStructure(self):
        """
        Compute the forces for the structure solver from the CFD solver resutls.
        """
        logger.debug("aeroToStructure")       
        # structure forces
        self.sfx = []
        self.sfy = []
        self.sfz = []
        self.smx = []
        self.smy = []
        self.smz = []
        
        # Aerodynamics forces
        self.afx = []
        self.afy = []
        self.afz = []

        # separates froces for each wings
        N = len(self.aircraftPartsNames)
        for i in range(N):
            self.afx.append(self.aircraftPointsAndForcesCFD[i][:,3])
            self.afy.append(self.aircraftPointsAndForcesCFD[i][:,4])
            self.afz.append(self.aircraftPointsAndForcesCFD[i][:,5])
        # logger.debug(self.afx[0].shape)
        # logger.debug(self.afx[1].shape)
        # logger.debug(self.afx[2].shape)
        # logger.debug(self.afx[3].shape)
        # Calls the function in order to compute the moment generated by each
        # force on the wing.
        self.computeMoments()

        # Computes the forces that act on the structure
        for i in range(N):
            # Computes the forces
            self.sfx.append(np.matmul(self.H[i].T,self.afx[i]))
            self.sfy.append(np.matmul(self.H[i].T,self.afy[i]))
            self.sfz.append(np.matmul(self.H[i].T,self.afz[i]))
            # Computes the moments
            self.smx.append(np.matmul(self.H[i].T,self.amx[i]))
            self.smy.append(np.matmul(self.H[i].T,self.amy[i]))
            self.smz.append(np.matmul(self.H[i].T,self.amz[i]))
        # Prints force felt by the structure
        logger.debug("sfx = \n"+str(self.sfx))
        logger.debug("sfy = \n"+str(self.sfy))
        logger.debug("sfz = \n"+str(self.sfz))
        # Prints moment felt by the structure
        logger.debug("smx = \n"+str(self.smx))
        logger.debug("smy = \n"+str(self.smy))
        logger.debug("smz = \n"+str(self.smz))

    def computeDistanceForMorments(self):
        """
        1) Retrieves the forces
        2) Retrieves the points
        3) Compute the distance between this point force location and all the
            lines nodes points location.
        4) Select the three closest points.
        5) Computes the point projection on both segments and on the closest
            point. Retrives: dx,dy,dz
        6) Select if the current point projection is one of the following
            three cases:
                1. On the real segement for both secments
                2. Only on one segment
                3. On no segment and hence exactly the distance is taken from
                  the exact projection
        7) Computes the moment
        """
        N = len(self.aircraftPointsAndForcesCFD)
        self.distanceMatrix = []

        for i in range(N):
            M = len(self.aircraftPointsAndForcesCFD[i])
            X = self.aircraftPointsAndForcesCFD[i][:,1:4]
            Y = self.geo.aircraftNodesPoints[i]
            
            # Computes the distance between each point of X and each point of
            # Y. This leads to an (NxM) matrix, M being the number of structure
            # nodes points.
            dist = sp.spatial.distance.cdist(X,Y,"euclidean")
            self.distanceMatrix.append(np.empty((M,3)))
            # Finds the minimal 3 values
            for j in range(M):
                point = self.aircraftPointsAndForcesCFD[i][j,1:4]
                point = point.astype('float64')
                # point[1] = np.around(point[1],6)
                # logger.debug("point = \n"+str(point))
                indexes = np.argsort(dist[j])[:3]
                # Stores the 3 points of interest
                p1 = self.geo.aircraftNodesPoints[i][indexes[0]]
                p2 = self.geo.aircraftNodesPoints[i][indexes[1]]
                p3 = self.geo.aircraftNodesPoints[i][indexes[2]]

                # Computes the two lines vectors
                v12 = sk.Vector(p2-p1)
                v23 = sk.Vector(p3-p2)

                line1 = sk.Line(point=p2, direction=v12)
                line2 = sk.Line(point=p2, direction=v23)
                # logger.debug(line1)
                # sys.exit()
                proj1 = line1.project_point(point)
                proj2 = line2.project_point(point)

                # Computes the distance between the projected point and the
                # most far way point. This permetis to test if the projection
                # is still on the structural mesh or not
                distP1Proj1 = np.linalg.norm(p1 - proj1)
                distP1P2 = np.linalg.norm(p1 - p2)
                distP3Proj2 = np.linalg.norm(p3 - proj2)
                distP2P3 = np.linalg.norm(p3 - p2)

                # the two selected segments are parallel.
                if proj1[0] == proj2[0] and \
                    proj1[1] == proj2[1] and \
                    proj1[2] == proj2[2]:
                    delta = -(point - proj1)

                # both projected points are on the mesh, we need to take the
                # the one that has the least distance to the line.
                elif distP1Proj1 < distP1P2 and distP3Proj2 < distP2P3:
                    norm1 = np.linalg.norm(point - proj1)
                    norm2 = np.linalg.norm(point - proj2)
                    norms = np.array([norm1,norm2])
                    projs = np.array([proj1,proj2])
                    delta = -(point - projs[np.argmin(norms)])

                # line 2 projected point is outside the mesh but not the line 1
                # projected point.
                elif distP1Proj1 > distP1P2 and distP3Proj2 < distP2P3:
                    delta = -(point - proj1)

                # line 1 projected point is outside the mesh but not the line 2
                # projected point.
                elif distP1Proj1 < distP1P2 and distP3Proj2 > distP2P3:
                    delta = -(point - proj2)

                # line 1 projected point is outside the mesh but not the line 2
                # projected point.
                elif distP1Proj1 > distP1P2 and distP3Proj2 > distP2P3:
                    delta = -(point - p2)
                self.distanceMatrix[i][j] = np.array([delta[0],delta[1],delta[2]])

    def computeMoments(self):
        """
        Computes the moments using the user input distance from the leading
        edge of the wing.
        """
        self.amx = []
        self.amy = []
        self.amz = []
        N = len(self.distanceMatrix)
        for i in range(N):
            # Computes the moment on the beam generated by all the forces.
            self.amx.append(self.distanceMatrix[i][:,2]*self.afy[i] +
                            self.distanceMatrix[i][:,1]*self.afz[i])
            self.amy.append(self.distanceMatrix[i][:,2]*self.afx[i] +
                            self.distanceMatrix[i][:,0]*self.afz[i])
            self.amz.append(self.distanceMatrix[i][:,0]*self.afy[i] +
                            self.distanceMatrix[i][:,1]*self.afx[i])

            # In order to limit the numerical error
            # self.amx[i][np.abs(self.amx[i])<1e-7] = 0.
            # self.amy[i][np.abs(self.amy[i])<1e-7] = 0.
            # self.amz[i][np.abs(self.amz[i])<1e-7] = 0.

    def structureToAero(self):
        """
        Converts the displacements from the structure mesh to the aerodynamic
        mesh.
        """
        # For debugging only
        plotting = False
        
        # structure displacements
        self.sux = []
        self.suy = []
        self.suz = []
        self.stx = []
        self.sty = []
        self.stz = []
        
        # aerodynamic displacements
        self.aux = []
        self.auy = []
        self.auz = []
        self.atx = []
        self.aty = []
        self.atz = []

        # number of beams
        N = len(self.geoP)
        # logger.debug("N = "+str(N))
        old = 0
        
        # Separates the results for each wing. This leads to an amazing quality
        # jump in the simulation results.
        for i in range(N):
            # Number of node for each beams
            # logger.debug()
            M = len(self.geoP[i])
            self.sux.append(self.csd.results.get('tensors').get('comp:U')["ux"][old:old+M])
            self.suy.append(self.csd.results.get('tensors').get('comp:U')["uy"][old:old+M])
            self.suz.append(self.csd.results.get('tensors').get('comp:U')["uz"][old:old+M])
            
            self.stx.append(self.csd.results.get('tensors').get('comp:U')["thx"][old:old+M])
            self.sty.append(self.csd.results.get('tensors').get('comp:U')["thy"][old:old+M])
            self.stz.append(self.csd.results.get('tensors').get('comp:U')["thz"][old:old+M])
            old += M
        np.set_printoptions(precision=3)
        logger.debug("\n"*20)
        logger.debug("==== sux ========================")
        logger.debug(self.sux)
        logger.debug("==== suy ========================")
        logger.debug(self.suy)
        logger.debug("==== suz ========================")
        logger.debug(self.suz)
        logger.debug("==== stx ========================")
        logger.debug(self.stx)
        logger.debug("==== sty ========================")
        logger.debug(self.sty)
        logger.debug("==== stz ========================")
        logger.debug(self.stz)
        
        logger.debug("\n"*20)
        # Computes the aerodynamic displacements and the aerodynamic angles
        for i in range(N):
            self.aux.append(np.matmul(self.H[i],self.sux[i]))
            self.auy.append(np.matmul(self.H[i],self.suy[i]))
            self.auz.append(np.matmul(self.H[i],self.suz[i]))
            
            self.atx.append(np.matmul(self.H[i],self.stx[i]))
            self.aty.append(np.matmul(self.H[i],self.sty[i]))
            self.atz.append(np.matmul(self.H[i],self.stz[i]))

        # logger.debug(self.atx[0])
        # logger.debug(self.atx[1])
        # logger.debug(self.atx[2])
        # logger.debug(self.atx[3])
        
        # Assembles the displacements into a big vector
        size = len(self.SU2)
        self.displacements = np.empty((size,3))
        logger.debug(self.displacements.shape)
        N = len(self.geoP) # in the optimale case N=4
        
        # Assembles the displacements
        ind = self.order[0] # CFD mesh index
        logger.debug(ind)
        dmx = self.aty[ind] * self.distanceMatrix[ind][:,2] + \
              self.atz[ind] * self.distanceMatrix[ind][:,1]
        dmy = self.atx[ind] * self.distanceMatrix[ind][:,2] + \
              self.atz[ind] * self.distanceMatrix[ind][:,0]
        dmz = self.atx[ind] * self.distanceMatrix[ind][:,1] + \
              self.aty[ind] * self.distanceMatrix[ind][:,0]
        dfx = self.aux[ind]
        dfy = self.auy[ind]
        dfz = self.auz[ind]
        for i in range(1,N):
            ind = self.order[i] # CFD mesh index
            logger.debug(ind)
            logger.debug(type(ind))
            dmxi = self.aty[ind] * self.distanceMatrix[ind][:,2] + \
                   self.atz[ind] * self.distanceMatrix[ind][:,1]
            dmx = np.concatenate((dmx,dmxi))
            dfx = np.concatenate((dfx,self.aux[ind]))

            dmyi = self.atx[ind] * self.distanceMatrix[ind][:,2] + \
                   self.atz[ind] * self.distanceMatrix[ind][:,0]
            dmy = np.concatenate((dmy,dmyi))
            dfy = np.concatenate((dfy,self.auy[ind]))

            dmzi = self.atx[ind] * self.distanceMatrix[ind][:,1] + \
                   self.aty[ind] * self.distanceMatrix[ind][:,0]
            dmz = np.concatenate((dmz,dmzi))
            dfz = np.concatenate((dfz,self.auz[ind]))
        self.displacements = np.array([dmx+dfx,dmy+dfy,dmz+dfz]).T
        logger.debug(self.displacements[0])
        logger.debug(self.displacements[1])
        logger.debug(self.displacements[1])
        logger.debug(np.max(self.displacements[0]))
        logger.debug(np.max(self.displacements[1]))
        logger.debug(np.max(self.displacements[2]))

        # Plots
        if self.plotting:
            fig = plt.figure("figure 2")
            ax = fig.add_subplot(111, projection='3d')
            old = 0
            # for j in range(len(self.aircraftPartsNames)):
            #     ax.scatter(self.aircraftPointsAndForcesCFD[j][:,1],
            #                self.aircraftPointsAndForcesCFD[j][:,2],
            #                self.aircraftPointsAndForcesCFD[j][:,3],
            #                label=self.geo.aircraftPartsUIDs[j])
            ax.scatter(self.SU2["x"],
                        self.SU2["y"],
                        self.SU2["z"],
                        label = "undeformed")

            ax.scatter(self.SU2["x"] + self.displacements[:,0],
                        self.SU2["y"] + self.displacements[:,1],
                        self.SU2["z"] + self.displacements[:,2],
                        label = "deformed")
            val = 15
            ax.set_xlim(-val,val)
            ax.set_ylim(-val,val)
            ax.set_zlim(-val,val)
            ax.legend()
            plt.show()

    # def deformVLM(self):
    #     pass
            