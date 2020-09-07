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
import pandas as pd
import time
import sys
import os
import fnmatch

logger = logging.getLogger(__name__)


class mapper:
    def __init__(self,forceFilePath,preMeshedStructre,csdSolverClassVar):
        """
        Initialises the class and separes the wing points for enhanced quality
        results.
        """
        # For debug purposes
        np.set_printoptions(precision=3)

        # Assembles matrices
        self.geo = preMeshedStructre
        self.geoP = preMeshedStructre.aircraftNodesPoints
        self.csd = csdSolverClassVar

        # Separates input points into each aircraft instance (fuselage, wings)
        SU2Data = pd.read_csv(forceFilePath)
        N = len(self.geo.aircraftPartsUIDs)
        logger.debug(N)
        aircraftPointsAndForcesCFD = []
        for i in range(N):
            logger.debug(i)
            aircraftPointsAndForcesCFD.append(SU2Data[SU2Data['marker'].str.contains(self.geo.aircraftPartsUIDs[i])].to_numpy())

        self.aircraftPartsNames = SU2Data["marker"].unique()
        del(SU2Data)
        logger.debug(self.aircraftPartsNames)
        logger.debug(self.geo.aircraftPartsUIDs)
        self.order = []

        # Constructs a dictionary for ease of use when regrouping displacements
        for i in self.aircraftPartsNames:
            index = self.geo.aircraftPartsUIDs.index(i)
            logger.debug(index)
            self.order.append(index)

    def computesTransformationsMatrices(self,forceFilePath):
        """
        Computes the transformation matrix in order for the mesh tranformation
        to follow the principle of virtual work.
        """

        # Computes transformations matrices
        self.iM = []
        self.A = []
        self.H = []

        # For debugging
        self.dzsGlob = []
        self.dzaGlob = []

        # Updates aircraftPointsAndForcesCFD
        SU2Data = pd.read_csv(forceFilePath)
        N = len(self.geo.aircraftPartsUIDs)
        logger.debug(N)
        aircraftPointsAndForcesCFD = []
        for i in range(N):
            logger.debug(i)
            aircraftPointsAndForcesCFD.append(SU2Data[SU2Data['marker'].str.contains(self.geo.aircraftPartsUIDs[i])].to_numpy())

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
            X = self.geoP[i]
            Y = aircraftPointsAndForcesCFD[i][:,1:4]
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
        # self.computeMoments(forceFilePath)
        # logger.info("Finised computing distance matrix")
        
        # Finds wings nodes close tu fuselage in order to correct the mesh
        # rotation
        # self.findCloseToWing(forceFilePath)
        
        # Frees memory. not sure it is still needed
        del(SU2Data)

    def phi(self,r,fun):
        """
        Set of radial basis functions that the user can choose of. After some
        test "Wendland C2" seems to be the better choice, but this is really
        up to user preference.
        """
        eps = 10
        # r = np.linalg.norm(x1-x2)
        if fun == "G":
            # Gaussian
            phi_x = np.exp(-eps*r**2)
        elif fun == "L":
            # Linear
            phi_x = r
        elif fun == 'cubic':
            # Cubic
            phi_x = r**3
        elif fun == "TPS":
            # Thin plate spline
            phi_x = r**2 * np.log(r)
        elif fun == "Po4":
            # Polyharmonic spline
            phi_x = r**3 * np.log(r**r)
        elif fun == "Po4":
            # Polyharmonic spline
            phi_x = r**5 * np.log(r**r)
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
        
    def aeroToStructure(self,args,forceFilePath,iteration):
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

        # Updates aircraftPointsAndForcesCFD
        SU2Data = pd.read_csv(forceFilePath)
        N = len(self.geo.aircraftPartsUIDs)
        logger.debug(N)
        aircraftPointsAndForcesCFD = []
        for i in range(N):
            logger.debug(i)
            aircraftPointsAndForcesCFD.append(SU2Data[SU2Data['marker'].str.contains(self.geo.aircraftPartsUIDs[i])].to_numpy())

        # separates froces for each wings
        N = len(self.aircraftPartsNames)
        for i in range(N):
            self.afx.append(aircraftPointsAndForcesCFD[i][:,4])
            self.afy.append(aircraftPointsAndForcesCFD[i][:,5])
            self.afz.append(aircraftPointsAndForcesCFD[i][:,6])
        # logger.debug(self.afx[0].shape)
        # logger.debug(self.afx[1].shape)
        # logger.debug(self.afx[2].shape)
        # logger.debug(self.afx[3].shape)
        # Calls the function in order to compute the moment generated by each
        # force on the wing.
        self.computeMoments(forceFilePath)

        # Computes the forces that act on the structure
        Fx = 0
        Fy = 0
        Fz = 0
        for i in range(N):
            # Computes the forces
            self.sfx.append(np.matmul(self.H[i].T,self.afx[i]))
            self.sfy.append(np.matmul(self.H[i].T,self.afy[i]))
            self.sfz.append(np.matmul(self.H[i].T,self.afz[i]))
            # Computes the moment part due to the aerodynamic force
            self.smx.append(np.matmul(self.H[i].T,self.amx[i]))
            self.smy.append(np.matmul(self.H[i].T,self.amy[i]))
            self.smz.append(np.matmul(self.H[i].T,self.amz[i]))
            
            # sums all the aerodynamic forces
            Fx += np.sum(self.afx[i])
            Fy += np.sum(self.afx[i])
            Fz += np.sum(self.afx[i])
            
            # Swept wing have a tendency to increase the central lift
            M = int(np.floor(len(self.sfx[i])/2))
            # Damps the inital and final jump
            self.sfx[i][0] = self.sfx[i][1]#*0
            self.sfx[i][-1] = self.sfx[i][-2]#*0
            self.sfx[i][M] = self.sfx[i][M-1]
            self.sfy[i][0] = self.sfy[i][1]#*0
            self.sfy[i][-1] = self.sfy[i][-2]#*0
            self.sfy[i][M] = self.sfy[i][M-1]
            self.sfz[i][0] = self.sfz[i][1]#*0
            self.sfz[i][-1] = self.sfz[i][-2]#*0
            self.sfz[i][M] = self.sfz[i][M-1]
            
            # Damps the inital and final jump
            self.smx[i][0] = self.smx[i][1]#*0
            self.smx[i][-1] = self.smx[i][-2]#*0
            self.smx[i][M] = self.smx[i][M-1]
            self.smy[i][0] = self.smy[i][1]#*0
            self.smy[i][-1] = self.smy[i][-2]#*0
            self.smy[i][M] = self.smy[i][M-1]
            self.smz[i][0] = self.smz[i][1]#*0
            self.smz[i][-1] = self.smz[i][-2]#*0
            self.smz[i][M] = self.smz[i][M-1]
        
        logger.debug(self.amy)
        logger.debug(self.smy)

        # Saves data for verificaiton
        df = pd.DataFrame()
        # Structure mesh node position
        df['x'] = self.geo.aircraftNodesPoints[0][:,0]
        df['y'] = self.geo.aircraftNodesPoints[0][:,1]
        df['z'] = self.geo.aircraftNodesPoints[0][:,2]
        # Forces
        df['Fx'] = pd.Series(self.sfx[0])
        df['Fy'] = pd.Series(self.sfy[0])
        df['Fz'] = pd.Series(self.sfz[0])
        # Moments/Torques
        df['Mx'] = pd.Series(self.smx[0])
        df['My'] = pd.Series(self.smy[0])
        df['Mz'] = pd.Series(self.smz[0])
        df.to_csv(args.cwd + '/CSD/results/FEM_frocesAndMoments'+str(iteration)+'.csv')
        
        if self.geo.settings['1G']:
            n = 1.0
            a_x = 0
            a_y = 0
            a_z = -9.81
        else:
            n = Fz/(self.geo.aircraftTotalMass * 9.81)
            a_x = Fx / self.geo.aircraftTotalMass
            a_y = Fy / self.geo.aircraftTotalMass
            a_z = n - 1
        # n = 0
        self.G = round(n,2)
        # logger.debug('a_x = ' + str(a_x))
        # logger.debug('a_y = ' + str(a_y))
        # logger.debug('a_z = ' + str(a_z))
        logger.debug('If G activated and 1G not activated G = '+str(self.G))
        logger.debug('Only used if G_load:true and 1G:false')
        # logger.debug('mass = ' + str(self.geo.aircraftTotalMass))
        # sys.exit()
        # Computes the force due to inertia on each strcutre node
        self.smf = []
        N = len(self.geo.aircraftSegementsMass)
        for i in range(N):
            # Since there is one more point then segment we need to add one at
            # the end.
            M = len(self.geo.aircraftSegementsMass[i]) + 1
            force = np.empty((M,3))
            for j in range(M):
                # loadFactor
                # Each node supports half of the weight to the left and half
                # to the right. This explains why the first and last nodes
                # support half the value of the mass.
                if j == 0:
                    massRight = self.geo.aircraftSegementsMass[i][j]
                    # Inertial force in the x direction, earth reference frame
                    force[j,0] = 1 * massRight * a_x * 0
                    # Inertial force in the x direction, earth reference frame
                    force[j,1] = 1 * massRight * a_y * 0
                    # Inertial force in the x direction, earth reference frame
                    force[j,2] = 0.5 * massRight * n * 9.81
                elif j == M-1:
                    massLeft = self.geo.aircraftSegementsMass[i][j-1]
                    # Inertial force in the x direction, earth reference frame
                    force[j,0] = 1 * massLeft * a_x * 0
                    # Inertial force in the x direction, earth reference frame
                    force[j,1] = 1 * massLeft * a_y * 0
                    # Inertial force in the x direction, earth reference frame
                    force[j,2] = 0.5 * massLeft * n * 9.81
                else:
                    massRight = self.geo.aircraftSegementsMass[i][j]
                    massLeft = self.geo.aircraftSegementsMass[i][j-1]
                    # Inertial force in the x direction, earth reference frame
                    force[j,0] = 0.5 * massRight * a_x * 0 + \
                                 0.5 * massLeft  * a_x * 0 
                    # Inertial force in the x direction, earth reference frame
                    force[j,1] = 0.5 * massRight * a_y * 0 + \
                                 0.5 * massLeft  * a_y * 0
                    # Inertial force in the x direction, earth reference frame
                    force[j,2] = 0.5 * massRight * n * 9.81 + \
                                 0.5 * massLeft  * n * 9.81
                        
            self.smf.append(force)
        # logger.debug(self.smf)
        # sys.exit()
        

        # Computes the moment due to inertia on each strcture node
        # for i i
        N = len(self.geo.aircraftMassDistances)
        self.smm = []
        for i in range(N):
            # Since there is one more point then segments we need to add one at
            # the end.
            M = len(self.geo.aircraftMassDistances[i])
            # logger.debug(M)
            # logger.debug(len(self.smf[i]))
            # sys.exit()
            moments = np.empty((M,3))
            for j in range(M):
                # WARNING only the vertical direction is implemented.
                dx = self.geo.aircraftMassDistances[i][j,0]
                dy = self.geo.aircraftMassDistances[i][j,1]
                dz = self.geo.aircraftMassDistances[i][j,2]
                moments[j,0] = 0
                moments[j,1] = np.sign(dx)*np.sqrt(dx**2 + dy**2) * self.smf[i][j,2]
                moments[j,2] = 0
            self.smm.append(moments)

        # logger.debug(moments)
        # sys.exit()
        # Computes the total of each force and moment in order to have an idea
        # of the information loss between two steps.
        self.totalAerodynamicFx = np.sum(self.afx)
        self.totalAerodynamicFy = np.sum(self.afy)
        self.totalAerodynamicFz = np.sum(self.afz)
        self.totalAerodynamicMx = np.sum(self.amx)
        self.totalAerodynamicMy = np.sum(self.amy)
        self.totalAerodynamicMz = np.sum(self.amz)
        self.totalStructureFx = np.sum(self.sfx)
        self.totalStructureFy = np.sum(self.sfy)
        self.totalStructureFz = np.sum(self.sfz)
        self.totalStructureMx = np.sum(self.smx)
        self.totalStructureMy = np.sum(self.smy)
        self.totalStructureMz = np.sum(self.smz)
        # logger.debug('Conservation of forces and moments')
        # logger.debug('\n'*5)
        # logger.debug(np.sum(self.afx))
        # logger.debug(np.sum(self.sfx))
        # logger.debug(np.sum(self.afy))
        # logger.debug(np.sum(self.sfy))
        # logger.debug(np.sum(self.afz))
        # logger.debug(np.sum(self.sfz))
        # logger.debug('\n'*5)
        # logger.debug(np.sum(self.amx))
        # logger.debug(np.sum(self.smx))
        # logger.debug(np.sum(self.amy))
        # logger.debug(np.sum(self.smy))
        # logger.debug(np.sum(self.amz))
        # logger.debug(np.sum(self.smz))

    def computeMoments(self,forceFilePath):
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
        # Updates aircraftPointsAndForcesCFD
        SU2Data = pd.read_csv(forceFilePath)
        N = len(self.geo.aircraftPartsUIDs)
        logger.debug(N)
        aircraftPointsAndForcesCFD = []
        for i in range(N):
            logger.debug(i)
            aircraftPointsAndForcesCFD.append(SU2Data[SU2Data['marker'].str.contains(self.geo.aircraftPartsUIDs[i])].to_numpy())

        N = len(aircraftPointsAndForcesCFD)
        self.distanceMatrix = []
        self.amx = []
        self.amy = []
        self.amz = []
        for i in range(N):
            M = len(aircraftPointsAndForcesCFD[i])
            X = aircraftPointsAndForcesCFD[i][:,1:4]
            Y = self.geo.aircraftNodesPoints[i]

            # Computes the distance between each point of X and each point of
            # Y. This leads to an (NxM) matrix, M being the number of structure
            # nodes points.
            dist = sp.spatial.distance.cdist(X,Y,"euclidean")
            # distances in the x,y,z coordinates in the airplane reference
            # frame.
            self.distanceMatrix.append(np.empty((M,3)))
            # Finds the minimal 3 values
            # logger.debug(dist)
            
            for j in range(M):
                point = aircraftPointsAndForcesCFD[i][j,1:4]
                point = point.astype('float64')
                # point[1] = np.around(point[1],6)
                # logger.debug("point = \n"+str(point))
                indexes = np.argsort(dist[j])[:3]
                # Stores the 3 points of interest
                p1 = self.geo.aircraftNodesPoints[i][indexes[0]]
                p2 = self.geo.aircraftNodesPoints[i][indexes[1]]
                p3 = self.geo.aircraftNodesPoints[i][indexes[2]]

                # Computes the two lines direction vectors
                # p1 will always be the closest structure point hence he will
                # be in the middle. Little drawing of what it looks like in
                # space below
                #
                # (P2 or P3)           (P1)
                #                             
                #                                   (P2 or P3)
                v12 = sk.Vector(p2-p1)
                v23 = sk.Vector(p3-p1)
                # logger.debug('Vectors')
                # logger.debug(v12)
                # logger.debug(v23)
                
                line1 = sk.Line(point=p2, direction=v12)
                line2 = sk.Line(point=p2, direction=v23)
                # logger.debug('lines')
                # logger.debug(line1)
                # logger.debug(line2)
                
                proj1 = line1.project_point(point)
                proj2 = line2.project_point(point)
                # logger.debug('Projected points+')
                # logger.debug(proj1)
                # logger.debug(proj2)
                
                # Computes the distance between the projected point and the
                # most far away point. This permetis to test if the projection
                # is still on the structural mesh or not
                distP1Proj1 = LA.norm(p2 - proj1)
                distP1P2 = LA.norm(p1 - p2)
                distP3Proj2 = LA.norm(p3 - proj2)
                distP2P3 = LA.norm(p1 - p3)
                # logger.debug('Projected points distance to center and to relative point')
                # logger.debug(distP1Proj1)
                # logger.debug(distP1P2)
                # logger.debug(distP3Proj2)
                # logger.debug(distP2P3)

                # the two selected segments are parallel.
                if proj1[0] == proj2[0] and \
                   proj1[1] == proj2[1] and \
                   proj1[2] == proj2[2]:
                    delta = -(point - proj1)

                # both projected points are on the mesh, we need to take the
                # the one that has the least distance to the line.
                elif distP1Proj1 < distP1P2 and distP3Proj2 < distP2P3:
                    norm1 = LA.norm(point - proj1)
                    norm2 = LA.norm(point - proj2)
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

                # line 1 projected point is outside the mesh and the line 2
                # projected point is also outside the structure mesh.
                elif distP1Proj1 > distP1P2 and distP3Proj2 > distP2P3:
                    delta = -(point - p1)
                self.distanceMatrix[i][j] = np.array([delta[0],delta[1],delta[2]])
            # logger.debug(self.distanceMatrix)
            # logger.debug(self.wingsPoints[i])
            
            # Computes the moment on the beam generated by all the forces.
            # logger.debug(self.afx[i])
            # logger.debug(self.afy[i])
            # logger.debug(self.afz[i])
            
            # In the airplane reference frame. in straight level flight
            # wind is positive in the x direction and the wings span
            # is in the y direction.
            self.amx.append(self.distanceMatrix[i][:,2]*self.afy[i] + # OK
                            self.distanceMatrix[i][:,1]*self.afz[i])  # OK
            self.amy.append(self.distanceMatrix[i][:,0]*self.afz[i]) # + # OK
                            # self.distanceMatrix[i][:,2]*self.afx[i])  # OK
            self.amz.append(self.distanceMatrix[i][:,0]*self.afy[i] + # OK
                            self.distanceMatrix[i][:,1]*self.afx[i])  # OK
        # for x, y, z in zip(self.distanceMatrix[i][:,0], self.afz[i], self.amy[i]):
        #     logger.debug(str(x) + ' ' + str(y) + ' ' + str(z))
        # logger.debug(np.sum(self.amy[i]))
        # sys.exit()

    def structureToAero(self,args,iteration,forceInitFilePath,forceFilePath):
        """
        Converts the displacements from the structure mesh to the aerodynamic
        mesh.
        """

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
            # if i == 1:
            #     corrDelta = self.csd.results.get('tensors').get('comp:U')["uz"][old+int(np.ceil(M/2))]
            #     corrTheta = self.csd.results.get('tensors').get('comp:U')["thy"][old+int(np.ceil(M/2))]
            # else:
            #     corrDelta = 0
            #     corrTheta = 0
            self.sux.append(self.csd.results.get('tensors').get('comp:U')["ux"][old:old+M])
            self.suy.append(self.csd.results.get('tensors').get('comp:U')["uy"][old:old+M])
            self.suz.append(self.csd.results.get('tensors').get('comp:U')["uz"][old:old+M]) 

            self.stx.append(self.csd.results.get('tensors').get('comp:U')["thx"][old:old+M])
            self.sty.append(self.csd.results.get('tensors').get('comp:U')["thy"][old:old+M])
            self.stz.append(self.csd.results.get('tensors').get('comp:U')["thz"][old:old+M])
            old += M
        
        # saves the results
        df = pd.DataFrame()
        df['x'] = self.geo.aircraftNodesPoints[0][:,0]
        df['y'] = self.geo.aircraftNodesPoints[0][:,1]
        df['z'] = self.geo.aircraftNodesPoints[0][:,2]
        # WARNING change of reference frame
        df['dx'] = self.suy[0]
        df['dy'] = self.sux[0]
        df['dz'] = self.suz[0]
        df['tx'] = self.sty[0]
        df['ty'] = self.stx[0]
        df['tz'] = self.stz[0]
        df.to_csv(args.cwd + '/CSD/results/FEM_displacementAndRotations'+str(iteration)+'.csv')
        
        # Computes the aerodynamic displacements and the aerodynamic angles
        for i in range(N):
            if i == 1:
                coef = 0.05
            else:
                coef = 0
            self.aux.append(np.matmul(self.H[i],self.sux[i]) - coef)
            self.auy.append(np.matmul(self.H[i],self.suy[i]) - coef)
            self.auz.append(np.matmul(self.H[i],self.suz[i]) - coef)

            self.atx.append(np.matmul(self.H[i],self.stx[i]) - coef)
            self.aty.append(np.matmul(self.H[i],self.sty[i]) - coef)
            self.atz.append(np.matmul(self.H[i],self.stz[i]) - coef)

        # Assembles the displacements into a big vector
        SU2Data = pd.read_csv(forceFilePath)
        size = len(SU2Data)
        self.displacements = np.empty((size,3))
        logger.debug(self.displacements.shape)
        N = len(self.geoP)  # in the optimale case N=4

        # Assembles the displacements
        ind = self.order[0]  # CFD mesh index
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
            ind = self.order[i]  # CFD mesh index
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
        # self.displacements = np.array([dfx,dfy,dfz]).T

        # Generates the deformation file
        # 'GlobalID_' +
        # New x,y,z positions
        N = len(self.displacements)
        # logger.debug(np.arange(N))
        # idx =  str()
        # idx = ['GlobalID_' + str(x) for x in np.arange(N)]
        idx = [str(x) for x in np.arange(N)]
        start = 0  # 9136
        idx = np.array(idx[start:])
        logger.debug(idx)
        const = 1
        # some_value = 'Wing'
        # var = SU2Data.loc[self.SU2_forces['marker'] == some_value]
        SU2DataInit = pd.read_csv(forceInitFilePath)
        x = SU2DataInit["x"] + const*self.displacements[start:,0]
        y = SU2DataInit["y"] + const*self.displacements[start:,1]
        z = SU2DataInit["z"] + const*self.displacements[start:,2]
        
        # X = np.array([x,y,z])
        # logger.debug(X)
        # s = X.shape
        # X = X.T.reshape(s[1],s[0])
        # X = X[self.wingCTFIndexes]
        # u,s,vh = np.linalg.svd(X)
        # logger.debug(vh)
        # sys.exit()
        
        path = os.path.join(self.geo.inputArgs.cwd,'CFD')
        caseName = os.listdir(path)
        caseName = fnmatch.filter(caseName, 'Case*')
        logger.debug(caseName)
        fname = path + '/' + caseName[0] + '/' + str(iteration) + '/disp.dat'
        logger.debug('filename: \n'+str(fname))
        # X = np.array()
        np.savetxt(fname,np.column_stack([idx,x,y,z]),delimiter='\t',fmt='%s')

        # Plots
        self.plotting = False
        if self.plotting:
            fig = plt.figure('figure 2')
            ax = fig.add_subplot(111, projection='3d')
            old = 0
            # for j in range(len(self.aircraftPartsNames)):
            #     ax.scatter(self.aircraftPointsAndForcesCFD[j][:,1],
            #                self.aircraftPointsAndForcesCFD[j][:,2],
            #                self.aircraftPointsAndForcesCFD[j][:,3],
            #                label=self.geo.aircraftPartsUIDs[j])
            ax.scatter(SU2Data[start:,'x'],
                       SU2Data[start:,'y'],
                       SU2Data[start:,'z'],
                       label='undeformed')

            ax.scatter(SU2Data[start:,'x'] + self.displacements[start:,0],
                       SU2Data[start:,'y'] + self.displacements[start:,1],
                       SU2Data[start:,'z'] + self.displacements[start:,2],
                       label='deformed')
            val = 15
            ax.set_xlim(-val,val)
            ax.set_ylim(-val,val)
            ax.set_zlim(-val,val)
            ax.legend()
            plt.show()
        del(SU2Data)

    def findCloseToWing(self,forceFilePath):
        """
        

        Parameters
        ----------
        forceFilePath : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # Separates the uids
        SU2Data = pd.read_csv(forceFilePath)
        N = len(self.geo.aircraftPartsUIDs)
        aircraftPointsAndForcesCFD = []

        for i in range(N):
            logger.debug(i)
            aircraftPointsAndForcesCFD.append(SU2Data[SU2Data['marker'].str.contains(self.geo.aircraftPartsUIDs[i])].to_numpy())
        logger.debug(aircraftPointsAndForcesCFD[0])
        
        # Stores fuselage points
        fuselagePoints = aircraftPointsAndForcesCFD[0][:,1:4]
        
        # Stores wing points
        wingPoints = aircraftPointsAndForcesCFD[1][:,1:4]
        
        # Compute the distance between each fuselage point to each wing point.
        # This is and NxM matrix with N number of fuselage points and M nbr.
        # of wing points.
        distance = sp.spatial.distance.cdist(fuselagePoints, wingPoints)
        # logger.debug(distance.shape)
        # logger.debug(distance)
        
        # Finds the the minimum distance between each fuselage point and the
        # the closest wing point.
        distanceFusePTWing = np.min(distance,axis=0)
        # Finds the the minimum distance between each wing point and the
        # the closest wing point.
        distanceWingPTFuse = np.min(distance,axis=1)
        # logger.debug(distanceFusePTWing.shape)
        # logger.debug(distanceWingPTFuse.shape)
        
        
        # Arbitrary value that selects the fuselage points that are new the
        # wing.
        distanceMax = 0.1
        
        # Wing point indexes selected with this criteras
        indexesWingTF = np.where(np.logical_and(distanceFusePTWing>=0, distanceFusePTWing<distanceMax))
        indexesWing = aircraftPointsAndForcesCFD[1][indexesWingTF,0]
        indexesWing = indexesWing[0]

        # Stores the wings point that satify the condition
        selectedWingPoints = aircraftPointsAndForcesCFD[1][indexesWingTF,1:4]
        s = selectedWingPoints.shape
        selectedWingPoints = selectedWingPoints.reshape(s[0]*s[1],s[2])
        selectedWingPoints = selectedWingPoints.astype(np.float32)
        u, s, vh = np.linalg.svd(selectedWingPoints)
        self.V = vh
        self.wingCTFIndexes = indexesWingTF
        logger.debug(vh)
        

        # sys.exit()
        
        # # Computes weighting coefficients
        # x = 1-localWingPointDist/distanceMax
        # # logger.debug(x)
        
        # coef = self.polynome(x)
        # # logger.debug(coef)

        # logger.debug('OG fuselage \n'+str(self.displacements[indexesFuselageTF]))
        # # logger.debug(type(globalWingPointIndex))
        # # logger.debug(type(self.displacements))
        # N = len(globalWingPointIndex)
        # temp = np.empty((N,3))
        # for i in range(N):
        #     # logger.debug()
        #     temp[i] = self.displacements[globalWingPointIndex[i]]
        # logger.debug('OG wing \n'+str(temp))
        # # logger.debug(self.displacements[globalWingPointIndex])
        # # # logger.debug((1 - coef))
        # # # logger.debug(coef)
        # # # logger.debug(( np.einsum('i,ij->ij', (1 - coef).T, self.displacements[indexesTF]) ))
        # # # logger.debug(np.einsum('i,ij->ij', coef.T, self.displacements[relatedFWingPointsIndex]))
        # fuselage = np.einsum('i,ij->ij', coef.T, self.displacements[indexesFuselageTF])
        # logger.debug('fuselage \n'+str(fuselage))
        
        # wing = np.einsum('i,ij->ij', (1 - coef).T, temp)
        # logger.debug('wing \n'+str(wing))
        # self.displacements[indexesFuselageTF] = np.zeros((N,3))  # fuselage + wing            
        # logger.debug(self.displacements[indexesFuselageTF])
        # # # logger.debug(self.displacements[indexesTF])
        # # # logger.debug(relatedFWingPointsIndex)
        # # sys.exit()
