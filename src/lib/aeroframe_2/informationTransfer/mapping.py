#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:56:01 2020

@author: Jean-Philippe Kuntzer

Mapping is a class that applies the principle of virtual work to aeroelaticity.
In this class the matrices linking the CFD mesh to the structure mesh are
computed.

OK TODO: For the computation of the moment, distance must be computed in the x,y
         plane.
OK TODO: verify the moment orientation
TODO: set to 0 small values to speed up the code and limit error.
"""

import logging
import numpy as np
import numpy.linalg as LA
import scipy as sp
import skspatial.objects as sk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

logger = logging.getLogger(__name__)


class mapper:
    def __init__(self,pytornadoVariables,preMeshedStructre,csdSolverClassVar):
        """
        Initialises the class and separes the wing points for enhanced quality
        results.
        """
        # For debug purposes
        plotting = False
        np.set_printoptions(precision=3)
        # Assembles matrices
        self.geo = preMeshedStructre
        self.lattice = pytornadoVariables[0]
        self.VLMdata = pytornadoVariables[1]
        self.geoP = preMeshedStructre.aircraftNodesPoints
        self.csd = csdSolverClassVar
        # Separates lattice.c into each wing instances
        self.wingsPoints = []
        self.limitsStart = []
        self.limitsEnd = []

        # Separates the lattice points understand the VLM mesh for each wing
        # hence leading to better results since each wing only contribute her
        # specific beam.
        number = len(self.lattice.bookkeeping_by_wing_uid)
        for i in self.lattice.bookkeeping_by_wing_uid:

            # Gets the data that needs to be sorted.
            listing = list(self.lattice.bookkeeping_by_wing_uid.get(i)[0][1])
            init = listing[0]
            logger.debug("init = "+str(init))

            listing = list(self.lattice.bookkeeping_by_wing_uid.get(i)[-1][1])
            panels = self.lattice.bookkeeping_by_wing_uid.get(i)[-1][2]

            # takes care of last segment hence the + 1 at the end
            end = listing[-1] + panels + 1
            logger.debug("number = "+str(number))
            logger.debug("end = "+str(end))
            self.limitsStart.append(init)
            self.limitsEnd.append(end)

            # Appends the separated points
            self.wingsPoints.append(self.lattice.c[init:end])
            number -= 1

        # Plot for debug purposes only
        if plotting:
            fig = plt.figure("figure 1")
            ax = fig.add_subplot(111, projection='3d')
            for i in range(len(self.wingsPoints)):
                ax.scatter(self.wingsPoints[i][:,0],
                           self.wingsPoints[i][:,1],
                           self.wingsPoints[i][:,2],
                           label='Wing '+str(i+1))
            val = 15
            ax.set_xlim(-val,val)
            ax.set_ylim(-val,val)
            ax.set_zlim(-val,val)
            ax.legend()
            plt.show()

    def computesTransformationsMatrices(self):
        # Computes transformations matrices
        self.iM = []
        self.A = []
        self.H = []
        # For debugging
        self.dzsGlob = []
        self.dzaGlob = []
        plotting = False

        for i in range(len(self.wingsPoints)):
            # Computes the matrix M and then invert it
            # permitted choices are: G,L,TPS,HMQ,HIMQ,C0,C2,C4,C6,EH see below
            # the definition
            fun = "L"
            n = self.geoP[i + self.geo.nFuselage].shape
            n = n[0]
            Mbeam = np.zeros((n,n))
            for k in range(n):
                for j in range(n):
                    x1 = self.geoP[i + self.geo.nFuselage][k]
                    x2 = self.geoP[i + self.geo.nFuselage][j]
                    Mbeam[k,j] = self.phi(x1,x2,fun)
            self.iM.append(LA.inv(Mbeam))

            # Computes the matrix Q which is also the matrix A transposed in
            # this specific case c.f. to the theory.
            m = self.wingsPoints[i].shape
            m = m[0]
            Q = np.zeros((n,m))
            logger.debug("n = "+str(n))
            logger.debug("m = "+str(m))
            for k in range(n):
                for j in range(m):
                    x1 = self.geoP[i + self.geo.nFuselage][k]
                    x2 = self.wingsPoints[i][j]
                    Q[k,j] = self.phi(x1,x2,fun)
            self.A.append(Q.T)
            self.H.append(np.matmul(self.A[i],self.iM[i]))
            logger.debug(self.lattice.c.shape)
            logger.debug("A "+str(self.A[0].shape))
            logger.debug("iM"+str(self.iM[0].shape))
            logger.debug("H "+str(self.H[0].shape))

            # tests the mapping:
            n = self.geoP[i + self.geo.nFuselage].shape
            n = n[0]
            dzs = np.zeros(n)
            for k in range(n):
                dzs[k] = 0.01 * self.geoP[i + self.geo.nFuselage][k,1]**2
            self.dzsGlob.append(dzs)
            dza = np.matmul(self.H[i],self.dzsGlob[i])
            self.dzaGlob.append(dza)

        # Plots line
        if plotting:
            fig = plt.figure("figure 2")
            ax = fig.add_subplot(111, projection='3d')
            for p in range(len(self.wingsPoints)):
                # ax.scatter(self.geoP[p + self.geo.nFuselage][:,0],
                #             self.geoP[p + self.geo.nFuselage][:,1],
                #             self.geoP[p + self.geo.nFuselage][:,2],
                #             label='beam wing '+str(p+1))
                # ax.scatter(self.geoP[p + self.geo.nFuselage][:,0],
                #             self.geoP[p + self.geo.nFuselage][:,1],
                #             self.geoP[p + self.geo.nFuselage][:,2]+self.dzsGlob[i],
                #             label='deformed beam wing '+str(p+1))
                ax.scatter(self.wingsPoints[p][:,0],
                           self.wingsPoints[p][:,1],
                           self.wingsPoints[p][:,2],
                           label='undeformed wing'+str(p+1))
                ax.scatter(self.wingsPoints[p][:,0],
                           self.wingsPoints[p][:,1],
                           self.wingsPoints[p][:,2]+self.dzaGlob[p],
                           label='deformed wing'+str(p+1))

            val = 15
            ax.set_xlim(-val,val)
            ax.set_ylim(-val,val)
            ax.set_zlim(-val,val)
            ax.legend()
            plt.show()

    def phi(self,x1,x2,fun):
        """
        Set of radial basis functions that the user can choose of. After some
        test "Linear" seems to be the better choice and only suitable choice.
        This is due to the fact that the aeroelasticity is done with a beam
        hence a line of point and not a surface. So when a RBF is defined with
        a mesh which is smaller than the chord length there is a problem since
        the RBF is zero at the leading a trailing edge.

        All the other functions are here in case someone finds a way to connect
        the solver to some 2D or 3D structure FEM solver.
        """
        eps = 1
        r = LA.norm(x1-x2)
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
        N = len(self.wingsPoints)
        for i in range(N):
            start = self.limitsStart[i]
            end = self.limitsEnd[i]
            self.afx.append(self.VLMdata.panelwise["fx"][start:end])
            self.afy.append(self.VLMdata.panelwise["fy"][start:end])
            self.afz.append(self.VLMdata.panelwise["fz"][start:end])

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

        # logger.debug(self.VLMdata.forces['x'])
        # logger.debug(self.VLMdata.forces['y'])
        # logger.debug(self.VLMdata.forces['z'])
        if not self.geo.settings['G_static']:
            a_x = self.VLMdata.forces['x'] / self.geo.aircraftTotalMass
            a_y = self.VLMdata.forces['y'] / self.geo.aircraftTotalMass
            a_z = (self.VLMdata.forces['z'] - self.geo.aircraftTotalMass*9.81)/self.geo.aircraftTotalMass
        else:
            a_x = 0
            a_y = 0
            a_z = -9.81
        self.G = self.VLMdata.forces['z']/(self.geo.aircraftTotalMass * 9.81)
        # logger.debug('a_x = ' + str(a_x))
        # logger.debug('a_y = ' + str(a_y))
        # logger.debug('a_z = ' + str(a_z))
        # logger.debug('G = '+str(self.G))
        # logger.debug('mass = ' + str(self.geo.aircraftTotalMass))

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
                n = a_z/9.81
                if j == 0:
                    massRight = self.geo.aircraftSegementsMass[i][j]
                    # force in the x direction, earth reference frame
                    force[j,0] = 0.5 * massRight * a_x * 0
                    # force in the x direction, earth reference frame
                    force[j,1] = 0.5 * massRight * a_y
                    # force in the x direction, earth reference frame
                    force[j,2] = 0.5 * massRight * n * 9.81
                elif j == M-1:
                    massLeft = self.geo.aircraftSegementsMass[i][j-1]
                    # force in the x direction, earth reference frame
                    force[j,0] = 0.5 * massLeft * a_x * 0
                    # force in the x direction, earth reference frame
                    force[j,1] = 0.5 * massLeft * a_y
                    # force in the x direction, earth reference frame
                    force[j,2] = 0.5 * massLeft * n * 9.81
                else:
                    massRight = self.geo.aircraftSegementsMass[i][j]
                    massLeft = self.geo.aircraftSegementsMass[i][j-1]
                    # force in the x direction, earth reference frame
                    force[j,0] = 0.5 * massRight * a_x * 0 + \
                                 0.5 * massLeft * a_x * 0
                    # force in the x direction, earth reference frame
                    force[j,1] = 0.5 * massRight * a_y + \
                                 0.5 * massLeft * a_y
                    # force in the x direction, earth reference frame
                    force[j,2] = 0.5 * massRight * n * 9.81 + \
                                 0.5 * massLeft * n * 9.81
            self.smf.append(force)

        # Computes the moment due to inertia on each strcture node
        # for i i
        N = len(self.geo.aircraftMassDistances)
        # logger.debug(N)
        # logger.debug(len(self.smf))
        # sys.exit()
        self.smm = []
        for i in range(N):
            # Since there is one more point then segment we need to ass one at
            # the end.
            M = len(self.geo.aircraftMassDistances[i])
            # logger.debug(M)
            # logger.debug(len(self.smf[i]))
            # sys.exit()
            moments = np.empty((M,3))
            for j in range(M):
                moments[j,0] = self.geo.aircraftMassDistances[i][j,0]
                moments[j,1] = self.geo.aircraftMassDistances[i][j,1]
                moments[j,2] = self.geo.aircraftMassDistances[i][j,2]
            self.smm.append(moments)

    def computeMoments(self):
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
        N = len(self.wingsPoints)
        self.distanceMatrix = []
        self.amx = []
        self.amy = []
        self.amz = []
        for i in range(N):
            M = len(self.wingsPoints[i])
            X = self.wingsPoints[i]
            Y = self.geo.aircraftNodesPoints[i+self.geo.nFuselage]

            # Computes the distance between each point of X and each point of
            # Y. This leads to an (NxM) matrix, M being the number of structure
            # nodes points.
            dist = sp.spatial.distance.cdist(X,Y,"euclidean")
            self.distanceMatrix.append(np.empty((M,3)))
            # Finds the minimal 3 values
            for j in range(M):
                point = self.wingsPoints[i][j]
                indexes = np.argsort(dist[j])[:3]
                # Stores the 3 points of interest
                p1 = self.geo.aircraftNodesPoints[i+self.geo.nFuselage][indexes[0]]
                p2 = self.geo.aircraftNodesPoints[i+self.geo.nFuselage][indexes[1]]
                p3 = self.geo.aircraftNodesPoints[i+self.geo.nFuselage][indexes[2]]

                # Computes the two lines vectors
                v12 = sk.Vector(p2-p1)
                v23 = sk.Vector(p3-p2)

                line1 = sk.Line(point=p2, direction=v12)
                line2 = sk.Line(point=p2, direction=v23)

                proj1 = line1.project_point(point)
                proj2 = line2.project_point(point)

                # Computes the distance between the projected point and the
                # most far way point. This permetis to test if the projection
                # is still on the structural mesh or not
                distP1Proj1 = LA.norm(p1 - proj1)
                distP1P2 = LA.norm(p1 - p2)
                distP3Proj2 = LA.norm(p3 - proj2)
                distP2P3 = LA.norm(p3 - p2)

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

                # line 1 projected point is outside the mesh but not the line 2
                # projected point.
                elif distP1Proj1 > distP1P2 and distP3Proj2 > distP2P3:
                    delta = -(point - p2)
                self.distanceMatrix[i][j] = np.array([delta[0],delta[1],delta[2]])

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
            M = len(self.geoP[i])
            self.sux.append(self.csd.results.get('tensors').get('comp:U')["ux"][old:old+M])
            self.suy.append(self.csd.results.get('tensors').get('comp:U')["uy"][old:old+M])
            self.suz.append(self.csd.results.get('tensors').get('comp:U')["uz"][old:old+M])

            self.stx.append(self.csd.results.get('tensors').get('comp:U')["thx"][old:old+M])
            self.sty.append(self.csd.results.get('tensors').get('comp:U')["thy"][old:old+M])
            self.stz.append(self.csd.results.get('tensors').get('comp:U')["thz"][old:old+M])
            old += M

        # Computes the aerodynamic displacements and the aerodynamic angles
        for i in range(N-self.geo.nFuselage):
            self.aux.append(np.matmul(self.H[i],self.sux[i+self.geo.nFuselage]))
            self.auy.append(np.matmul(self.H[i],self.suy[i+self.geo.nFuselage]))
            self.auz.append(np.matmul(self.H[i],self.suz[i+self.geo.nFuselage]))

            self.atx.append(np.matmul(self.H[i],self.stx[i+self.geo.nFuselage]))
            self.aty.append(np.matmul(self.H[i],self.sty[i+self.geo.nFuselage]))
            self.atz.append(np.matmul(self.H[i],self.stz[i+self.geo.nFuselage]))

        # Assembles the displacements into a big vector
        N = len(self.lattice.c)
        self.displacements = np.empty((N,3))
        N = len(self.geoP)
        old = 0
        coef1 = 1  # for debugging
        coef2 = 1  # for debugging
        for i in range(N-self.geo.nFuselage):
            M = len(self.aux[i])
            for j in range(M):
                # Displacements due to moments:
                dmx = self.aty[i][j] * self.distanceMatrix[i][:,2] + \
                      self.atz[i][j] * self.distanceMatrix[i][:,1]
                dmy = self.atx[i][j] * self.distanceMatrix[i][:,2] + \
                      self.atz[i][j] * self.distanceMatrix[i][:,0]
                dmz = self.atx[i][j] * self.distanceMatrix[i][:,1] + \
                      self.aty[i][j] * self.distanceMatrix[i][:,0]
                self.displacements[old+j][0] = coef1*self.aux[i][j] + coef2*dmx[j]
                self.displacements[old+j][1] = coef1*self.auy[i][j] + coef2*dmy[j]
                self.displacements[old+j][2] = coef1*self.auz[i][j] + coef2*dmz[j]
            old = old + M

        # For debugging only
        if plotting:
            fig = plt.figure("figure 2")
            ax = fig.add_subplot(111, projection='3d')
            # for p in range(len(self.wingsPoints)):
            #     ax.scatter(self.wingsPoints[p][:,0],
            #                self.wingsPoints[p][:,1],
            #                self.wingsPoints[p][:,2],
            #                label='undeformed wing'+str(p+1))
            #     ax.scatter(self.wingsPoints[p][:,0]+self.aux[p],
            #                self.wingsPoints[p][:,1]+self.auy[p],
            #                self.wingsPoints[p][:,2]+self.auz[p],
            #                label='deformed wing'+str(p+1))
            ax.scatter(self.lattice.c[:,0],
                       self.lattice.c[:,1],
                       self.lattice.c[:,2],
                       label='wing')
            ax.scatter(self.lattice.c[:,0] + self.displacements[:,0],
                       self.lattice.c[:,1] + self.displacements[:,1],
                       self.lattice.c[:,2] + self.displacements[:,2],
                       label='deformed wing')
            val = 15
            ax.set_xlim(-val,val)
            ax.set_ylim(-val,val)
            ax.set_zlim(-val,val)
            ax.legend()
            plt.show()
