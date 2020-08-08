#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:56:01 2020

@author: Jean-Philippe Kuntzer

TODO: For the computation of the moment, distance must be computed in the x,y
      plane.
"""

import logging
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

logger = logging.getLogger(__name__)


class mapper:
    def __init__(self,pytornadoVariables,preMeshedStructre,csdSolverClassVar):
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
        number = len(self.lattice.bookkeeping_by_wing_uid)
        for i in self.lattice.bookkeeping_by_wing_uid:
            # Gets the data for separating the wing points
            listing = list(self.lattice.bookkeeping_by_wing_uid.get(i)[0][1])
            init = listing[0]
            logger.debug("init = "+str(init))
            N = len(list(self.lattice.bookkeeping_by_wing_uid.get(i)))-1
            # logger.debug(N)
            listing = list(self.lattice.bookkeeping_by_wing_uid.get(i)[-1][1])
            panels = self.lattice.bookkeeping_by_wing_uid.get(i)[-1][2]
            # takes care of last segment
            # if number == 1:
            #     end = listing[-1]
            # else:
            end = listing[-1] + panels +1
            logger.debug("number = "+str(number))
            logger.debug("end = "+str(end))
            self.limitsStart.append(init)
            self.limitsEnd.append(end)
            # Appends the separated points
            self.wingsPoints.append(self.lattice.c[init:end])
            # logger.debug("Initial position"+str(init))
            # logger.debug("Final position"+str(listing[-1]))
            # logger.debug("Final position"+str(end))
            # logger.debug("\n")
            number -= 1

        # # Plot for debug purposes
        # if plotting:
        #     fig = plt.figure("figure 1")
        #     ax = fig.add_subplot(111, projection='3d')
        #     for i in range(len(self.wingsPoints)):
        #         ax.scatter(self.wingsPoints[i][:,0],
        #                     self.wingsPoints[i][:,1],
        #                     self.wingsPoints[i][:,2],
        #                     label='Wing '+str(i+1))
        #     val = 15
        #     ax.set_xlim(-val,val)
        #     ax.set_ylim(-val,val)
        #     ax.set_zlim(-val,val)
        #     ax.legend()
        #     plt.show()

        # Computes transformations matrices
        self.aPoints = self.lattice.c
        self.iM= []
        self.A = []
        self.H = []
        # For debugging
        self.dzsGlob = []
        self.dzaGlob = []
        for i in range(len(self.wingsPoints)):
            # Computes the matrix M and then invert it
            # permitted choices are: G,TPS,HMQ,HIMQ,C0,C2,C4,C6,EH see below
            # the definition
            fun = "C2"
            n = self.geoP[i + self.geo.nFuselage].shape
            n = n[0]
            Mbeam = np.zeros((n,n))
            for k in range(n):
                for j in range(n):
                    x1 = self.geoP[i + self.geo.nFuselage][k]
                    x2 = self.geoP[i + self.geo.nFuselage][j]
                    Mbeam[k,j] = self.phi(x1,x2,fun)
            self.iM.append(np.linalg.inv(Mbeam))

            # Computes the matrix Q
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

        # # Plots line
        # if plotting:
        #     fig = plt.figure("figure 2")
        #     ax = fig.add_subplot(111, projection='3d')
        #     for p in range(len(self.wingsPoints)):
        #         # ax.scatter(self.geoP[p + self.geo.nFuselage][:,0],
        #         #             self.geoP[p + self.geo.nFuselage][:,1],
        #         #             self.geoP[p + self.geo.nFuselage][:,2],
        #         #             label='beam wing '+str(p+1))
        #         # ax.scatter(self.geoP[p + self.geo.nFuselage][:,0],
        #         #             self.geoP[p + self.geo.nFuselage][:,1],
        #         #             self.geoP[p + self.geo.nFuselage][:,2]+self.dzsGlob[i],
        #                     # label='deformed beam wing '+str(p+1))
        #         ax.scatter(self.wingsPoints[p][:,0],
        #                    self.wingsPoints[p][:,1],
        #                    self.wingsPoints[p][:,2],
        #                    label='undeformed wing'+str(p+1))
        #         ax.scatter(self.wingsPoints[p][:,0],
        #                    self.wingsPoints[p][:,1],
        #                    self.wingsPoints[p][:,2]+self.dzaGlob[p],
        #                    label='deformed wing'+str(p+1))
        #     val = 15
        #     ax.set_xlim(-val,val)
        #     ax.set_ylim(-val,val)
        #     ax.set_zlim(-val,val)
        #     ax.legend()
        #     plt.show()

    def phi(self,x1,x2,fun):
        """
        set of radial basis functions
        """
        eps = 1
        r = np.linalg.norm(x1-x2)
        if fun == "G":
            # Gaussian
            phi_x = np.exp(-eps*r**2)
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
        Compute the forces for the structure solver from the CFD solver.
        """
        logger.debug("aeroToStructure")       
        self.sfx = []
        self.sfy = []
        self.sfz = []
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

        # Computes the forces that act on the structure
        for i in range(N):
            self.sfx.append(np.matmul(self.H[i].T,self.afx[i]))
            self.sfy.append(np.matmul(self.H[i].T,self.afy[i]))
            self.sfz.append(np.matmul(self.H[i].T,self.afz[i]))
        logger.debug("sfx = \n"+str(self.sfx))

    def structureToAero(self):
        """
        """
        plotting = False
        self.sux = []
        self.suy = []
        self.suz = []
        self.stx = []
        self.sty = []
        self.stz = []
        
        self.aux = []
        self.auy = []
        self.auz = []
        self.atx = []
        self.aty = []
        self.atz = []
        # separates froces for each wings
        self.Us = []
        self.Ts = []
        self.Ua = []
        self.Ta = []
        # number of beams
        N = len(self.geoP)
        logger.debug("N = "+str(N))
        old = 0
        
        for i in range(N):
            # Number of node for each beams
            # logger.debug()
            M = len(self.geoP[i])
            self.sux.append(self.csd.results.get('tensors').get('comp:U')["ux"][old:old+M])
            self.suy.append(self.csd.results.get('tensors').get('comp:U')["uy"][old:old+M])
            self.suz.append(self.csd.results.get('tensors').get('comp:U')["uz"][old:old+M])
            old += M
        for i in range(N-self.geo.nFuselage):
            self.aux.append(np.matmul(self.H[i],self.sux[i+self.geo.nFuselage]))
            self.auy.append(np.matmul(self.H[i],self.suy[i+self.geo.nFuselage]))
            self.auz.append(np.matmul(self.H[i],self.suz[i+self.geo.nFuselage]))
        logger.debug(self.suz[0])
        logger.debug(self.auz[0])
        # sys.exit()
        # Assembles the displacements into a big vector
        N = len(self.lattice.c)
        self.displacements = np.empty((N,3))
        N = len(self.geoP)
        logger.debug("N = "+str(N))
        old = 0
        for i in range(N-self.geo.nFuselage):
            logger.debug(i+self.geo.nFuselage)
            M = len(self.aux[i])
            for j in range(M):
                self.displacements[old+j][0] = self.aux[i][j]
                self.displacements[old+j][1] = self.auy[i][j]
                self.displacements[old+j][2] = self.auz[i][j]
            old = old + M
        # Plots

        if plotting:
            fig = plt.figure("figure 2")
            ax = fig.add_subplot(111, projection='3d')
            # for p in range(len(self.wingsPoints)):
            #     ax.scatter(self.wingsPoints[p][:,0],
            #                 self.wingsPoints[p][:,1],
            #                 self.wingsPoints[p][:,2],
            #                 label='undeformed wing'+str(p+1))
            #     ax.scatter(self.wingsPoints[p][:,0]+self.aux[p],
            #                 self.wingsPoints[p][:,1]+self.auy[p],
            #                 self.wingsPoints[p][:,2]+self.auz[p],
            #                 label='deformed wing'+str(p+1))
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
    def deformVLM(self):
        pass
            