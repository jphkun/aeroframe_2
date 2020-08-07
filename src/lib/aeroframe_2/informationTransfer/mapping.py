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
    def __init__(self,lattice,vlmdata,csdGeometry):
        # For debug purposes
        plotting = False
        np.set_printoptions(precision=3)
        # Assembles matrices
        self.geo = csdGeometry
        self.vlm = lattice
        self.data = vlmdata
        self.geoP = csdGeometry.aircraftNodesPoints
        # Separates lattice.c into each wing instances
        self.wingsPoints = []
        self.limitsStart = []
        self.limitsEnd = []
        number = len(lattice.bookkeeping_by_wing_uid)
        for i in lattice.bookkeeping_by_wing_uid:
            # Gets the data for separating the wing points
            listing = list(lattice.bookkeeping_by_wing_uid.get(i)[0][1])
            init = listing[0]
            N = len(list(lattice.bookkeeping_by_wing_uid.get(i)))-1
            # logger.debug(N)
            listing = list(lattice.bookkeeping_by_wing_uid.get(i)[N][1])
            panels = lattice.bookkeeping_by_wing_uid.get(i)[N][2]
            # takes care of last segment
            if number == 1:
                end = listing[-1]
            else:
                end = listing[-1] + panels
            self.limitsStart.append(init)
            self.limitsEnd.append(end)
            # Appends the separated points
            self.wingsPoints.append(lattice.c[init:end])
            # logger.debug("Initial position"+str(init))
            # logger.debug("Final position"+str(listing[-1]))
            # logger.debug("Final position"+str(end))
            # logger.debug("\n")
            number -= 1

        # Plot for debug purposes
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

        # Computes transformations matrices
        self.aPoints = self.vlm.c
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
            for k in range(n):
                for j in range(m):
                    x1 = self.geoP[i + self.geo.nFuselage][k]
                    x2 = self.wingsPoints[i][j]
                    Q[k,j] = self.phi(x1,x2,fun)
            self.A.append(Q.T)
            self.H.append(np.matmul(self.A[i],self.iM[i]))

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
                            # label='deformed beam wing '+str(p+1))
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
       
        self.sfx = []
        self.sfy = []
        self.sfz = []
        self.afx = []  # self.data.panelwise["fx"]
        self.afy = []  # self.data.panelwise["fy"]
        self.afz = []  # self.data.panelwise["fz"]
        # separates froces for each wings
        N = len(self.wingsPoints)
        for i in range(N):
            start = self.limitsStart[i]
            end = self.limitsEnd[i]
            self.afx.append(self.data.panelwise["fx"][start:end])
            self.afy.append(self.data.panelwise["fy"][start:end])
            self.afz.append(self.data.panelwise["fz"][start:end])
        
        for i in range(N):
            logger.debug("aeroToStructure")
            self.sfx.append(np.matmul(self.H[i].T,self.afx[i]))
            self.sfy.append(np.matmul(self.H[i].T,self.afy[i]))
            self.sfz.append(np.matmul(self.H[i].T,self.afz[i]))
        logger.debug("sfx = \n"+str(self.sfx))
            # pass
# xs = np.array([[0,0,0],
#                [0,1,0],
#                [0,2,0],
#                [0,3,0]])
# n = xs.shape

# # Plots line
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(xs[:,0],xs[:,1],xs[:,2], label='beam')


# def phi(x1,x2,eps):
#     norm = np.linalg.norm(x1-x2)
#     r = norm
#     # Gaussian: np.exp(eps*norm)**2
#     # tin plate: r**2 * np.log(r)
#     # multiquadratic: (1+(eps*r)**2)**0.5
#     return (1+(eps*r)**2)**0.5

# # print(xs)
# # print("\n")
# # P = np.ones((1,n[0]))
# # P = np.concatenate((P,xs.T))
# # # 4xN
# # print("P")
# # print(P)
# # print("\n")
# # p_shape = P.shape

# eps = 1
# M = np.empty((n[0],n[0]))
# for i in range(n[0]):
#     for j in range(n[0]):
#         M[i,j] = phi(xs[i],xs[j],eps)        
# print("M")
# print(M)
# print("\n")

# # zeros = np.zeros((p_shape[0],p_shape[0]))
# # Css1 = np.concatenate((zeros, P.T),axis=0)
# # print(Css1)
# # Css2 = np.concatenate((P, M),axis=0)
# # Css = np.concatenate((Css1,Css2),axis=1)

# # np.set_printoptions(precision=1)
# # print("Css")
# # print(Css)
# # print("\n")

# deltaZs = np.array([0.0,
#                     0.1,
#                     0.2,
#                     0.4])
# deltaZs = 0.5*deltaZs
# deltaTxs = np.array([0.0,
#                      0.1,
#                      0.2,
#                      0.4]) 
# deltaTxs = 0.5*deltaTxs

# eigs = np.linalg.eig(M)
# print("eigs")
# print(eigs[0])
# print("\n")
# invM = np.linalg.inv(M)
# print("inv")
# print(invM)
# print("\n")

# # _lambda = np.matmul(invM, deltaXs)
# # print("Lambda")
# # print(_lambda)
# # print(_lambda.shape)
# # print("\n")

# # Plots surface points
# xa = np.array([[-1,0,0],
#                [-1,1,0],
#                [-1,2,0],
#                [-1,3,0],
#                [-1,0.5,0],
#                [-1,1.5,0],
#                [-1,2.5,0],
#                [-1,3.5,0],
#                [ 1,0,0],
#                [ 1,1,0],
#                [ 1,2,0],
#                [ 1,3,0],
#                [ 1,0.5,0],
#                [ 1,1.5,0],
#                [ 1,2.5,0],
#                [ 1,3.5,0],
#                [ 0,0,0],
#                [ 0,1,0],
#                [ 0,2,0],
#                [ 0,3,0],
#                [ 0,0.5,0],
#                [ 0,1.5,0],
#                [ 0,2.5,0],
#                [ 0,3.5,0]])
# ax.scatter(xa[:,0],xa[:,1],xa[:,2], label='surface')
# m = xa.shape

# # print("xa")
# # print(xa)
# # print("\n")
# # Q = np.ones((1,m[0]))
# # Q = np.concatenate((Q,xa.T))
# # # 4xN
# # print("Q")
# # print(Q)
# # q_shape = Q.shape
# # print(q_shape)
# # print("\n")
# eps = 1
# As = np.empty((n[0],m[0]))
# for i in range(n[0]):
#     for j in range(m[0]):
#         As[i,j] = phi(xa[i],xa[j],eps)

# print("A")
# print(As)

# As_shape = As.shape
# # print(k_shape)
# # print("\n")
# # As = np.concatenate((Q.T, K.T),axis=1)
# # print("As")
# # print(As.shape)
# # print(As)
# H_s = np.matmul(As.T,invM)
# print("H_s")
# print(H_s)
# print("\n")
# deltaZa = np.matmul(H_s,deltaZs)
# deltaTxa = np.matmul(H_s,deltaTxs)
# print("Delta X")
# print(deltaZa)
# print("Delta T")
# print(deltaTxa)
# print("\n")

# def tranferRotation(p,b,deltaTxa):
#     # Finds the two closest points
#     # Computes the vector of the beam
#     # Computes the distance
#     # Multiply the distance by deltaTxa
#     # WARNING: Il faut savoir si il est à gauche ou à droite de la ligne
    
#     # Finds the two point
#     N = len(b)
#     dist = np.empty(N)
#     for i in range(N):
#         dist[i] = np.linalg.norm(p-b[i])
#     index1 = np.argmin(dist)
#     print("index1",index1)
#     print("dist 1 = ",dist)
#     dist[index1] = 1e15
#     index2 = np.argmin(dist)
#     print("index2",index2)
#     print("dist 2 = ",dist)
    
#     # Computes the line director vector
#     u = b[index1]-b[index2]
#     AB = p-b[index1]
#     crossProduct = np.cross(AB,u)

#     if LA.norm(u) == 0: u = 1e10
#     d = LA.norm(crossProduct)/LA.norm(u)
#     print("p = ",p)
#     print("b1 = ",b[index1])
#     print("b2 = ",b[index2])
#     print("d = ",d)
#     print("cross",crossProduct)
#     print("\n")
#     if p[0] < b[index1,0]:
#         dz = d*deltaTxa
#     else:
#         dz = -d*deltaTxa
#     return dz
# N = len(xa)
# dz = np.empty(N)
# for i in range(N):
#     dz[i] = tranferRotation(xa[i], xs, deltaTxa[i])
# ax.scatter(xa[:,0],
#            xa[:,1],
#            xa[:,2]+deltaZa[:]+dz[:], label='surface deformed')
# val = 3
# ax.set_xlim(-val,val)
# ax.set_ylim(-val,val)
# ax.set_zlim(-val,val)
# ax.legend()
# plt.show()