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
        self.geo = csdGeometry
        self.vlm = lattice
        self.data = vlmdata
        self.geoP = csdGeometry.aircraftNodesPoints
        logger.debug("Lattice shape = "+str(self.vlm.p.shape))
        logger.debug("Lattice shape = "+str(self.vlm.v.shape))
        logger.debug("Lattice shape = "+str(self.vlm.c.shape))
        sys.exit()
        # assembles all structural points into one vector for ease of use and
        # computation speed.
        # if self.geo.nFuselage > 0:
        #     self.sPoints = self.geoP[1]
        #     Np = len(self.geoP)-2
        #     for i in range(Np):
        #          self.sPoints = np.concatenate((self.sPoints, self.geoP[i+2]))
        # else:
        #     self.sPoints = self.geoP[0]
        #     Np = len(self.geoP)-1
        #     for i in range(Np):
        #          self.sPoints = np.concatenate((self.sPoints, self.geoP[i+1]))
        
        self.aPoints = self.vlm.c
        plotting = False
        
        self.nBeams = len(self.geoP)
        self.M = []
        self.iM= []
        self.A = []
        self.H = []
        for i in range(self.nBeams):

            # logger.debug("sPoint = \n" +str(self.sPoints))
            # logger.debug("vlm.c = \n" + str(self.vlm.c))
            # G,TPS,HMQ,HIMQ,C0,C2,C4,C6,EH
            fun = "C2"
            logger.debug(type(self.geoP[i]))
            # Computes the matrix M
            n = self.geoP[i].shape
            n = n[0]
            Mbeam = np.empty((n,n))
            for k in range(n):
                for j in range(n):
                    x1 = self.geoP[i][k]
                    x2 = self.geoP[i][j]
                    Mbeam[k,j] = self.phi(x1,x2,fun)
            self.M.append(Mbeam)
            
            # logger.debug("M = "+str(self.M))
            # sys.exit()
            # np.set_printoptions(precision=1)
            # logger.debug(M)
            self.iM.append(np.linalg.inv(self.M[i]))
            # Computes the matrix Q
            m = self.aPoints.shape
            m = m[0]
            Q = np.empty((n,m))
            for k in range(n):
                for j in range(m):
                    x1 = self.geoP[i][k]
                    x2 = self.aPoints[j]
                    Q[k,j] = self.phi(x1,x2,fun)
            
            self.A.append(Q.T)
            # logger.debug("Q shape:"+str(Q.shape))
            # logger.debug("A shape:"+str(A.shape))
            np.set_printoptions(precision=1)
            # logger.debug("M = \n"+str(M))
            # logger.debug("inv(M) = \n"+str(iM))
            # logger.debug("A = \n"+str(A))
        
            # Computes coupling matrix H
            self.H.append(np.matmul(self.A[i],self.iM[i]))
            # logger.debug("H = \n"+str(self.H))
            # tests the mapping:
            # dzs = np.empty(n)
            # for i in range(n):
            #     dzs[i] = 0.01*self.sPoints[i,1]**2
            # dza = np.matmul(self.H,dzs)
        
        # # # Plots line
        # # if plotting:
        # #     fig = plt.figure()
        # #     ax = fig.add_subplot(111, projection='3d')
        # #     ax.scatter(self.sPoints[:,0],
        # #                self.sPoints[:,1],
        # #                self.sPoints[:,2],
        # #                label='beam')
        # #     ax.scatter(self.sPoints[:,0],
        # #                self.sPoints[:,1],
        # #                self.sPoints[:,2]+dzs,
        # #                label='deformed beam')
        # #     ax.scatter(self.aPoints[:,0],
        # #                self.aPoints[:,1],
        # #                self.aPoints[:,2],
        # #                label='VLM')
        # #     ax.scatter(self.aPoints[:,0],
        # #                self.aPoints[:,1],
        # #                self.aPoints[:,2]+dza,
        # #                label='deformed VLM')
        # #     val = 15
        # #     ax.set_xlim(-val,val)
        # #     ax.set_ylim(-val,val)
        # #     ax.set_zlim(-val,val)
        # #     ax.legend()
        # #     plt.show()

        # sys.exit()


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
        for i in range(self.nBeams):
            logger.debug("aeroToStructure")
            afx = self.data.panelwise["fx"]
            afy = self.data.panelwise["fy"]
            afz = self.data.panelwise["fz"]
            # nx = self.data.panelwise["fx"].shape
            # ny = self.data.panelwise["fy"].shape
            # nz = self.data.panelwise["fz"].shape
            # logger.debug("nx = "+str(nx))
            # logger.debug("ny = "+str(ny))
            # logger.debug("nz = "+str(nz))
            # N = self.vlm.c.shape
            # logger.debug("N = "+str(N))
            self.sfx.append(np.matmul(self.H[i].T,afx))
            self.sfy.append(np.matmul(self.H[i].T,afy))
            self.sfz.append(np.matmul(self.H[i].T,afz))
            # snx = self.sfx[i].shape
            # sny = self.sfy[i].shape
            # snz = self.sfz[i].shape
            # logger.debug("afx = \n"+str(afx))
            # logger.debug("afy = \n"+str(afy))
            # logger.debug("afz = \n"+str(afz))
            # logger.debug("sfx = \n"+str(self.sfx))
            # logger.debug("sfy = \n"+str(self.sfy))
            # logger.debug("sfz = \n"+str(self.sfz))
            # logger.debug(snx)
            # logger.debug(sny)
            # logger.debug(snz)
            # sys.exit()
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