#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:56:01 2020

@author: Jean-Philippe Kuntzer
"""

import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

xs = np.array([[0,0,0],
               [0,1,0],
               [0,2,0],
               [0,3,0]])
n = xs.shape

# Plots line
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs[:,0],xs[:,1],xs[:,2], label='beam')


def phi(x1,x2,eps):
    norm = np.linalg.norm(x1-x2)
    r = norm
    # Gaussian: np.exp(eps*norm)**2
    # tin plate: r**2 * np.log(r)
    # multiquadratic: (1+(eps*r)**2)**0.5
    return (1+(eps*r)**2)**0.5

# print(xs)
# print("\n")
# P = np.ones((1,n[0]))
# P = np.concatenate((P,xs.T))
# # 4xN
# print("P")
# print(P)
# print("\n")
# p_shape = P.shape

eps = 1
M = np.empty((n[0],n[0]))
for i in range(n[0]):
    for j in range(n[0]):
        M[i,j] = phi(xs[i],xs[j],eps)        
print("M")
print(M)
print("\n")

# zeros = np.zeros((p_shape[0],p_shape[0]))
# Css1 = np.concatenate((zeros, P.T),axis=0)
# print(Css1)
# Css2 = np.concatenate((P, M),axis=0)
# Css = np.concatenate((Css1,Css2),axis=1)

# np.set_printoptions(precision=1)
# print("Css")
# print(Css)
# print("\n")

deltaZs = np.array([0.0,
                    0.1,
                    0.2,
                    0.4])
deltaZs = 0.5*deltaZs
deltaTxs = np.array([0.0,
                     0.1,
                     0.2,
                     0.4]) 
deltaTxs = 0.5*deltaTxs

eigs = np.linalg.eig(M)
print("eigs")
print(eigs[0])
print("\n")
invM = np.linalg.inv(M)
print("inv")
print(invM)
print("\n")

# _lambda = np.matmul(invM, deltaXs)
# print("Lambda")
# print(_lambda)
# print(_lambda.shape)
# print("\n")

# Plots surface points
xa = np.array([[-1,0,0],
               [-1,1,0],
               [-1,2,0],
               [-1,3,0],
               [-1,0.5,0],
               [-1,1.5,0],
               [-1,2.5,0],
               [-1,3.5,0],
               [ 1,0,0],
               [ 1,1,0],
               [ 1,2,0],
               [ 1,3,0],
               [ 1,0.5,0],
               [ 1,1.5,0],
               [ 1,2.5,0],
               [ 1,3.5,0],
               [ 0,0,0],
               [ 0,1,0],
               [ 0,2,0],
               [ 0,3,0],
               [ 0,0.5,0],
               [ 0,1.5,0],
               [ 0,2.5,0],
               [ 0,3.5,0]])
ax.scatter(xa[:,0],xa[:,1],xa[:,2], label='surface')
m = xa.shape

# print("xa")
# print(xa)
# print("\n")
# Q = np.ones((1,m[0]))
# Q = np.concatenate((Q,xa.T))
# # 4xN
# print("Q")
# print(Q)
# q_shape = Q.shape
# print(q_shape)
# print("\n")
eps = 1
As = np.empty((n[0],m[0]))
for i in range(n[0]):
    for j in range(m[0]):
        As[i,j] = phi(xa[i],xa[j],eps)

print("A")
print(As)

As_shape = As.shape
# print(k_shape)
# print("\n")
# As = np.concatenate((Q.T, K.T),axis=1)
# print("As")
# print(As.shape)
# print(As)
H_s = np.matmul(As.T,invM)
print("H_s")
print(H_s)
print("\n")
deltaZa = np.matmul(H_s,deltaZs)
deltaTxa = np.matmul(H_s,deltaTxs)
print("Delta X")
print(deltaZa)
print("Delta T")
print(deltaTxa)
print("\n")

def tranferRotation(p,b,deltaTxa):
    # Finds the two closest points
    # Computes the vector of the beam
    # Computes the distance
    # Multiply the distance by deltaTxa
    # WARNING: Il faut savoir si il est à gauche ou à droite de la ligne
    
    # Finds the two point
    N = len(b)
    dist = np.empty(N)
    for i in range(N):
        dist[i] = np.linalg.norm(p-b[i])
    index1 = np.argmin(dist)
    print("index1",index1)
    print("dist 1 = ",dist)
    dist[index1] = 1e15
    index2 = np.argmin(dist)
    print("index2",index2)
    print("dist 2 = ",dist)
    
    # Computes the line director vector
    u = b[index1]-b[index2]
    AB = p-b[index1]
    crossProduct = np.cross(AB,u)

    if LA.norm(u) == 0: u = 1e10
    d = LA.norm(crossProduct)/LA.norm(u)
    print("p = ",p)
    print("b1 = ",b[index1])
    print("b2 = ",b[index2])
    print("d = ",d)
    print("cross",crossProduct)
    print("\n")
    if p[0] < b[index1,0]:
        dz = d*deltaTxa
    else:
        dz = -d*deltaTxa
    return dz
N = len(xa)
dz = np.empty(N)
for i in range(N):
    dz[i] = tranferRotation(xa[i], xs, deltaTxa[i])
ax.scatter(xa[:,0],
           xa[:,1],
           xa[:,2]+deltaZa[:]+dz[:], label='surface deformed')
val = 3
ax.set_xlim(-val,val)
ax.set_ylim(-val,val)
ax.set_zlim(-val,val)
ax.legend()
plt.show()