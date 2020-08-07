#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 18:19:27 2020

@author: Jean-Philippe Kuntzer
"""

import logging
import numpy as np
import pandas as pd
from numpy import linalg as LA
from numpy.core.umath_tests import inner1d
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


logger = logging.getLogger(__name__)


class Mesh_Def:

    def __init__(self, lattice,transformCurrent):
        """
        *_p : lattice points
        *r_ p: lattice points reshaped array for ease of use
        *_v : horseshoe vortex points
        *r_v : horseshoe vortex points reshaped for ease of use
        *_c : cell collocation points
        *_n : normal vector directions
        *_b : bound leg midpoints
        f_a : final area of each individual panels
        u_* : defomation of the given type of point. follows the same naming
              pattern as for the previous 5 items
        ----------
        lattice : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.newLattice = transformCurrent
        # stores lattice shapes (only theses two are needed, the others are
        # of identical shape)
        self.s_p = lattice.p.shape
        self.s_v = lattice.v.shape
        self.s_c = lattice.c.shape
        self.s_b = lattice.bound_leg_midpoints.shape
        self.s_n = lattice.n.shape

        # stores lattice intitial (i_) data
        self.i_p = np.copy(lattice.p)
        self.ir_p = self.i_p.reshape((self.s_p[0] * self.s_p[1], self.s_p[2]))
        self.i_v = np.copy(lattice.v)
        self.ir_v = self.i_v.reshape((self.s_v[0] * self.s_v[1], self.s_v[2]))
        self.i_c = np.copy(lattice.c)
        self.i_n = np.copy(lattice.n)
        self.i_b = np.copy(lattice.bound_leg_midpoints)

        # stores lattice final (f_) data
        self.f_p = np.zeros([self.s_p[0], self.s_p[1], self.s_p[2]])
        self.fr_p = np.zeros([self.s_p[0] * self.s_p[1], self.s_p[2]])
        self.f_v = np.zeros([self.s_p[0], self.s_p[1], self.s_p[2]])
        self.fr_v = np.zeros([self.s_p[0] * self.s_p[1], self.s_p[2]])
        self.f_c = np.zeros([self.s_c[0] * self.s_c[1]])
        self.f_n = np.zeros([self.s_c[0] * self.s_c[1]])
        self.f_b = np.zeros([self.s_c[0] * self.s_c[1]])
        self.f_a = np.zeros([self.s_c[0] * self.s_c[1]])

        # Cells absolute y corrdinates (needed for testing and debug)
        self.y_p = np.abs(self.ir_p[:,1])
        self.y_v = np.abs(self.ir_v[:,1])
        self.y_c = np.abs(self.i_c[:,1])
        self.y_b = np.abs(self.i_b[:,1])
        self.x_p = np.abs(self.ir_p[:,0])
        self.x_v = np.abs(self.ir_v[:,0])
        self.x_c = np.abs(self.i_c[:,0])
        self.x_b = np.abs(self.i_b[:,0])

        # Mesh displacement
        self.u_p = np.zeros((self.s_p[0] * self.s_p[1], self.s_p[2]))
        self.u_v = np.zeros((self.s_p[0] * self.s_p[1], self.s_p[2]))
        self.u_c = np.zeros((self.s_c[0], self.s_c[1]))
        self.u_b = np.zeros((self.s_c[0], self.s_c[1]))
        # logger.debug(self.u_p)

    def newPoints(self):
        """
        Loads a displacement file of format .csv and up

        Returns
        -------
        RBF explained
        https://www.youtube.com/watch?v=OOpfU3CvUkM

        None.
        TODO: take into accound the potential rotations! if file is constructed
              with beams.
        """
        logger.debug("=== deformation function called ===")

        x = self.i_c[:,0]
        y = self.i_c[:,1]
        z = self.i_c[:,2]
        d = self.newLattice.displacements
        logger.debug(self.newLattice.displacements)
        logger.debug(self.newLattice.displacements.shape)
        # sys.exit()
        # s = dataset.shape


        # h = list(disp.columns.values)
        # N_headers = len(h)

        # Sorts out which type of FEM simulation was done (beam or shell)
        # TODO: separate the airplane if half using the x axis. At the moment
        #       there is an issue with the center of the airplane.
        # if s[1] == 6:
        logger.info("Input deformation data is of type surface")
        # interpolates the points (lattice.p)
        rbfi = Rbf(x,y,z,d,function='linear',mode="N-D")
        self.u_p = rbfi(self.ir_p[:,0],self.ir_p[:,1],self.ir_p[:,2])

        # interpolates the vortex horseshoe points (lattice.v)
        for i in range(len(self.ir_v)):
            if (i % 4) == 1:
                self.u_v[i] = rbfi(self.ir_v[i,0],
                                   self.ir_v[i,1],
                                   self.ir_v[i,2])
                self.u_v[i-1] = self.u_v[i]
            elif (i % 4) == 2:
                self.u_v[i] = rbfi(self.ir_v[i,0],
                                   self.ir_v[i,1],
                                   self.ir_v[i,2])
                self.u_v[i+1] = self.u_v[i]
        # interpolates the collocation points (lattice.c)
        self.u_c = rbfi(self.i_c[:,0],self.i_c[:,1],self.i_c[:,2])

        # interpolates the bound leg mid-points (lattice.blm)
        self.u_b = rbfi(self.i_b[:,0],self.i_b[:,1],self.i_b[:,2])
        
        # Feed values to the deformed points (f for final).
        self.fr_p = self.ir_p + self.u_p
        self.f_p = self.i_p + self.u_p.reshape(self.s_p[0],
                                               self.s_p[1],
                                               self.s_p[2])
        self.fr_v = self.ir_v + self.u_v
        self.f_v = self.i_v + self.u_v.reshape(self.s_v[0],
                                               self.s_v[1],
                                               self.s_v[2])
        self.f_c = self.i_c + self.u_c
        self.f_b = self.i_b + self.u_b

    def plotNewMesh(self):
        fig = plt.figure("figure 2")
        ax = fig.add_subplot(111, projection='3d')
        # for p in range(len(self.wingsPoints)):
        ax.scatter(self.i_c[:,0],
                   self.i_c[:,1],
                   self.i_c[:,2],
                   label='undeformed wing')
        ax.scatter(self.f_c[:,0],
                   self.f_c[:,1],
                   self.f_c[:,2],
                   label='deformed wing')

        ax.scatter(self.i_b[:,0],
                   self.i_b[:,1],
                   self.i_b[:,2],
                   label='undeformed wing')
        ax.scatter(self.f_b[:,0],
                   self.f_b[:,1],
                   self.f_b[:,2],
                   label='undeformed wing')

        val = 15
        ax.set_xlim(-val,val)
        ax.set_ylim(-val,val)
        ax.set_zlim(-val,val)
        ax.legend()
        plt.show()
        
    def deformation(self):
        """
        This function deforms the mesh, computes the new parameters p, v, c, n
        a. The newly computed parameters are then fed back into the lattice
        class variable in the stdrun.run function.

        The stdrun.run function will then continue to compute the simulation
        with the deformed mesh.

        Parameters
        ----------
        settings : class variable
            Variable of class settings. This variable is used for checking
            which simulation should be done especially during the debug and
            testing phase. It also provides the path of the current simulation

        Returns
        -------
        None.
        """

        logger.info("=== Starts deformation function ===")

        # Computes the initial normal vector for each panel. SVD has a
        # proprety that in the vh matrix, all the vectors are unitary
        # orthonormal vectors. This allows to have the reference frame and
        # compute the angles between old undeformed mesh and the new
        # deformed reference frame for the panel.
        G = np.concatenate((self.i_c, self.i_c, self.i_c, self.i_c), axis=1)
        mat = self.i_p - G.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        u, s, vh_i = LA.svd(mat)

        # user input choice
        # TODO make it simpler for future use
        
        # path = str(settings.paths('f_deformation'))
        framat = ["framat","Framat","FramAT","framAT"]
        logger.debug("Proceed to shape function selection")

        # Computes the new points
        self.newPoints()

        # Computes the deformed reference frame by using the same SVD proprety
        # as before.
        G = np.concatenate((self.f_c, self.f_c, self.f_c, self.f_c), axis=1)
        mat = self.f_p - G.reshape(self.s_p[0],self.s_p[1],self.s_p[2])
        u, s, vh_f = LA.svd(mat)

        # Computes the roation pivot vector. Equivalent of a hinge axis for
        # rotation. This is useful for the quaternion description
        rot_g = np.cross(vh_f[:,2,:],vh_i[:,2,:])
        rot = rot_g / np.linalg.norm(rot_g,axis=1,)[:,np.newaxis]
        rot[np.isnan(rot)] = 0.0

        # Computes the angle between the intial and deformed normal vector
        # dot product of vector "a" (initial state) and "b" (deformed state).
        ab = inner1d(vh_f[:,2,:],vh_i[:,2,:])
        a = LA.norm(vh_f[:,2,:], axis=1)
        b = LA.norm(vh_i[:,2,:], axis=1)
        angle = np.arccos(ab / (a*b))
        angle[np.isnan(angle)] = 0.0

        # Some angles might be computed in the opposite direction, hence being
        # greater than pi/2. A correction is done just below
        corrector = np.zeros(angle.shape)
        corrector[angle > np.pi/2] = np.pi
        angle = angle - corrector

        # Rotates the vector using a "quaternion". It was thought of this way
        # but scipy permits to describe it this way.
        quat = np.einsum('i,ij->ij',-angle,rot)
        r = R.from_rotvec(quat)
        self.f_n = r.apply(self.i_n)

        # Computes the new surface area by using a first order method. Could
        # be impoved but for consistency reasions it is done in the exact same
        # way as how it's computed in the "c_lattice.cpp".
        s = 0.5 * ((self.f_p[:,1] - self.f_p[:,0])
                 + (self.f_p[:,2] - self.f_p[:,3]))
        c = 0.5 * ((self.f_p[:,3] - self.f_p[:,0])
                 + (self.f_p[:,2] - self.f_p[:,1]))
        s = LA.norm(s,axis=1)
        c = LA.norm(c,axis=1)
        # New surface area
        self.f_a = s * c
        # self.plotNewMesh()
        # Checks if small values are present and changes them by a small
        # number
        # self.is_zero()
