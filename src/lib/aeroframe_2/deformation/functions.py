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


logger = logging.getLogger(__name__)


class Mesh_Def:

    def __init__(self, lattice,file=None):
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
        self.file = file
        
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

    def mesh_update(self):
        """
        Feeds deformed values back to the mesh

        Returns
        -------
        None.

        """
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

    def cantilever(self, y, q, L, E, Ix):
        return q*y**2 * (6*L**2 - 4*L*y + y**2) / (24*E*Ix)

    def shape_1(self,settings):
        """
        Shape function 1. This functions respresents line of slope 1. The only
        purpouse of this shape function is to test if the deformation is done
        correctly.
        """
        csv_save = True
        const = 0.0
        logger.info("Shape function 1 is selected")
        case = settings.settings["aircraft"]
        if case == "1_flat_funcActivated.json":
            m = 0.0
        elif case == "2_dih_funcActivated.json":
            m = 0.1
        elif case == "3_anh_funcActivated.json":
            m = -0.1
        elif case == "4_flat_funcActivated.json":
            m = 0
        elif case == "5_dih_funcActivated.json":
            m = 0.1
        elif case == "6_anh_funcActivated.json":
            m = -0.1
        elif case == "7_flat_funcActivated.json":
            m = 0
        elif case == "8_dih_funcActivated.json":
            m = 0.1
        elif case == "9_anh_funcActivated.json":
            m = -0.1
        elif case == "10_flat_funcActivated.json":
            m = 0
        elif case == "11_dih_funcActivated.json":
            m = 0.1
        elif case == "12_anh_funcActivated.json":
            m = -0.1
        elif case == "AircraftMalo-std_funActivated.xml":
            m = const
        elif case == "B7772VSP_v3.1_funActivated.xml":
            m = const
        elif case == "BWB_102_VTP1_v3.1_funActivated.xml":
            m = const
        elif case == "BWB_ACFA_cpacs_v3.1_funActivated.xml":
            m = const
        elif case == "Circlewing_Test.v_3.1_funActivated.xml":
            m = const
        elif case == "D150_AGILE_Hangar_funActivated.xml":
            m = const
        elif case == "Boxwing_AGILE_Hangar_funActivated_v3.1.xml":
            m = const
        elif case == "Optimale_Tornado_SU2_funActivated.xml":
            m = const
        elif case == "Optimale.xml":
            m = 1e-6
        elif case == "EbeeX_d0_q0.xml":
            m = 0.1975
        elif "Wing_" in case:
            logger.debug(int(case[15:-4]))
            m = np.linspace(0,0.1,13)
            m = m[int(case[15:-4])-1]
            
        else:
            logger.warning("Deformation input UNEXPECTED")

        h = 0

        self.u_p[:,2] = m * self.y_p + h
        self.u_v[:,2] = m * self.y_v + h
        self.u_c[:,2] = m * self.y_c + h
        self.u_b[:,2] = m * self.y_b + h
        self.mesh_update()
        # logger.debug(self.u_p.shape)

        # Saves data to csv. This part is put here ease of use during debug
        # phase. "csv_save" will be set to false when this phase is done
        # TODO set "csv_save" to false when debug is finished
        if csv_save:
            headers = ["x","y","z","dx","dy","dz"]
            points = np.concatenate((self.i_c,self.u_c),axis=1)
            filepath = str(settings.paths('f_deformation'))
            name = "deformation_data.csv"
            dataset = pd.DataFrame(points,columns=headers)
            dataset.to_csv(filepath[:-4]+name,index=False,float_format='%.18E')
            logger.info("csv file saved")

    def shape_2(self,settings):
        """
        Shape function2. This function computes the slope at each y location
        for a cantilever beam of length "L", with a distributed load of "q".
        The Young modulus is imposed for steel and "I" the second moment of
        inertia is also imposed.
        """
        csv_save = True
        case = settings.settings["aircraft"]

        logger.info("Shape function 2 is selected")
        if case == "EbeeX_d0_q0.xml":
            m = 0.0
        if case == "Optimale.xml":
            m = 0.0
        else:
            logger.warning("Deformation input UNEXPECTED")

        h = 0
        self.u_p[:,2] = m * self.y_p + h
        self.u_v[:,2] = m * self.y_v + h
        self.u_c[:,2] = m * self.y_c + h
        self.u_b[:,2] = m * self.y_b + h
        cst2 = 1e-1
        self.u_p[:,2] = self.u_p[:,2] + cst2*self.y_p**2
        self.u_v[:,2] = self.u_v[:,2] + cst2*self.y_v**2
        self.u_c[:,2] = self.u_c[:,2] + cst2*self.y_c**2
        self.u_b[:,2] = self.u_b[:,2] + cst2*self.y_b**2
        self.mesh_update()
        if csv_save:
            headers = ["x","y","z","dx","dy","dz"]
            points = np.concatenate((self.i_c,self.u_c),axis=1)
            filepath = str(settings.paths('f_deformation'))
            name = "deformation_data.csv"
            dataset = pd.DataFrame(points,columns=headers)
            dataset.to_csv(filepath[:-4]+name,index=False,float_format='%.18E')
            logger.info("csv file saved")

        # # [N/m] Distributed load
        # q = 2000000
        # # [m] Wing span
        # L = 1
        # # logger.debug("L = " + str(L) )
        # # [Pa] Elasticity modulus
        # E = 210e9
        # # [m**4] Second moment of area
        # Ix = 1.330e-6
        # # Computes beam deformation
        # self.u_p[:,2] = self.cantilever(self.y_p, q, L, E, Ix)
        # self.u_v[:,2] = self.cantilever(self.y_v, q, L, E, Ix)
        # self.u_c[:,2] = self.cantilever(self.y_c, q, L, E, Ix)
        # self.u_b[:,2] = self.cantilever(self.y_b, q, L, E, Ix)
        # # logger.debug(self.i_v)
        # self.mesh_update()
        # s = self.f_v.shape
        # var = self.f_v.reshape(s[0]*s[1],s[2]) - self.i_v.reshape(s[0]*s[1],s[2])
        # logger.debug(np.max(var[:,0]))

    def csv_deformation(self,settings):
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
        logger.debug("=== csv deformation function called ===")
        path = self.file
        if path is None:
            logger.error("NO DEFORMATION FILE")
        
        logger.debug(path)
        try:
            dataset = pd.read_csv(path)
            dataset = dataset.to_numpy()
            x = dataset[:,0]
            y = dataset[:,1]
            z = dataset[:,2]
            d = dataset[:,3:]
            s = dataset.shape
            # logger.debug(y)

            # angle = angle - corrector
        except FileNotFoundError:
            logger.error("No such deformation file or directiory" + str(path))

        # h = list(disp.columns.values)
        # N_headers = len(h)

        # Sorts out which type of FEM simulation was done (beam or shell)
        # TODO: separate the airplane if half using the x axis. At the moment
        #       there is an issue with the center of the airplane.
        if s[1] == 6:
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

    def is_zero(self):
        pass
        val = 1e-4
        # self.f_p[np.abs(self.f_p) < val] = 0
        # # logger.debug(self.f_p)
        # self.f_v[np.abs(self.f_v) < val] = 0
        # # self.f_c[np.abs(self.f_c) < val] = val
        # # self.f_n[np.abs(self.f_n) < val] = val
        # self.f_b[np.abs(self.f_b) < val] = 0
        # # self.f_a[np.abs(self.f_a) < val] = val

    def deformation(self,settings):
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
        
        path = str(settings.paths('f_deformation'))
        logger.debug("Proceed to shape function selection")
        if settings.settings["deformation_method"] == "shape_1":
            logger.debug(settings.settings["deformation_method"])
            self.shape_1(settings)
        elif settings.settings["deformation_method"] == "shape_2":
            logger.debug(settings.settings["deformation_method"])
            self.shape_2(settings)
        elif settings.settings["deformation_method"] == "load_from_csv":
            self.csv_deformation(settings)
        else:
            logger.error("No shape function selected")
        # else:
        #     # TODO implement the deformation via CSV file here
        #     logger.error("Deformation from file")
        #     self.csv_deformation(settings)

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

        # Checks if small values are present and changes them by a small
        # number
        self.is_zero()
