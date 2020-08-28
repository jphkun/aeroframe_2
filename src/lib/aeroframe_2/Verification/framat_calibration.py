#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 19:47:17 2020

@author: cfse2
"""

# The Model object is used to set up the entire structure model, and to run a
# beam analysis
from framat import Model
import numpy as np

# Create a new instance of the Model object
model = Model()

# ===== MATERIAL =====
# Create a material definition which can be referenced when creating beams.
# Note that you can add as many materials as you want. Just provide a different
# UID (unique identifier) for each new material. Below we define the Young's
# modulus, the shear modulus and the density.
mat = model.add_feature('material', uid='dummy')
mat.set('E', 70e9)
mat.set('G', 27e9)
mat.set('rho', 2100)

# ===== CROSS SECTION =====
# Besides material data, we also need cross section geometry, or more
# specifically, the cross section area, the second moments of area, and the
# torsional constant.
cs = model.add_feature('cross_section', uid='dummy')
cs.set('A', 1e-2)
cs.set('Iy',8e-5)
cs.set('Iz',1e-3)
cs.set('J', 1e-3)

# ===== BEAM =====
# Next, let's add a beam! We define the geometry using "named nodes", that is,
# we provide the coordinates of some "support nodes" which can be referred to
# with their UIDs.
beam = model.add_feature('beam')
beam.add('node', [-5.0,  0, 0], uid='a')
beam.add('node', [-3.75, 0, 0], uid='b')
beam.add('node', [-2.5,  0, 0], uid='c')
beam.add('node', [-1.25, 0, 0], uid='d')
beam.add('node', [ 0.0,  0, 0], uid='e')
beam.add('node', [ 1.25, 0, 0], uid='f')
beam.add('node', [ 2.5,  0, 0], uid='g')
beam.add('node', [ 3.75, 0, 0], uid='h')
beam.add('node', [ 5.0,  0, 0], uid='i')


# Set the number of elements for the beam.
beam.set('nelem', 4)
# Set the material, cross section and cross section orientation
v1 = 'a'
v2 = 'i'
beam.add('material', {'from': v1, 'to': v2, 'uid': 'dummy'})
beam.add('cross_section', {'from': v1, 'to': v2, 'uid': 'dummy'})
beam.add('orientation', {'from': v1, 'to': v2, 'up': [0, 0, 1]})
# Add some line loads [N/m] and point loads [N]
# beam.add('point_load', {'at': 'a', 'load': [12.812, 0.005,    1763.31, 7108.881,  24.792,  -54.809]})
# beam.add('point_load', {'at': 'b', 'load': [15.575, 0.011,    3821.05, 13430.958, 511.137, -60.007]})
# beam.add('point_load', {'at': 'c', 'load': [9.402,  0.014,    4744.77, 11413.5,   845.287, -26.826]})
# beam.add('point_load', {'at': 'd', 'load': [4.908,  0.014,    5100.26, 6182.208,  1022.159, -7.982]})
# beam.add('point_load', {'at': 'e', 'load': [3.447,  0.0 ,     5198.69, 0.001,     1076.735, 0.0]})
# beam.add('point_load', {'at': 'f', 'load': [4.908,  0.014,    5100.26, 6182.208,  1022.159, -7.982]})
# beam.add('point_load', {'at': 'g', 'load': [9.402,  0.014,    4744.77, 11413.5,   845.287, -26.826]})
# beam.add('point_load', {'at': 'h', 'load': [15.575, 0.011,    3821.05, 13430.958, 511.137, -60.007]})
# beam.add('point_load', {'at': 'i', 'load': [12.812, 0.005,    1763.31, 7108.881,  24.792,  -54.809]})

# beam.add('point_load', {'at': 'a', 'load': [0, 0, 1e4, 0, 0, 0]})
# beam.add('point_load', {'at': 'b', 'load': [0, 0, 1e5, 0, 0, 0]})
# beam.add('point_load', {'at': 'c', 'load': [0, 0, 1e4, 0, 0, 0]})
# beam.add('point_load', {'at': 'd', 'load': [0, 0, 1e3, 0, 0, 0]})
# beam.add('point_load', {'at': 'e', 'load': [0, 0, 1e2, 0, 0, 0]})
# beam.add('point_load', {'at': 'f', 'load': [0, 0, 1e3, 0, 0, 0]})
# beam.add('point_load', {'at': 'g', 'load': [0, 0, 1e4, 0, 0, 0]})
# beam.add('point_load', {'at': 'h', 'load': [0, 0, 1e5, 0, 0, 0]})
# beam.add('point_load', {'at': 'i', 'load': [0, 0, 1e4, 0, 0, 0]})

# beam.add('point_load', {'at': 'a', 'load': [12.812, 0.005, 1763.31,0*  7108.881, 0*  24.792,0* -54.809]})
# beam.add('point_load', {'at': 'b', 'load': [15.575, 0.011, 3821.05,0* 13430.958, 0* 511.137,0* -60.007]})
# beam.add('point_load', {'at': 'c', 'load': [9.402,  0.014, 4744.77,0*   11413.5, 0* 845.287,0* -26.826]})
# beam.add('point_load', {'at': 'd', 'load': [4.908,  0.014, 5100.26,0*  6182.208, 0*1022.159,0*  -7.982]})
# beam.add('point_load', {'at': 'e', 'load': [3.447,  0.0 ,  5198.69,0*     0.001, 0*1076.735,0*     0.0]})
# beam.add('point_load', {'at': 'f', 'load': [4.908,  0.014, 5100.26,0*  6182.208, 0*1022.159,0*  -7.982]})
# beam.add('point_load', {'at': 'g', 'load': [9.402,  0.014, 4744.77,0*   11413.5, 0* 845.287,0* -26.826]})
# beam.add('point_load', {'at': 'h', 'load': [15.575, 0.011, 3821.05,0* 13430.958, 0* 511.137,0* -60.007]})
# beam.add('point_load', {'at': 'i', 'load': [12.812, 0.005, 1763.31,0*  7108.881, 0*  24.792,0* -54.809]})


beam.add('point_load', {'at': 'a', 'load': [0*12.812,0* 0.005,0* 1763.31,  7108.881,   24.792, -54.809]})
beam.add('point_load', {'at': 'b', 'load': [0*15.575,0* 0.011,0* 3821.05, 13430.958,  511.137, -60.007]})
beam.add('point_load', {'at': 'c', 'load': [0*9.402, 0* 0.014,0* 4744.77,   11413.5,  845.287, -26.826]})
beam.add('point_load', {'at': 'd', 'load': [0*4.908, 0* 0.014,0* 5100.26,  6182.208, 1022.159,  -7.982]})
beam.add('point_load', {'at': 'e', 'load': [0*3.447, 0* 0.0 , 0* 5198.69,     0.001, 1076.735,     0.0]})
beam.add('point_load', {'at': 'f', 'load': [0*4.908, 0* 0.014,0* 5100.26,  6182.208, 1022.159,  -7.982]})
beam.add('point_load', {'at': 'g', 'load': [0*9.402, 0* 0.014,0* 4744.77,   11413.5,  845.287, -26.826]})
beam.add('point_load', {'at': 'h', 'load': [0*15.575,0* 0.011,0* 3821.05, 13430.958,  511.137, -60.007]})
beam.add('point_load', {'at': 'i', 'load': [0*12.812,0* 0.005,0* 1763.31,  7108.881,   24.792, -54.809]})


# ===== BOUNDARY CONDITIONS =====
# We also must constrain our model. Below, we fix the nodes 'a' and 'd'
bc = model.set_feature('bc')
bc.add('fix', {'node': 'e', 'fix': ['all']})
# bc.add('fix', {'node': 'd', 'fix': ['all']})

# ===== POST-PROCESSING =====
# By default the analysis is run without any GUI, but to get a visual
# representation of the results we can create a plot
# pp = model.set_feature('post_proc')
# pp.set('plot_settings', {'show': True})
# pp.add('plot', ['undeformed', 'deformed', 'node_uids', 'nodes', 'forces'])

# Run the beam analysis
results = model.run()

# ===== RESULTS =====
# The result object contains all relevant results. For instance, we may fetch
# the global load vector.
load_vector = results.get('tensors').get('F')[2::6]
# load_vector = results.get('tensors').get('comp:U')["uz"]
# print(load_vector)
np.set_printoptions(3)
# print(results.get('tensors').get('comp:U')["ux"])
# print(results.get('tensors').get('comp:U')["uy"])
print(results.get('tensors').get('comp:U')["uz"])
