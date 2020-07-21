#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:48:39 2020

@author: Jean-Philippe Kuntzer
"""

##############################################################################
#
#   CASE 1
#
##############################################################################

# from framat import Model                  # Import the FramAT model object
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# model = Model.example()  # Load the cantilever example model
# model.set_feature('post_proc')
# model.add('plot', ['undeformed', 'deformed', 'nodes'])
# result = model.run()                               # Run the pre-defined analysis

# # Mesh data
# # print(model.get("beam"))

# def add_deformed_undeformed(m):
#     # to_show = m.get('post_proc').get('plot')[plot_num]
#     abm = m.results.get('mesh').get('abm')

#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.gca(projection='3d')

#     for beam_id in abm.beams.keys():
#         xyz = abm.get_sup_points(beam_id)
#         x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
#         ax.plot(x, y, z)

#         # ==================
#         d = m.results.get('beam')[0].get('deformation')
#         xd = x + d['ux']
#         yd = y + d['uy']
#         zd = z + d['uz']
#         ax.plot(xd, yd, zd)

# # Deformationd data
# print(result.get("beam")[0].get('deformation')["ux"])
# print(result.get("beam")[0].get('deformation')["uy"])
# print(result.get("beam")[0].get('deformation')["uz"])
# print(result.get("beam")[0].get('deformation')["thx"])
# print(result.get("beam")[0].get('deformation')["thy"])
# print(result.get("beam")[0].get('deformation')["thz"])
# add_deformed_undeformed(model)

##############################################################################
#
#   CASE 2
#
##############################################################################

# from framat import Model


# model = Model()

# mat = model.add_feature('material', uid='dummy')
# mat.set('E', 1)
# mat.set('G', 1)
# mat.set('rho', 1)

# cs = model.add_feature('cross_section', uid='dummy')
# cs.set('A', 1)
# cs.set('Iy', 1)
# cs.set('Iz', 1)
# cs.set('J', 1)

# beam = model.add_feature('beam')
# beam.add('node', [0, 0, 0], uid='root')
# beam.add('node', [1, 0, 0], uid='corner')
# beam.add('node', [1, 1, 0], uid='tip')
# beam.set('nelem', 10)
# beam.add('material', {'from': 'root', 'to': 'tip', 'uid': 'dummy'})
# beam.add('cross_section', {'from': 'root', 'to': 'tip', 'uid': 'dummy'})
# beam.add('orientation', {'from': 'root', 'to': 'tip', 'up': [0, 0, 1]})
# beam.add('point_load', {'at': 'corner', 'load': [0, 0, -1, 0, 0, 0]})
# 5
# bc = model.set_feature('bc')
# bc.add('fix', {'node': 'root', 'fix': ['all']})

# # pp = model.set_feature('post_proc')
# # pp.set('plot_settings', {'Show': True})# Not needed anymore in this version
# # pp.add('plot', ('undeformed'))

# results = model.run()
# print(results.get("matrices").get("F_react"))


##############################################################################
#
#   CASE 3
#
##############################################################################

# The Model object is used to set up the entire structure model, and to run a
# beam analysis
from framat import Model

# Create a new instance of the Model object
model = Model()

# ===== MATERIAL =====
# Create a material definition which can be referenced when creating beams.
# Note that you can add as many materials as you want. Just provide a different
# UID (unique identifier) for each new material. Below we define the Young's
# modulus, the shear modulus and the density.
mat = model.add_feature('material', uid='dummy')
mat.set('E', 1)
mat.set('G', 1)
mat.set('rho', 1)

# ===== CROSS SECTION =====
# Besides material data, we also need cross section geometry, or more
# specifically, the cross section area, the second moments of area, and the
# torsional constant.
cs = model.add_feature('cross_section', uid='dummy')
cs.set('A', 1)
cs.set('Iy', 1)
cs.set('Iz', 1)
cs.set('J', 1)

# ===== BEAM =====
# Next, let's add a beam! We define the geometry using "named nodes", that is,
# we provide the coordinates of some "support nodes" which can be referred to
# with their UIDs.
beam = model.add_feature('beam')
beam.add('node', [0.0, 0, 0], uid='a')
beam.add('node', [1.5, 0, 0], uid='b')
beam.add('node', [1.5, 3, 0], uid='c')
beam.add('node', [0.0, 3, 0], uid='d')
# Set the number of elements for the beam.
beam.set('nelem', 40)
# Set the material, cross section and cross section orientation
beam.add('material', {'from': 'a', 'to': 'd', 'uid': 'dummy'})
beam.add('cross_section', {'from': 'a', 'to': 'd', 'uid': 'dummy'})
beam.add('orientation', {'from': 'a', 'to': 'd', 'up': [0, 0, 1]})
# Add some line loads [N/m] and point loads [N]
beam.add('distr_load', {'from': 'a', 'to': 'b', 'load': [0, 0, -2, 0, 0, 0]})
# beam.add('distr_load', {'from': 'c', 'to': 'd', 'load': [0, 0, -2, 0, 0, 0]})
# beam.add('distr_load', {'from': 'b', 'to': 'c', 'load': [0, 0, 1, 0, 0, 0]})
beam.add('point_load', {'at': 'b', 'load': [+0.1, +0.2, +0.3, 0, 0, 0]})
beam.add('point_load', {'at': 'c', 'load': [-0.1, -0.2, -0.3, 0, 0, 0]})

# ===== BOUNDARY CONDITIONS =====
# We also must constrain our model. Below, we fix the nodes 'a' and 'd'
bc = model.set_feature('bc')
bc.add('fix', {'node': 'a', 'fix': ['all']})
bc.add('fix', {'node': 'd', 'fix': ['all']})

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
load_vector = results.get('tensors').get('F')
print(load_vector)