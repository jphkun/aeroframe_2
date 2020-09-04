#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 10:08:44 2020

@author: Jean-Philippe Kuntzer
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

cwd1 = '/home/cfse2/Documents/aeroframe_2/test/static/0_WingValidationPytornadoFramAT_case0/CFD/_results/'
cwd2 = '/home/cfse2/Documents/aeroframe_2/test/static/0_WingValidationPytornadoFramAT_case0/'
# Unloaded beam
filename_1 = '10_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_2 = '20_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_3 = '40_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_4 = '80_NoSL_NoAL_FEM_displacementAndRotations0.csv'

# Weighted beam with no aerodynamic forces
filename_6 = '10_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_7 = '20_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_8 = '40_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_9 = '80_SL_NoAL_FEM_displacementAndRotations0.csv'

# Influence of aerodynamic mesh
filename_10 = '10_NoSL_05x20AL_results.csv'
filename_11 = '10_NoSL_05x40AL_results.csv'
filename_12 = '10_NoSL_10x40AL_results.csv'
filename_13 = '10_NoSL_10x80AL_results.csv'

# Influence of aerodynamic mesh
filename_14 = '20_NoSL_05x20AL_results.csv'
filename_15 = '20_NoSL_05x40AL_results.csv'
filename_16 = '20_NoSL_10x40AL_results.csv'
filename_17 = '20_NoSL_10x80AL_results.csv'

# Influence of aerodynamic mesh
filename_18 = '40_NoSL_05x20AL_results.csv'
filename_19 = '40_NoSL_05x40AL_results.csv'
filename_20 = '40_NoSL_10x40AL_results.csv'
filename_21 = '40_NoSL_10x80AL_results.csv'

# Influence of aerodynamic mesh
filename_22 = '80_NoSL_05x20AL_results.csv'
filename_23 = '80_NoSL_05x40AL_results.csv'
filename_24 = '80_NoSL_10x40AL_results.csv'
filename_25 = '80_NoSL_10x80AL_results.csv'

filenames_NoSL_NoAL = [filename_4, filename_3, filename_2, filename_1]
dataFrames_NoSL_NoAL = []
for filename in filenames_NoSL_NoAL:
    df = pd.read_csv(cwd1 + filename)
    dataFrames_NoSL_NoAL.append(df)

filenames_SL_NoAL = [filename_6, filename_7, filename_8, filename_9]
dataFrames_SL_NoAL = []
for filename in filenames_SL_NoAL:
    df = pd.read_csv(cwd1 + filename)
    dataFrames_SL_NoAL.append(df)

# Influence of aerodynamic mesh
filenames_10_NoSL_AL = [filename_10, filename_11, filename_12, filename_13]
dataFrames_10NoSL_AL = []
for filename in filenames_10_NoSL_AL:
    df = pd.read_csv(cwd2 + filename,sep=';')
    dataFrames_10NoSL_AL.append(df)

# Influence of aerodynamic mesh
filenames_20_NoSL_AL = [filename_14, filename_15, filename_16, filename_17]
dataFrames_20NoSL_AL = []
for filename in filenames_20_NoSL_AL:
    df = pd.read_csv(cwd2 + filename,sep=';')
    dataFrames_20NoSL_AL.append(df)

# Influence of aerodynamic mesh
filenames_40NoSL_AL = [filename_18, filename_19, filename_20, filename_21]
dataFrames_40NoSL_AL = []
for filename in filenames_40NoSL_AL:
    df = pd.read_csv(cwd2 + filename,sep=';')
    dataFrames_40NoSL_AL.append(df)

# Influence of aerodynamic mesh
filenames_80NoSL_AL = [filename_22, filename_23, filename_24, filename_25]
dataFrames_80NoSL_AL = []
for filename in filenames_80NoSL_AL:
    df = pd.read_csv(cwd2 + filename,sep=';')
    dataFrames_80NoSL_AL.append(df)

# Influence of strcutre mesh
filenames_NoSL_10x80AL = [filename_13, filename_17, filename_21, filename_25]
dataFrames_NoSL_10x80AL = []
for filename in filenames_NoSL_10x80AL:
    df = pd.read_csv(cwd2 + filename,sep=';')
    dataFrames_NoSL_10x80AL.append(df)

# Graph properties
titleSize = 30
textSize = 20
width = 3


# ##############################################################################
# # Convergeance of unloaded beam
# ##############################################################################
# plt.figure(1)
# plt.title('Convergeance of unloaded beam',fontsize=titleSize)
# for dataFrame, filename in zip(dataFramesNoSL_NoAL, filenamesNoSL_NoAL):
#     plt.plot(dataFrame['y'],dataFrame['dz'],'-o',
#              linewidth=width,label='Nodes: '+filename[:2])
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Wingspan position [m]',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.show()


# ##############################################################################
# # Convergeance of weighted beam
# ##############################################################################
# def w(q,y,L,E,I):
#     """
#     Theoretical result of a cantilever beam
#     """
#     y = np.abs(y)
#     w = q*y**2 * (6*L**2 - 4*L*y + y**2)/(24*E*I)
#     return w

# L = 10
# rho = 2100
# A = 3.591149e-01
# I = 4E-05
# E = 70e9
# y = np.linspace(-L/2,L/2,1000)
# m = L*A*rho/L
# mg = -9.81*m

# plt.figure(1)
# plt.title('Convergeance of weighted beam',fontsize=titleSize)
# for dataFrame, filename in zip(dataFramesSL_NoAL, filenamesSL_NoAL):
#     err = (np.min(dataFrame['dz']) - np.min(w(mg,y,L/2,E,I)))
#     err = (np.min(w(mg,y,L/2,E,I)) - np.min(dataFrame['dz']))
#     # print(np.max(w(mg,y,L/2,E,I)))
#     # print(np.max(dataFrame['dz']))
#     err = 100*err / np.min(w(mg,y,L/2,E,I))
#     err = np.abs(err.round(decimals=3))
#     plt.plot(dataFrame['y'],dataFrame['dz'],'-o',
#              linewidth=width,label='Nodes: '+filename[:2]+', error:'+str(err)+'%')
# plt.plot(y,w(mg,y,L/2,E,I),linewidth=width,label='Analytical solution',)
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Wingspan position [m]',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.show()

##############################################################################
# Convergeance of aerodynamic mesh
##############################################################################

# plt.figure(1)
# plt.title('Aeroelastic convergence with 10 structure nodes',fontsize=titleSize)
# for dataFrame, filename in zip(dataFrames_10NoSL_AL, filenames_10_NoSL_AL):
#     plt.plot(dataFrame['Relative error'].values,'-o',
#               linewidth=width,
#               label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Aeroelastic iteration',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.grid()
# plt.show()


# plt.figure(2)
# plt.title('Aeroelastic convergence with 20 structure nodes',fontsize=titleSize)
# for dataFrame, filename in zip(dataFrames_20NoSL_AL, filenames_20_NoSL_AL):
#     plt.plot(dataFrame['Relative error'].values,'-o',
#               linewidth=width,
#               label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Aeroelastic iteration',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.grid()
# plt.show()


# plt.figure(3)
# plt.title('Aeroelastic convergence with 40 structure nodes',fontsize=titleSize)
# for dataFrame, filename in zip(dataFrames_40NoSL_AL, filenames_40NoSL_AL):
#     plt.plot(dataFrame['Relative error'].values,'-o',
#               linewidth=width,
#               label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Aeroelastic iteration',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.grid()
# plt.show()

# plt.figure(3)
# plt.title('Aeroelastic convergence with 80 structure nodes',fontsize=titleSize)
# for dataFrame, filename in zip(dataFrames_80NoSL_AL, filenames_80NoSL_AL):
#     plt.plot(dataFrame['Relative error'].values,'-o',
#               linewidth=width,
#               label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Aeroelastic iteration',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.grid()
# plt.show()

plt.figure(4)
plt.title('Aeroelastic varying strcuture nodes from 10 to 80',fontsize=titleSize)
for dataFrame, filename in zip(dataFrames_NoSL_10x80AL, filenames_NoSL_10x80AL):
    plt.plot(dataFrame['Relative error'].values,'-o',
              linewidth=width,
              label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.xlabel('Aeroelastic iteration',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.grid()
plt.show()
