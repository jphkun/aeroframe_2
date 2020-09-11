#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 10:08:44 2020

@author: Jean-Philippe Kuntzer
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

cwd1 = '../../../../test/static/0_WingValidationPytornadoFramAT_case0/CFD/_results/'
cwd2 = '../../../../test/static/0_WingValidationPytornadoFramAT_case0/'


# Unloaded beam
filename_1 = '5_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_2 = '10_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_3 = '20_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_4 = '40_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_5 = '80_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_6 = '160_NoSL_NoAL_FEM_displacementAndRotations0.csv'
filename_7 = '320_NoSL_NoAL_FEM_displacementAndRotations0.csv'

# Weighted beam with no aerodynamic forces
filename_8 = '5_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_9 = '10_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_10 = '20_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_11 = '40_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_12 = '80_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_13 = '160_SL_NoAL_FEM_displacementAndRotations0.csv'
filename_14 = '320_SL_NoAL_FEM_displacementAndRotations0.csv'

# Influence of aerodynamic mesh
filename_15 = '10_NoSL_05x20AL_results.csv'
filename_16 = '10_NoSL_05x40AL_results.csv'
filename_17 = '10_NoSL_10x40AL_results.csv'
filename_18 = '10_NoSL_10x80AL_results.csv'

# Influence of aerodynamic mesh
filename_19 = '20_NoSL_05x20AL_results.csv'
filename_20 = '20_NoSL_05x40AL_results.csv'
filename_21 = '20_NoSL_10x40AL_results.csv'
filename_22 = '20_NoSL_10x80AL_results.csv'

# Influence of aerodynamic mesh
filename_23 = '40_NoSL_05x20AL_results.csv'
filename_24 = '40_NoSL_05x40AL_results.csv'
filename_25 = '40_NoSL_10x40AL_results.csv'
filename_26 = '40_NoSL_10x80AL_results.csv'

# Influence of aerodynamic mesh
filename_27 = '80_NoSL_05x20AL_results.csv'
filename_28 = '80_NoSL_05x40AL_results.csv'
filename_29 = '80_NoSL_10x40AL_results.csv'
filename_30 = '80_NoSL_10x80AL_results.csv'


# Beam mesh convergeance
# filename_26 = 's005_results.csv'
# filename_27 = 's010_results.csv'
# filename_28 = 's020_results.csv'
# filename_29 = 's040_results.csv'
# filename_30 = 's080_results.csv'
# filename_31 = 's160_results.csv'
# filename_32 = 's320_results.csv'

# No loads test
filenames_NoSL_NoAL = [filename_1, filename_2, filename_3, filename_4, filename_5, filename_6]
dataFrames_NoSL_NoAL = []
for filename in filenames_NoSL_NoAL:
    df = pd.read_csv(cwd1 + filename)
    dataFrames_NoSL_NoAL.append(df)

# Constant load test
filenames_SL_NoAL = [filename_8, filename_9, filename_10, filename_11, filename_12, filename_13]
dataFrames_SL_NoAL = []
for filename in filenames_SL_NoAL:
    df = pd.read_csv(cwd1 + filename)
    dataFrames_SL_NoAL.append(df)

# Influence of aerodynamic mesh
filenames_10_NoSL_AL = [filename_15, filename_16, filename_17, filename_18]
dataFrames_10NoSL_AL = []
for filename in filenames_10_NoSL_AL:
    df1 = pd.read_csv(cwd2 + filename,sep=';')
    df2 = df1.iloc[0:1,:]
    df2['Relative error'] = 0
    df = pd.concat([df2,df1])
    dataFrames_10NoSL_AL.append(df)

# Influence of aerodynamic mesh
filenames_20_NoSL_AL = [filename_19, filename_20, filename_21, filename_22]
dataFrames_20NoSL_AL = []
for filename in filenames_20_NoSL_AL:
    df1 = pd.read_csv(cwd2 + filename,sep=';')
    df2 = df1.iloc[0:1,:]
    df2['Relative error'] = 0
    df = pd.concat([df2,df1])
    dataFrames_20NoSL_AL.append(df)

# Influence of aerodynamic mesh
filenames_40NoSL_AL = [filename_23, filename_24, filename_25, filename_26]
dataFrames_40NoSL_AL = []
for filename in filenames_40NoSL_AL:
    df1 = pd.read_csv(cwd2 + filename,sep=';')
    df2 = df1.iloc[0:1,:]
    df2['Relative error'] = 0
    df = pd.concat([df2,df1])
    dataFrames_40NoSL_AL.append(df)

# Influence of aerodynamic mesh
filenames_80NoSL_AL = [filename_27, filename_28, filename_29, filename_30]
dataFrames_80NoSL_AL = []
for filename in filenames_80NoSL_AL:
    df1 = pd.read_csv(cwd2 + filename,sep=';')
    df2 = df1.iloc[0:1,:]
    df2['Relative error'] = 0
    df = pd.concat([df2,df1])
    dataFrames_80NoSL_AL.append(df)

# Influence of strcutre mesh
filenames_NoSL_10x80AL = [filename_18, filename_22, filename_26, filename_30]
dataFrames_NoSL_10x80AL = []
for filename in filenames_NoSL_10x80AL:
    df1 = pd.read_csv(cwd2 + filename,sep=';')
    df2 = df1.iloc[0:1,:]
    df2['Relative error'] = 0
    df = pd.concat([df2,df1])
    dataFrames_NoSL_10x80AL.append(df)


# # Influence of strcutre mesh
# filenames_NoSL_10x80AL = [filename_26, filename_27, filename_28, filename_29, filename_30, filename_31, filename_32]
# dataFrames_NoSL_10x80AL = []
# for filename in filenames_NoSL_10x80AL:
#     df1 = pd.read_csv(cwd2 + filename,sep=';')
#     # print(list(df1.columns))
#     df2 = df1.iloc[0:1,:]
#     df2['Relative error'] = 0
#     df = pd.concat([df2,df1])
#     # print(df1)
#     # print(df2)
#     # print(df)
#     # sys.exit()
#     dataFrames_NoSL_10x80AL.append(df)

# Graph properties
titleSize = 35
textSize = 25
width = 5


# ##############################################################################
# # Convergeance of unloaded beam
# ##############################################################################
# i = 0
# plt.figure(1)
# plt.set_cmap('Set2')
# plt.title('Convergeance of unloaded beam',fontsize=titleSize)
# for dataFrame, filename in zip(dataFrames_NoSL_NoAL, filenames_NoSL_NoAL):
#     plt.plot(dataFrame['y'],dataFrame['dz'],'-o',
#               linewidth=width,label='Nodes: '+str(5*2**i))
#     i += 1
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Wingspan position [m]',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.xticks(fontsize=textSize)
# plt.yticks(fontsize=textSize)


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
# i = 0

# plt.figure(1)
# plt.title('Convergeance of weighted beam',fontsize=titleSize)
# for dataFrame, filename in zip(dataFrames_SL_NoAL, filenames_SL_NoAL):
#     err = (np.min(dataFrame['dz']) - np.min(w(mg,y,L/2,E,I)))
#     err = (np.min(w(mg,y,L/2,E,I)) - np.min(dataFrame['dz']))
#     # print(np.max(w(mg,y,L/2,E,I)))
#     # print(np.max(dataFrame['dz']))
#     err = 100*err / np.min(w(mg,y,L/2,E,I))
#     err = np.abs(err.round(decimals=3))
#     plt.plot(dataFrame['y'],dataFrame['dz'],'-o',
#               linewidth=width,label='Nodes: '+str(5*2**i)+', error:'+str(err)+'%')
#     i += 1
# plt.plot(y,w(mg,y,L/2,E,I),'r',linewidth=width,label='Analytical solution')
# plt.ylabel('Displacement [m]',fontsize=textSize)
# plt.xlabel('Wingspan position [m]',fontsize=textSize)
# plt.legend(fontsize=textSize)
# plt.show()

# #############################################################################
# Convergeance of aerodynamic mesh
# #############################################################################


plt.figure(1)
plt.title('Aeroelastic convergence with 10 structure nodes',fontsize=titleSize)
for dataFrame, filename in zip(dataFrames_10NoSL_AL, filenames_10_NoSL_AL):
    plt.plot(dataFrame['Relative error'].values,'-o',
              linewidth=width,
              label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.xlabel('Aeroelastic iteration',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()



plt.figure(2)
plt.title('Aeroelastic convergence with 20 structure nodes',fontsize=titleSize)
for dataFrame, filename in zip(dataFrames_20NoSL_AL, filenames_20_NoSL_AL):
    plt.plot(dataFrame['Relative error'].values,'-o',
              linewidth=width,
              label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.xlabel('Aeroelastic iteration',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()



plt.figure(3)
plt.title('Aeroelastic convergence with 40 structure nodes',fontsize=titleSize)
for dataFrame, filename in zip(dataFrames_40NoSL_AL, filenames_40NoSL_AL):
    plt.plot(dataFrame['Relative error'].values,'-o',
              linewidth=width,
              label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.xlabel('Aeroelastic iteration',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()


plt.figure(4)
plt.title('Aeroelastic convergence with 80 structure nodes',fontsize=titleSize)
for dataFrame, filename in zip(dataFrames_80NoSL_AL, filenames_80NoSL_AL):
    plt.plot(dataFrame['Relative error'].values,'-o',
              linewidth=width,
              label='Structure nodes: '+filename[:2]+'\n VLM nodes:'+filename[8:13])
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.xlabel('Aeroelastic iteration',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()


i = 0
plt.figure(5)
for dataFrame, filename in zip(dataFrames_NoSL_10x80AL, filenames_NoSL_10x80AL):
    plt.plot(dataFrame['Relative error'].values,'-o',
              linewidth=width,
              label='Structure nodes: '+str(10*2**i))
    i += 1
plt.title('Aeroelastic final value analysis',fontsize=titleSize)
plt.xlabel('Iteration',fontsize=textSize)
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)


# error =   np.array([44.3,10.8,0.7,-2.7,-4.0,-4.5,-4.7])
# Npoints = np.array([5,   10,  20, 40,  80,  160, 320])
# plt.figure('First step error with mesh size')
# plt.semilogx(Npoints,error,'-o',
#           linewidth=width,
#           label='Relative error with analytical solution')
# plt.title('Comparing first step error',fontsize=titleSize)
# plt.xlabel('Number of nodes',fontsize=textSize)
# plt.ylabel('Error [%]',fontsize=textSize)
# plt.grid()
# plt.legend(fontsize=textSize)
# plt.xticks(fontsize=textSize)
# plt.yticks(fontsize=textSize)



plt.show()
