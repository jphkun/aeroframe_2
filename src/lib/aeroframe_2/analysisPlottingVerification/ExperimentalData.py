#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 13:17:55 2020

@author: Jean-Philippe Kuntzer

TODO: Balance norm
TODO: Simulate with pytornado
TODO: Compute CL, pytornado and SU2
TODO: 

"""

import numpy as np
import numpy.linalg as LA
import pandas as pd
import matplotlib.pyplot as plt
# import skspatial.objects as sk
from mpl_toolkits.mplot3d import Axes3D
import sys


pathBalance = '../../../../test/ExperimentalData/Balance/'
nameBalance = 's100a10p.xlsx'

pathOptiTrack = '../../../../test/ExperimentalData/OptiTrack/'
nameOptiTrack = 's100a10p.csv'

# Reads balance data
b_header = ['Time','ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8','Fx','Fy','Fz','Mx','My','Mz']
balance = pd.read_excel(pathBalance + nameBalance, skiprows=49)
balance.columns = b_header

# Reads OptiTrack data
oT_header = ['Frame','Time',
             'Rx','Ry','Rz','Px','Py','Pz',
             'Px11','Py11','Pz11',
             'Px12','Py12','Pz12',
             'Px13','Py13','Pz13',
             'Px14','Py14','Pz14',
             'Px21','Py21','Pz21',
             'Px22','Py22','Pz22',
             'Px23','Py23','Pz23',
             'Px24','Py24','Pz24']
optiTrack = pd.read_csv(pathOptiTrack + nameOptiTrack,skiprows=7,usecols=range(32))
optiTrack.columns = oT_header
optiTrack.dropna

# Processses OptiTrackData
vectorsFuselage = []
vectorsTip = []
angle = []



for x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4 in zip(optiTrack['Px11'],optiTrack['Py11'],optiTrack['Pz11'],
                    optiTrack['Px12'],optiTrack['Py12'],optiTrack['Pz12'],
                    optiTrack['Px13'],optiTrack['Py13'],optiTrack['Pz13'],
                    optiTrack['Px14'],optiTrack['Py14'],optiTrack['Pz14']):
    # Fuselage nose
    point2 = np.array([x2,y2,z2])#.round(4)
    # Fuselage tail
    point1 = np.array([x1,y1,z1])#.round(4)
    # wing tip leading edge
    point3 = np.array([x3,y3,z3])#.round(4)
    # wing tip trailing edge
    point4 = np.array([x4,y4,z4])#.round(4)
    # Fuselage vector
    v1 = point1-point2
    # Wing tip vector
    v2 = point4-point3
    
    vectorsFuselage.append(v1)
    vectorsTip.append(v2)
    
    dot = np.dot(v1,v2)
    n = np.linalg.norm(v1)
    v = np.linalg.norm(v2)
    angle.append(np.arccos(dot/(n*v)))
avg = 10   
optiTrack['angle'] = angle
optiTrack['angle'] = optiTrack['angle'].rolling(window=avg).mean() - optiTrack['angle'][0:avg*1].sum()/avg*1
# print(optiTrack['angle'])


norm1 = optiTrack[['Px11','Py11','Pz11']].rolling(window=avg).mean() - optiTrack[['Px11','Py11','Pz11']][0:avg].sum()/avg
norm2 = optiTrack[['Px12','Py12','Pz12']].rolling(window=avg).mean() - optiTrack[['Px12','Py12','Pz12']][0:avg].sum()/avg
norm3 = optiTrack[['Px13','Py13','Pz13']].rolling(window=avg).mean() - optiTrack[['Px13','Py13','Pz13']][0:avg].sum()/avg
norm4 = optiTrack[['Px14','Py14','Pz14']].rolling(window=avg).mean() - optiTrack[['Px14','Py14','Pz14']][0:avg].sum()/avg
optiTrack['normP1'] = LA.norm(norm1,axis=1)
optiTrack['normP2'] = LA.norm(norm2,axis=1)
optiTrack['normP3'] = LA.norm(norm3,axis=1)
optiTrack['normP4'] = LA.norm(norm4,axis=1)

# Computes wing properties
rootChord = np.linalg.norm(point2-point1)
wingtipChord = np.linalg.norm(point4-point3)
semiSpan = point3[0] - point1[0]
semiSpan = point4[0] - point1[0]
fusellageCenter = 0.5*(point1 + point2)
wingTipCenter = 0.5*(point3 + point4)
print('root chord:      ' + str(rootChord))
print('wing tip chord:  ' + str(wingtipChord))
print('semi-span:       ' + str(semiSpan))
print('fuselage center: ' + str(fusellageCenter))
print('wing tip center: ' + str(wingTipCenter))

# Graph properties
avg = avg  # graph moving average
titleSize = 35
textSize = 25
width = 5




# Balance data
# print(balance)
plt.figure('Balance forces')
plt.title('Balance forces',fontsize=titleSize)
plt.plot(balance['Time'].rolling(window=avg).mean(), balance['Fx'].rolling(window=avg).mean(), '-' , label='Fx', linewidth=width)
plt.plot(balance['Time'].rolling(window=avg).mean(), balance['Fy'].rolling(window=avg).mean(), 'r-', label='Fy', linewidth=width)
plt.plot(balance['Time'].rolling(window=avg).mean(), balance['Fz'].rolling(window=avg).mean(), '-' , label='Fz', linewidth=width)
plt.legend()
plt.xlabel('Time [s]',fontsize=textSize)
plt.ylabel('Force [N]',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()
# plt.show()

plt.figure('Balance moments')
plt.title('Balance moments',fontsize=titleSize)
plt.plot(balance['Time'].rolling(window=avg).mean(), balance['Mx'].rolling(window=avg).mean(), '-' , label='Mx', linewidth=width)
plt.plot(balance['Time'].rolling(window=avg).mean(), balance['My'].rolling(window=avg).mean(), 'r-', label='My', linewidth=width)
plt.plot(balance['Time'].rolling(window=avg).mean(), balance['Mz'].rolling(window=avg).mean(), '-' , label='Mz', linewidth=width)
plt.legend()
plt.xlabel('Time [s]',fontsize=textSize)
plt.ylabel('Force [Nm]',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()
# plt.show()


# Actual Z axis in aeroframe_2
plt.figure('OptiTrack displacements norm')
plt.title('Vertical displacement',fontsize=titleSize)
plt.plot(optiTrack['Time'].rolling(window=avg).mean(), optiTrack['normP1'], label='P1', linewidth=width)
plt.plot(optiTrack['Time'].rolling(window=avg).mean(), optiTrack['normP2'], label='P2', linewidth=width)
plt.plot(optiTrack['Time'].rolling(window=avg).mean(), optiTrack['normP3'], label='P3', linewidth=width)
plt.plot(optiTrack['Time'].rolling(window=avg).mean(), optiTrack['normP4'], label='P4', linewidth=width)
plt.legend()
plt.xlabel('Time [s]',fontsize=textSize)
plt.ylabel('Dismplacement [m]',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.grid()


# Optitrack Data
angleRelative = optiTrack['angle'].rolling(window=avg).mean() - optiTrack['Px11'][0:avg].sum()/avg
angleRelative = np.rad2deg(angleRelative)

plt.figure('OptiTrack twist angle')
plt.title('Angle between fuselage and wingtip',fontsize=titleSize)
plt.plot(optiTrack['Time'],optiTrack['angle'],label='angle', linewidth=width)
plt.xlabel('Time [s]',fontsize=textSize)
plt.ylabel('angle [deg]',fontsize=textSize)
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.ylim(-1e-1,1e-1)
plt.grid()



# var = 10
# lim = 1
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(optiTrack['Px11'][var],optiTrack['Py11'][var],optiTrack['Pz11'][var])
# ax.scatter(optiTrack['Px12'][var],optiTrack['Py12'][var],optiTrack['Pz12'][var])
# ax.scatter(optiTrack['Px13'][var],optiTrack['Py13'][var],optiTrack['Pz13'][var])
# ax.scatter(optiTrack['Px14'][var],optiTrack['Py14'][var],optiTrack['Pz14'][var])
# ax.quiver(optiTrack['Px12'][var],optiTrack['Py12'][var],optiTrack['Pz12'][var],
#           vectorsFuselage[var][0],vectorsFuselage[var][1],vectorsFuselage[var][2])
# ax.quiver(optiTrack['Px13'][var],optiTrack['Py13'][var],optiTrack['Pz13'][var],
#           vectorsTip[var][0],vectorsTip[var][1],vectorsTip[var][2])
# ax.set_xlabel('X Label')
# ax.set_ylabel('Y Label')
# ax.set_zlabel('Z Label')
# ax.set_xlim(0,lim)
# ax.set_ylim(0,lim)
# ax.set_zlim(0,lim)

plt.show()