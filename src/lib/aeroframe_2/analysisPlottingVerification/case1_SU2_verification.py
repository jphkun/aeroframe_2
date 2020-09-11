#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 10:26:54 2020

@author: Jean-Philippe Kuntzer
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os

# Choices: SU2, Pytornado
# ssh cfs11
solver = 'SU2'
dy = 0.0625*10
if solver == 'SU2':
    cwd1 = '../../../../test/static/4_WingValidationSU2FramAT_case1/CSD/results/'
    cwd2 = '../../../../test/static/4_WingValidationSU2FramAT_case1/'
    filename4 = 'case1.json'
    l = 5

else:

    # cwd1 = '../../../../test/static/0_WingValidationPytornadoFramAT_case0/CFD/_results/'
    # cwd2 = '../../../../test/static/0_WingValidationPytornadoFramAT_case0/'
    # filename4 = 'case0.json'
    # l = 5

    # cwd1 = '../../../../test/static/1_WingValidationPytornadoFramAT_case1/CFD/_results/'
    # cwd2 = '../../../../test/static/1_WingValidationPytornadoFramAT_case1/'
    # filename4 = 'case1.json'
    # l = 5

    cwd1 = '../../../../test/static/2_WingValidationPytornadoFramAT_case2/CFD/_results/'
    cwd2 = '../../../../test/static/2_WingValidationPytornadoFramAT_case2/'
    filename4 = 'case2.json'
    l = 5

    # cwd1 = '../../../../test/static/7_wingWindTunnelPytornado/CFD/_results/'
    # cwd2 = '../../../../test/static/7_wingWindTunnelPytornado/'
    # filename4 = 'case5.json'
    # l = 0.5

# /home/cfse2/Documents/aeroframe_2/test/static/1_WingValidationPytornado_case1/
if solver == 'SU2':
    filename1 = 'force0.csv'
else:
    filename1 = 'forces0.csv'

filename2 = 'FEM_frocesAndMoments0.csv'
filename3 = 'FEM_displacementAndRotations0.csv'

if solver == 'SU2':
    filename5 = 'force1.csv'
    fx = 'fx'
    fy = 'fy'
    fz = 'fz'
else:
    filename5 = 'forces1.csv'
    fx = 'Fx'
    fy = 'Fy'
    fz = 'Fz'


# distance between the elastic center and the leading edge
print(os.getcwd())
try:
    with open(cwd2+filename4) as f:
        settings = json.load(f)
    f.close()
except FileNotFoundError:
    print('error : \n'+cwd2+filename4)

e = settings['wing1']['elasticAxis']


# Computes the analytical response

E = settings['wing1']['materialProperties']['E']
G = settings['wing1']['materialProperties']['G']
rho = settings['wing1']['materialProperties']['rho']
A = settings['wing1']['mechanicalProperties']['A']
Iy = settings['wing1']['mechanicalProperties']['Iy']
Iz = settings['wing1']['mechanicalProperties']['Iz']
J = settings['wing1']['mechanicalProperties']['J']
d = 2 * settings['wing1']['massAxis']

mass = 2*l*A*rho
linearMass = 1*A*rho
print('A:    ' + str(A))
print('E:    ' + str(E))
print('G:    ' + str(G))
print('Iy:   ' + str(Iy))
print('Iz:   ' + str(Iz))
print('J:    ' + str(J))
print('mass: ' + str(mass))
print('lmass:' + str(linearMass))

# Imports data
df1 = pd.read_csv(cwd1 + filename1)
# Iteration 0 CSD forces
df2 = pd.read_csv(cwd1 + filename2)
# Iteration 0 CSD displacements
df3 = pd.read_csv(cwd1 + filename3)
# Iteration 1 CFD forces
df5 = pd.read_csv(cwd1 + filename5)


# Lumps all the froces in that have the same y coordinate
# Iteration 0 CFD forces

dy2 = np.abs(df2['y'][2] - df2['y'][3])
N = int(2*l/dy)
upperLimit = l

# CFDmoment = 12127.657819960401
# MAPPEDmoment = np.sum(df2['My'])
# print(round(100*(CFDmoment - MAPPEDmoment)/CFDmoment ,1))
# sys.exit()

# Computes the aerodynamic moment for each point
print(df1)

df1['My'] = (-df1['x'] + (e-0.5)*2) * df1[fz]

print((-df1['x'] + (e-0.5)*2))
print(df1['My'])


df_1 = pd.DataFrame()
for i in range(N):
    df = df1[(df1['y'] <= upperLimit-i*dy) & (df1['y'] > upperLimit-(i+1)*dy)].sum(axis=0)
    n = len(df)
    df['x'] = 0
    df['y'] = 0.5 * ((upperLimit-i*dy) + upperLimit-(i+1)*dy)
    df['z'] = 0
    df[fx] = df[fx] / dy
    df[fy] = df[fy] / dy
    df[fz] = df[fz] / dy
    df['My'] = df['My']
    # print(df['My'])
    if solver == 'SU2':
        del df['ids']
        del df['marker']
    df_1 = df_1.append(df, ignore_index=True)

# Iteration 1 CFD forces
df_5 = pd.DataFrame()
for i in range(N):
    df = df1[(df1['y'] <= upperLimit-i*dy) & (df1['y'] > upperLimit-(i+1)*dy)].sum(axis=0)
    n = len(df)
    df['x'] = 0
    df['y'] = 0.5 * ((upperLimit-i*dy) + upperLimit-(i+1)*dy)
    df['z'] = 0
    df[fx] = df[fx] / dy
    df[fy] = df[fy] / dy
    df[fz] = df[fz] / dy
    df['My'] = df['My']
    if solver == 'SU2':
        del df['ids']
        del df['marker']
    df_5 = df_1.append(df, ignore_index=True)

print('\n')
print(np.sum(df2['My']))
print(np.sum(df1['My']))
print(np.sum(df_1['My']))
# sys.exit()
# Converts data to the correct format for polynomial fit and plotting
Fz_cfd0 = df_1[fz].values
Fz_cfd1 = df_5[fz].values
My_cfd0 = df_1['My'].values

# beam moment
Mb_cfd0 = My_cfd0

# Polynomial fitting of both curves
cFz, rFz, _, _, _ = np.polyfit(df_1['y'], df_1[fz], 7, full=True)
cMb, rMb, _, _, _ = np.polyfit(df_1['y'], df_1['My']/dy, 7, full=True)
pFz = np.poly1d(cFz)
print('Lift distribution function:')
print(pFz)
pMb = np.poly1d(cMb)
print('\n')
print('Moment distribution function:')
print(pMb)
cFz = np.flip(cFz)

y_pol = np.linspace(-l,l,100)
Fz_pol = pFz(y_pol)
Mb_pol = pMb(y_pol)


n = 0

c11 = -((cFz[6]*l**7)/7 +  (cFz[4]*l**5)/5  + (cFz[2]*l**3)/3 +   (cFz[0] - n*linearMass*9.81)*l**1)
c12 = -((cFz[6]*l**8)/56 + (cFz[4]*l**6)/30 + (cFz[2]*l**4)/12 + ((cFz[0] - n*linearMass*9.81)*l**2)/2 +c11*l)
def w(y):
    w = (1/(E*Iy)) * \
        (((cFz[0] - n*linearMass*9.81)*y**4) / 24 + \
        (cFz[2]*y**6) / 360 + \
        (cFz[4]*y**8) / 1680 + \
        (cFz[6]*y**10)/ 5040 + \
        (c11*y**3)     / 6 +\
        (c12*y**2)     / 2)
    return w


# def f(x):
#     w = cFz[7]*x**7 + cFz[6]*x**6 + cFz[5]*x*5 + cFz[4]*x**4 + cFz[3]*x**3 + cFz[2]*x**2 +cFz[1]*x + cFz[0]
#     return w


def t(y):
    yMax = l
    t = (1/(G*J)) * \
        ((((cMb[7] + d*n*linearMass*9.81)*yMax**1) / 1 +  \
          (cMb[5]*yMax**3) / 3 + \
          (cMb[3]*yMax**5) / 5 + \
          (cMb[1]*yMax**7) / 7 )*y -
        (((cMb[7] + d*n*linearMass*9.81)*y**2) / 2 +  \
          (cMb[5]*y**4) / 12 + \
          (cMb[3]*y**6) / 30 + \
          (cMb[1]*y**8) / 56 ))
    return t


# Graph properties
titleSize = 35
textSize = 25
width = 5

y_an = np.linspace(-l,l,100)


# Lift distribtion
plt.figure(7)
plt.plot(df_1['y'],df_1[fz],'-',label='CFD value',      linewidth=width)
plt.plot(y_pol,    Fz_pol,    'r-',label='Polynomial fit', linewidth=width)
plt.plot(df2['y'],df2['Fz']/dy2,'-',label='FEM value', linewidth=width)
plt.title('Lift distribution',fontsize=titleSize)
plt.xlabel('position [m]',fontsize=textSize)
plt.ylabel('Lift [N/m]',fontsize=textSize)
plt.xlim(-l,l)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)

# Vertical displacement
plt.figure(8)
plt.plot(df3['y'],df3['dz'],'-r',label='FEM', linewidth=width)
plt.plot(y_an,w(np.abs(y_an)),'-',label='Analytical solution', linewidth=width)
errorMax = 100*(np.max(np.abs(df3['dz'])) - np.max(np.abs(w(np.abs(y_an))))) / np.max(np.abs(w(np.abs(y_an))))
errorMax = round(errorMax,1)
# print(df3['dz'])
plt.title('Displacement \n Error max: '+str(errorMax)+'%',fontsize=titleSize)
plt.xlabel('position [m]',fontsize=textSize)
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)


# Torque distribution
plt.figure(9)
# Plots
plt.plot(df_1['y'],df_1['My']/dy,'-',label='CFD', linewidth=width)
plt.plot(y_pol,Mb_pol,'r-',label='Polynomial fitting', linewidth=width)
plt.plot(df2['y'],df2['My']/dy2,'-',label='FEM value', linewidth=width)
plt.xlabel('position [m]',fontsize=textSize)
plt.ylabel('Torque [Nm/m]',fontsize=textSize)
plt.title('Torque distribution',fontsize=titleSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)

# Rotation angle
plt.figure(10)
plt.plot(df3['y'],np.rad2deg(np.abs(df3['tx'])),'-r',label='FEM', linewidth=width)
plt.plot(y_an,np.rad2deg(t(np.abs(y_an))),'-',label='Analytical solution', linewidth=width)
errorMax2 = 100*(np.max(np.abs(df3['tx'])) - np.max(np.abs(t(np.abs(y_an))))) / np.max(np.abs(t(np.abs(y_an))))
errorMax2 = round(errorMax2, 1)
# print(errorMax2)
# np.set_printoptions(precision=3)
plt.title('Rotation angle \n Max error: '+str(errorMax2)+'%',fontsize=titleSize)
plt.xlabel('position [m]',fontsize=textSize)
plt.ylabel('Rotation angle [deg]',fontsize=textSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.show()

# # a_0
# plt.figure(11)
# # plt.plot(df_1['y'],Fz_cfd0,'r-',label='L cfd', linewidth=width)
# plt.plot(df_1['y'],Fz_cfd0,'o-',
#           label='iteration 0 $c_{l}$ cfd',
#           linewidth=width)
# theta = t(np.abs(df_1['y']))
# plt.plot(df_1['y'], Fz_cfd0*(1+theta/np.deg2rad(10)),'r-',
#           label='Iteration 1 $c_{l}$ analytical',
#           linewidth=width)

# # plt.plot(df_1['y'],Fz_cfd1,'go-',
# #           label='Iteration 1 $c_{l}$ cfd',
# #           linewidth=width)
# # errMax3 = 100*(Fz_cfd0*(1+theta/np.deg2rad(10)) - Fz_cfd1)/(Fz_cfd0*(1+theta/np.deg2rad(10)))
# # errMax3 = errMax3.round(decimals=1)
# # errMax3 = np.max(errMax3)
# # plt.title('Lift distribution at step 1\n Max error: '+str(errMax3)+'%',fontsize=titleSize)
# plt.xlabel('y position [m]',fontsize=textSize)
# plt.ylabel('Lift [N/m]',fontsize=textSize)
# plt.xlim(-5,5)
# plt.grid()
# plt.legend(fontsize=textSize)
# plt.xticks(fontsize=textSize)
# plt.yticks(fontsize=textSize)





plt.show()


error =   np.array([44.3,10.8,0.7,-2.7,-4.0,-4.5,-4.7])
Npoints = np.array([5,   10,  20, 40,  80,  160, 320])