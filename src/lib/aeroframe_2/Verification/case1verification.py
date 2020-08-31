#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 08:15:28 2020

@author: Jean-Philippe Kuntzer
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import json

cwd1 = '/home/cfse2/Documents/aeroframe_2/test/static/1_WingValidationPytornado_case1//CFD/_results/'
cwd2 = '/home/cfse2/Documents/aeroframe_2/test/static/1_WingValidationPytornado_case1//'
# /home/cfse2/Documents/aeroframe_2/test/static/1_WingValidationPytornado_case1/
filename1 = 'forces0.csv'
filename2 = 'FEM_frocesAndMoments0.csv'
filename3 = 'FEM_displacementAndRotations0.csv'
filename4 = 'case1.json'
filename5 = 'forces1.csv'
# distance between the elastic center and the leading edge

try:
    with open(cwd2+filename4) as f:
        settings = json.load(f)
    f.close()
except FileNotFoundError:
    print('error : \n'+cwd2+filename4)

a = settings['wing1']['elasticAxis']

# Imports data
df1 = pd.read_csv(cwd1 + filename1)
df2 = pd.read_csv(cwd1 + filename2)
df3 = pd.read_csv(cwd1 + filename3)
df5 = pd.read_csv(cwd1 + filename5)
dfF0 = df1.groupby(['y'])['y','Fz'].agg('sum')
aggregation_functions = {'y': 'first', 'Fy': 'sum'}
df5['y'] = df5['y'].round(decimals=2)
df_new = df5.groupby(df5['y']).aggregate(aggregation_functions)
dfF1 = df5.groupby(['y'])['y','Fz'].sum()

# print(dfF0)
# print(dfF1)
# sys.exit()

# Computes the moment
df1['x_my'] = df1['x']
# print(df1['x'])
# print(df1['x_my'])
df1['My'] = (-df1['x']+(a-0.5)*2)*df1['Fz']
# print(-df1['x'])
# print(-df1['x']+(a-0.5)*2)
# sys.exit()
dfM = df1.groupby(['y'])['y','My'].agg('sum')
# print(df2)

# Converts data to the correct format for polynomial fit and plotting
Fz_cfd0 = dfF0.values
Fz_cfd1 = dfF1.values
My_cfd = dfM.values
y_cfd = Fz_cfd0[:,0]/5

# Computes the dy to make the force and moment idependent of the mesh
dy1 = np.abs(y_cfd[1]-y_cfd[0])
dy2 = np.abs(df2['y'][1]-df2['y'][0])
# print(dy1)

Fz_cfd0 = Fz_cfd0[:,1]/dy1
Fz_cfd1 = Fz_cfd1[:,1]/dy1
# print(Fz_cfd0)
# print(Fz_cfd1)
# sys.exit()
My_cfd = My_cfd[:,1]/dy1
# beam moment
Mb_cfd = My_cfd

# Polynomial fitting of both curves
cFz, rFz, _, _, _ = np.polyfit(y_cfd, Fz_cfd0, 7, full=True)
cMb, rMb, _, _, _ = np.polyfit(y_cfd, Mb_cfd, 7, full=True)
pFz = np.poly1d(cFz)
print(pFz)
pMb = np.poly1d(cMb)
print(pMb)
cFz = np.flip(cFz)


y_pol = np.linspace(-5,5,100)
Fz_pol = pFz(y_pol)
Mb_pol = pMb(y_pol)

# Computes the analytical response
E = settings['wing1']['materialProperties']['E']
G = settings['wing1']['materialProperties']['G']
Iy = settings['wing1']['mechanicalProperties']['Iy']
Iz = settings['wing1']['mechanicalProperties']['Iz']
l = 5

c1 = -((cFz[6]*l**7)/7 +  (cFz[4]*l**5)/5  + (cFz[2]*l**3)/3 +  (cFz[0]*l**1)/1)
c2 = -((cFz[6]*l**8)/56 + (cFz[4]*l**6)/30 + (cFz[2]*l**4)/12 + (cFz[0]*l**2)/2 +c1*l)
# print(cFz)

def w(y):
    w = (1/(E*Iy)) * \
        ((cFz[0]*y**4) / 24 + \
          (cFz[2]*y**6) / 360 + \
          (cFz[4]*y**8) / 1680 + \
          (cFz[6]*y**10)/ 5040 + \
          (c1*y**3)     / 6 +\
          (c2*y**2)     / 2)
    return w

def f(x):
    # w = cFz[7]*y**0 + cFz[6]*y**1 + cFz[5]*y*2 + cFz[4]*y**3 + cFz[3]*y**4 + cFz[2]*y**5 + cFz[1]*y**6 + cFz[0]*y**7
    w = cFz[7]*x**7 + cFz[6]*x**6 + cFz[5]*x*5 + cFz[4]*x**4 + cFz[3]*x**3 + cFz[2]*x**2 +cFz[1]*x + cFz[0]
    return w

# print(cMb[1])
# print(cMb[3])
# print(cMb[5])
# print(cMb[6])
def t(y):
    yMax = 5
    t = (1/(G*Iz)) * \
        (((cMb[7]*yMax**1) / 1 +  \
          (cMb[5]*yMax**3) / 3 + \
          (cMb[3]*yMax**5) / 5 + \
          (cMb[1]*yMax**7) / 7 )*y -
        ((cMb[7]*y**2) / 2 +  \
          (cMb[5]*y**4) / 12 + \
          (cMb[3]*y**6) / 30 + \
          (cMb[1]*y**8) / 56 ))
    # t = (1/(G*Iz)) * \
    #     ((cMb[0]*y**2) / 2 +  \
    #      (cMb[2]*y**4) / 12 + \
    #      (cMb[4]*y**6) / 30 + \
    #      (cMb[6]*y**8) / 56 )
    return t

def a0(L):
    a = L/213.802
    return a

def a1(L,theta):
    a = L/(6125*(np.deg2rad(2) + theta))
    return a

# Graph properties
titleSize = 30
textSize = 20
width = 3

# Lift distribtion
plt.figure(2)
plt.plot(y_pol,Fz_pol,'r-',label='Polynomial fit', linewidth=width)
plt.plot(y_cfd,Fz_cfd0,'o-',label='CFD value', linewidth=width)
plt.plot(df2['y'],df2['Fz']/dy2,'o-',label='FEM value', linewidth=width)
# plt.plot(y_cfd,f(y_cfd),'-',label='FEM value', linewidth=width)
plt.title('Lift distribution',fontsize=titleSize)
plt.xlabel('y position [m]',fontsize=textSize)
plt.ylabel('Lift [N/m]',fontsize=textSize)
plt.xlim(-5,5)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)

# Vertical displacement
plt.figure(3)
plt.plot(df3['y'],df3['dz'],'-ro',label='FEM', linewidth=width)
plt.plot(y_cfd,w(np.abs(y_cfd)),'-',label='Analytical solution', linewidth=width)
errorMax = 100*(np.max(df3['dz']) - np.max(w(np.abs(y_cfd)))) / np.max(w(np.abs(y_cfd)))
errorMax = errorMax.round(decimals=1)
plt.title('Displacement \n Error max: '+str(errorMax)+'%',fontsize=titleSize)
plt.xlabel('postion [m]',fontsize=textSize)
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)


# print(df2)
# print(df3)

# Torque distribution
plt.figure(4)
# Plots
plt.plot(y_cfd,My_cfd,'o',label='CFD', linewidth=width)
plt.plot(y_pol,Mb_pol,'r-',label='Polynomial fitting', linewidth=width)
plt.plot(df2['y'],df2['My']/dy2,'o-',label='FEM value', linewidth=width)
plt.xlabel('postion [m]',fontsize=textSize)
plt.ylabel('Troque [Nm/m]',fontsize=textSize)
plt.title('Torque distribution',fontsize=titleSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)


# Rotation angle
plt.figure(5)
plt.plot(df3['y'],df3['tx'],'-ro',label='FEM', linewidth=width)
plt.plot(y_cfd,t(np.abs(y_cfd)),'-',label='Analytical solution', linewidth=width)
errorMax2 = 100*(np.max(df3['tx']) - np.max(t(np.abs(y_cfd)))) / np.max(t(np.abs(y_cfd)))
errorMax2 = errorMax2.round(decimals=1)
print(errorMax2)
# np.set_printoptions(precision=3)
plt.title('Rotation angle \n Max error: '+str(errorMax2)+'%',fontsize=titleSize)
plt.xlabel('positon [m]',fontsize=textSize)
plt.ylabel('Rotation angle [rad]',fontsize=textSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
# plt.show()

# a_0
plt.figure(6)
# plt.plot(y_cfd,Fz_cfd0,'r-',label='L cfd', linewidth=width)
plt.plot(y_cfd,Fz_cfd0,'o-',
         label='iteration 0 $c_{l}$ cfd',
         linewidth=width)
theta = t(np.abs(y_cfd))
print()
# 
plt.plot(y_cfd,6125*(np.deg2rad(2.0) + theta) * a0(Fz_cfd0),'r-',
         label='Iteration 1 $c_{l}$ analytical',
         linewidth=width)

plt.plot(y_cfd,Fz_cfd1,'go-',
         label='Iteration 1 $c_{l}$ cfd',
         linewidth=width)
errMax3 = (np.max(6125*(np.deg2rad(2.0) + theta) * a0(Fz_cfd0)) - np.max(Fz_cfd1))/np.max(6125*(np.deg2rad(2.0) + theta) * a0(Fz_cfd0))
errMax3 = 100*errMax3.round(decimals=3)
plt.title('Lift distribution at step 1\n Max error: '+str(errMax3)+'%',fontsize=titleSize)
plt.xlabel('y position [m]',fontsize=textSize)
plt.ylabel('$C_{l}$ [1/(m)]',fontsize=textSize)
plt.xlim(-5,5)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)
plt.show()


