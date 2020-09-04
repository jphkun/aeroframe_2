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

cwd1 = '/home/cfse2/Documents/aeroframe_2/test/static/2_WingValidationPytornadoFramAT_case2/CFD/_results/'
cwd2 = '/home/cfse2/Documents/aeroframe_2/test/static/2_WingValidationPytornadoFramAT_case2/'
# /home/cfse2/Documents/aeroframe_2/test/static/1_WingValidationPytornado_case1/
filename1 = 'forces0.csv'
filename2 = 'FEM_frocesAndMoments0.csv'
filename3 = 'FEM_displacementAndRotations0.csv'
filename4 = 'case2.json'
filename5 = 'forces1.csv'
# distance between the elastic center and the leading edge

try:
    with open(cwd2+filename4) as f:
        settings = json.load(f)
    f.close()
except FileNotFoundError:
    print('error : \n'+cwd2+filename4)

a = settings['wing1']['elasticAxis']


# Computes the analytical response

E = settings['wing1']['materialProperties']['E']
G = settings['wing1']['materialProperties']['G']
rho = settings['wing1']['materialProperties']['rho']
A = settings['wing1']['mechanicalProperties']['A']
Iy = settings['wing1']['mechanicalProperties']['Iy']
Iz = settings['wing1']['mechanicalProperties']['Iz']
J = settings['wing1']['mechanicalProperties']['J']
d = 2 * settings['wing1']['massAxis']
l = 5
mass = 2*l*A*rho
linearMass = A*rho
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
dfF0 = df1.groupby(['y'])['y','Fz'].agg('sum')
aggregation_functions = {'y': 'first', 'Fy': 'sum'}

# Lumps all the froces in that have the same y coordinate
# Iteration 0 CFD forces
dfF2 = df2.groupby(['y'])['y','Fz'].agg('sum')
aggregation_functions = {'y': 'first', 'Fy': 'sum'}

# Iteration 1 CFD forces
df5['y'] = df1['y']
df_new = df5.groupby(df5['y']).aggregate(aggregation_functions)
dfF1 = df5.groupby(['y'])['y','Fz'].sum()



# print(dfF0)
# print(df5)
# sys.exit()

# Computes the moment
df1['x_my'] = df1['x']
# aerodynamic moment
df1['My'] = (-df1['x'] + (a-0.5)*2) * df1['Fz']
dfM = df1.groupby(['y'])['y','My'].agg('sum')
# print(dfM)

# Converts data to the correct format for polynomial fit and plotting
Fz_cfd0 = dfF0.values
Fz_cfd1 = dfF1.values
Fz_csd2 = dfF2.values
My_cfd = dfM.values
# Number of chordwise nodes
N = 10
y_cfd = Fz_cfd0[:,0]/N
y_csd = Fz_csd2[:,0]

# Computes the dy to make the force and moment idependent of the mesh
dy1 = np.abs(y_cfd[1]-y_cfd[0])
dy2 = np.abs(df2['y'][1]-df2['y'][0])
dy3 = np.abs(y_csd[1]-y_csd[0])
# print(dy1)

Fz_cfd0 = Fz_cfd0[:,1]/dy1
Fz_cfd1 = Fz_cfd1[:,1]/dy1
Fz_csd2 = Fz_csd2[:,1]/dy3
My_cfd = My_cfd[:,1]/dy1
# beam moment
Mb_cfd = My_cfd
# print(Fz_cfd0)
# print(Fz_csd2)
# sys.exit()
# Polynomial fitting of both curves
cFz,  rFz,  _, _, _ = np.polyfit(y_cfd, Fz_cfd0, 7, full=True)
cFz2, rFz2, _, _, _ = np.polyfit(y_csd, Fz_csd2, 7, full=True)
cMb,  rMb,  _, _, _ = np.polyfit(y_cfd, Mb_cfd, 7, full=True)
pFz = np.poly1d(cFz)
pFz2 = np.poly1d(cFz2)
print('Lift distribution function:')
print(pFz)
pMb = np.poly1d(cMb)
print('\n')
print('Moment distribution function:')
print(pMb)
cFz = np.flip(cFz)
cFz2 = np.flip(cFz2)


y_pol = np.linspace(-5,5,100)
Fz_pol = pFz(y_pol)
Fz2_pol = pFz2(y_pol)
Mb_pol = pMb(y_pol)


c11 = -((cFz[6]*l**7)/7 +  (cFz[4]*l**5)/5  + (cFz[2]*l**3)/3 +   (cFz[0] - linearMass*9.81)*l**1)
c12 = -((cFz[6]*l**8)/56 + (cFz[4]*l**6)/30 + (cFz[2]*l**4)/12 + ((cFz[0] - linearMass*9.81)*l**2)/2 +c11*l)
def w1(y):
    w = (1/(E*Iy)) * \
        (((cFz[0] - linearMass*9.81)*y**4) / 24 + \
        (cFz[2]*y**6) / 360 + \
        (cFz[4]*y**8) / 1680 + \
        (cFz[6]*y**10)/ 5040 + \
        (c11*y**3)     / 6 +\
        (c12*y**2)     / 2)
    return w


c21 = -((cFz2[6]*l**7)/7 + (cFz2[4]*l**5)/5  +(cFz2[2]*l**3)/3 +  (cFz2[0] - linearMass*9.81)*l**1)
c22 = -((cFz2[6]*l**8)/56 +(cFz2[4]*l**6)/30 +(cFz2[2]*l**4)/12 +((cFz2[0] - linearMass*9.81)*l**2)/2 +c21*l)
def w2(y):
    w = (1/(E*Iy)) * \
        (((cFz2[0] - linearMass*9.81)*y**4) / 24 + \
        (cFz2[2]*y**6) / 360 + \
        (cFz2[4]*y**8) / 1680 + \
        (cFz2[6]*y**10)/ 5040 + \
        (c21*y**3)     / 6 +\
        (c22*y**2)     / 2)
    return w

def t(y):
    yMax = 5
    t = (1/(G*J)) * \
        ((((cMb[7] + d*linearMass*9.81)*yMax**1) / 1 +  \
          (cMb[5]*yMax**3) / 3 + \
          (cMb[3]*yMax**5) / 5 + \
          (cMb[1]*yMax**7) / 7 )*y -
        (((cMb[7] + d*linearMass*9.81)*y**2) / 2 +  \
          (cMb[5]*y**4) / 12 + \
          (cMb[3]*y**6) / 30 + \
          (cMb[1]*y**8) / 56 ))
    return t

# Graph properties
titleSize = 30
textSize = 20
width = 3

# Lift distribtion
plt.figure(7)

plt.plot(y_cfd,Fz_cfd0,'o-',label='CFD value', linewidth=width)
plt.plot(y_pol,Fz_pol,'r-',label='Polynomial fit', linewidth=width)
plt.plot(df2['y'],df2['Fz']/dy2,'o-',label='FEM value', linewidth=width)
# plt.plot(y_pol,Fz2_pol,'-',label='Polynomial fit', linewidth=width)
plt.title('Lift distribution',fontsize=titleSize)
plt.xlabel('y position [m]',fontsize=textSize)
plt.ylabel('Lift [N/m]',fontsize=textSize)
plt.xlim(-5,5)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)

# Vertical displacement
plt.figure(8)
plt.plot(df3['y'],df3['dz'],'-ro',label='FEM', linewidth=width)
plt.plot(y_cfd,w1(np.abs(y_cfd)),'-',label='Analytical solution', linewidth=width)
# plt.plot(y_cfd,w2(np.abs(y_cfd)),'-',label='Analytical solution', linewidth=width)
errorMax = 100*(np.max(np.abs(df3['dz'])) - np.max(np.abs(w1(np.abs(y_cfd))))) / np.max(np.abs(w1(np.abs(y_cfd))))
errorMax = errorMax.round(decimals=1)
plt.title('Displacement \n Error max: '+str(errorMax)+'%',fontsize=titleSize)
plt.xlabel('postion [m]',fontsize=textSize)
plt.ylabel('Displacement [m]',fontsize=textSize)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)


# Torque distribution
plt.figure(9)
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

# # Since FEM are corrected this needs to be corrected also
# N = len(df3['tx'])
# for i in range(N):
#     # print()
    
#     # print()
#     if df3['y'][i] < 0:
#         df3.at[str(i+1), 'tx'] = df3['tx'][i]

# Rotation angle
plt.figure(10)
plt.plot(df3['y'],np.abs(df3['tx']),'-ro',label='FEM', linewidth=width)
plt.plot(y_cfd,t(np.abs(y_cfd)),'-',label='Analytical solution', linewidth=width)
errorMax2 = 100*(np.max(np.abs(df3['tx'])) - np.max(t(np.abs(y_cfd)))) / np.max(t(np.abs(y_cfd)))
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
plt.figure(11)
# plt.plot(y_cfd,Fz_cfd0,'r-',label='L cfd', linewidth=width)
plt.plot(y_cfd,Fz_cfd0,'o-',
          label='iteration 0 $c_{l}$ cfd',
          linewidth=width)
theta = t(np.abs(y_cfd))
plt.plot(y_cfd, Fz_cfd0*(1+theta/np.deg2rad(10)),'r-',
          label='Iteration 1 $c_{l}$ analytical',
          linewidth=width)

plt.plot(y_cfd,Fz_cfd1,'go-',
          label='Iteration 1 $c_{l}$ cfd',
          linewidth=width)
errMax3 = 100*(Fz_cfd0*(1+theta/np.deg2rad(10)) - Fz_cfd1)/(Fz_cfd0*(1+theta/np.deg2rad(10)))
errMax3 = errMax3.round(decimals=1)
errMax3 = np.max(errMax3)
plt.title('Lift distribution at step 1\n Max error: '+str(errMax3)+'%',fontsize=titleSize)
plt.xlabel('y position [m]',fontsize=textSize)
plt.ylabel('Lift [N/m]',fontsize=textSize)
plt.xlim(-5,5)
plt.grid()
plt.legend(fontsize=textSize)
plt.xticks(fontsize=textSize)
plt.yticks(fontsize=textSize)





plt.show()


