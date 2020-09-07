#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 10:59:01 2020

@author: Jean-Philippe Kuntzer
"""


import numpy as np
import matplotlib.pyplot as plt


def main():
    x = np.linspace(-1,1,10000)
    xp = np.abs(x)
    y1 = phi(xp,'L')
    y2 = phi(xp,'TPS')
    y3 = phi(xp,'G')
    y4 = phi(xp,'HMQ')
    y5 = phi(xp,'HIMQ')
    y6 = phi(xp,'C0')
    y7 = phi(xp,'C2')
    y8 = phi(xp,'C4')
    y9 = phi(xp,'C6')    
    y10 = phi(xp,'EH')
    
    # linewidth
    lwdth = 10
    titleSize = 28
    labelSize = 20
    
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(x, y1, label='$Linear: \phi = r $', linewidth=lwdth)
    # ax.plot(x, y2, label='$\phi = r^2 ln(r) $', linewidth=lwdth)
    ax.plot(x, y3, 'r', label='Gaussian: $\phi = e^{-\epsilon r^2}$', linewidth=lwdth)
    # ax.plot(x, y4, label='$\phi = \sqrt{1+(\epsilon r)^2}$', linewidth=lwdth)
    # ax.plot(x, y5, label='$\phi = (\sqrt{1+(\epsilon r)^2})^{-1}$', linewidth=lwdth)
    ax.plot(x, y6, label='Wendland $C_0: \phi = (1−r)^2$', linewidth=lwdth)
    # ax.plot(x, y7, label='$\phi = (1−r)^4 (4r+ 1)$', linewidth=lwdth)
    # ax.plot(x, y8, label='$\phi = (1−r)^6 (35r^2+ 18r+ 3)$', linewidth=lwdth)
    # ax.plot(x, y9, label='$\phi = (1−r)^8 (32r^3+ 25r^2+ 8r+ 1)$', linewidth=lwdth)
    # ax.plot(x, y10, label='Euclid hat $ \phi = \pi*((1/12 r^3) - r \epsilon^2 + 4/3 \epsilon^3)$', linewidth=lwdth)
    
    plt.title('Common radial basis functions', fontsize=titleSize)
    ax.legend(loc='upper right',fontsize=labelSize)
    plt.xlabel('x',fontsize=labelSize)
    plt.ylabel('$\phi(x)$',fontsize=labelSize)
    plt.xticks(fontsize=labelSize)
    plt.yticks(fontsize=labelSize)
    plt.show()

# φ(r) = 
# φ(r) = (1−r)4(4r+ 1)
# φ(r) = (1−r)6(35r2+ 18r+ 3)
# φ(r) = (1−r)8(32r3+ 25r2+ 8r+ 1)

def phi(r,fun):
    """
    Set of radial basis functions that the user can choose of. After some
    test "Wendland C2" seems to be the better choice, but this is really
    up to user preference.
    """
    eps = 1
    # r = np.linalg.norm(x1-x2)
    if fun == "G":
        # Gaussian
        phi_x = np.exp(-eps*r**2)
    elif fun == "L":
        # Linear
        phi_x = r
    elif fun == 'cubic':
        # Cubic
        phi_x = r**3
    elif fun == "TPS":
        # Thin plate spline
        phi_x = r**2 * np.log(r)
    elif fun == "Po4":
        # Polyharmonic spline
        phi_x = r**3 * np.log(r**r)
    elif fun == "Po4":
        # Polyharmonic spline
        phi_x = r**5 * np.log(r**r)
    elif fun == "HMQ":
        # Hardy's multiquadratic
        phi_x = (eps**2 + r**2)**0.5
    elif fun == "HIMQ":
        # Hardy's inverse multiquadratic
        phi_x = 1/(eps**2 + r**2)**0.5
    elif fun == "C0":
        # Wendland's C0
        phi_x = (1-r)**2
    elif fun == "C2":
        # Wendland's C2
        phi_x = (1-r)**4 * (4*r + 1)
    elif fun == "C4":
        # Wendland's C4
        phi_x = (1-r)**6 * (35*r**2 + 18*r + 3)
        # print('C4')
    elif fun == "C6":
        # Wendland's C6
        phi_x = (1-r)**8 * (32*r**3 + 25*r**2 + 8*r + 1)
    elif fun == "EH":
        # Euclid's hat
        phi_x = np.pi*((1/12*r**3) - r*eps**2 + 4/3*eps**3)
    return phi_x

if __name__ == "__main__":
    main()
