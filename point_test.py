#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:33:50 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt

# x = np.linspace(-0.5,0.5,10, endpoint=True)
# y = np.linspace(-0.5,0.5,10, endpoint=True)
# x, y = np.meshgrid(x, y)

r = np.linspace(0,0.5,10, endpoint=True)
t = np.linspace(0,2*np.pi,30, endpoint=True)
r, t = np.meshgrid(r, t)

x = r*np.cos(t)
y = r*np.sin(t)


a, b, c = 1, 1, 2

z = c*np.sqrt(1-x**2/a**2-y**2/b**2)

phi = np.arctan2(b*x, (a*y))
theta = np.sqrt(x**2+y**2)
theta = theta/np.max(theta)*np.pi/2

x = a*np.sin(theta)*np.cos(phi)
y = b*np.sin(theta)*np.sin(phi)
z = c*np.cos(theta)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(x, y, z, linewidth=2, edgecolors='k', antialiased=False)