#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 13:18:22 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt

R = 2.0
c = 10
r = np.linspace(0, R, 10)
theta = 0

rho = np.sqrt(r)
# aux =
z = c*np.sqrt(1 - (R-rho)**2/R**2)

plt.figure(1, clear=True)
plt.plot([0,R],[0,c], 'k')
plt.plot(r, z, '--o', lw=3)
plt.plot(r, rho, '--o', lw=3)
plt.plot(rho, z, '--', lw=3)