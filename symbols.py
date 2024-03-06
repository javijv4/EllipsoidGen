#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 12:54:14 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

r, R, theta, c, k = sp.symbols('r, R, theta, c, k')
x, y = sp.symbols('x, y')

x = r*sp.cos(theta)
y = r*sp.sin(theta)
z = c*sp.sqrt(1-x**2/R**2-y**2/R**2)
dzdr = sp.diff(z, r)

eq = dzdr - k
sol = sp.solve(eq, r)

R = 2.0
