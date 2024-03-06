#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 18:36:44 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
import meshio as io
from EllipsoidGen import gen_ellipsoid_mesh


ndiv_r = 15
ndiv_t = 3
order = 2

radius1 = 20
height1 = 46
radius2 = 32
height2 = 53
theta_max = np.pi/2*1.3

mesh = gen_ellipsoid_mesh(radius1, radius2, height1, height2, theta_max, ndiv_t, ndiv_r, order)

io.write('ellipsoid.vtu', mesh)