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
import cheartio as chio


ndiv_r = 20
ndiv_t = 6
order = 2

radius1 = 20
height1 = 46
radius2 = 32
height2 = 53
theta_max = np.pi/2*1.3

mesh, bdata = gen_ellipsoid_mesh(radius1, radius2, height1, height2, theta_max, ndiv_t, ndiv_r, order)

chio.write_mesh('mesh3/ellipsoid', mesh.points, mesh.cells[0].data)
chio.write_bfile('mesh3/ellipsoid', bdata)

io.write('mesh3/ellipsoid.vtu', mesh)