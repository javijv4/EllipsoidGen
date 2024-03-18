#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 18:36:44 2024

@author: Javiera Jilberto Vallejos
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import meshio as io
from EllipsoidGen import gen_ellipsoid_mesh
import cheartio as chio


out_path = 'mesh2'
if not os.path.exists(out_path): os.mkdir(out_path)

ndiv_r = 14
ndiv_t = 6
order = 2

radius1 = 20
height1 = 46
radius2 = 32
height2 = 53
theta_max = np.pi/2*1.3

mesh, bdata = gen_ellipsoid_mesh(radius1, radius2, height1, height2, theta_max, ndiv_t, ndiv_r, order)

chio.write_mesh(out_path + '/ellipsoid_quad', mesh.points, mesh.cells[0].data)
chio.write_bfile(out_path + '/ellipsoid_quad', bdata)

apex_id = np.where(np.isclose(mesh.points[:,0],0)*np.isclose(mesh.points[:,1],0))[0]
chio.write_specific(out_path + '/apex', np.array([apex_id]), np.zeros(len(apex_id)))
print(list(apex_id))
io.write(out_path + '/ellipsoid_quad.vtu', mesh)