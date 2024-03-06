#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:51:17 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
import pygmsh
import meshio as io

r = 1.0
side = 0.6
ndiv = 3
ndiv_r = 3
with pygmsh.geo.Geometry() as geom:
    lcar = 0.1
    p0 = geom.add_point([0.0, 0.0], lcar)
    p1 = geom.add_point([-side/2, -side/2], lcar)
    p2 = geom.add_point([side/2, -side/2], lcar)
    p3 = geom.add_point([side/2, side/2], lcar)
    p4 = geom.add_point([-side/2, side/2], lcar)

    p5 = geom.add_point([r*np.cos(5*np.pi/4), r*np.sin(5*np.pi/4)], lcar)
    p6 = geom.add_point([r*np.cos(7*np.pi/4), r*np.sin(7*np.pi/4)], lcar)
    p7 = geom.add_point([r*np.cos(np.pi/4), r*np.sin(np.pi/4)], lcar)
    p8 = geom.add_point([r*np.cos(3*np.pi/4), r*np.sin(3*np.pi/4)], lcar)

    l1 = geom.add_line(p1, p2)
    l2 = geom.add_line(p2, p3)
    l3 = geom.add_line(p3, p4)
    l4 = geom.add_line(p4, p1)

    l5 = geom.add_line(p1, p5)
    l6 = geom.add_line(p2, p6)
    l7 = geom.add_line(p3, p7)
    l8 = geom.add_line(p4, p8)

    c1 = geom.add_circle_arc(p5, p0, p6)
    c2 = geom.add_circle_arc(p6, p0, p7)
    c3 = geom.add_circle_arc(p7, p0, p8)
    c4 = geom.add_circle_arc(p8, p0, p5)


    geom.set_transfinite_curve(l1, ndiv, "Progression", 1.0)
    geom.set_transfinite_curve(l2, ndiv, "Progression", 1.0)
    geom.set_transfinite_curve(l3, ndiv, "Progression", 1.0)
    geom.set_transfinite_curve(l4, ndiv, "Progression", 1.0)

    geom.set_transfinite_curve(l5, ndiv_r, "Progression", 1.0)
    geom.set_transfinite_curve(l6, ndiv_r, "Progression", 1.0)
    geom.set_transfinite_curve(l7, ndiv_r, "Progression", 1.0)
    geom.set_transfinite_curve(l8, ndiv_r, "Progression", 1.0)

    ll = geom.add_curve_loop([l1,l2,l3,l4])

    rect = geom.add_plane_surface(ll)

    geom.set_transfinite_curve(c1, ndiv, "Progression", 1.0)
    geom.set_transfinite_curve(c2, ndiv, "Progression", 1.0)
    geom.set_transfinite_curve(c3, ndiv, "Progression", 1.0)
    geom.set_transfinite_curve(c4, ndiv, "Progression", 1.0)

    cc1 = geom.add_curve_loop([l1,l6,-c1,-l5])
    cc2 = geom.add_curve_loop([l2,l7,-c2,-l6])
    cc3 = geom.add_curve_loop([l3,l8,-c3,-l7])
    cc4 = geom.add_curve_loop([l4,l5,-c4,-l8])

    circ1 = geom.add_plane_surface(cc1)
    circ2 = geom.add_plane_surface(cc2)
    circ3 = geom.add_plane_surface(cc3)
    circ4 = geom.add_plane_surface(cc4)

    geom.set_transfinite_surface(rect, "Left", [p1,p2,p3,p4])
    geom.set_transfinite_surface(circ1, "Left", [p1,p2,p6,p5])
    geom.set_transfinite_surface(circ2, "Left", [p2,p3,p7,p6])
    geom.set_transfinite_surface(circ3, "Left", [p3,p4,p8,p7])
    geom.set_transfinite_surface(circ4, "Left", [p4,p1,p5,p8])
    geom.set_recombined_surfaces([rect, circ1, circ2, circ3, circ4])

    mesh = geom.generate_mesh(order=2)

mesh = io.Mesh(mesh.points, {'quad9': mesh.cells_dict['quad9']})

ien = mesh.cells[0].data
arr = np.array([0,3,2,1,7,6,5,4,8])
ien[0:(ndiv-1)*(ndiv-1)] = ien[0:(ndiv-1)*(ndiv-1), arr] 

io.write('circle.vtu', mesh)

