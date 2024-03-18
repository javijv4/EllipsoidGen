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


def fix_orientation(mesh, order, ndiv):
    if order == 1:
        ien = mesh.cells_dict['quad']
        arr = np.array([0,3,2,1])
        ien[0:(ndiv-1)*(ndiv-1)] = ien[0:(ndiv-1)*(ndiv-1), arr]
        mesh = io.Mesh(mesh.points, {'quad': ien})
    elif order == 2:
        ien = mesh.cells_dict['quad9']
        arr = np.array([0,3,2,1,7,6,5,4,8])
        ien[0:(ndiv-1)*(ndiv-1)] = ien[0:(ndiv-1)*(ndiv-1), arr]
        mesh = io.Mesh(mesh.points, {'quad9': ien})

    return mesh


def fix_circular_mesh(mesh, order, ndiv, r):
    ien = mesh.cells[0].data

    xyz = mesh.points
    cnodes = np.unique(ien[(ndiv-1)*(ndiv-1):])

    if order == 1:
        ncirc = ndiv*4-4
    elif order == 2:
        ncirc = (ndiv*2-1)*4-4

    angles = np.linspace(np.pi, -np.pi, ncirc+1, endpoint=True)[:-1]

    node_angles = np.arctan2(xyz[:,1], xyz[:,0])
    node_angles[np.isclose(node_angles, -np.pi)] = np.pi

    sort = np.argsort(node_angles[cnodes])[::-1]
    cnodes = cnodes[sort].reshape([-1, len(cnodes)//len(angles)])

    xs = np.concatenate([np.full(ncirc//8,-side/2),
                         np.linspace(-side/2, side/2, ncirc//4+1, endpoint=True)[:-1],
                         np.full(ncirc//4,side/2),
                         np.linspace(side/2, -side/2, ncirc//4+1, endpoint=True)[:-1],
                         np.full(ncirc//8,-side/2)])
    ys = np.concatenate([np.linspace(0, side/2, ncirc//8+1, endpoint=True)[:-1],
                         np.full(ncirc//4,side/2),
                         np.linspace(side/2, -side/2, ncirc//4+1, endpoint=True)[:-1],
                         np.full(ncirc//4,-side/2),
                         np.linspace(-side/2, 0, ncirc//8+1, endpoint=True)[:-1]])
    angles_s = np.arctan2(ys, xs)

    for i in range(len(angles)):
        theta = angles[i]
        theta_s = angles_s[i]
        nodes = cnodes[i]
        x, y = xyz[nodes,0:2].T

        radius = np.sqrt(x**2+y**2)
        sort = np.argsort(radius)
        nodes = nodes[sort]

        min_r = np.min(radius)
        max_r = r

        alpha = np.linspace(0,1,len(nodes), endpoint=True)

        radius = np.linspace(min_r, max_r, len(nodes), endpoint=1)
        xyz[nodes[:],0] = radius[:]*(np.cos(theta)*alpha + np.cos(theta_s)*(1-alpha))
        xyz[nodes[:],1] = radius[:]*(np.sin(theta)*alpha + np.sin(theta_s)*(1-alpha))

    mesh.points = xyz
    return mesh

r = 1.0
side = 0.6
ndiv = 5
ndiv_r = 3
order=2

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

    mesh = geom.generate_mesh(order=order)

if order == 1:
    mesh = io.Mesh(mesh.points, {'quad': mesh.cells_dict['quad']})
elif order == 2:
    mesh = io.Mesh(mesh.points, {'quad9': mesh.cells_dict['quad9']})




# Need to get rid of node 0 (does not belong to any element)
mesh.points = mesh.points[1:]
mesh.cells[0].data = mesh.cells[0].data - 1

mesh = fix_orientation(mesh, order, ndiv)
mesh = fix_circular_mesh(mesh, order, ndiv, r)


io.write('circle.vtu', mesh)

