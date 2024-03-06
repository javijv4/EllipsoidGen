#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 18:34:37 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
import meshio as io
import pygmsh

def gen_circle(ndiv, ndiv_r, side=0.6,  r = 1.0, order=2):
    assert (order == 1) or (order == 2)

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
        ien = mesh.cells_dict['quad9']
        arr = np.array([0,3,2,1,7,6,5,4,8])
        ien[0:(ndiv-1)*(ndiv-1)] = ien[0:(ndiv-1)*(ndiv-1), arr] 
        mesh = io.Mesh(mesh.points, {'quad9': ien})

    return mesh


def warp_to_ellipsoid(circle_mesh, a, b, c, theta_max=None, zmax=None):
    if theta_max is None and zmax is None: return
    elif theta_max is None:
        theta_max = np.arccos(zmax/c)

    x, y = circle_mesh.points[:,0:2].T

    phi = np.arctan2(b*x, a*y)
    theta = np.sqrt(x**2+y**2)
    theta = theta/np.max(theta)*theta_max

    x = a*np.sin(theta)*np.cos(phi)
    y = b*np.sin(theta)*np.sin(phi)
    z = c*np.cos(theta)

    xyz = np.array([x, y, z]).T
    cells = circle_mesh.cells[0]
    mesh = io.Mesh(xyz, {cells.type: cells.data})

    return mesh


def build_hexahedral_mesh(ellipsoid_mesh1, ellipsoid_mesh2, ndiv_t, order=1):
    xyz_endo = ellipsoid_mesh1.points
    ien_quad = ellipsoid_mesh1.cells[0].data
    xyz_epi = ellipsoid_mesh2.points

    if order == 1:
        ien = []
        xyz = xyz_endo
        ien_past = ien_quad

        for i in range(ndiv_t-1):
            xyz_mid = (xyz_epi - xyz_endo)*(i+1)/ndiv_t + xyz_endo
            ien_mid = ien_quad + len(xyz)

            xyz = np.vstack([xyz, xyz_mid])
            ien_hex = np.hstack([ien_past, ien_mid])
            ien.append(ien_hex)
            ien_past = ien_mid

        # Adding epi
        ien_epi = ien_quad + len(xyz)
        ien_hex = np.hstack([ien_past, ien_epi])

        xyz = np.vstack([xyz, xyz_epi])
        ien.append(ien_hex)

        ien = np.vstack(ien)
        mesh = io.Mesh(xyz, {'hexahedron': ien})

    elif order == 2:
        ien = []
        xyz = xyz_endo
        xyz_nodes = len(xyz_endo)
        ien0 = ien_quad
        xyz0 = xyz_endo
        for i in range(ndiv_t-1):
            xyz2 = (xyz_epi - xyz_endo)*(i+1)/ndiv_t + xyz_endo
            xyz1 = (xyz2 - xyz0)*0.5 + xyz0

            ien1 = ien_quad + len(xyz)
            ien2 = ien_quad + len(xyz) + xyz_nodes

            xyz = np.vstack([xyz, xyz1, xyz2])
            ien_hex = np.hstack([ien0[:,:4], ien2[:,:4], ien0[:,4:8], ien2[:,4:8], 
                            ien1[:,0:4], ien1[:,[7,5,4,6]], 
                            ien0[:,-1,None], ien2[:,-1,None], ien1[:,-1,None]])
            ien.append(ien_hex)
            ien0 = ien2
            xyz0 = xyz2


        # Adding epi
        ien2 = ien_quad + len(xyz)
        xyz2 = xyz_epi
        xyz1 = (xyz2 - xyz0)*0.5 + xyz0

        ien1 = ien_quad + len(xyz)
        ien2 = ien_quad + len(xyz) + xyz_nodes
        ien_hex = np.hstack([ien0[:,:4], ien2[:,:4], ien0[:,4:8], ien2[:,4:8], 
                        ien1[:,0:4], ien1[:,[7,5,4,6]], 
                        ien0[:,-1,None], ien2[:,-1,None], ien1[:,-1,None]])

        xyz = np.vstack([xyz, xyz1, xyz2])
        ien.append(ien_hex)
        ien = np.vstack(ien)                       

        mesh = io.Mesh(xyz, {'hexahedron27': ien})
        
    return mesh


def gen_ellipsoid_mesh(radius1, radius2, height1, height2, theta_max, ndiv_t, ndiv_r, order):
    ndiv = ndiv_r//2

    a1, b1 = radius1, radius1
    c1 = height1
    a2, b2 = radius2, radius2
    c2 = height2

    circle_mesh1 = gen_circle(ndiv, ndiv_r, order=order, side=0.4)
    ellipsoid_mesh1 = warp_to_ellipsoid(circle_mesh1, a1, b1, c1, theta_max=theta_max)

    circle_mesh2 = gen_circle(ndiv, ndiv_r, order=order, side=0.4)
    ellipsoid_mesh2 = warp_to_ellipsoid(circle_mesh2, a2, b2, c2, zmax=np.min(ellipsoid_mesh1.points[:,2]))

    mesh = build_hexahedral_mesh(ellipsoid_mesh1, ellipsoid_mesh2, ndiv_t, order=order)

    return mesh
