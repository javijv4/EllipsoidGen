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
from scipy.spatial import KDTree

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


    mesh = fix_orientation(mesh, order, ndiv)
    mesh = fix_circular_mesh(mesh, order, ndiv, r)

    return mesh


def fix_orientation(mesh, order, ndiv):
    if order == 1:
        mesh = io.Mesh(mesh.points, {'quad': mesh.cells_dict['quad']})
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

    for i in range(len(angles)):
        theta = angles[i]
        nodes = cnodes[i]
        x, y = xyz[nodes,0:2].T

        radius = np.sqrt(x**2+y**2)
        sort = np.argsort(radius)
        nodes = nodes[sort]

        min_r = np.min(radius)
        max_r = r

        radius = np.linspace(min_r, max_r, len(nodes), endpoint=1)
        xyz[nodes[1:],0] = radius[1:]*np.cos(theta)
        xyz[nodes[1:],1] = radius[1:]*np.sin(theta)

    mesh.points = xyz
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


def get_surface_mesh(mesh):
    ien = mesh.cells[0].data

    if ien.shape[1] == 8:   # Assuming hex
        array = np.array([[0,1,5,4],
                          [1,2,6,5],
                          [2,3,7,6],
                          [3,0,4,7],
                          [0,1,2,3],
                          [4,5,6,7]])
        nelems = np.repeat(np.arange(ien.shape[0]),6)
        faces = np.vstack(ien[:,array])
        sort_faces = np.sort(faces,axis=1)

        f, i, c = np.unique(sort_faces, axis=0, return_counts=True, return_index=True)
        ind = i[np.where(c==1)[0]]
        bfaces = faces[ind]
        belem = nelems[ind]

    elif ien.shape[1] == 27:   # Assuming hex27
        array = np.array([[0,1,5,4,8,17,12,16,22],
                          [1,2,6,5,9,18,13,17,21],
                          [2,3,7,6,10,19,14,18,23],
                          [3,0,4,7,11,16,15,19,20],
                          [0,1,2,3,8,9,10,11,24],
                          [4,5,6,7,12,13,14,15,25]])
        nelems = np.repeat(np.arange(ien.shape[0]),6)
        faces = np.vstack(ien[:,array])
        sort_faces = np.sort(faces,axis=1)

        f, i, c = np.unique(sort_faces, axis=0, return_counts=True, return_index=True)
        ind = i[np.where(c==1)[0]]
        bfaces = faces[ind]
        belem = nelems[ind]

    return belem, bfaces


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
        nodes_minus_epi = len(xyz)
        ien2 = ien_quad + nodes_minus_epi
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

    bdata = get_boundary_data(mesh, ellipsoid_mesh1, ellipsoid_mesh2, order)


    return mesh, bdata

def get_boundary_data(mesh, ellipsoid_mesh1, ellipsoid_mesh2, order):

    if order == 1:
        xyz = mesh.points
        belems, bfaces = get_surface_mesh(mesh)

    elif order == 2:
        xyz = mesh.points
        belems, bfaces = get_surface_mesh(mesh)

    midpoints = np.mean(xyz[bfaces], axis=1)

    xyz_elems1 = ellipsoid_mesh1.points[ellipsoid_mesh1.cells[0].data]
    xyz_elems2 = ellipsoid_mesh2.points[ellipsoid_mesh2.cells[0].data]

    midpoints1 = np.mean(xyz_elems1, axis=1)
    midpoints2 = np.mean(xyz_elems2, axis=1)

    tree = KDTree(midpoints)
    _, corr1 = tree.query(midpoints1)
    _, corr2 = tree.query(midpoints2)

    labels = np.zeros(len(bfaces), dtype=int) + 3
    labels[corr1] = 1
    labels[corr2] = 2

    bdata = np.vstack([belems, bfaces.T, labels]).T

    return bdata



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

    mesh, bdata = build_hexahedral_mesh(ellipsoid_mesh1, ellipsoid_mesh2, ndiv_t, order=order)

    return mesh, bdata
