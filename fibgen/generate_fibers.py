import meshio as io
import cheartio as chio
import numpy as np
import meshio as io


@staticmethod
def rotate_basis(eC, eL, eT, alpha):
    eC = eC/np.linalg.norm(eC, axis=1)[:,None]
    eT = eT/np.linalg.norm(eT, axis=1)[:,None]
    eL = eL/np.linalg.norm(eL, axis=1)[:,None]

    # Matrix of directional vectors
    Q = np.stack([eC, eL, eT], axis=-1)
    Q = np.transpose(Q, (2, 1, 0))

    # Create rotation matrix - from Doste code
    axis = eT
    R = np.array([[np.cos(alpha) + (axis[:, 0]**2)*(1 - np.cos(alpha)), axis[:,0] * axis[:,1]*(1 - np.cos(alpha)) - axis[:,2]*np.sin(alpha), axis[:,0]*axis[:,2]*(1 - np.cos(alpha)) + axis[:,1]*np.sin(alpha)],
                            [axis[:,1]*axis[:,0]*(1 - np.cos(alpha)) + axis[:,2]*np.sin(alpha), np.cos(alpha) + (axis[:,1]**2)*(1 - np.cos(alpha)), axis[:,1]*axis[:, 2]*(1 - np.cos(alpha)) - axis[:, 0]*np.sin(alpha)],
                            [axis[:,2]*axis[:,0]*(1 - np.cos(alpha)) - axis[:,1]*np.sin(alpha), axis[:,2]*axis[:,1]*(1 - np.cos(alpha)) + axis[:, 0]*np.sin(alpha), np.cos(alpha)+(axis[:, 2]**2)*(1 - np.cos(alpha))]])

    # Rotate the circumferential direction around the transmural direction
    QX = np.zeros_like(R)
    for i in range(len(eC)):
        QX[:, :, i] = np.matmul(Q[:, :, i], R[:, :, i])

    return QX

def get_discontinouts_mesh(mesh, gl):

    # Generate a discontinuous mesh
    new_ien = mesh.cells[0].data
    apex_nodes = np.where(np.isclose(mesh.points[:,0],0)*np.isclose(mesh.points[:,1],0))[0]
    inter_nodes = apex_nodes
    new_xyz =np.copy(mesh.points)

    inter_elems = np.where(np.any(np.isin(new_ien, inter_nodes), axis=1))[0]
    midpoints = np.mean(new_xyz[new_ien[inter_elems]], axis=1)
    angles = np.arctan2(midpoints[:,1], midpoints[:,0])
    order = np.argsort(angles)
    ind = inter_elems[order.reshape([4,-1])]
    midpoints = midpoints[order].reshape([4,-1,3])

    new_gl = gl.copy()
    map_new_nodes = np.arange(len(new_xyz))
    map_node_faces = []
    for i in range(4):
        faces = new_ien[ind[i]][:,0:4]
        midpoints = np.mean(new_xyz[faces], axis=1)
        order = np.argsort(midpoints[:,-1])
        faces = np.vstack([faces[order], new_ien[ind[i][order[-1]],4:]])
        midpoints = np.mean(new_xyz[faces], axis=1)
        vectors = midpoints  - new_xyz[inter_nodes]
        faces = faces[~np.isin(faces, inter_nodes)].reshape([-1,3])


        if i == 0:  # No add nodes, but the gl needs to be fix
            new_gl[inter_nodes] = vectors/np.linalg.norm(vectors, axis=1)[:,None]
            node_face = np.vstack([inter_nodes, faces.T]).T

        else:
            new_nodes = np.arange(len(inter_nodes)) + len(new_xyz)
            new_gl = np.vstack([new_gl, vectors/np.linalg.norm(vectors, axis=1)[:,None]])

            map_new_nodes[inter_nodes] = np.arange(len(inter_nodes)) + len(new_xyz)

            new_ien[ind[i]] = map_new_nodes[new_ien[ind[i]]]
            new_xyz = np.vstack([new_xyz, new_xyz[inter_nodes]])
            map_new_nodes = np.append(map_new_nodes, inter_nodes)

            map_new_nodes[0:len(mesh.points)] = np.arange(len(mesh.points))

            node_face = np.vstack([new_nodes, faces.T]).T
        
        map_node_faces.append(node_face)

    map_node_faces = np.vstack(map_node_faces)
    mesh = io.Mesh(new_xyz, {'hexahedron': new_ien})
    io.write(path + 'fiber.vtu', mesh)

    return mesh, map_new_nodes, new_gl, map_node_faces



def get_fibers(gl, gt, trans):

    eL_lv = gl/np.linalg.norm(gl, axis=1)[:,None]
    eT_lv = gt/np.linalg.norm(gt, axis=1)[:,None]

    # circumferential
    eC_lv = np.cross(eL_lv, eT_lv, axisa=1, axisb=1)
    eC_lv = eC_lv/np.linalg.norm(eC_lv, axis=1)[:,None]

    # Ensuring orthogonality
    eL_lv = np.cross(eT_lv, eC_lv, axisa=1, axisb=1)
    eL_lv = eL_lv/np.linalg.norm(eL_lv, axis=1)[:,None]


    alpha = (trans*2-1)*np.pi/3
    Qlv = rotate_basis(eC_lv, eL_lv, eT_lv, alpha)

    f, n, s = Qlv
    f = f.T
    n = n.T
    s = s.T

    return f, s, n


path = '../mesh2/'
mesh_quad = chio.read_mesh(path + 'ellipsoid_quad', meshio=True)
mesh = chio.read_mesh(path + 'ellipsoid', meshio=True)

map_quad_lin = chio.map_between_meshes(mesh_quad, mesh)

long = chio.read_dfile(path + 'Long-1.D')[map_quad_lin]
trans = chio.read_dfile(path + 'Trans-1.D')[map_quad_lin]
gl = chio.read_dfile(path + 'GradLong-1.D')
gt = chio.read_dfile(path + 'GradTrans-1.D')

lin_fibers = get_fibers(gl, gt, trans)

mesh.point_data['f'] = lin_fibers[0]
mesh.point_data['s'] = lin_fibers[1]
mesh.point_data['n'] = lin_fibers[2]
mesh.point_data['trans'] = trans
mesh.point_data['long'] = long
mesh.point_data['gl'] = gl

save = np.hstack(lin_fibers)
chio.write_dfile(path + 'fiber_lin.field', save)
io.write(path + 'fiber_lin.vtu', mesh)


disc_mesh, map_new_nodes, disc_gl, map_node_faces = get_discontinouts_mesh(mesh, gl)
disc_fibers = get_fibers(disc_gl, gt[map_new_nodes], trans[map_new_nodes])

# f = disc_fibers[0]
# s = disc_fibers[1]
# n = disc_fibers[2]

# for i in range(len(map_node_faces)):
#     node = map_node_faces[i,0]
#     face = map_node_faces[i]
#     midpoint = np.mean(disc_mesh.points[face], axis=0)
#     vector = midpoint - disc_mesh.points[node]
#     aux = np.copy(vector)
#     vector[0] = aux[1]
#     vector[1] = -aux[0]

#     f[node] = vector/np.linalg.norm(vector)
#     vector = np.cross(f[node], s[node])
#     n[node] = vector/np.linalg.norm(vector)

# arr = np.array([[0,1,0],[2/np.sqrt(2),2/np.sqrt(2),0],[1,0,0]])
# arr = np.mean(arr, axis=0)
# print(arr/np.linalg.norm(arr))

# print(list(map_node_faces))
# disc_fibers = [f,s,n]

disc_mesh.point_data['f'] = disc_fibers[0]
disc_mesh.point_data['s'] = disc_fibers[1]
disc_mesh.point_data['n'] = disc_fibers[2]
disc_mesh.point_data['trans'] = trans[map_new_nodes]
disc_mesh.point_data['long'] = long[map_new_nodes]
disc_mesh.point_data['gl'] = disc_gl
disc_mesh.point_data['gt'] = gt[map_new_nodes]



save = np.hstack(disc_fibers)
chio.write_dfile(path + 'fiber.field', save)
io.write(path + 'fiber.vtu', disc_mesh)
chio.write_mesh(path + 'fiber', disc_mesh.points, disc_mesh.cells[0].data)

x = np.mean([0.845723, 0.937863, 0.5])
y = np.mean([0.5, -0.230757, -0.845723])
z = np.mean([0.186419, 0.259159, 0.186419])
vec = np.array([x,y,z])
vec = vec/np.linalg.norm(vec)
print(vec)