import polyscope as ps
import numpy as np

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/MESH_READER.m
# https://en.wikipedia.org/wiki/OFF_(file_format)
# https://stackoverflow.com/questions/31129968/off-files-on-python
def read_off(file):
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    n_vertices, n_faces, n_edges = tuple([int(s) for s in file.readline().strip().split(' ')])
    vertices = np.asarray([[float(s) for s in file.readline().strip().split(' ')] for _ in range(n_vertices)])
    faces = np.asarray([[int(s) for s in file.readline().strip().split(' ')][1:] for _ in range(n_faces)])
    return vertices, faces



def normv(v):
    """L2 norm of vector"""
    s = 0
    for x in v: s += x
    return np.sqrt(s)


# https://raw.githubusercontent.com/azencot/gp.fns/master/code/MESH.m
class Mesh:
    def __init__(name, vertices, faces):
        self.name = name
        self.vertices = vertices
        self.triangles = triangles # only triangles.
        for t in self.triangles:
            assert len(t) == 3

        self.nv = len(self.vertices)
        self.nf = len(self.faces)

        X = self.vertices; T = self.triangles
        NN1 = X[T[:, 0], :] - X[T[:, 1], :]
        NN2 = X[T[:, 0], :] - X[T[:, 2], :]
        NN = np.cross(NN1, NN2)
        II = np.concatenate((T[:, 0], T[:, 1], T[:, 2]))
        JJ = np.concatenate((T[:, 1], T[:, 2], T[:, 0]))




# https://raw.githubusercontent.com/azencot/gp.fns/master/code/test_sphere_jet.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/test_sphere_rot.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/vortices_func.m

if __name__ == "__main__":
    ps.init()
    with open("./gp.fns/data/sphere_s3.off") as f:
      (vertices, faces) = read_off(f)
      ps.register_surface_mesh("s3", vertices, faces, smooth_shade=True)
    ps.show()

