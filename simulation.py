import polyscope as ps

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/MESH_READER.m
# https://en.wikipedia.org/wiki/OFF_(file_format)
# https://stackoverflow.com/questions/31129968/off-files-on-python
def read_off(file):
    if 'OFF' != file.readline().strip():
        raise('Not a valid OFF header')
    n_verts, n_faces, n_edges = tuple([int(s) for s in file.readline().strip().split(' ')])
    verts = [[float(s) for s in file.readline().strip().split(' ')] for _ in range(n_verts)]
    faces = [[int(s) for s in file.readline().strip().split(' ')][1:] for _ in range(n_faces)]
    return verts, faces


# https://raw.githubusercontent.com/azencot/gp.fns/master/code/MESH.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/test_sphere_jet.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/test_sphere_rot.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/vortices_func.m
