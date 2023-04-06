import polyscope as ps
import numpy as np
# from scipy.sparse import sparse
import scipy.sparse.linalg as sla
from scipy.sparse import coo_matrix
from scipy.sparse import spdiags
from scipy.sparse import csr_matrix




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

        self.X = self.vertices; self.T = self.triangles

        NN1 = X[T[:, 0], :] - X[T[:, 1], :]
        NN2 = X[T[:, 0], :] - X[T[:, 2], :]
        self.NN = np.cross(NN1, NN2)

        self.ta = normV(NN)/2;
        # self.N = NN./repmat(MESH.normv(NN),1,3);
        self.N = NN / np.tile(MESH.normv(NN), (1, 3))


        self.II = np.concatenate((T[:, 0], T[:, 1], T[:, 2]))
        self.JJ = np.concatenate((T[:, 1], T[:, 2], T[:, 0]))

        self.IIn = np.concatenate((self.II, self.JJ, self.II, self.JJ))
        self.JJn = np.concatenate((self.JJ, self.II, self.II, self.JJ))


        # Use the Barycentric areas
        self.M = self.mass_matrix_barycentric()
        AA = np.sum(self.M, axis=1)
        self.va = np.array(AA).flatten() 
        self.va_inv = np.ones(self.va.shape) / self.va
        self.A = spdiags(self.va, 0, self.nv, self.nv, format='csr')
        self.A_INV = spdiags(self.va_inv, 0, self.nv, self.nv, format='csr')

        # compute cotan Laplacian
        self.W = self.cotLaplacian()
        self.L = self.A_INV @ self.W
        self.LL = self.L @ self.L

        # compute some eigenvectors of L
        neigs = min(50, self.nv)
        evals, evecs = eigs(self.W, M=self.A, k=neigs, sigma=0, which='LM', tol=1e-5)
        ii = np.argsort(np.diag(evals))
        evals = evals[ii, ii]
        evecs = evecs[:, ii]
        self.LB_basis = evecs
        self.LB_basisi = np.dot(evecs.T, self.A.todense())
        self.LB_evals = np.diag(evals)

        # create preconditioner for W
        opts = {'fill_factor': 1E-5, 'drop_tol': 1E-5, 'type': 'ict'}
        ilu = spilu(self.W.tocsc(), drop_tol=opts['drop_tol'])
        self.W_precon = LinearOperator(self.W.shape, ilu.solve)

        # initialize data needed for vf2op
        self.JJc1, self.JJc2, self.JJc3 = self.vf2op_initialization()

        # remember for iterative solve
        self.last_laplace_inverse = None

    def dot(self, f, g):
        # l2 = np.dot(f, self.va @ g) # should this be * or @? probably @
        l2 = np.dot(f, self.va * g) # should this be * or @? probably @
        return l2


    def vort2energy(mesh, w):
        s = -mesh.laplaceInverseF(w)
        ke = -np.dot(s, w * mesh.va) # weight by barycentric area (mesh.va)
        return ke
    
    def dotV(mesh, U, V):
        a = U * V
        sa = np.sum(a, axis=1)
        l2 = np.dot(mesh.ta, sa)
        return l2

    def norm(mesh, f):
        n2 = mesh.dot(f, f)
        nrm = np.sqrt(n2)
        return nrm

    def normV(mesh, V):
        n2 = mesh.dotV(V, V)
        nrm = np.sqrt(n2)
        return nrm

    def trace_dot(mesh, U, V):
        dg = np.dot(U.T, mesh.A.dot(V))
        dp = 0.5 * np.sum(dg)
        dp2 = 0.5 * mesh.dot(mesh.va, dg)
        return dp, dp2

    def J(mesh, V):
        JV = np.cross(mesh.N, V)
        return JV

    def grad(mesh, f):
        # input function is defined on vertices (nv x 1)
        T = mesh.triangles
        Ar = np.tile(mesh.ta, (1, 3))
        G = (np.tile(f[T[:, 0]], (1, 3)) * mesh.JJc1
            + np.tile(f[T[:, 1]], (1, 3)) * mesh.JJc2
            + np.tile(f[T[:, 2]], (1, 3)) * mesh.JJc3)
        gradF = G / (2 * Ar)
        return gradF

    def curlF(mesh, s):
        grad_s = mesh.grad(s)
        V = -mesh.J(np.cross(mesh.N, grad_s))
        return V

    def curlV(mesh, V):
        X = mesh.vertices
        T = mesh.triangles

        # Compute the edges of each triangle
        C1 = X[T[:, 2]] - X[T[:, 1]]
        C2 = X[T[:, 3]] - X[T[:, 1]]
        C3 = X[T[:, 1]] - X[T[:, 2]]

        # Compute the surface normals
        N = np.cross(C1, C2)
        area = np.sqrt(np.sum(N**2, axis=1))
        N = N / np.tile(area, (3, 1)).T

        # Compute the curl of the vector field
        I = np.concatenate((T[:, 0], T[:, 1], T[:, 2]))
        J = np.ones_like(I)
        S = np.concatenate((np.sum(np.cross(V[T[:, 1]], C1), axis=1),
                            np.sum(np.cross(V[T[:, 2]], C2), axis=1),
                            np.sum(np.cross(V[T[:, 0]], C3), axis=1)))
        curlw = sparse.coo_matrix((S / (2 * np.tile(area, 3))), (I, J), shape=(mesh.nv, 1)).todense()
        curlw = np.array(curlw).ravel()

        return curlw / mesh.va

    def project(mesh, f):
        fp = f - (mesh.LB_basisi[0] @ f) * mesh.LB_basis[:,0]
        return fp

    def laplaceF(mesh, f):
        lf = mesh.L @ f
        return lf

    def laplaceInverseF(mesh, f, guess=None, prec=1e-12):
        if guess is None:
            guess = mesh.last_laplace_inverse
        if mesh.nv < 200:
            lif = sla.spsolve(mesh.L, f)
        else:
            WL = mesh.W
            PC = mesh.W_precon
            b = mesh.project(mesh.A @ f)
            lif, flag = sla.pcg(WL, b, tol=prec, maxiter=500, M=PC, callback=None, atol=None, x0=guess)
            if flag:
                lif = sla.spsolve(mesh.L, f)
        lif = mesh.project(lif)
        mesh.last_laplace_inverse = lif
        return lif

    def mass_matrix_barycentric(mesh):
        Ar = mesh.ta
        Mij = 1/12 * np.concatenate([Ar, Ar, Ar])
        Mji = Mij
        Mii = 1/6 * np.concatenate([Ar, Ar, Ar])
        Mn = np.concatenate([Mij, Mji, Mii])
        In_MM = np.concatenate([mesh.II, mesh.JJ, mesh.II])
        Jn_MM = np.concatenate([mesh.JJ, mesh.II, mesh.II])
        # convert to CSR format to be compatible with most functions in module.
        M = coo_matrix((Mn, (In_MM, Jn_MM)), shape=(mesh.nv, mesh.nv)).tocsr()
        return M

    # Note: The MESH.normv function in Matlab computes the Euclidean norm of
    # each row of a matrix. In the Python code, this is achieved using
    # np.linalg.norm with axis=1.
    def cotLaplacian(mesh):
        X = mesh.vertices
        T = mesh.triangles

        # Find orig edge lengths and angles
        L1 = np.linalg.norm(X[T[:, 1], :] - X[T[:, 2], :], axis=1)
        L2 = np.linalg.norm(X[T[:, 2], :] - X[T[:, 0], :], axis=1)
        L3 = np.linalg.norm(X[T[:, 0], :] - X[T[:, 1], :], axis=1)

        A1 = (L2**2 + L3**2 - L1**2) / (2 * L2 * L3)
        A2 = (L3**2 + L1**2 - L2**2) / (2 * L3 * L1)
        A3 = (L1**2 + L2**2 - L3**2) / (2 * L1 * L2)
        AA = np.vstack((A1, A2, A3)).T
        AA = np.arccos(AA)

        # The Cot Laplacian
        S = 0.5 * np.hstack((1 / np.tan(AA[:, 2]), 1 / np.tan(AA[:, 0]), 1 / np.tan(AA[:, 1]), 1 / np.tan(AA[:, 2])))
        Sn = np.hstack((-S, -S, S, S))
        W = csr_matrix((Sn, (mesh.IIn, mesh.JJn)), shape=(mesh.nv, mesh.nv))

        return W

    def vf2op_initialization(mesh):
        X = mesh.vertices
        T = mesh.triangles
        Nf = mesh.N
        v1 = X[T[:,0],:]
        v2 = X[T[:,1],:]
        v3 = X[T[:,2],:]
        C1 = v3 - v2
        C2 = v1 - v3
        C3 = v2 - v1
        JC1 = np.cross(Nf, C1)
        JC2 = np.cross(Nf, C2)
        JC3 = np.cross(Nf, C3)
        return JC1, JC2, JC3


    def DV(mesh, Vf):
        WW = vf2opWW(mesh, Vf)
        op = mesh.A_INV @ WW
        return op

    def adDV(mesh, Vf):
        WW = vf2opWW(mesh, Vf)
        op = mesh.A_INV @ WW.T
        return op

    # In MATLAB, spdiags is a function that creates a sparse matrix from
    # diagonals provided as input arguments. It is used to construct
    # sparse matrices with known structure or to extract diagonals from
    # an existing matrix.
    def vf2opWW(mesh, Vf):
        Sij = 1/6 * np.vstack((np.sum(mesh.JJc2 * Vf, axis=1),
                            np.sum(mesh.JJc3 * Vf, axis=1),
                            np.sum(mesh.JJc1 * Vf, axis=1)))
        Sji = 1/6 * np.vstack((np.sum(mesh.JJc1 * Vf, axis=1),
                            np.sum(mesh.JJc2 * Vf, axis=1),
                            np.sum(mesh.JJc3 * Vf, axis=1)))
        Sn = np.vstack((Sij, Sji, -Sij, -Sji))
        WW = spdiags(Sn.T, [0, 1, 2], mesh.nv, mesh.nv)
        return WW





# https://raw.githubusercontent.com/azencot/gp.fns/master/code/test_sphere_jet.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/test_sphere_rot.m

# https://raw.githubusercontent.com/azencot/gp.fns/master/code/vortices_func.m

if __name__ == "__main__":
    ps.init()
    with open("./gp.fns/data/sphere_s3.off") as f:
      (vertices, faces) = read_off(f)
      ps.register_surface_mesh("s3", vertices, faces, smooth_shade=True)
    ps.show()

