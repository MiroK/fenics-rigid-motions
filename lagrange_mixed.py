from dolfin import *
from block import block_mat, block_vec, block_transpose
from block.iterative import MinRes
from block.algebraic.petsc import AMG
import rigid_motions

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"


def lagrange_mixed(lmbda, mu, f, h, mesh, Z=None):
    '''
    Solves 

        -div(sigma) = f in  Omega
            sigma.n = h on  boundary

    where sigma(u) = 2*mu*eps(u) + lambda*div(u)*I. The problem is reformulated by
    Lagrange multiplier nu to inforce orthogonality with the space of rigid
    motions. To get robustnes in lmbda solid pressure p = lambda*div u is introduced. 
    The system to be solved with MinRes is 

        P*[A C  B; *[u, = P*[L,
           C' D 0;   p,      0,
           B' 0 0]   nu]     0]

    with P a precondtioner. We run on series of meshes to show mesh independence
    of the solver.
    '''
    if not isinstance(mesh, Mesh):
        # NOTE: You can precompute the 'symbolic' basis and pass it here
        return [lagrange_mixed(lmbda, mu, f, h, mesh_) for mesh_ in mesh]
    
    # For cube
    V = VectorFunctionSpace(mesh, 'CG', 2)
    Q = FunctionSpace(mesh, 'CG', 1)

    u, v = TrialFunction(V), TestFunction(V)
    p, q = TrialFunction(Q), TestFunction(Q)

    # Strain
    epsilon = lambda u: sym(grad(u))
    # Stress
    gdim = mesh.geometry().dim()
    sigma = lambda u: 2*mu*epsilon(u) + lmbda*tr(epsilon(u))*Identity(gdim)
    a = 2*mu*inner(sym(grad(u)), sym(grad(v)))*dx

    A = assemble(a)

    c = inner(div(v), p)*dx
    C = assemble(c)

    d = -(inner(p, q)/Constant(lmbda))*dx
    D = assemble(d)

    m = inner(u, v)*dx
    M = assemble(m)

    # NOTE: Avoiding use of Q space in the assembly - dense blocks!
    X = VectorFunctionSpace(mesh, 'R', 0, dim=6)
    Zh = rigid_motions.RMBasis(V, X, Z)  # L^2 orthogonal
    B = M*Zh

    # System operator
    AA = block_mat([[A,                  C, B], 
                    [block_transpose(C), D, 0],
                    [block_transpose(B), 0, 0]])

    # Right hand side
    L = inner(f, v)*dx + inner(h, v)*ds
    b0 = assemble(L)
    b1 = assemble(inner(Constant(0), q)*dx)
    # Equivalent to assemble(inner(Constant((0, )*6), q)*dx) but cheaper
    b2 = Function(X).vector()
    bb = block_vec([b0, b1, b2])

    # Block diagonal preconditioner
    IV = assemble(a + m)
    IQ = assemble(inner(p, q)*dx)
    IX = rigid_motions.identity_matrix(X)
    BB = block_mat([[AMG(IV), 0,              0], 
                    [0,       AMG(IQ),        0],
                    [0,       0,            IX]])

    # Solve, using random initial guess
    x0 = AA.create_vec()
    [as_backend_type(xi).vec().setRandom() for xi in x0]

    AAinv = MinRes(AA, precond=BB, initial_guess=x0, maxiter=120, tolerance=1E-8,
                   show=2, relativeconv=True)

    x = AAinv*bb

    # # Functions from coefficients
    # # uh = Function(V, x[0])     # Displacement
    # # ph = Function(Q, x[1])     # Solid pressure
    # # nuh = Zh.rigid_motion(x[2])  # Function in V

    niters = len(AAinv.residuals) - 1
    assert niters < 120

    P = rigid_motions.Projector(Zh)

    P*x0[0]  # to get orthogonality

    if MPI.rank(mesh.mpi_comm()) == 0:
        print '\033[1;37;31m%s\033[0m' % ('Orthogonality %g' % max(P.alphas))
        pass

    return V.dim() + Q.dim() + 6, niters


def test_lagrange_mixed():
    '''Number of iterations should not blow up'''
    from numpy import savetxt, array, c_

    mu = Constant(1)
    f = Expression(('A*sin(2*x[0])', 'A*cos(3*(x[0]+x[1]+x[2]))', 'A*sin(x[2])'),
                    degree=3, A=0.01)
    h = Constant((0, 0, 0))

    comm = mpi_comm_world().tompi4py()

    Ns = [2, 4, 8]
    if comm.size > 3:
        Ns.extend([16, 32, 64])

    data = []
    lmbdas = 1E12, 1E8, 1E4, 1E0
    for lmbda in lmbdas:
        meshes = (BoxMesh(Point(1, 1, 1), Point(2, 1.5, 1.25), N, N, N) for N in Ns)
        converged = lagrange_mixed(Constant(lmbda), mu, f, h, meshes)
        assert all(converged)

        if len(data) == 0:
            data = array(converged)
        else:
            data = c_[data, array(converged)[:, -1]]

    # Dump data for plotting
    if comm.rank == 0:
        header = ' '.join(map(str, lmbdas))
        savetxt('./.lagrange_mixed.txt', data, fmt=['%d']*data.shape[1],
                header=header)

    return True

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    set_log_level(PROGRESS)
    assert test_lagrange_mixed()
