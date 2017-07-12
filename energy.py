from dolfin import *
from block import block_mat, block_vec, block_transpose
from block.iterative import ConjGrad
from block.algebraic.petsc import AMG
from rigid_motions import first
import rigid_motions


def energy(lmbda, mu, f, h, mesh, Z=None):
    '''
    Solves 

        -div(sigma) = f in  Omega
            sigma.n = h on  boundary

    where sigma(u) = 2*mu*eps(u) + lambda*div(u)*I. The problem is reformulated by
    considering a complemented strain energy (for which rigid motions are not
    tranparent). The system to be solved with CG as

        P*[A+E]*[u] = P'*b

    with P a precondtioner. We run on series of meshes to show mesh independence
    of the solver.
    '''
    if not isinstance(mesh, Mesh):
        # Precompute the 'symbolic' basis
        mesh0 = first(mesh)
        Z = rigid_motions.rm_basis(mesh0)
        return [energy(lmbda, mu, f, h, mesh_, Z) for mesh_ in mesh]
    
    # For cube
    V = VectorFunctionSpace(mesh, 'CG', 1)
    u, v = TrialFunction(V), TestFunction(V)
    # Strain
    epsilon = lambda u: sym(grad(u))
    # Stress
    gdim = mesh.geometry().dim()
    sigma = lambda u: 2*mu*epsilon(u) + lmbda*tr(epsilon(u))*Identity(gdim)
    # Energy of elastic deformation
    a = inner(sigma(u), epsilon(v))*dx
    A = assemble(a)
    # Mass matrix for B
    m = inner(u, v)*dx
    M = assemble(m)

    # NOTE: Avoiding use of Q space in the assembly - dense blocks!
    Q = VectorFunctionSpace(mesh, 'R', 0, dim=6)
    Zh = rigid_motions.RMBasis(V, Q, Z)  # L^2 orthogonal
    B = M*Zh

    # System operator
    AA = A + B*block_transpose(B)

    # Right hand side
    L = inner(f, v)*dx + inner(h, v)*ds
    # Orthogonalize
    P = rigid_motions.Projector(Zh)
    b = assemble(L)
    b0 = block_transpose(P)*b

    # Preconditioner
    AM = assemble(a + m)
    BB = AMG(AM)

    # Solve, using random initial guess
    x0 = AA.create_vec()
    as_backend_type(x0).vec().setRandom()
    
    AAinv = ConjGrad(AA, precond=BB, initial_guess=x0, maxiter=100, tolerance=1E-8,
                     show=2, relativeconv=True)

    x = AAinv*b0
    
    # Functions from coefficients
    # uh = Function(V, x)  # Displacement

    niters = len(AAinv.residuals) - 1
    assert niters < 100

    P*x  # to get orthogonality
    if MPI.rank(mesh.mpi_comm()) == 0:
        print '\033[1;37;31m%s\033[0m' % ('Orthogonality %g' % max(P.alphas))

    return V.dim(), niters


def test_energy():
    '''Number of iterations should not blow up'''
    lmbda = Constant(1)
    mu = Constant(1)
    f = Expression(('A*sin(2*x[0])', 'A*cos(3*(x[0]+x[1]+x[2]))', 'A*sin(x[2])'),
                    degree=3, A=0.01)
    h = Constant((0, 0, 0))

    comm = mpi_comm_world().tompi4py()

    Ns = [2, 4, 8, 16, 32]
    if comm.size > 2:
        Ns.extend([64, 128])

    meshes = (BoxMesh(Point(1, 1, 1), Point(2, 1.5, 1.25), N, N, N) for N in Ns)

    converged = energy(lmbda, mu, f, h, meshes)
    assert all(converged)

    # Dump data for plotting
    if comm.rank == 0:
        from numpy import savetxt, array
        savetxt('./.energy.txt', array(converged), fmt=['%d', '%d'])

    return True

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    set_log_level(PROGRESS)
    assert test_energy()
