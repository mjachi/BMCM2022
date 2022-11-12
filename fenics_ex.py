"""
Solve the constant velocity scalar wave equation in an arbitrary number of dimensions.
It injects a point source with a time-dependent source time function.
"""
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

set_log_level(30)

def ricker_wavelet(t, sigma=40): 
    return 2 / (np.sqrt(3 * sigma) * np.pi ** (1/4)) * (1 - (t / sigma)**2) * exp(- (t / (2 * sigma **2))**2)

def ricker_source(t, f=40):
    t -= 2 / f
    return (1 - 2 * (np.pi*f*t)**2) * np.exp(-(np.pi*f*t)**2)

def sine_source(t, f=40):
    return np.sin(2 * np.pi*f*t)


def awefem(mesh, t, source_loc=None):

    # Function space
    V = FunctionSpace(mesh, "Lagrange", 1)

    # Boundary condition
    bc = DirichletBC(V, Constant(0), "on_boundary")

    # Trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Discretization
    c = 6
    dt = t[1] - t[0]
    u0 = Function(V)  # u0 = uN-1
    u1 = Function(V)  # u1 = uN1

    # Variational formulation
    F = (u - 2 * u1 + u0) * v * dx + (dt * c) ** 2 * dot(
        grad(u + 2 * u1 + u0) / 4, grad(v) ) * dx
    a, L = lhs(F), rhs(F)

    # Solver
    A, b = assemble_system(a, L)
    solver = LUSolver(A, "mumps")
    solver.parameters["symmetric"] = True
    bc.apply(A, b)

    # Solution
    u = Function(V)  # uN+1

    # Source
    if source_loc is None:
        mesh_center = np.mean(mesh.coordinates(), axis=0)
        source_loc = Point(mesh_center)
    else:
        source_loc = Point(source_loc)

    # Time stepping
    for i, t_ in enumerate(t[1:]):
        b = assemble(L)
        delta = PointSource(V, source_loc, sine_source(t_) * dt**2)
        delta.apply(b)
        solver.solve(u.vector(), b)

        u0.assign(u1)
        u1.assign(u)

        if t_ in [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0]:
            plot(
                u,
                vmin=.0,     # sets a minimum to the color scale
                vmax=0.003,
                cmap='rainbow', # the color map style
                alpha=1,        # transparency of the mesh
                title="Approximation at time step {:03f}".format(t_)
            )  # continue execution
        plt.show()

if __name__ == "__main__":

    ot, dt, nt = 0.0, 1e-3, 1000
    t = ot + np.arange(nt) * dt

    mesh = UnitSquareMesh(200, 200, "crossed")
    #plot(mesh, title="Finite Element Mesh")

    awefem(mesh, t, source_loc=(0.5, 0.1))
