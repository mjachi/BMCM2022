from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

def ricker_wavelet(t, sigma=40): 
    return 2 / (np.sqrt(3 * sigma) * np.pi ** (1/4)) * (1 - (t / sigma)**2) * exp(- (t / (2 * sigma **2))**2)

def sine_source(t, f=40):
    return np.sin(2 * np.pi*f*t)

def fourier(term_count, t, L = 1/15): 
    modes = list(range(1, 2*term_count, 2))
    expansion = 0 
    for k in modes: 
        expansion += 1 / k * np.sin(k * np.pi * t / L)
    return 4 / np.pi * expansion 

def wefem_multisource(mesh, t, source_locs): 
    # Function space
    V = FunctionSpace(mesh, "Lagrange", 1)

    # NOTE: FeniCS automatically specifies a 0 Neumann boundary condition
    bc = DirichletBC(V, Constant(0), "on_boundary")

    # Trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Discretization
    c = 3
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

    # Solution
    u = Function(V)  # uN+1

    source_points = []

    for source_loc in source_locs: 
        source_points.append(Point(source_loc))

    for _, t_ in enumerate(t[1:]):
        b = assemble(L)
        if t_ < 0.1: 
            for source_point in source_points: 
                delta = PointSource(V, source_point, sine_source(t_) * dt**2)
            #delta = PointSource(V, source_loc, fourier(3, t_) * dt**2)
            
                delta.apply(b)
        
        solver.solve(u.vector(), b)

        u0.assign(u1)
        u1.assign(u)

        if t_ in [0.05, 0.1, 0.15, 0.2, 0.25, 0.33, 0.5, 0.75]:
            plot(
                u,
                vmin=.0,     
                vmax=0.003,
                cmap='rainbow',
                alpha=0.5,        
                title="Approximation without Friction at Time Step {:03f}".format(t_)
            )  
        plt.show()


def awefem(mesh, t, source_loc):
    # Function space
    V = FunctionSpace(mesh, "Lagrange", 1)

    # NOTE: FeniCS automatically specifies a 0 Neumann boundary condition
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

    # Solution
    u = Function(V)  # uN+1

    source_loc = Point(source_loc)

    # Time stepping
    for _, t_ in enumerate(t[1:]):
        b = assemble(L)
        if t_ < 0.1: 
            delta = PointSource(V, source_loc, sine_source(t_) * dt**2)
            #delta = PointSource(V, source_loc, fourier(3, t_) * dt**2)
            
            delta.apply(b)
        solver.solve(u.vector(), b)

        u0.assign(u1)
        u1.assign(u)

        if t_ in [0.05, 0.1, 0.15, 0.2, 0.25, 0.33, 0.5]:
            plot(
                u,
                vmin=.0,     
                vmax=0.003,
                cmap='rainbow',
                alpha=0.5,        
                title="Approximation without Friction at Time Step {:03f}".format(t_)
            )  
        plt.show()

if __name__ == "__main__":
    ot, dt, nt = 0.0, 1e-3, 1000
    t = ot + np.arange(nt) * dt

    #mesh = UnitSquareMesh(200, 200, "crossed")
    #mesh = Mesh("circle_mesh.xml")
    mesh = Mesh("amphi2D.xml")

    #awefem(mesh, t, source_loc=(0., 0.1))
    wefem_multisource(mesh, t, source_locs=[(0.1, 0.1), (0, 0.1), (-0.1, 0.1)])
