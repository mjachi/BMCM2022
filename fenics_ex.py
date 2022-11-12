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


def awefem(mesh, t, V, source_loc=None):

    # Function space
    CGFS = FunctionSpace(mesh, 'CG', 1)
    DGFS = FunctionSpace(mesh, 'DG', 0)

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

    resarray = []

    # Time stepping
    for i, t_ in enumerate(t[1:]):
        b = assemble(L)
        delta = PointSource(V, source_loc, sine_source(t_) * dt**2)
        delta.apply(b)
        solver.solve(u.vector(), b)

        u0.assign(u1)
        u1.assign(u)

        resarray.append(u.vector().get_local())

        #resarray.append(interpolate(project(sqrt(inner(u,u)), V), DGFS).vector()[:])

        #if t_ in [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0]:
        #    plot(
        #        u,
        #        vmin=.0,     # sets a minimum to the color scale
        #        vmax=0.003,
        #        cmap='rainbow', # the color map style
        #        alpha=1,        # transparency of the mesh
        #        title="Approximation at time step {:03f}".format(t_)
        #    )  # continue execution
        #plt.show()

    return np.array(resarray)

if __name__ == "__main__":

    ot, dt, nt = 0.0, 1e-3, 100
    t = ot + np.arange(nt) * dt

    mesh = UnitSquareMesh(50, 50, "crossed")

    V = FunctionSpace(mesh, "Lagrange", 1)

    res = awefem(mesh, t, V, source_loc=(0.5, 0.1))
    #element_averages = []
#
    #for element in range(len(res[0])): 
    #    elt_avg = 0 
    #    for time_step in range(len(t) - 1): 
    #        elt_avg += res[time_step][element]
    #    
    #    elt_avg = elt_avg / len(t)
    #    element_averages.append(elt_avg)

    dof_averages = np.mean(res, axis=0)
    n = V.dim() 
    d = mesh.geometry().dim() 
    dof_coordinates = V.tabulate_dof_coordinates() 
    dof_coordinates.resize((n, d)) 
    dof_x = dof_coordinates[:, 0]
    dof_y = dof_coordinates[:, 1]

    print(res.shape)
    print(len(dof_coordinates))
    print(len(dof_averages))

    fig = plt.figure() 
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(dof_x, dof_y, dof_averages, cmap='virdis', marker='.', alpha=0.5, )
    ax.plot_trisurf(dof_x, dof_y, dof_averages, triangles=dof_coordinates)
    plt.title("Unit Volume Average Energies")

    plt.show() 

    
