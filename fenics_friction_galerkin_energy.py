from dolfin import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

set_log_level(30)

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

def awefem(mesh, t, V, source_loc=None):
    # Trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Discretization
    c = 343.9
    k = 0.1
    dt = t[1] - t[0]
    u0 = Function(V)  # u0 = uN-1
    u1 = Function(V)  # u1 = uN1

    # Variational formulation
    #F = (u - 2 * u1 + u0) * v * dx + (dt * c) ** 2 * dot(
    #    grad(u + 2 * u1 + u0) / 4, grad(v) ) * dx
    F = (u  - (2 + k * dt) / (k + dt) * u1 - u0) * v * dx + (c * dt) ** 2 * dot(grad(u1), grad(v)) * dx 
    a, L = lhs(F), rhs(F)

    # Solver
    A, b = assemble_system(a, L)
    solver = LUSolver(A, "mumps")
    solver.parameters["symmetric"] = True

    # Solution
    u = Function(V)  # uN+1

    resarray = []

    # Time stepping
    for i, t_ in enumerate(t[1:]):
        b = assemble(L)
        #delta = PointSource(V, source_loc, fourier(1, t_) * dt**2)
        delta = PointSource(V, source_loc, sine_source(t_) * dt **2)
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

    ot, dt, nt = 0.0, 1e-3, 1000
    t = ot + np.arange(nt) * dt

    mesh = Mesh("amphi2D.xml")
    V = FunctionSpace(mesh, "Lagrange", 1)

    res = awefem(mesh, t, V, source_loc=(0.0, 0.1))

    #res = awefem(mesh, t, V, source_loc=(0.5, 0.1))
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

    phys_const = 98.21 ** 2 / (2 * 1.163 * 343.9)

    fig = plt.figure() 
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(dof_x, dof_y, dof_averages, cmap='virdis', marker='.', alpha=0.5, )
    trisurf = ax.plot_trisurf(dof_x, dof_y, phys_const + 0.5 * 1.163 * dof_averages ** 2 / (dt * nt), triangles=dof_coordinates, cmap=cm.coolwarm)
    #x, Vel = CV.legend_elements(str_format=lambda x: f'{x:.1f}')
    #plt.legend(x, Vel)
    #fig.colorbar(trisurf, ax =ax, shrink=0.5, aspect=5, format="%0.2f")
    plt.title(r"Averages $\mathcal{E}_{avg}[u](x, t)$ for Full Time")

    plt.show() 