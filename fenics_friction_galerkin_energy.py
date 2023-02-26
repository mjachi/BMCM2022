from dolfin import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

def ricker_wavelet(t, sigma=40): 
    return 2 / (np.sqrt(3 * sigma) * np.pi ** (1/4)) * (1 - (t / sigma)**2) * exp(- (t / (2 * sigma **2))**2)

def sine_source(t, f=40):
    return 10 * np.sin(2 * np.pi*f*t)

def fourier(term_count, t, L = 1/15): 
    modes = list(range(1, 2*term_count, 2))
    expansion = 0 
    for k in modes: 
        expansion += 1 / k * np.sin(k * np.pi * t / L)
    return 10 * 4 / np.pi * expansion 

def wefem_multisource(mesh, t, source_locs): 
    # Function space
    V = FunctionSpace(mesh, "Lagrange", 1)

    # NOTE: FeniCS automatically specifies a 0 Neumann boundary condition
    #bc = DirichletBC(V, Constant(0), "on_boundary")

    # Trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Discretization
    c = 1000
    k = 0.99995
    dt = t[1] - t[0]
    u0 = Function(V)  # u0 = uN-1
    u1 = Function(V)  # u1 = uN1

    # Variational formulation
    #F = (u - 2 * u1 + u0) * v * dx + (dt * c) ** 2 * dot(
    #    grad(u + 2 * u1 + u0) / 4, grad(v) ) * dx
    #F = (u  - (2 + k * dt) / (k + dt) * u1 + u0) * v * dx + (c * dt) ** 2 * dot(grad(u - (2 + k * dt) / (k + dt) * u1 + u0) / 4, grad(v)) * dx 
    F = (u - ( 2 + k * dt) / (k + dt) * u1 + u0) * v * dx + (c * dt) ** 2 * dot(grad(u + 2 * u1 + u0) / 4, grad(v)) * dx 
    a, L = lhs(F), rhs(F)

    # Solver
    A, b = assemble_system(a, L)
    solver = LUSolver(A, "mumps")
    solver.parameters["symmetric"] = True
    #bc.apply(A, b)

    # Solution
    u = Function(V)  # uN+1

    source_points = []

    for source_loc in source_locs: 
        source_points.append(Point(source_loc))

    res_array = []

    for _, t_ in enumerate(t[1:]):
        b = assemble(L)
        if t_ < 0.1: 
            for source_point in source_points: 
                delta = PointSource(V, source_point, 10000 * sine_source(t_) * dt **2)
            #delta = PointSource(V, source_point, fourier(4, t_) * dt**2)
            #delta = PointSource(V, source_loc, fourier(3, t_) * dt**2)
            
                delta.apply(b)
        
        solver.solve(u.vector(), b)

        u0.assign(u1)
        u1.assign(u)

        res_array.append(u.vector().get_local())

        #if t_ in [0.05, 0.1, 0.15, 0.25, 0.35, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 2.75, 3.0, 4.75]:
        #    c = plot(
        #        u,
        #        vmin=.0,     
        #        vmax=0.003,
        #        cmap='rainbow',
        #        alpha=0.5,        
        #        title="Approximation with Friction at Time Step {:03f}".format(t_)
        #    )  
        #    plt.colorbar(c) 
        #plt.show()

    return np.array(res_array)

if __name__ == "__main__":
    ot, dt, nt = 0.0, 1e-3, 150
    t = ot + np.arange(nt) * dt

    #mesh = Mesh("triangle.xml")
    #mesh = UnitSquareMesh(100, 100, "crossed")
    #mesh = Mesh("circle_mesh.xml")
    #mesh = Mesh("amphi2D.xml")
    #mesh = Mesh("ellipse.xml")
    #mesh = Mesh("rect.xml")
    mesh = Mesh("penta2.xml")
    V = FunctionSpace(mesh, "Lagrange", 1)

    #awefem(mesh, t, source_loc=(0., 0.1))
    #wefem_multisource(mesh, t, source_locs=[ (0.25, 0.25), (-0.25, 0.25), (0, -0.25)])

    #source_points = [(0.52345127713513175, 0.2869721966085012), (0.541656158894040896, 0.4591626583411702), (0.251656158894040896, 0.1991626583411702), (0.1, 0.1), (0.99, 0)]
    source_points = []
    for idx in np.arange(3): 
        source_points.append((np.random.uniform(0, 1/2), np.random.uniform(0,1)))

    print(source_points)
    #source_points = [(0., 0.), (0.1, 0.1), (0.25, 0.1), (1/2, 1/3), ]

    avg_collection = []
    for source_point in source_points:
        res = wefem_multisource(mesh,t, source_locs=[source_point])

        avg_collection.append(np.mean(res, axis=0) ** 2)

    avg_collection = np.array(avg_collection)
    print(avg_collection.shape)

    E0 = 100 * sine_source(0) ** 2 * np.ones(avg_collection[0].shape)

    E_mean = np.mean(avg_collection, axis=0)


    print(E_mean)
    print(np.mean(E_mean))
    print(np.linalg.norm(E_mean - E0))
    

    #n = V.dim() 
    #d = mesh.geometry().dim() 
    #dof_coordinates = V.tabulate_dof_coordinates() 
    #dof_coordinates.resize((n, d)) 
    #dof_x = dof_coordinates[:, 0]
    #dof_y = dof_coordinates[:, 1]



    #fig = plt.figure() 

    #ax = fig.add_subplot(111, projection='3d')

    #trisurf = ax.plot_trisurf(dof_x, dof_y, 10000 * dof_averages **2, triangles=dof_coordinates, cmap=cm.coolwarm)

    #plt.title(r"Averages $\mathcal{E}_{avg}[u](x, t)$ for Full Time")
    #plt.show() 
