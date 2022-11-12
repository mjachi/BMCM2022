from dolfin import *
import numpy as np

import matplotlib.pyplot as plt


c=5000
mesh = UnitSquareMesh(8, 8)
#mesh = RectangleMesh(-2, -2, 2, 2,80,80)
V=FunctionSpace(mesh, "Lagrange", 1)

# Time variables
dt = 0.00004; t = 0; T = 0.004

# Previous and current solution
u1= interpolate(Constant(0.0), V)
u0= interpolate(Constant(0.0), V)

# Variational problem at each time
u = TrialFunction(V)
v = TestFunction(V)

a = u*v*dx + dt*dt*c*c*inner(grad(u), grad(v))*dx
L = 2*u1*v*dx-u0*v*dx

bc = DirichletBC(V, 0, "on_boundary")
A, b = assemble_system(a, L, bc)

u=Function(V)
while t <= T:
    A, b = assemble_system(a, L, bc)
    #delta = PointSource(V, Point(0.5, 0.5), sin(c * 10 * t))
    delta = PointSource(V, Point(0.5, 0.5), (1 - 2 * (np.pi * 40 * t) ** 2 * np.exp(-(np.pi * 40 * t) ** 2)))
    delta.apply(b)
    solve(A, u.vector(), b)
    u0.assign(u1)
    u1.assign(u)
    t += dt

    # Reduce the range of the solution so that we can see the waves
    j = 0
    for i in u.vector():
        i = min(.01, i)
        i = max(-.01, i)
        u.vector()[j] = i;
        j += 1

    plot(u, interactive=False)

plot(u, interactive=True)
plt.show() 