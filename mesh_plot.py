from dolfin import *
import numpy as np
import matplotlib.pyplot as plt


mesh = UnitSquareMesh(10, 10, "crossed")

plot(mesh, title="Finite Element Mesh")
plt.show()
