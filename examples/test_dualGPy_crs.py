import meshio
from dualGPy.Mesh import Mesh2D
from dualGPy.Graph import Graph2D
import numpy as np

n=5
Mesh1 = Mesh2D(n,False)
Mesh1.get_boundary_faces()
Mesh1.ComputeGeometry()
Mesh1.boundary_detection()
# graph
Graph1 = Graph2D(Mesh1)
Graph1.get_CSR()
print(Mesh1.connectivity)
