import meshio
from dualGPy.Mesh import Mesh2D
from dualGPy.Graph import Graph2D
import numpy as np

points = []
cells_test = []
n=5
for i in range(n):
 for j in range(n):
  points.append([j,i])
for k,element in enumerate(points):
  if (k+1)%n!=0 and (k < (len(points)-n)) :
   cells_test.append([k,k+1,n+k+1,n+k])
print(len(cells_test))
cells =[("quad",cells_test)]
mesh = meshio.Mesh(
    points,
    cells,
)
Mesh1 = Mesh2D(mesh)
Mesh1.get_boundary_faces()
Mesh1.ComputeVolume()
Mesh1.ComputeArea()
Mesh1.boundary_detection()
# graph
Graph1 = Graph2D(Mesh1)
Graph1.get_CSR()
print(Mesh1.connectivity)
