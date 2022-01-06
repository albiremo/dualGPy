from CoMMA import *
import meshio
from dualGPy.CGNS_Adapter import CGNS_Adapter
from dualGPy.Graph import Graph2D
from dualGPy.Mesh import Mesh2D
import pandas as pd

name = "rae.cgns"
adapter = CGNS_Adapter(name)
points,cells =  adapter.prepare_cgns_for_meshio(0)


mesh = meshio.Mesh(
    points,
    cells,
)

mesh.write("input.vtk")
Mesh1 = Mesh2D(mesh)
print("setup done")
Mesh1.get_boundary_faces()
Dfaces = pd.DataFrame.from_dict(Mesh1.Dfaces, orient="index")
Dfaces.to_csv("Dfaces.csv")
faces = pd.DataFrame.from_dict(Mesh1.faces, orient="index")
faces.tp_csv("faces.csv")
print("area computed")
Mesh1.ComputeVolume()
print("volume comuted")
Mesh1.ComputeArea()
rint("bf computed")
Mesh1.boundary_detection()
print("boundary detection performed")
# graph
Graph1 = Graph2D(Mesh1)
Graph1.get_adj_matrix()
print("adjacency matrix written")
adjacency = np.asarrray(Graph1.adj)
np.savetxt('adj.csv', adjacency, delimiter=',')
Graph1.adj_to_csr()
edges = np.asarrray(Graph1.edges)
vertex = np.asarrray(Graph1.vertex)
np.savetxt('edges.csv', edges, delimiter=',')
np.savetxt('vertex.csv', vertex, delimiter=',')
print("shift to CSR done")

