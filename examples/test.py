import meshio
from dualGPy import Mesh as dm
from dualGPy import Graph as graph
import numpy as np
import networkx as nx

# EXAMPLE 1: CREATE A MESH
# 1. Set up points and cells
# 2. Use them to get a meshio mesh
# 3. Use meshio mesh to build a dualGPy mesh

points = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [2.0, 1.0],
    [3.0, 0.0],
    [3.0, 1.0],
    [3.0, 2.0],
    [2.0, 2.0],
    [1.0, 2.0],
    [0.0, 2.0],
    [0.0, -1.0],
    [1.0, -1.0],
    [2.0, -1.0],
    [3.0, -1.0],
    [4.0, -1.0],
    [4.0,0.0],
    [4.0,1.0],
    [4.0,2.0],
    [4.0,3.0],
    [3.0,3.0],
    [2.0,3.0],
    [1.0,3.0],
    [0.0,3.0],
    [0.0,-2.0],
    [1.0,-2.0],
    [2.0,-2.0],
    [3.0,-2.0],
    [4.0,-2.0]
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3],[4,6,7,5],[5,7,8,9],[3,5,9,10],[3,10,11,2],[14,15,6,4],[13,14,4,1],[12,13,1,0],[15,16,17,6],[6,17,18,7],[7,18,19,8],[8,19,20,21],[9,8,21,22],[10,9,22,23],[11,10,23,24],[25,26,13,12],[26,27,14,13],[27,28,15,14],[28,29,16,15]]),
]

mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    # Each item in cell data must match the cells array
)

Mesh1 = dm.Mesh2D(mesh)
Mesh1.get_boundary_faces()
Mesh1.ComputeGeometry()
Mesh1.boundary_detection()

# EXAMPLE 2: CREATE AND EXPORT A GRAPH
# To create a graph, just pass a dualGPy mesh.
# MIND: one should have called get_boundary_faces beforehand in order to set up the connectivity
Graph1 = graph.Graph2D(Mesh1)
Graph1.get_CSR()
# Graph1.adj_to_csr()
# Graph inherits also adj_to_csr, but it basically does the same thing as get_CSR

Graph1.draw_graph("finale.png")
