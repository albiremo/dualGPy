import meshio
import dualGPy as dm
# For plotting both the mesh and dual mesh
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
# two triangles and one quad
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
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3],[4,6,7,5],[5,7,8,9],[3,5,9,10],[3,10,11,2],[14,15,6,4],[13,14,4,1],[12,13,1,0]]),
]

mesh = meshio.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    # Each item in cell data must match the cells array
    cell_sets = {
    "Boundary": {"line": [[0, 2],[2,11],[11,10],[10,9],[9,8],[8,7],[7,6],[15,6],[14,15],[13,14],[12,13],[12,0]]},
}
)

Mesh1 = dm.Mesh(mesh)
Mesh1.get_dual()
graph , numarray = Mesh1.generate_graph()
print(numarray)
## initialization of the figure
plt.figure()
## cycle on the elements
for elemento in Mesh1.cells:
## cycle on the point of the elements
 for i in range(len(elemento)):
# # plot the grid
  x_value = [mesh.points[elemento[i-1],0],mesh.points[elemento[i],0]]
  y_value = [mesh.points[elemento[i-1],1],mesh.points[elemento[i],1]]
  plt.plot(x_value,y_value)  
#draw the adjacency graph
nx.draw(graph,pos=Mesh1.centers, with_labels=True)
plt.savefig('foo.png')
Mesh1.adj_to_csr(numarray)
