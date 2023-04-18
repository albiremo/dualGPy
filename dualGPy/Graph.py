import abc
import mypy
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from dualGPy import Utils as ut
class Graph(abc.ABC):
    """Base class for the Graph2D and Graph3D representation.
     WARNING! 2D and 3D are considered from a geometric point of view:

     mesh : dualGPy mesh representation
    """
    def __init__(self, mesh):
     self.mesh = mesh
     self.graph={}
     self.nx_graph = nx.Graph(self.mesh.connectivity)
     self.edges = []
     self.vertex = []
     self.weight = []

    @abc.abstractmethod
    def get_CSR(self):
     raise NotImplementedError

    def adj_to_csr(self):
     """Build the adjacency matrix via `networkx` and parse it in a CSR (Compressed Sparse Row) representation"""
     # the CSR representation is built starting from the definition, hence
     # we cycle over the adjacency matrix, we identify the non zero entry
     # and finally we fill the vertex vector with the position where we find the
     # non zero elements in the edges vector. We hence fill the edges vector at the same
     # way

     # numpy adjacency matrix
     adj = nx.to_numpy_array(self.nx_graph, nodelist=range(len(self.mesh.cells)))
     non_zero = np.count_nonzero(adj,axis=1)
     for i in range(len(adj[:,1])):
       self.vertex.append(np.sum(non_zero[:i]))
       for j in range(len(adj[i,:])):
           if adj[i,j]!=0:
            self.edges.append(j)
     self.vertex.append(np.sum(np.count_nonzero(adj)))


class Graph2D(Graph):
    def __init__(self, mesh):
        super().__init__(mesh)

    def get_CSR(self):
     """Build a CSR representation of the graph"""
     self.vertex.append(0)
     somma = 0
     # Initialize the keys of the graph dictionary
     # cycle on the points
     for key in self.mesh.connectivity:
        self.edges.extend(self.mesh.connectivity[key])
        somma += len(self.mesh.connectivity[key])
        self.vertex.append(somma)



    def draw_graph(self, string, mesh = None):
     """Draw the mesh and the graph.
Parameters:
* string: str
  Name of the file in which the graph will be printed
* mesh: dualGPy mesh, default: None
  mesh whose graph will be drawn. If None, the mesh used to set up the graph is used
"""
     m = self.mesh if mesh is None else mesh
     plt.figure()
     ### cycle on the elements
     for elemento in m.cells:
      ### cycle on the point of the elements
         for i in range(len(elemento)):
       ## # plot the grid
          x_value = [m.mesh.points[elemento[i-1],0],m.mesh.points[elemento[i],0]]
          y_value = [m.mesh.points[elemento[i-1],1],m.mesh.points[elemento[i],1]]
          plt.plot(x_value,y_value,c='r')
     ##draw the adjacency graph
     nx.draw(self.nx_graph,pos=m.centers,with_labels = True)
     plt.savefig(string)

