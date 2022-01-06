import abc
import mypy
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from dualGPy import Utils as ut
class Graph(abc.ABC): 
    """Base class for the Graph2D and Graph3D representation.
     WARNING! 2D and 3D are considered from a geometric point of
     view:

     mesh : Mesh representation from meshio library
    """
    def __init__(self, mesh):
     self.mesh = mesh
     self.graph={}
     self.nx_graph = nx.Graph(self.graph)
     self.adj = []
     self.edges = []
     self.vertex = []
     self.weight = []
    
    @abc.abstractmethod
    def get_adj_matrix(self):
     raise NotImplementedError

    def adj_to_csr(self):
     """Parsing the adjacency matrix in numpy in a CSR (Compressed Sparse Row) representation
     Parameters:

     adj : Adjacency matrix

     """
     # the CSR representation is built starting from the definition, hence
     # we cycle over the adjacency matrix, we identify the non zero entry
     # and finally we fill the vertex vector with the position where we find the 
     # non zero elements in the edges vector. We hence fill the edges vector at the same
     # way
#     for i in range(len(self.adj[:,1])): 
#       self.vertex.append(np.sum(np.count_nonzero(self.adj,axis=1)[:i]))
     non_zero = np.count_nonzero(self.adj,axis=1)
     for i,e in enumerate(self.adj[:,1]):
       print(i) 
       self.vertex.append(np.sum(non_zero[:i]))
#       edges = [ j for j,f in enumerate(self.adj[i,:]) if self.adj[i,j]!=0] 
       for j,f in enumerate(self.adj[i,:]):
           if self.adj[i,j]!=0:
            self.edges.append(j)
     self.vertex.append(np.sum(np.count_nonzero(self.adj)))


class Graph2D(Graph):
    def __init__(self, mesh):
        super().__init__(mesh)

    def get_adj_matrix(self):
     # Initialize the keys of the graph dictionary
     for i,e in enumerate(self.mesh.cells):
        self.graph.update({i :[]})
     # cycle on the points
     for idx,e in enumerate(self.mesh.mesh.points):
        print(idx)
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        compliant_cells = ut.get_dual_points(self.mesh.cells, idx)
        # in this part we build the graph: for each point of the mesh we have the compliant cells
        # and we cycle over the compliant cells (two nested loop, with an if that avoids to inspect the same cell)
        # me create the list inter that check the common point between two vectors (that can have also different 
        # dimension, considering that they can represent cells of completely different shape. 
        # checked that we have more than two vertex in common (WE ARE IN 2D HERE), and that the node is not already
        # connected with the analysed cell, we add it to the respective dictionary key.
        for i in compliant_cells:
          for j in compliant_cells:
             if i!=j:
        # we use the intersect method of the lists
               inter = list(set(self.mesh.cells[i]).intersection(self.mesh.cells[j]))
        # we define a bidirected graph
               if ((len(inter)>=2) and (j not in self.graph[i])):
                 self.graph[i].append(j)     

     # we use a directed graph because we do not want the bi-directed references 
     # in the adjacency matrix
     self.nx_graph = nx.Graph(self.graph)
     # numpy adjacency matrix
     self.adj = nx.to_numpy_array(self.nx_graph, nodelist=range(len(self.mesh.cells)))
    def draw_graph(self,Mesh1,string):
     """Draw the mesh and the graph"""
     plt.figure()
     ### cycle on the elements
     for elemento in Mesh1.cells:
      ### cycle on the point of the elements
         for i in range(len(elemento)):
       ## # plot the grid
          x_value = [Mesh1.mesh.points[elemento[i-1],0],Mesh1.mesh.points[elemento[i],0]]
          y_value = [Mesh1.mesh.points[elemento[i-1],1],Mesh1.mesh.points[elemento[i],1]]
          plt.plot(x_value,y_value)  
     ##draw the adjacency graph
     nx.draw(self.nx_graph,pos=Mesh1.centers, with_labels=True)
     plt.savefig(string)

