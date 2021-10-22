import abc
import mypy
import numpy as np
import networkx as nx
from Utils import get_dual_points
class Graph(abc.ABC): 
    def __init__(self, mesh):
        self.mesh = mesh
        self.graph={}
        self.adj = []
        self.edges = []
        self.vertex = []
        self.weight = []
    
    @abc.abstractmethod
    def get_adj_matrix(self):
        raise NotImplementedError

    def generate_graph(self) -> nx.Graph:
     """ Generating a networkx graph element starting from the graph dictionary previously generated
     and the numpy respective adjacency matrix"""
     
     # we use a directed graph because we do not want the bi-directed references 
     # in the adjacency matrix
     g = nx.DiGraph(self.graph)
     # numpy adjacency matrix
     numg = nx.to_numpy_array(g, nodelist=range(len(self.mesh.cells)))
     return g , numg
    
    def adj_to_csr(self):
     """Parsing the adjacency matrix in numpy in a CSR (Compressed Sparse Row) representation
     Parameters:

     adj : Adjacency matrix

     """
     # initialization vertex and edges
     # TODO : weight
     self.vertex = np.zeros(len(adj))

     # the CSR representation is built starting from the definition, hence
     # we cycle over the adjacency matrix, we identify the non zero entry
     # and finally we fill the vertex vector with the position where we find the 
     # non zero elements in the edges vector. We hence fill the edges vector at the same
     # way

     for i in range(len(self.adj[:,1])):         
       self.vertex[i] =  np.sum(np.count_nonzero(self.adj,axis=1)[:i])
       for j in range(len(self.adj[1,:])):
           if adj[i,j]!=0:
            self.edges.append(j)


class Graph2D(Graph):
    def __init__(self, mesh):
        super().__init__(mesh)

    def get_adj_matrix(self):
     # Get the first set of points of the dual mesh
     for i in range(len(self.mesh.cells)):
        self.graph.update({i :[]})
     print(self.graph)
     # cycle on the points
     for idx in range(len(self.mesh.mesh.points)):
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        compliant_cells = get_dual_points(self.mesh.cells, idx)
        # in this part we build the graph: for each point of the mesh we have the compliant cells
        # and we cycle over the compliant cells (two nested loop, with an if that avoids to inspect the same cell)
        # me create the list inter that check the common point between two vectors (that can have also different 
        # dimension, considering that they can represent cells of completely different shape. 
        # checked that we have more than two vertex in common (WE ARE IN 2D HERE), and that the node is not already
        # connected with the analysed cell, we add it to the respective dictionary key.
        for i in compliant_cells:
          for j in compliant_cells:
             if i!=j:
               inter = list(set(self.mesh.cells[i]).intersection(self.mesh.cells[j]))
        # the last statement avoid the bi-directed graph
               if ((len(inter)>=2) and (j not in self.graph[i]) and (i not in self.graph[j])):
                 self.graph[i].append(j)     


