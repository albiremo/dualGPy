import meshio
import mypy
import numpy as np
import networkx as nx
import itertools

class Mesh() :
    def __init__(self, mesh):
        self.mesh = mesh
        self.graph = {}
        self.cells = []  
        self.boundary_cells = []
        self.boundary_faces = []
        self.centers = []
        #TODO: implement in setup_mesh the method to fill it
        self.volume = []
        self.area = []

    def setup_mesh(self):
        """ Method that setup the elements to elaborate the graph and hence casts all the cells in a vector starting from a meshio.Mesh object. It fills the cells list of array and the boundary list. 
        The method fills also the vector of the center points of the mesh.
        Results:
            - class.boundary faces
            - class.cells
            - class.centers """
        for i in range(len(self.mesh.cells)):
          for j in range(len(self.mesh.cells[i][1])):
           self.cells.append(self.mesh.cells[i][1][j])
        # to parse the boundaries we make use of the dictionary iterators
        for key, value in self.mesh.cell_sets.items():
          for key1,value1 in self.mesh.cell_sets[key].items():
             for i in range(len(value1)):
               self.boundary_faces.append(value1[i])
        x_centerpoint=0
        y_centerpoint=0
        # initialization of the vector of all the centerpoints of the cells
        for elemento in self.cells:
           # cycle on the point of the elements
           for i in range(len(elemento)):
           # compute the centerpoint (accumulating)
             x_centerpoint += self.mesh.points[elemento[i],0]
             y_centerpoint += self.mesh.points[elemento[i],1]
             # define the center point
           x_centerpoint/=len(elemento)
           y_centerpoint/=len(elemento)
        #  append to the centerpoints vector the element computed
           self.centers.append([x_centerpoint,y_centerpoint])
        #  re-initialize the accumulation vectors
           x_centerpoint=0
           y_centerpoint=0
        # fill the area vector by cycling on the cells and use the 
        # chain rule by calling the static method ger_area
        # TODO: area computation by chain rule not working properly
#        points_cells = []
#        for cellula in self.cells:
#            for j in range(len(cellula)):
#                points_cells.append(self.mesh.points[cellula[j]])
#            print(points_cells)
#            self.area.append(self.get_area(points_cells))

            # re-initialize point cells
            
         #   point_cells = []

    @staticmethod
    def get_area(points) -> np.float64:
      """ Static method that returns the area of a polygon given the vertices.
         Parameters:
         points:     numpy.ndarray
            Vertices of the polygon
         Warning: the points need to be ordered clockwise or anticlockwise."""

      # shift all the points by one
        
      shifted = np.roll(points, 1, axis=0)
      print(points[:,0])
        # Use the shoelace formula

      area = 0.5 * np.sum((shifted[:, 0] + points[:, 0])*(shifted[:, 1] - points[:, 1]))
      return np.abs(area)

    @staticmethod    
    def get_dual_points(compliant_cells : int, index : int) -> int:
      """ Static method that returns the points of the dual mesh nearest to the point in the mesh given by the index.
          Parameters:
          mesh:       meshio.Mesh object
              Input mesh.
          index:      int
              Index of the point in the input mesh for which to calculate the nearest points of the dual mesh.
          Return:
              compliant 
        # Find the cells where the given index appears, REMEMBER the where statement gives you immediately back the indexof the compliant cell
        # building the compliant cells list
      """       
      compliant=[]
      for i in range(len(compliant_cells)):
      # compress with the use of any to determine the compliant cells
         if any(compliant_cells[i]==index):
           compliant.append(i) 
      # Find the centers of all the cells
      return compliant
    
    def get_dual(self):
     """ Returns the dual mesh held in a dictionary Graph with dual["points"] giving the coordinates and
     dual["cells"] giving the indicies of all the cells of the dual mesh.
     """
     self.setup_mesh()
     # Get the first set of points of the dual mesh
     compliant_cells = Mesh.get_dual_points(self.cells, 0)    
     # define a new key of the dictionary for each cell of the mesh
     faces = {}
     for i in range(len(self.cells)):
        self.graph.update({i :[]})
        faces.update({i :[]})
     # cycle on the points
     for idx in range(1, len(self.mesh.points)):
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        compliant_cells = Mesh.get_dual_points(self.cells, idx)
        # in this part we build the graph: for each point of the mesh we have the compliant cells
        # and we cycle over the compliant cells (two nested loop, with an if that avoids to inspect the same cell)
        # me create the list inter that check the common point between two vectors (that can have also different 
        # dimension, considering that they can represent cells of completely different shape. 
        # checked that we have more than two vertex in common (WE ARE IN 2D HERE), and that the node is not already
        # connected with the analysed cell, we add it to the respective dictionary key.
        for i in compliant_cells:
          for j in compliant_cells:
             if i!=j:
               inter = list(set(self.cells[i]).intersection(self.cells[j]))
               # the last statement avoid the bi-directed graph
               if ((len(inter)>=2) and (j not in self.graph[i]) and (i not in self.graph[j])):
                 self.graph[i].append(j)     
               # in the faces part we have to associate with each cell all the faces
               # like in a bi-directed graph
               if ((len(inter)>=2) and (inter not in faces[i])):
                 faces[i].append(inter)
     boundary_test = []  
# TODO Experimental part to automatically determine the boundary condition
# Right now we generate all the combination of possible faces and we see if they
# are in the neightborhood. If they are not we print them out. Two issues identified
# with this strategy:
# - The diagonal are not considered as a neightborhood and they would be considered as boundary
#   for this reason we need a way to exclude them

     for i in range(len(self.cells)):
     # solution from https://stackoverflow.com/questions/69618239/missing-couple-of-elements-in-a-vector
        combination = itertools.combinations(self.cells[i],2)
        inter_boundary = set(combination).difference(map(tuple,faces[i]))
        print(i) 
        print("faces:", faces[i])
        list_inter_boundary = list(inter_boundary)
        loop_boundary = list_inter_boundary.copy()
        # We create a copy of the list with [:] because we cannot remove elements from a list we are looping
        # https://stackoverflow.com/questions/14126726/python-throws-valueerror-list-removex-x-not-in-list   
        # We check the presence of the iverted faces and we free the list of the boundaries.
        for element in loop_boundary:
           for face_cell in faces[i]:
             if all(np.flip(element) == face_cell):
                 list_inter_boundary.remove(element)
        print("cleaned list:", list_inter_boundary)

     # Define the boundary cells in a similar way we did for the graph
     self.boundary_cells = []       
     for i in range(len(self.cells)):
     # We cycle on the boundary cells marked
         for j in self.boundary_faces:
            inter = list(set(self.cells[i]).intersection(j))
            if ((len(inter)>=2) and (i not in self.boundary_cells)):
              self.boundary_cells.append(i) 
     print(self.boundary_cells)     

    def generate_graph(self) -> nx.Graph:
     """ Generating a networkx graph element starting from the graph dictionary previously generated
     and the numpy respective adjacency matrix"""
     
     # we use a directed graph because we do not want the bi-directed references 
     # in the adjacency matrix
     g = nx.DiGraph(self.graph)
     # numpy adjacency matrix
     numg = nx.to_numpy_array(g, nodelist=range(len(self.cells)))
     return g , numg
    
    @staticmethod
    def adj_to_csr(adj : int):
     """Parsing the adjacency matrix in numpy in a CSR (Compressed Sparse Row) representation
     Parameters:

     adj : Adjacency matrix

     """
     # initialization vertex and edges
     # TODO : weight
     vertex = np.zeros(len(adj))
     edges = []
     weight = []

     # the CSR representation is built starting from the definition, hence
     # we cycle over the adjacency matrix, we identify the non zero entry
     # and finally we fill the vertex vector with the position where we find the 
     # non zero elements in the edges vector. We hence fill the edges vector at the same
     # way

     for i in range(len(adj[:,1])):         
       vertex[i] =  np.sum(np.count_nonzero(adj,axis=1)[:i])
       for j in range(len(adj[1,:])):
           if adj[i,j]!=0:
            edges.append(j)
     return vertex, edges



