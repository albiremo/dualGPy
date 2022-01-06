import abc 
import meshio
import mypy
import numpy as np
import itertools
from dualGPy import Utils as ut 

class Mesh(abc.ABC):
    """ Interface class to compute all the geometrical characteristics of the mesh """
    def __init__(self, mesh):
        self.mesh = mesh
        self.cells = []  
        # non directional faces (in graph sense)
        self.faces = {}
        # Directional faces in graph sense
        self.Dfaces ={}
        self.centers = []
        self.volume = []
        self.area = []
        self.onValley  = []
        self.onRidge  = []
        self.onCorner  = []
        self.boundary_cells = []

    @abc.abstractmethod
    def ComputeVolume(self,points):
        
        """ compute the volume of the cells with respect to the dimensionality (2D area, 3D volume) """
        raise NotImplementedError
    @abc.abstractmethod
    def ComputeArea(self,points):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError
    @abc.abstractmethod
    def get_boundary_faces(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError

    @abc.abstractmethod
    def boundary_detection(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError


    def setup_mesh(self):
        """ Method that setup the elements to elaborate the graph and hence casts all the cells in a vector starting from a meshio.Mesh object. It fills the cells list of array and the boundary list. 
        The method fills also the vector of the center points of the mesh.
        Results:
            - class.boundary faces
            - class.cells
            - class.centers """
        for i,e in enumerate(self.mesh.cells):
          for j,f in enumerate(self.mesh.cells[i][1]):
           self.cells.append(self.mesh.cells[i][1][j])
     # to activate only if automatic detection is not active
     # to parse the boundaries we make use of the dictionary iterators
     #   for key, value in self.mesh.cell_sets.items():
     #     for key1,value1 in self.mesh.cell_sets[key].items():
     #        for i in range(len(value1)):
     #          self.boundary_faces.append(value1[i])
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

class Mesh2D(Mesh) :
    
    def __init__(self, mesh):
        super().__init__(mesh)
        self.setup_mesh()

    def ComputeVolume(self):
        """ In the case of the 2D class it will be an Area """
        # points of the specific cell
        cell_points=[]
        # cycle on the cells
        for cell in self.cells:
        # cycle on the indexes
            for index in cell:
        # accumulating the cell points
                cell_points.extend(self.mesh.points[index])
        # applying shoelace formula
        # 1. reshape the cell points vector to operate directly with vectors, avoiding unnecessary loops
        # 2. apply the shoelace
            cell_points_reshaped=np.reshape(cell_points,(len(cell),2))
            shifted = np.roll(cell_points_reshaped, 1, axis=0)
            volume = 0.5 * np.sum((shifted[:, 0] + cell_points_reshaped[:, 0])*(shifted[:, 1] - cell_points_reshaped[:, 1]))
            self.volume.append(abs(volume))
            cell_points = []

    def ComputeArea(self):
     """ Compute the area that in case of the 2D is the length of the segment associated with the faces of the cell. As in the graph the segmen    t are not repeated for different cells but are considered 
     """ 
     for key, value in self.Dfaces.items():
         for segment in value:
             leng = np.sqrt((self.mesh.points[segment[1]][1]-self.mesh.points[segment[0]][1])**2+(self.mesh.points[segment[1]][0]-self.mesh.points[segment[0]][0])**2)
             self.area.append(leng)
 
    def get_boundary_faces(self):
     """ Returns the dual mesh held in a dictionary Graph with dual["points"] giving the coordinates and
     dual["cells"] giving the indicies of all the cells of the dual mesh.
     """
     # Get the first set of points of the dual mesh
     compliant_cells = ut.get_dual_points(self.cells, 0)    
     for i in range(len(self.cells)):
        self.faces.update({i :[]}) 
        self.Dfaces.update({i :[]}) 
     # cycle on the points
     for idx in range(1, len(self.mesh.points)):
        print(idx)
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        compliant_cells = ut.get_dual_points(self.cells, idx)
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
        # in the faces part we have to associate with each cell all the faces
        # like in a bi-directed graph
               if ((len(inter)>=2) and (inter not in self.faces[i]) and (inter not in self.faces[j])):
                 self.Dfaces[i].append(inter) 
               if ((len(inter)>=2) and (inter not in self.faces[i])):
                 self.faces[i].append(inter)


    def boundary_detection(self):
     """ automatically determine the boundary condition
         Right now we generate all the combination of possible faces and we see if they
         are in the neightborhood. If they are not we print them out."""
     # define the dictionary of boundaries to determine later exactly the boundaty based on number of diagonals
     boundary_dict = {}
     for i in range(len(self.cells)):
        # initialize boundary dict
        boundary_dict.update({i :[]})
        # solution from https://stackoverflow.com/questions/69618239/missing-couple-of-elements-in-a-vector
        combination = itertools.combinations(self.cells[i],2)
        inter_boundary = set(combination).difference(map(tuple,self.faces[i]))
        list_inter_boundary = list(map(list,inter_boundary))
        loop_boundary = list_inter_boundary.copy()
        # We create a copy of the list with copy method because we cannot remove elements from a list we are looping
        # https://stackoverflow.com/questions/14126726/python-throws-valueerror-list-removex-x-not-in-list   
        # We check the presence of the iverted faces and we free the list of the boundaries.
        for element in loop_boundary:
           for face_cell in self.faces[i]:
             if all(np.flip(element) == face_cell):
                 list_inter_boundary.remove(element)
        boundary_dict[i].extend(list_inter_boundary)
        # if is more than number of diagonal i should add it to the boundary cells, because
        # one face is the boundary (I am not currently interested in which are the boundary faces)
        # and build the on boundary vector
        num_diag = len(self.cells[i])*(len(self.cells[i])-3)/2
        num_boundaries = len(boundary_dict[i]) - num_diag
        if (len(boundary_dict[i]) > num_diag) : 
            self.boundary_cells.append(np.int(num_boundaries))
            if (num_boundaries == 1):
                self.onValley.append(i)
            elif (num_boundaries == 2):
                self.onRidge.append(i)
            elif (num_boundaries>=3):
                self.onCorner.append(i)
        else:
            self.boundary_cells.append(np.int(0))
        # Alternative way to define the boundary cells    
#----------------------------------------------------------------
#     self.boundary_cells = []       
#     for i in range(len(self.cells)):
#     # We cycle on the boundary cells marked
#         for j in self.boundary_faces:
#            inter = list(set(self.cells[i]).intersection(j))
#            if ((len(inter)>=2) and (i not in self.boundary_cells)):
#              self.boundary_cells.append(i) gg
#---------------------------------------------------------------

