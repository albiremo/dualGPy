import abc
import meshio
import mypy
import matplotlib.pyplot as plt
import numpy as np
import itertools
import time
from dualGPy import Utils as ut
from dualGPy.Geometry import Face2D,Tetra,Hexa
from dualGPy.Graph import Graph2D
from collections import defaultdict
from itertools import product
from enum import IntEnum

class CellType(IntEnum):
    """Enum for cell types according to how many faces are on the boundary: interior, 0 faces; valley, 1; ridge, 2;
    corner, 3."""
    INTERIOR = 0
    VALLEY   = 1
    RIDGE    = 2
    CORNER   = 3

class Mesh(abc.ABC):
    """ This is a class interface representing a generic Mesh, 2D or 3D. The :class:`Mesh` is a 
        wrapper of the meshio object Mesh, containing methods to compute the area and volume (depending 
        of the dimensionality of the problem (2D or 3D).
        :param mesh: a :mod:`meshio` object  """
    def __init__(self, mesh):
        self.mesh = mesh
        # the Cells are parsed from the meshio object in a simpler structure (TODO: maybe create a dictionary between
        # cells and cell type
        self.cells = []
        self.cell_type = []
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
        self.connectivity = {}



    @abc.abstractmethod
    def ComputeGeometry(self,points):
        """ compute the volume of the cells with respect to the dimensionality (2D area, 3D volume) """
        raise NotImplementedError
    @abc.abstractmethod
    def get_boundary_faces(self):
        """ compute the area of the faces of the cells with respect to the dimensionality (2D length of the segment, 3D area of the face) """
        raise NotImplementedError

    @abc.abstractmethod
    def boundary_detection(self):
        """ Run the detection of the connectivity (it verifies the neightborhood for each cell) and fills the dictionaries Dfaces and faces"""
        raise NotImplementedError

    @abc.abstractmethod
    def setup_mesh(self):
        """ Method that setup the elements to elaborate the graph and hence casts all the cells in a vector starting from a meshio.Mesh object. It fills the cells list of array and the boundary list.
        The method fills also the vector of the center points of the mesh.
        Results:
            - class.cells : cells parsed in a list of list
            - class.cell_type : number of point constituting the cell
            - class.centers : centerpoint of each cell (useful for the post processing and visualizations) """
        raise NotImplementedError

class Mesh2D(Mesh) :
    """ This is a concrete class representing a 2D mesh. The :class:`Mesh2D` derives from the :class:`Mesh` 
      and contains  concrete methods coming from :class:`Mesh`.
       
      The constructor of the class takes a variadic list of arguments. It is possible or to give a :mod:`meshio` object (an external mesh
      given in the datacard or 2 parameters that allows to the constructor to build a :mod:`meshio` object that will represent a 2D mesh
      of squares (number of element in x will be the same of number of elements in y): 

      :param n: number of elements in the 2D mesh 
      :param anisotropic: or true or false it allows to introduce the anisotropicity in one direction
    """ 
    def __init__(self, *args):
       if len(args)==1:
          mesh=args[0]
       if len(args)>1:
           points =[]
           cells_test = []
           n = args[0]
           anisotropic = args[1]
           points = [ [j,i] for i,j in product(range(n),range(n)) ]
           for k,element in enumerate(points):
                if (k+1)%n!=0 and (k < (len(points)-n)) :
                  cells_test.append([k,k+1,n+k+1,n+k])
           cells =[("quad",cells_test)]
           mesh = meshio.Mesh(points,cells)
       super().__init__(mesh)
       # look at this to understand what has been done https://stackoverflow.com/questions/2728346/passing-parameter-to-base-class-constructor-or-using-instance-variable
       self.setup_mesh()
    def setup_mesh(self):
        for i,e in enumerate(self.mesh.cells):
          for j,f in enumerate(self.mesh.cells[i][1]):
           self.cells.append(self.mesh.cells[i][1][j])
           self.cell_type.append(len(self.mesh.cells[i][1][j]))
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




    def draw_graph(self,string):
     """Draw the mesh and the undirected graph
        :param string: string of the name of the figure produced"""
     plt.figure()
    ### cycle on the elements
     for i,elemento in enumerate(self.cells):
      ### cycle on the point of the elements
         print(i)
     #to print numero cell
     #    if i<100:
     #     print("true")
     #     plt.text(self.centers[i][0],self.centers[i][1], i, fontsize = 10)
         for i in range(len(elemento)):
       ## # plot the grid
           x_value = [self.mesh.points[elemento[i-1],0],self.mesh.points[elemento[i],0]]
           y_value = [self.mesh.points[elemento[i-1],1],self.mesh.points[elemento[i],1]]
           plt.plot(x_value,y_value,c='k')
         #else:
         # break
     x1,x2,y1,y2 = plt.axis()
     plt.axis((0.0,0.25,-0.1,0.1))
     ##draw the adjacency graph
     plt.savefig(string)

    def draw_aggl_lines_full(self,string,lines_dict):
     """Draw the mesh and the agglomeration lines that are defined in a dictionary given as an
        input. This version is the one that plot the whole mesh.
        :param string: string of the name of the figure produced
        :param lines_dict: dictionary where the key is the number of the line and the value is a list of the cells representing this line"""
     plt.figure()
     plt.axis('equal')
    ### cycle on the elements
     for i,elemento in enumerate(self.cells):
      ### cycle on the point of the elements
         print(i)
     #to print numero cell
     #    if i<100:
     #     print("true")
     #     plt.text(self.centers[i][0],self.centers[i][1], i, fontsize = 10)
         for i in range(len(elemento)):
       ## # plot the grid
           x_value = [self.mesh.points[elemento[i-1],0],self.mesh.points[elemento[i],0]]
           y_value = [self.mesh.points[elemento[i-1],1],self.mesh.points[elemento[i],1]]
           plt.plot(x_value,y_value,c='b')
         #else:
         # break
     x_line = []
     y_line = []
     for key,value in lines_dict.items():
         for cell in value:
             x_line.append(self.centers[cell][0])
             y_line.append(self.centers[cell][1])
         plt.plot(x_line,y_line,linewidth=2.0,c='r')
         x_line = []
         y_line = []
     ##draw the adjacency graph
     plt.savefig(string)

    def draw_bnd(self,string,vector):
     """Draw the mesh and the priority boundaries (for the agglomeration) 
        :param string: string of the name of the figure produced
        :param vector: vector of the boundaries to draw
        """
     plt.figure()
     plt.axis('equal')
    ### cycle on the elements
     for i,elemento in enumerate(self.cells):
      ### cycle on the point of the elements
         print(i)
     #to print numero cell
     #    if i<100:
     #     print("true")
     #     plt.text(self.centers[i][0],self.centers[i][1], i, fontsize = 10)
         for i in range(len(elemento)):
       ## # plot the grid
           x_value = [self.mesh.points[elemento[i-1],0],self.mesh.points[elemento[i],0]]
           y_value = [self.mesh.points[elemento[i-1],1],self.mesh.points[elemento[i],1]]
           plt.plot(x_value,y_value,c='b')
         #else:
         # break
     for key,value in enumerate(vector):
         x_line = self.centers[value][0]
         y_line = self.centers[value][1]
         plt.plot(x_line,y_line,c='r',marker="x")
     ##draw the adjacency graph
     plt.savefig(string)



    def draw_aggl_lines(self,string,lines_dict,xlim_min,xlim_max,ylim_min,ylim_max):
     """Draw the mesh and the agglomeration lines that are defined in a dictionary given as an
        input. This version is the one that plot a window on the mesh.
        :param string: string of the name of the figure produced
        :param lines_dict: dictionary where the key is the number of the line and the value is a list of the cells representing this line
        :param xlim_min: minimum x to plot
        :param xlim_max: maximum x to plot
        :param ylim_min: minimum y to plot
        :param ylim_max: maximum y to plot"""
     plt.figure()
    ### cycle on the elements
     for i,elemento in enumerate(self.cells):
      ### cycle on the point of the elements
         print(i)
     #to print numero cell
     #    if i<100:
     #     print("true")
     #     plt.text(self.centers[i][0],self.centers[i][1], i, fontsize = 10)
         for i in range(len(elemento)):
       ## # plot the grid
           x_value = [self.mesh.points[elemento[i-1],0],self.mesh.points[elemento[i],0]]
           y_value = [self.mesh.points[elemento[i-1],1],self.mesh.points[elemento[i],1]]
           plt.plot(x_value,y_value,c='b')
         #else:
         # break
     x_line = []
     y_line = []
     for key,value in lines_dict.items():
         for cell in value:
             x_line.append(self.centers[cell][0])
             y_line.append(self.centers[cell][1])
         plt.plot(x_line,y_line,linewidth=2.0,c='r')
         x_line = []
         y_line = []
     x1,x2,y1,y2 = plt.axis()
     plt.axis((xlim_min,xlim_max,ylim_min,ylim_max))
     ##draw the adjacency graph
     plt.savefig(string)


    def ComputeGeometry(self):
        """ In the case of the 2D class it will be an Area """
        # points of the specific cell
        cell_points=[]
        # cycle on the cells
        for i,cell in enumerate(self.cells):
        # cycle on the indexes
            for index in cell:
        # accumulating the cell points
                cell_points.extend(self.mesh.points[index])
        # applying shoelace formula
        # 1. reshape the cell points vector to operate directly with vectors, avoiding unnecessary loops
        # 2. apply the shoelace
            cella = Face2D(self.faces[i],cell_points)
            cella.ComputeArea()
            cella.ComputeLength(self.mesh.points)
            self.volume.append(cella.area)
            self.area.extend(cella.leng_segments)
            cell_points =[]

    def boundary_detection_Easy(self):
      """ automatically determine the boundary condition, starting from the given graph
         Right now we generate all the combination of possible faces and we see if they
         are in the neighborhood. If they are not we print them out."""
      # define the dictionary of boundaries to determine later exactly the boundary based on number of diagonals
      for k,v in self.connectivity.items():
         connections = len(v)

         if (self.cell_type[k]==3):
            num_boundaries = 3-connections
         else: 
            num_boundaries = 4-connections
         if (num_boundaries==CellType.VALLEY): 
            self.onValley.append(k)
            self.boundary_cells.append(np.int(num_boundaries))
         elif ((num_boundaries) == CellType.RIDGE):
            self.onRidge.append(k)
            self.boundary_cells.append(np.int(num_boundaries))
         elif ((num_boundaries) >= CellType.CORNER):
            self.onCorner.append(k)
            self.boundary_cells.append(np.int(num_boundaries))

    def get_boundary_faces(self):
     """ Returns the dual mesh held in a dictionary Graph with dual["points"] giving the coordinates and
     dual["cells"] giving the indicies of all the cells of the dual mesh.
     """
     # Get the first set of points of the dual mesh
     d = ut.dict_of_indices(self.cells)
     # Initializing
     self.faces = {i:[] for i in range(len(self.cells))}
     self.Dfaces = {i:[] for i in range(len(self.cells))}
     self.connectivity = {i:[] for i in range(len(self.cells))}
     # cycle on the points
     for idx in range(len(self.mesh.points)):
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
       # compliant_cells = ut.get_dual_points(self.cells, idx)
        compliant_cells = d[idx]
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
                 self.connectivity[i].append(j)
        print(idx)


    def boundary_detection(self):
     """ automatically determine the boundary condition
         Right now we generate all the combination of possible faces and we see if they
         are in the neighborhood. If they are not we print them out."""
     # define the dictionary of boundaries to determine later exactly the boundary based on number of diagonals
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
        # We check the presence of the inverted faces and we free the list of the boundaries.
        for element in loop_boundary:
           for face_cell in self.faces[i]:
              if sorted(element) == sorted(face_cell):
                 list_inter_boundary.remove(element)
        boundary_dict[i].extend(list_inter_boundary)
        # if is more than number of diagonal i should add it to the boundary cells, because
        # one face is the boundary (I am not currently interested in which are the boundary faces)
        # and build the on boundary vector
        num_diag = len(self.cells[i])*(len(self.cells[i])-3)/2
        num_boundaries = len(boundary_dict[i]) - num_diag
        if (len(boundary_dict[i]) > num_diag) :
            self.boundary_cells.append(np.int(num_boundaries))
            if (num_boundaries == CellType.VALLEY):
                self.onValley.append(i)
            elif (num_boundaries == CellType.RIDGE):
                self.onRidge.append(i)
            elif (num_boundaries >= CellType.CORNER):
                self.onCorner.append(i)
        else:
            self.boundary_cells.append(np.int(0))


class Mesh3D(Mesh):
    """ Implements the 3D mesh: ATTENTION! Right now only tetra supported, but flexible to implement also
        hexa and pyramids""" 
    def __init__(self, *args):
       if len(args)==1:
          mesh=args[0]
       if len(args)>1:
           points = []
           cells_test = []
           n = args[0]
           anisotropic = args[1]
           print(n,anisotropic)
           if anisotropic == False :
            for i in range(n):
             for j in range(n):
              for k in range(n):
                points.append([k,j,i])
            for fila in range(n-1):
             for k in range(n*n):
               if (k+1)%n!=0 and (k < (n*(n-1))) :
                 cells_test.append([k+(fila*n*n),k+1+(fila*n*n),n+k+1+(fila*n*n),n+k+(fila*n*n),k+(n*n)+(fila*n*n),k+1+(n*n)+(fila*n*n),n+k+1+(n*n)+(fila*n*n),n+k+(n*n)+(fila*n*n)])
            cells =[("hexahedron",cells_test)]
            mesh = meshio.Mesh(points,cells)
       super().__init__(mesh)
       # look at this to understand what has been done https://stackoverflow.com/questions/2728346/passing-parameter-to-base-class-constructor-or-using-instance-variable
       self.setup_mesh()

    def setup_mesh(self):
        for i,e in enumerate(self.mesh.cells):
          for j,f in enumerate(self.mesh.cells[i][1]):
           self.cells.append(self.mesh.cells[i][1][j])
        x_centerpoint=0
        y_centerpoint=0
        z_centerpoint=0
        # initialization of the vector of all the centerpoints of the cells
        for elemento in self.cells:
        # cycle on the point of the elements
           for i in range(len(elemento)):
        # compute the centerpoint (accumulating)
             x_centerpoint += self.mesh.points[elemento[i],0]
             y_centerpoint += self.mesh.points[elemento[i],1]
             z_centerpoint += self.mesh.points[elemento[i],2]
        # define the center point
           x_centerpoint/=len(elemento)
           y_centerpoint/=len(elemento)
           z_centerpoint/=len(elemento)
        #  append to the centerpoints vector the element computed
           self.centers.append([x_centerpoint,y_centerpoint,z_centerpoint])
        #  re-initialize the accumulation vectors
           x_centerpoint=0
           y_centerpoint=0
           z_centerpoint=0


    def get_boundary_faces(self):
     """ Returns the dual mesh held in a dictionary Graph with dual["points"] giving the coordinates and
     dual["cells"] giving the indices of all the cells of the dual mesh.
     """
     # Get the first set of points of the dual mesh
     d = ut.dict_of_indices(self.cells)
     # Initializing
     self.faces = {i:[] for i in range(len(self.cells))}
     self.Dfaces = {i:[] for i in range(len(self.cells))}
     self.connectivity = {i:[] for i in range(len(self.cells))}
     # cycle on the points
     for idx in range(len(self.mesh.points)):
        # Get the dual mesh points for a given mesh vertex and the compliant cells to be analysed
        compliant_cells = d[idx]
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
               if ((len(inter)>=3) and (inter not in self.faces[i]) and (inter not in self.faces[j])):
                 self.Dfaces[i].append(inter)
               if ((len(inter)>=3) and (inter not in self.faces[i])):
                 self.faces[i].append(inter)
                 self.connectivity[i].append(j)
        print(idx)


    def ComputeGeometry(self):
        """ In the case of the 2D class it will be an Area """
        # points of the specific cell
        cell_points=[]
        # cycle on the cells
        for i,cell in enumerate(self.cells):
        # cycle on the indexes
            for index in cell:
        # accumulating the cell points
                cell_points.extend(self.mesh.points[index])
        # applying shoelace formula
        # 1. reshape the cell points vector to operate directly with vectors, avoiding unnecessary loops
        # 2. apply the shoelace
            if len(cell_points)==4:
               cella = Tetra(cell_points,self.faces[i])
            else:
               cella= Hexa(cell_points,self.faces[i])
            cella.ComputeArea(self.mesh.points)
            cella.ComputeVolume()
            self.volume.append(cella.volume)
            self.area.extend(cella.AreaFaces)
            cell_points =[]

    def boundary_detection(self):
     """ automatically determine the boundary condition
         Right now we generate all the combination of possible faces and we see if they
         are in the neighborhood. If they are not we print them out."""
     # define the dictionary of boundaries to determine later exactly the boundary based on number of diagonals
     boundary_dict = {}
     for i in range(len(self.cells)):
        # initialize boundary dict
        boundary_dict.update({i :[]})
        # solution from https://stackoverflow.com/questions/69618239/missing-couple-of-elements-in-a-vector
        combination = itertools.combinations(self.cells[i],3)
        inter_boundary = set(combination).difference(map(tuple,self.faces[i]))
        list_inter_boundary = list(map(list,inter_boundary))
        loop_boundary = list_inter_boundary.copy()
        # We create a copy of the list with copy method because we cannot remove elements from a list we are looping
        # https://stackoverflow.com/questions/14126726/python-throws-valueerror-list-removex-x-not-in-list
        # We check the presence of the inverted faces and we free the list of the boundaries.
        for element in loop_boundary:
           for face_cell in self.faces[i]:
              if sorted(element) == sorted(face_cell):
                 list_inter_boundary.remove(element)
        boundary_dict[i].extend(list_inter_boundary)
        # In case of tetra we do not have diagonals, hence given 3 points
        # we define exactly the boundary
        num_boundaries = len(boundary_dict[i])
        self.boundary_cells.append(np.int(num_boundaries))
        if (num_boundaries == CellType.VALLEY):
            self.onValley.append(i)
        elif (num_boundaries == CellType.RIDGE):
            self.onRidge.append(i)
        elif (num_boundaries >= CellType.CORNER):
            self.onCorner.append(i)


