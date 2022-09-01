import abc
import numpy as np
from numba import njit, prange
import dualGPy.Utils as ut
import itertools

class Face(abc.ABC):
    def __init__(self,points):
        self.points = points
    @property
    def area(self):
        return self._area

    @area.setter
    def area(self, new_area):
        if (new_area>0):
            self._area = new_area
        else:
            print("Error: the area of this solid is 0")
    @area.deleter
    def area(self):
        del self._area
    @property
    def n_vertices(self):
        return self._n_vertices
    @n_vertices.setter
    def n_vertices(self, new_n_vertices):
        if (new_n_vertices > 0):
           self._n_vertices = new_n_vertices
        else:
           print("Error: the n_segments of a face is always > 0")
    @n_vertices.deleter
    def n_vertices(self):
        del self._n_vertices
    @abc.abstractmethod
    def ComputeArea(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError

class Face2D(Face) :
    def __init__(self, segments,*args, **kwargs):
        self.segments = segments
        self.leng_segments = []
        super().__init__(*args, **kwargs)
        self.n_vertices = self.points.shape[0]

    @staticmethod
    @njit
    def algebric_area(points,shifted):
        # TODO look if we should shift all to the father: alias to have a global method ComputeArea and re-structure the concept of interface base-class to a simple abstract class. The static method has been implemented to try to use Numba.
        area_0 = 0.5 * np.sum((shifted[:, 0] + points[:, 0])*(shifted[:, 1] - points[:, 1]))
        return area_0

    def ComputeArea(self):
        shifted = np.roll(self.points, 1, axis=0)
        area_0 = self.algebric_area(self.points,shifted)
        self.area = abs(area_0)

    def ComputeLength(self,global_points):
        for segment in self.segments:
             leng = np.sqrt((global_points[segment[1]][1]-global_points[segment[0]][1])**2+(global_points[segment[1]][0]-global_points[segment[0]][0])**2)
             self.leng_segments.append(leng)

class Face3D(Face) :
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_vertices = self.points.shape[0]
        self.a = self.points[0,:]
        self.b = self.points[1,:]
        self.c = self.points[2,:]

    def unit_normal(self):
        x = np.linalg.det([[1,self.a[1],self.a[2]],
             [1,self.b[1],self.b[2]],
             [1,self.c[1],self.c[2]]])
        y = np.linalg.det([[self.a[0],1,self.a[2]],
             [self.b[0],1,self.b[2]],
             [self.c[0],1,self.c[2]]])
        z = np.linalg.det([[self.a[0],self.a[1],1],
             [self.b[0],self.b[1],1],
             [self.c[0],self.c[1],1]])
        magnitude = (x**2 + y**2 + z**2)**.5
        #if (magnitude==0):
        #   magnitude = 1e-3
        return ([x/magnitude, y/magnitude, z/magnitude])

    def ComputeArea(self):
        #https://stackoverflow.com/questions/12642256/find-area-of-polygon-from-xyz-coordinates
        #shape (N, 3)
        poly = self.points
        #all edges
        edges = poly[1:] - poly[0:1]
        # row wise cross product
        cross_product = np.cross(edges[:-1],edges[1:], axis=1)
        #area of all triangles
        area = np.linalg.norm(cross_product, axis=1)/2
        self.area = sum(area)

class Solid(abc.ABC):
    """ Interface class to compute the different characteristics of an element of
#        a 3D mesh. It takes as an input the global points and build in the constructor the 
         cell points. """
    def __init__(self,c_points,g_points,Faces):
        self.Faces = Faces 
        self.global_points = g_points
        self.cell_points = c_points
    # TODO: - reconstruct the cell_points from global points deleting an argument            
    # -we can think to set it as a property
    #  look at https://stackoverflow.com/questions/37564798/python-property-on-a-list
        self.AreaFaces = []
        self.n_faces = int(len(self.Faces))
    @property
    def n_vertices(self):
        return self._n_vertices
    @n_vertices.setter
    def n_vertices(self, new_n_vertices):
        if (new_n_vertices > 0):
           self._n_vertices = new_n_vertices
        else:
           print("Error: the n_segments of a face is always > 0")
    @n_vertices.deleter
    def n_vertices(self):
        del self._n_vertices

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, new_volume):
        if new_volume > 0:
            self._volume = new_volume
        else:
            print("Error: the volume of this solid is 0")
    @volume.deleter
    def volume(self):
        del self._volume
    @property
    def n_faces(self):
        return self._n_faces
    @n_faces.setter
    def n_faces(self, new_n_faces):
        if (new_n_faces >= 0):
           self._n_vertices = new_n_faces
        else:
           print("Error: the n_faces of a face is always > 0")
    @n_faces.deleter
    def n_faces(self):
        del self._n_faces

    @abc.abstractmethod
    def ComputeArea(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError

    @abc.abstractmethod
    def ComputeVolume(self):
        """ compute the Volume of the faces of the cell with respect to the dimensionality """
        raise NotImplementedError


class Tetra(Solid):
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_vertices = self.cell_points.shape[0]
        assert(self.n_vertices==4)
    def ComputeArea(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        if self.Faces:
         for faccia in self.Faces:
             points = self.global_points[faccia, :]
             faccia_el = Face3D(points)
             faccia_el.ComputeArea()
             self.AreaFaces.append(faccia_el.area)
    def ComputeVolume(self):
        """ compute the Volume of  of the cells with respect to the dimensionality """
        # https://stackoverflow.com/questions/9866452/calculate-volume-of-any-tetrahedron-given-4-points
        mat_1 = self.cell_points.transpose()
        mat_2 = np.vstack([mat_1,np.ones((1,4))])
        self.volume= 1/6*abs(np.linalg.det(mat_2))

class Hexa(Solid):
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.edges = []
        self.n_vertices = self.cell_points.shape[0]
        assert(self.n_vertices==8)
    def ComputeArea(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        if self.Faces:
         for faccia in self.Faces:
             points = self.global_points[faccia, :]
             faccia_el = Face3D(points)
             faccia_el.ComputeArea()
             self.AreaFaces.append(faccia_el.area)
    def RetriveEdges(self):       
        """ Retrive de edges of the Hexa cell analyzed"""
        d = ut.dict_of_indices(self.Faces)
        for key,value in d.items():
             for i,j in itertools.combinations(value,2):
#            for i in value:
#                for j in value:
                    # if i!=j: 
                       inter = list(set(self.Faces[i]).intersection(self.Faces[j]))
                       if ((len(inter)>=2) and (inter not in self.edges)):
                          self.edges.append(inter)
    def ComputeVolume(self):
        """ compute the Volume of  of the cells with respect to the dimensionality """
        prodotto = []
        self.RetriveEdges()
        d = ut.dict_of_indices(self.edges)
        for key,value in d.items():
            if len(value)==3:
               for it,index in enumerate(value):
                   segment = self.edges[index]
                   p1 = self.global_points[segment[0]]
                   p2 = self.global_points[segment[1]]
                   prodotto.append(p2-p1)
               self.volume = abs(np.dot(prodotto[0],np.cross(prodotto[1],prodotto[2])))
               break

class Wedge(Solid):
    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self.edges = []
        self.n_vertices = self.cell_points.shape[0]
        assert(self.n_vertices==6)
    def ComputeArea(self):
        """ compute the Area of the faces of the cell with respect to the dimensionality """
        if self.Faces:
         for faccia in self.Faces:
             points = self.global_points[faccia, :]
             faccia_el = Face3D(points)
             faccia_el.ComputeArea()
             self.AreaFaces.append(faccia_el.area)
    def RetriveEdges(self):       
        """ Retrive de edges of the Hexa cell analyzed"""
        d = ut.dict_of_indices(self.Faces)
        for key,value in d.items():
             for i,j in itertools.combinations(value,2):
#            for i in value:
#                for j in value:
                    # if i!=j: 
                       inter = list(set(self.Faces[i]).intersection(self.Faces[j]))
                       if ((len(inter)>=2) and (inter not in self.edges)):
                          self.edges.append(inter)
    def ComputeVolume(self):
        """ compute the Volume of  of the cells with respect to the dimensionality. The volume of the
            wedge is considered to be the half of the volume of the corresponding Hexa """
        prodotto = []
        self.RetriveEdges()
        d = ut.dict_of_indices(self.edges)
        for key,value in d.items():
            if len(value)==3:
               for it,index in enumerate(value):
                   segment = self.edges[index]
                   p1 = self.global_points[segment[0]]
                   p2 = self.global_points[segment[1]]
                   prodotto.append(p2-p1)
               self.volume = abs(np.dot(prodotto[0],np.cross(prodotto[1],prodotto[2])))/2
               break

